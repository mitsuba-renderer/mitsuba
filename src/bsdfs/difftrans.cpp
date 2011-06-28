/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/consttexture.h>
#include <mitsuba/hw/renderer.h>

MTS_NAMESPACE_BEGIN


/*! \plugin{difftrans}{Diffuse transmitter}
 *
 * \parameters{
 *     \lastparameter{transmittance}{\Spectrum\Or\Texture}{
 *       Specifies the diffuse transmittance / albedo of the material
 *       \default{0.5}
 *     }
 * }
 *
 * \renderings{
 *     \rendering{Homogeneous reflectance, see \lstref{lambertian-uniform}}{bsdf_lambertian_plain}
 *     \rendering{Textured reflectance, see \lstref{lambertian-textured}}{bsdf_lambertian_textured}
 * }
 *
 * The Lambertian material represents an ideally diffuse material
 * with a user-specified amount of reflectance. Apart from a 
 * homogeneous reflectance value, the plugin can also accept a nested 
 * or referenced texture map to be used as the source of reflectance 
 * information, which is then mapped onto the shape based on its UV
 * parameterization.
 * When no parameters are specified, the model uses the default 
 * of 50% reflectance.
 *
 * \begin{xml}[caption=Reflectance specified as an sRGB color, label=lst:lambertian-uniform]
 * <bsdf type="lambertian">
 *     <srgb name="reflectance" value="#6d7185"/>
 * </bsdf>
 * \end{xml}
 * \begin{xml}[caption=Lambertian material with a texture map, label=lst:lambertian-textured]
 * <bsdf type="lambertian">
 *     <texture type="ldrtexture" name="reflectance">
 *       <string name="filename" value="wood.jpg"/>
 *     </texture>
 * </bsdf>
 * \end{xml}
 */

class DiffuseTransmitter : public BSDF {
public:
	DiffuseTransmitter(const Properties &props) 
		: BSDF(props) {
		m_transmittance = new ConstantTexture(
			props.getSpectrum("transmittance", Spectrum(.5f)));
		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_combinedType = m_type[0] = EDiffuseTransmission | EFrontSide | EBackSide;
		m_usesRayDifferentials = false;
	}

	DiffuseTransmitter(Stream *stream, InstanceManager *manager) 
		: BSDF(stream, manager) {
		m_transmittance = static_cast<Texture *>(manager->getInstance(stream));
		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_combinedType = m_type[0] = EDiffuseTransmission | EFrontSide | EBackSide;
		m_usesRayDifferentials = m_transmittance->usesRayDifferentials();
	}

	virtual ~DiffuseTransmitter() {
		delete[] m_type;
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return Spectrum(0.0f);
	}

	Spectrum f(const BSDFQueryRecord &bRec) const {
		if (!(bRec.typeMask & m_combinedType)
			|| bRec.wi.z*bRec.wo.z >= 0)
			return Spectrum(0.0f);

		return m_transmittance->getValue(bRec.its) * INV_PI;
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		if (bRec.wi.z*bRec.wo.z >= 0)
			return 0.0f;
		return std::abs(Frame::cosTheta(bRec.wo)) * INV_PI;
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		if (!(bRec.typeMask & m_combinedType))
			return Spectrum(0.0f);
		bRec.wo = squareToHemispherePSA(sample);
		if (bRec.wi.z > 0)
			bRec.wo.z *= -1;
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseTransmission;
		if (Frame::cosTheta(bRec.wo) == 0)
			return Spectrum(0.0f);
		return m_transmittance->getValue(bRec.its) / std::abs(Frame::cosTheta(bRec.wo));
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (!(bRec.typeMask & m_combinedType)) 
			return Spectrum(0.0f);
		bRec.wo = squareToHemispherePSA(sample);
		if (bRec.wi.z > 0)
			bRec.wo.z *= -1;
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseTransmission;
		pdf = std::abs(Frame::cosTheta(bRec.wo)) * INV_PI;
		if (Frame::cosTheta(bRec.wo) == 0)
			return Spectrum(0.0f);
		return m_transmittance->getValue(bRec.its) * INV_PI;
	}
		
	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "transmittance") {
			m_transmittance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_transmittance->usesRayDifferentials();
		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_transmittance.get());
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "DiffuseTransmitter[" << endl
			<< "  transmittance = " << indent(m_transmittance->toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_transmittance;
};

// ================ Hardware shader implementation ================ 

class DiffuseTransmitterShader : public Shader {
public:
	DiffuseTransmitterShader(Renderer *renderer, const Texture *reflectance) 
		: Shader(renderer, EBSDFShader), m_transmittance(reflectance) {
		m_transmittanceShader = renderer->registerShaderForResource(m_transmittance.get());
	}

	bool isComplete() const {
		return m_transmittanceShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_transmittance.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_transmittanceShader.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (wi.z*wo.z >= 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return " << depNames[0] << "(uv) * 0.31831;" << endl
			<< "}" << endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << evalName << "(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_transmittance;
	ref<Shader> m_transmittanceShader;
};

Shader *DiffuseTransmitter::createShader(Renderer *renderer) const { 
	return new DiffuseTransmitterShader(renderer, m_transmittance.get());
}

MTS_IMPLEMENT_CLASS(DiffuseTransmitterShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(DiffuseTransmitter, false, BSDF)
MTS_EXPORT_PLUGIN(DiffuseTransmitter, "Diffuse transmitter")
MTS_NAMESPACE_END

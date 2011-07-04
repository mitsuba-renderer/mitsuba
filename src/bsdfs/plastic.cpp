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

MTS_NAMESPACE_BEGIN

/*! \plugin{plastic}{Smooth plastic material}
 *
 * \parameters{
 *     \parameter{intIOR}{\Float}{Interior index of refraction \default{1.5046}}
 *     \parameter{extIOR}{\Float}{Exterior index of refraction \default{1.0}}
 *     \parameter{diffuse\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor used to modulate the diffuse reflectance component\default{1.0}}
 *     \parameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor used to modulate the specular reflectance component\default{1.0}}
 * }
 */
class SmoothPlastic : public BSDF {
public:
	SmoothPlastic(const Properties &props) 
			: BSDF(props) {
		/* Specifies the internal index of refraction at the interface */
		m_intIOR = props.getFloat("intIOR", 1.5046f);
		/* Specifies the external index of refraction at the interface */
		m_extIOR = props.getFloat("extIOR", 1);

		m_diffuseReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("diffuseReflectance", Spectrum(0.5f)));
		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDiffuseReflection | EFrontSide;
		m_type[1] = EDeltaTransmission | EFrontSide;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = false;
	}

	SmoothPlastic(Stream *stream, InstanceManager *manager) 
			: BSDF(stream, manager) {
		m_intIOR = stream->readFloat();
		m_extIOR = stream->readFloat();
		m_diffuseReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDiffuseReflection | EFrontSide;
		m_type[1] = EDeltaTransmission | EFrontSide;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = 
			m_diffuseReflectance->usesRayDifferentials() ||
			m_specularReflectance->usesRayDifferentials();
	}

	virtual ~SmoothPlastic() {
		delete[] m_type;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeFloat(m_intIOR);
		stream->writeFloat(m_extIOR);
		manager->serialize(stream, m_diffuseReflectance.get());
		manager->serialize(stream, m_specularReflectance.get());
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "diffuseReflectance") {
			m_diffuseReflectance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_diffuseReflectance->usesRayDifferentials();
		} else if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "specularReflectance") {
			m_specularReflectance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_specularReflectance->usesRayDifferentials();
		} else {
			BSDF::addChild(name, child);
		}
	}

	inline Vector reflect(const Vector &wi) const {
		return Vector(-wi.x, -wi.y, wi.z);
	}

	Spectrum f(const BSDFQueryRecord &bRec) const {
		bool sampleDiffuse = (bRec.typeMask & EDiffuseReflection)
			&& (bRec.component == -1 || bRec.component == 0);

		if (Frame::cosTheta(bRec.wi) <= 0 || 
			Frame::cosTheta(bRec.wo) <= 0 || !sampleDiffuse)
			return Spectrum(0.0f);

		return m_diffuseReflectance->getValue(bRec.its) * INV_PI * 
			(1 - fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR));
	}

	Spectrum fDelta(const BSDFQueryRecord &bRec) const {
		bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
			&& (bRec.component == -1 || bRec.component == 1);
		if (Frame::cosTheta(bRec.wi) <= 0 || 
			Frame::cosTheta(bRec.wo) <= 0 || !sampleSpecular)
			return Spectrum(0.0f);

		return m_specularReflectance->getValue(bRec.its) * 
			fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		bool sampleDiffuse = (bRec.typeMask & EDiffuseReflection)
			&& (bRec.component == -1 || bRec.component == 0);

		if (Frame::cosTheta(bRec.wi) <= 0 || 
			Frame::cosTheta(bRec.wo) <= 0 || !sampleDiffuse)
			return 0.0f;

		bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
			&& (bRec.component == -1 || bRec.component == 1);

		Float pdf = Frame::cosTheta(bRec.wo) * INV_PI;
		if (sampleSpecular)
			pdf *= 1 - fresnel(Frame::cosTheta(bRec.wi),
					m_extIOR, m_intIOR);

		return pdf;
	}

	Float pdfDelta(const BSDFQueryRecord &bRec) const {
		bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
			&& (bRec.component == -1 || bRec.component == 1);

		if (Frame::cosTheta(bRec.wi) <= 0 || 
			Frame::cosTheta(bRec.wo) <= 0 || !sampleSpecular)
			return 0.0f;

		bool sampleDiffuse = (bRec.typeMask & EDiffuseReflection)
			&& (bRec.component == -1 || bRec.component == 0);

		Float pdf = std::abs(Frame::cosTheta(bRec.wo));
		if (sampleDiffuse)
			pdf *= fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);

		return pdf;
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &_sample) const {
		bool sampleDiffuse = (bRec.typeMask & EDiffuseReflection)
			&& (bRec.component == -1 || bRec.component == 0);
		bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
			&& (bRec.component == -1 || bRec.component == 1);

		if ((!sampleDiffuse && !sampleSpecular) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);
			
		Float Fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);
		Point2 sample(_sample);

		if (sampleDiffuse && sampleSpecular) {
			if (sample.x > Fr) {
				sample.x = (sample.x - Fr) / (1 - Fr);
				bRec.wo = squareToHemispherePSA(sample);
				bRec.sampledComponent = 0;
				bRec.sampledType = EDiffuseReflection;
				pdf = Frame::cosTheta(bRec.wo) * INV_PI * (1-Fr);
				return m_diffuseReflectance->getValue(bRec.its) * INV_PI * (1-Fr);
			} else {
				bRec.sampledComponent = 1;
				bRec.sampledType = EDeltaReflection;
				bRec.wo = reflect(bRec.wi);
				pdf = std::abs(Frame::cosTheta(bRec.wo)) * Fr;
				return m_specularReflectance->getValue(bRec.its) * Fr;
			}
		} else if (sampleDiffuse) {
			bRec.wo = squareToHemispherePSA(sample);
			bRec.sampledComponent = 0;
			bRec.sampledType = EDiffuseReflection;
			pdf = Frame::cosTheta(bRec.wo) * INV_PI;
			return m_diffuseReflectance->getValue(bRec.its) * INV_PI * (1-Fr);
		} else {
			bRec.sampledComponent = 1;
			bRec.sampledType = EDeltaReflection;
			bRec.wo = reflect(bRec.wi);
			pdf = std::abs(Frame::cosTheta(bRec.wo));
			return m_specularReflectance->getValue(bRec.its) * Fr;
		}
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &_sample) const {
		bool sampleDiffuse = (bRec.typeMask & EDiffuseReflection)
			&& (bRec.component == -1 || bRec.component == 0);
		bool sampleSpecular = (bRec.typeMask & EDeltaReflection)
			&& (bRec.component == -1 || bRec.component == 1);

		if ((!sampleDiffuse && !sampleSpecular) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);
			
		Float Fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);
		Point2 sample(_sample);

		if (sampleDiffuse && sampleSpecular) {
			if (sample.x > Fr) {
				sample.x = (sample.x - Fr) / (1 - Fr);
				bRec.wo = squareToHemispherePSA(sample);
				bRec.sampledComponent = 0;
				bRec.sampledType = EDiffuseReflection;
				return m_diffuseReflectance->getValue(bRec.its)
				   / Frame::cosTheta(bRec.wo);
			} else {
				bRec.sampledComponent = 1;
				bRec.sampledType = EDeltaReflection;
				bRec.wo = reflect(bRec.wi);
				return m_specularReflectance->getValue(bRec.its)
					/ std::abs(Frame::cosTheta(bRec.wo));
			}
		} else if (sampleDiffuse) {
			bRec.wo = squareToHemispherePSA(sample);
			bRec.sampledComponent = 0;
			bRec.sampledType = EDiffuseReflection;
			return m_diffuseReflectance->getValue(bRec.its) * (1-Fr)
				/ Frame::cosTheta(bRec.wo);
		} else {
			bRec.sampledComponent = 1;
			bRec.sampledType = EDeltaReflection;
			bRec.wo = reflect(bRec.wi);
			return m_specularReflectance->getValue(bRec.its) * Fr
				/ std::abs(Frame::cosTheta(bRec.wo));
		}
	}


	std::string toString() const {
		std::ostringstream oss;
		oss << "SmoothPlastic[" << endl
			<< "  intIOR = " << m_intIOR << "," << endl 
			<< "  extIOR = " << m_extIOR << "," << endl
			<< "  diffuseReflectance = " << indent(m_diffuseReflectance->toString()) << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	Float m_intIOR, m_extIOR;
	ref<Texture> m_specularReflectance;
	ref<Texture> m_diffuseReflectance;
};

/* Fake glass shader -- it is really hopeless to visualize
   this material in the VPL renderer, so let's try to do at least 
   something that suggests the presence of a transparent boundary */
class SmoothPlasticShader : public Shader {
public:
	SmoothPlasticShader(Renderer *renderer) :
		Shader(renderer, EBSDFShader) {
		m_flags = ETransparent;
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return vec3(0.08);" << endl
			<< "}" << endl;
		oss << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return vec3(0.08);" << endl
			<< "}" << endl;
	}
	MTS_DECLARE_CLASS()
};

Shader *SmoothPlastic::createShader(Renderer *renderer) const { 
	return new SmoothPlasticShader(renderer);
}

MTS_IMPLEMENT_CLASS(SmoothPlasticShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(SmoothPlastic, false, BSDF)
MTS_EXPORT_PLUGIN(SmoothPlastic, "Smooth plastic BSDF");
MTS_NAMESPACE_END

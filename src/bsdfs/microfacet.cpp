/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/**
 * Microfacet BRDF model based on
 * "Microfacet Models for Refraction through Rough Surfaces"
 * by Bruce Walter, Stephen R. Marschner, Hongsong Li
 * and Kenneth E. Torrance.
 *
 * Can be used to simulate a diffuse material coated with
 * a rough dielectric layer -- the diffuse reflectance is
 * modulated by the Fresnel transmittance.
 */
class Microfacet : public BSDF {
public:
	Microfacet(const Properties &props) 
		: BSDF(props) {
		m_diffuseReflectance = new ConstantTexture(
			props.getSpectrum("diffuseReflectance", Spectrum(0.0f)));
		m_specularReflectance = new ConstantTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));

		m_kd = props.getFloat("diffuseAmount", 1.0f);
		m_ks = props.getFloat("specularAmount", 1.0f);

		bool verifyEnergyConservation = props.getBoolean("verifyEnergyConservation", true);

		if (verifyEnergyConservation && (m_kd * m_diffuseReflectance->getMaximum().max() 
				+ m_ks * m_specularReflectance->getMaximum().max() > 1.0f)) {
			Log(EWarn, "%s: Energy conservation is potentially violated!", props.getID().c_str());
			Log(EWarn, "Max. diffuse reflectance = %f * %f = %f", m_kd, m_diffuseReflectance->getMaximum().max(), m_kd*m_diffuseReflectance->getMaximum().max());
			Log(EWarn, "Max. specular reflectance = %f * %f = %f", m_ks, m_specularReflectance->getMaximum().max(), m_ks*m_specularReflectance->getMaximum().max());
			Float normalization = 1/(m_kd * m_diffuseReflectance->getMaximum().max() + m_ks * m_specularReflectance->getMaximum().max());
			Log(EWarn, "Reducing the albedo to %.1f%% of the original value to be on the safe side. "
				"Specify verifyEnergyConservation=false to prevent this.", normalization * 100);
			m_kd *= normalization; m_ks *= normalization;
		}

		Float avgDiffReflectance = m_diffuseReflectance->getAverage().average() * m_kd;
		Float avgSpecularReflectance = m_specularReflectance->getAverage().average() * m_ks;

		m_specularSamplingWeight = props.getFloat("specularSamplingWeight", 
			avgSpecularReflectance / (avgDiffReflectance + avgSpecularReflectance));
		m_diffuseSamplingWeight = 1.0f - m_specularSamplingWeight;

		m_alphaB = props.getFloat("alphaB", .1f);
		m_intIOR = props.getFloat("intIOR", 1.5f);
		m_extIOR = props.getFloat("extIOR", 1.0f);

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDiffuseReflection;
		m_type[1] = EGlossyReflection;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = 
			m_diffuseReflectance->usesRayDifferentials() ||
			m_specularReflectance->usesRayDifferentials();
	}

	Microfacet(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager) {
		m_diffuseReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaB = stream->readFloat();
		m_kd = stream->readFloat();
		m_ks = stream->readFloat();
		m_intIOR = stream->readFloat();
		m_extIOR = stream->readFloat();
		m_specularSamplingWeight = stream->readFloat();
		m_diffuseSamplingWeight = stream->readFloat();

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDiffuseReflection;
		m_type[1] = EGlossyReflection;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = 
			m_diffuseReflectance->usesRayDifferentials() ||
			m_specularReflectance->usesRayDifferentials();
	}

	virtual ~Microfacet() {
		delete[] m_type;
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return m_diffuseReflectance->getValue(its) * m_kd;
	}

	/**
	 * Beckmann distribution function for gaussian random surfaces
	 * @param thetaM Tangent of the angle between M and N.
	 */
	Float beckmannD(const Vector &m) const {
		Float ex = Frame::tanTheta(m) / m_alphaB;
		return std::exp(-(ex*ex)) / (M_PI * m_alphaB*m_alphaB * 
			std::pow(Frame::cosTheta(m), (Float) 4.0f));
	}

	/**
	 * Sample microsurface normals according to 
	 * the Beckmann distribution
	 */
	Normal sampleBeckmannD(Point2 sample) const {
		Float thetaM = std::atan(std::sqrt(-m_alphaB*m_alphaB 
			* std::log(1.0f - sample.x)));
		Float phiM = (2.0f * M_PI) * sample.y;
		return Normal(sphericalDirection(thetaM, phiM));
	}

	/**
	 * Smith's shadow-masking function G1 for the Beckmann distribution
	 * @param m The microsurface normal
	 * @param v An arbitrary direction
	 */
	Float smithBeckmannG1(const Vector &v, const Vector &m) const {
		if (dot(v, m)*Frame::cosTheta(v) <= 0)
			return 0.0;

		const Float tanTheta = Frame::tanTheta(v);

		if (tanTheta == 0.0f)
			return 1.0f;

		const Float a = 1.0f / (m_alphaB * tanTheta);
		const Float aSqr = a * a;

		if (a >= 1.6f)
			return 1.0f;

		return (3.535f * a + 2.181f * aSqr)/(1.0f + 2.276f * a + 2.577f * aSqr);
	}

	inline Vector reflect(const Vector &wi, const Normal &n) const {
		return Vector(n*(2.0f*dot(n, wi))) - wi;
	}

	inline Spectrum fSpec(const BSDFQueryRecord &bRec, const Vector &Hr) const {
		/* Microsurface normal distribution */
		Float D = beckmannD(Hr);
		/* Smith's shadow-masking function for the Beckmann distribution */
		Float G = smithBeckmannG1(bRec.wi, Hr) * smithBeckmannG1(bRec.wo, Hr);
		/* Calculate the total amount of specular reflection */
		Float specRef = D * G / 
			(4.0f * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo));

		return m_specularReflectance->getValue(bRec.its) * specRef; 
	}
	
	inline Spectrum f(const BSDFQueryRecord &bRec) const {
		Spectrum result(0.0f);

		if (!(bRec.typeMask & m_combinedType)
			|| bRec.wi.z <= 0 || bRec.wo.z <= 0)
			return result;

		bool hasDiffuse = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasGlossy   = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		/* Fresnel factor */
		Vector Hr = normalize(bRec.wi+bRec.wo);
		Float F = fresnel(dot(bRec.wi, Hr), m_extIOR, m_intIOR);

		if (hasGlossy)
			result += fSpec(bRec, Hr) * (F * m_ks);
		if (hasDiffuse)
			result += m_diffuseReflectance->getValue(bRec.its) * (INV_PI * (1-F) * m_kd);

		return result;
	}

	inline Float pdfSpec(const BSDFQueryRecord &bRec) const {
		Vector Hr = normalize(bRec.wi + bRec.wo);
		/* Jacobian of the half-direction transform. */
		Float dwhr_dwo = 1.0f / (4.0f * absDot(bRec.wo, Hr));
		return beckmannD(Hr) * Frame::cosTheta(Hr) * dwhr_dwo;
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		bool hasDiffuse = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasGlossy   = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		if (bRec.wi.z <= 0 || bRec.wo.z <= 0)
			return 0.0f;
		
		if (hasDiffuse && hasGlossy) {
			return m_specularSamplingWeight * pdfSpec(bRec) +
				   m_diffuseSamplingWeight * pdfDiffuse(bRec);
		} else if (hasDiffuse) {
			return pdfDiffuse(bRec);
		} else if (hasGlossy) {
			return pdfSpec(bRec);
		}

		return 0.0f;
	}
	
	inline Spectrum sampleSpecular(BSDFQueryRecord &bRec) const {
		/* Sample M, the microsurface normal */
		Normal m = sampleBeckmannD(bRec.sample);
		/* Perfect specular reflection along the microsurface normal */
		bRec.wo = reflect(bRec.wi, m);

		bRec.sampledComponent = 1;
		bRec.sampledType = EGlossyReflection;

		if (bRec.wo.z <= 0)
			return Spectrum(0.0f);

		return f(bRec) / pdf(bRec);
	}

	inline Float pdfDiffuse(const BSDFQueryRecord &bRec) const {
		return Frame::cosTheta(bRec.wo) * INV_PI;
	}

	inline Spectrum sampleLambertian(BSDFQueryRecord &bRec) const {
		bRec.wo = squareToHemispherePSA(bRec.sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		return f(bRec) / pdf(bRec);
	}

	Spectrum sample(BSDFQueryRecord &bRec) const {
		if (bRec.wi.z <= 0)
			return Spectrum(0.0f);

		bool sampleDiffuse = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleGlossy   = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		if (sampleDiffuse && sampleGlossy) {
			if (bRec.sample.x <= m_specularSamplingWeight) {
				bRec.sample.x = bRec.sample.x / m_specularSamplingWeight;
				return sampleSpecular(bRec);
			} else {
				bRec.sample.x = (bRec.sample.x - m_specularSamplingWeight)
					/ m_diffuseSamplingWeight;
				return sampleLambertian(bRec);
			}
		} else if (sampleDiffuse) {
			return sampleLambertian(bRec);
		} else if (sampleGlossy) {
			return sampleSpecular(bRec);
		}

		return Spectrum(0.0f);
	}
	
	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(Texture::m_theClass) && name == "diffuseReflectance") {
			m_diffuseReflectance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_diffuseReflectance->usesRayDifferentials();
		} else if (child->getClass()->derivesFrom(Texture::m_theClass) && name == "specularReflectance") {
			m_specularReflectance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_specularReflectance->usesRayDifferentials();
		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_diffuseReflectance.get());
		manager->serialize(stream, m_specularReflectance.get());
		stream->writeFloat(m_alphaB);
		stream->writeFloat(m_kd);
		stream->writeFloat(m_ks);
		stream->writeFloat(m_intIOR);
		stream->writeFloat(m_extIOR);
		stream->writeFloat(m_specularSamplingWeight);
		stream->writeFloat(m_diffuseSamplingWeight);
	}

	Shader *createShader(Renderer *renderer) const; 

	std::string toString() const {
		std::ostringstream oss;
		oss << "Microfacet[" << endl
			<< "  diffuseReflectance = " << indent(m_diffuseReflectance->toString()) << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  intIOR = " << m_intIOR << "," << endl
			<< "  extIOR = " << m_extIOR << "," << endl
			<< "  alphaB = " << m_alphaB << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_diffuseReflectance;
	ref<Texture> m_specularReflectance;
	Float m_alphaB, m_intIOR, m_extIOR;
	Float m_specularSamplingWeight;
	Float m_diffuseSamplingWeight;
	Float m_ks, m_kd;
};

// ================ Hardware shader implementation ================ 

class MicrofacetShader : public Shader {
public:
	MicrofacetShader(Renderer *renderer, 
			const Texture *diffuseReflectance,
			const Texture *specularReflectance,
			Float alphaB, Float etaInt, Float etaExt,
			Float ks, Float kd) : Shader(renderer, EBSDFShader), 
			m_diffuseReflectance(diffuseReflectance),
			m_specularReflectance(specularReflectance),
			m_alphaB(alphaB), m_etaInt(etaInt), m_etaExt(etaExt),
			m_ks(ks), m_kd(kd) {
		m_diffuseReflectanceShader = renderer->registerShaderForResource(m_diffuseReflectance.get());
		m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());
	}

	bool isComplete() const {
		return m_diffuseReflectanceShader.get() != NULL &&
			   m_specularReflectanceShader.get() != NULL;
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_diffuseReflectanceShader.get());
		deps.push_back(m_specularReflectanceShader.get());
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_diffuseReflectance.get());
		renderer->unregisterShaderForResource(m_specularReflectance.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform float " << evalName << "_alphaB;" << endl
			<< "uniform float " << evalName << "_etaInt;" << endl
			<< "uniform float " << evalName << "_etaExt;" << endl
			<< "uniform float " << evalName << "_ks;" << endl
			<< "uniform float " << evalName << "_kd;" << endl
			<< endl
			<< "float " << evalName << "_beckmannD(vec3 m) {" << endl
			<< "   float ex = tanTheta(m) / " << evalName << "_alphaB;" << endl
			<< "   return exp(-(ex*ex))/(3.14159 * " << evalName << "_alphaB*" << evalName << "_alphaB*" << endl
			<< "          pow(cosTheta(m), 4.0));" << endl
			<< "}" << endl
			<< endl
			<< "float " << evalName << "_smithBeckmannG1(vec3 v, vec3 m) {" << endl
			<< "   if (dot(v, m)*cosTheta(v) <= 0)" << endl
			<< "      return 0.0;" << endl
			<< "   float tt = tanTheta(v);" << endl
			<< "   float a = 1.0/(" << evalName << "_alphaB*tt);" << endl
			<< "   if (a > 1.6)" << endl
			<< "     return 1.0;" << endl
			<< "   float aSqr = a*a;" << endl
			<< "   return (3.535 * a + 2.181 * aSqr)/(1.0 + 2.276 * a + 2.577 * aSqr);" << endl
			<< "}" << endl
			<< endl
			<< "float " << evalName << "_fSpec(vec3 wi, vec3 wo, vec3 hr) {" << endl
			<< "   float D = " << evalName << "_beckmannD(hr)" << ";" << endl
			<< "   float G = " << evalName << "_smithBeckmannG1(wi, hr) * " << endl
			<< "             " << evalName << "_smithBeckmannG1(wo, hr);"  << endl
			<< "   return D * G / (4*cosTheta(wi) * cosTheta(wo));" << endl
			<< "}" << endl
			<< endl
			<< "float " << evalName << "_fresnel(vec3 wi, vec3 n, float etaExt, float etaInt) {" << endl
			<< "    float eta = etaExt / etaInt; " << endl
			<< "    float cosTheta1 = dot(wi, n); " << endl
			<< "    float cosTheta2 = sqrt(1.0 - eta*eta*" << endl
			<< "                    (1.0 - cosTheta1 * cosTheta1));" << endl
			<< "    float Rs = (etaExt * cosTheta1 - etaInt * cosTheta2)" << endl
			<< "             / (etaExt * cosTheta1 + etaInt * cosTheta2);" << endl
			<< "    float Rp = (etaInt * cosTheta1 - etaExt * cosTheta2)" << endl
			<< "             / (etaInt * cosTheta1 + etaExt * cosTheta2);" << endl
			<< "    return (Rs * Rs + Rp * Rp) / 2.0;" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (wi.z <= 0 || wo.z <= 0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    vec3 hr = normalize(wi + wo);" << endl
			<< "    float Fr = " << evalName << "_fresnel(wi, hr, " << evalName << "_etaExt, " << evalName << "_etaInt);"<< endl 
			<< "    float Ft = 1-Fr;"<< endl 
			<< "    return " << depNames[0] << "(uv) * (0.31831 * Ft) * " << evalName << "_kd" << endl
			<< "           + " << depNames[1] << "(uv) * (" << evalName << "_fSpec(wi, wo, hr) * Fr) * " << evalName << "_ks" << endl
			<< "}" << endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (wi.z <= 0 || wo.z <= 0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return " << depNames[0] << "(uv) * 0.31831;" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_alphaB", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_etaInt", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_etaExt", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_ks", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_kd", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_alphaB);
		program->setParameter(parameterIDs[1], m_etaInt);
		program->setParameter(parameterIDs[2], m_etaExt);
		program->setParameter(parameterIDs[3], m_ks);
		program->setParameter(parameterIDs[4], m_kd);
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_diffuseReflectance;
	ref<const Texture> m_specularReflectance;
	ref<Shader> m_diffuseReflectanceShader;
	ref<Shader> m_specularReflectanceShader;
	Float m_alphaB, m_etaInt, m_etaExt, m_ks, m_kd;
};

Shader *Microfacet::createShader(Renderer *renderer) const { 
	return new MicrofacetShader(renderer, m_diffuseReflectance.get(),
		m_specularReflectance.get(), m_alphaB, m_intIOR, m_extIOR,
		m_ks, m_kd);
}

MTS_IMPLEMENT_CLASS(MicrofacetShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Microfacet, false, BSDF)
MTS_EXPORT_PLUGIN(Microfacet, "Microfacet BRDF");
MTS_NAMESPACE_END

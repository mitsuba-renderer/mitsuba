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
#include <mitsuba/hw/basicshader.h>
#include "../medium/materials.h"
#include "ior.h"

MTS_NAMESPACE_BEGIN

/**
 * Integrated Jensen-style dipole BRDF -- not to be
 * used just by itself. Please refer to the sssbrdf plugin
 */
class DipoleBRDF : public BSDF {
public:
	DipoleBRDF(const Properties &props) 
		: BSDF(props) {

		Spectrum sigmaS, sigmaA;
		Float eta = 1.33f;
		lookupMaterial(props, sigmaS, sigmaA, &eta, false);

		/* Specifies the internal index of refraction at the interface */
		m_intIOR = lookupIOR(props, "intIOR", eta);

		/* Specifies the external index of refraction at the interface */
		m_extIOR = lookupIOR(props, "extIOR", "air");

		/* Mean cosine angle of the phase function */
		m_g = props.getFloat("g", 0.0f); 

		/* Scattering coefficient of the layer */
		m_sigmaS = new ConstantSpectrumTexture(
			props.getSpectrum("sigmaS", sigmaS));

		/* Absorption coefficient of the layer */
		m_sigmaA = new ConstantSpectrumTexture(
			props.getSpectrum("sigmaA", sigmaA));

		if (props.hasProperty("sigmaT"))
			m_sigmaT = new ConstantSpectrumTexture(
				props.getSpectrum("sigmaT"));
		if (props.hasProperty("albedo"))
			m_albedo = new ConstantSpectrumTexture(
				props.getSpectrum("albedo"));
	}

	DipoleBRDF(Stream *stream, InstanceManager *manager) 
		: BSDF(stream, manager) {

		m_intIOR = stream->readFloat();
		m_extIOR = stream->readFloat();
		m_g	  = stream->readFloat();

		m_sigmaS = static_cast<Texture *>(manager->getInstance(stream));
		m_sigmaA = static_cast<Texture *>(manager->getInstance(stream));

		configure();
	}

	Float evalFdr(Float eta) const {
		if (eta > 1) {
			return  -1.440f / (eta * eta) + 0.710f / eta + 0.668f + 0.0636f * eta;
		} else if (eta < 1) {
			return  -0.4399f + 0.7099f / eta - 0.3319f / (eta * eta)
					+ 0.0636f / (eta * eta * eta);
		} else {
			/* eta == 1 */
			return 0.0f;
		}
	}

	void configure() {
		if (m_sigmaT != NULL || m_albedo != NULL) {
			/* Support for the alternative scattering/absorption
			 * coefficient parameter passing convention */
			if (m_sigmaT == NULL || m_albedo == NULL)
				SLog(EError, "Please provide *both* sigmaT & albedo!");

			m_sigmaS = new SpectrumProductTexture(m_sigmaT, m_albedo);
			m_sigmaA = new SpectrumSubtractionTexture(m_sigmaT, m_sigmaS);
			m_sigmaT = NULL;
			m_albedo = NULL;
		}

		m_components.clear();
		m_components.push_back(EDiffuseReflection | EFrontSide
			| (m_sigmaS->isConstant() && m_sigmaA->isConstant() ? 0 : ESpatiallyVarying));
		m_usesRayDifferentials = m_sigmaS->usesRayDifferentials()
			|| m_sigmaA->usesRayDifferentials();

		/* Numerically approximate the diffuse Fresnel reflectance */
		const Float Fdr = fresnelDiffuseReflectance(m_extIOR / m_intIOR, false);

		/* Compute the extrapolation distance */
		m_A = (1 + Fdr) / (1 - Fdr);

		BSDF::configure();
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		Spectrum sigmaA = m_sigmaA->getValue(its),
				 sigmaS = m_sigmaS->getValue(its),
				 reducedAlbedo;

		/* ==================================================================== */
		/*            Diffuse Reflectance due to Multiple Scattering            */
		/* ==================================================================== */
	
		/* Reduced scattering albedo */
		const Spectrum reducedSigmaS = sigmaS * (1.0f - m_g),
		               reducedSigmaT = reducedSigmaS + sigmaA;

		/* Avoid divisions by 0 */
		for (int i = 0; i < SPECTRUM_SAMPLES; i++)
			reducedAlbedo[i] = reducedSigmaT[i] > 0.0f ?
				(reducedSigmaS[i]/reducedSigmaT[i]) : 0.0f;
		
		/* Diffuse Reflectance */
		const Spectrum rootExp = ((Spectrum(1.0f) - reducedAlbedo) * 3.0f).sqrt();
		return (reducedAlbedo * 0.5f) * (Spectrum(1.0f) +
				(-rootExp*(4.0f/3.0f*m_A)).exp()) * (-rootExp).exp();
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0 
			|| Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		return getDiffuseReflectance(bRec.its) * (INV_PI * Frame::cosTheta(bRec.wo));
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0 
			|| Frame::cosTheta(bRec.wo) <= 0)
			return 0.0f;

		return Frame::cosTheta(bRec.wo) * INV_PI;
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0) 
			return Spectrum(0.0f);
		
		bRec.wo = squareToHemispherePSA(sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;

		return getDiffuseReflectance(bRec.its);
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0) 
			return Spectrum(0.0f);
		
		bRec.wo = squareToHemispherePSA(sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;

		pdf = DipoleBRDF::pdf(bRec, ESolidAngle);

		return getDiffuseReflectance(bRec.its);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "sigmaS")
				m_sigmaS = static_cast<Texture *>(child);
			else if (name == "sigmaA")
				m_sigmaA = static_cast<Texture *>(child);
			else if (name == "sigmaT")
				m_sigmaT = static_cast<Texture *>(child);
			else if (name == "albedo")
				m_albedo = static_cast<Texture *>(child);
			else
				BSDF::addChild(name, child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeFloat(m_intIOR);
		stream->writeFloat(m_extIOR);
		stream->writeFloat(m_g);

		manager->serialize(stream, m_sigmaS.get());
		manager->serialize(stream, m_sigmaA.get());
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "DipoleBRDF[" << endl
			<< "  intIOR = " << m_intIOR  << "," << endl 
			<< "  extIOR = " << m_extIOR  << "," << endl
   			<< "  g = " << m_g << "," << endl
   			<< "  sigmaS = " << indent(m_sigmaS.toString()) << "," << endl
   			<< "  sigmaA = " << indent(m_sigmaA.toString()) << "," << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	Float m_intIOR, m_extIOR, m_g, m_A;
	ref<Texture> m_sigmaS, m_sigmaA;
	/* Temporary fields */
	ref<Texture> m_sigmaT;
	ref<Texture> m_albedo;
};

// ================ Hardware shader implementation ================ 

/**
 * This is a relatively approximate GLSL shader for the Dipole BRDF model.
 * It assumes that the layer is infinitely thick (i.e. there is no
 * transmission) and that the phase function is isotropic
 */
class DipoleBRDFShader : public Shader {
public:
	DipoleBRDFShader(Renderer *renderer, const Texture *sigmaS, const Texture *sigmaA, 
		Float A, Float g) : Shader(renderer, EBSDFShader), m_sigmaS(sigmaS), 
		  m_sigmaA(sigmaA), m_A(A), m_g(g) {
		m_sigmaSShader = renderer->registerShaderForResource(m_sigmaS.get());
		m_sigmaAShader = renderer->registerShaderForResource(m_sigmaA.get());
	}

	bool isComplete() const {
		return m_sigmaSShader.get() != NULL
			&& m_sigmaAShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_sigmaS.get());
		renderer->unregisterShaderForResource(m_sigmaA.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_sigmaSShader.get());
		deps.push_back(m_sigmaAShader.get());
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_A", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_g", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_A);
		program->setParameter(parameterIDs[1], m_g);
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform float " << evalName << "_A;" << endl
			<< "uniform float " << evalName << "_g;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    float cosThetaI = cosTheta(wi);" << endl
			<< "    float cosThetaO = cosTheta(wo);" << endl
			<< "    if (cosThetaI < 0.0 || cosThetaO < 0.0)" << endl
			<< "        return vec3(0.0);" << endl
			<< "    vec3 sigmaS = " << depNames[0] << "(uv);" << endl
			<< "    vec3 sigmaA = " << depNames[1] << "(uv);" << endl
			<< "    vec3 reducedSigmaS = (1.0-" << evalName << "_g)*sigmaS;" << endl
			<< "    vec3 reducedSigmaT = reducedSigmaS + sigmaA, reducedAlbedo;" << endl
			<< "    for (int i=0; i<3; ++i)" << endl
			<< "        reducedAlbedo[i] = reducedSigmaT[i] > 0.0 ? reducedSigmaS[i]/reducedSigmaT[i] : 0.0;" << endl
			<< "    vec3 rootExp = sqrt((1.0 - reducedAlbedo) * 3.0);" << endl
			<< "    return (reducedAlbedo * 0.5) * (1+exp(-rootExp*(4.0/3.0 * " << endl
			<< "      " << evalName << "_A" << "))) * exp(-rootExp) * inv_pi * cosThetaO;" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << evalName << "(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_sigmaS;
	ref<const Texture> m_sigmaA;
	ref<Shader> m_sigmaSShader;
	ref<Shader> m_sigmaAShader;

	Float m_A, m_g;
};

Shader *DipoleBRDF::createShader(Renderer *renderer) const { 
	return new DipoleBRDFShader(renderer, m_sigmaS.get(), 
			m_sigmaA.get(), m_A, m_g);
}

MTS_IMPLEMENT_CLASS(DipoleBRDFShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(DipoleBRDF, false, BSDF)
MTS_EXPORT_PLUGIN(DipoleBRDF, "Dipole BRDF")
MTS_NAMESPACE_END

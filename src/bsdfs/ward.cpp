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
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

/**
 * Anisotropic Ward BRDF model based on the following four papers
 *
 *   "Measuring and Modeling Anisotropic Reflection" by 
 *   Gregory J. Ward, SIGGRAPH 1992
 *
 *   "Notes on the Ward BRDF" by Bruce Walter, Technical Report 
 *   PCG-05-06, Cornell University
 *
 *   "An Improved Normalization for the Ward Reflectance Model"
 *   by Arne Duer, Journal of Graphics Tools 11, 1 (2006), 51–59
 *
 *   "A New Ward BRDF Model with Bounded Albedo" by
 *   by David Geisler-Moroder and Arne Dur, 
 *   Computer Graphics Forum, Volume 29, Issue 4
 */
class Ward : public BSDF {
public:
	/// Supported model types
	enum EModelType {
		/// The original Ward model
		EWard = 0,
		/// Ward model with correction by Arne Duer
		EWardDuer = 1,
		/// Energy-balanced Ward model
		EBalanced = 2
	};

	Ward(const Properties &props) 
		: BSDF(props) {
		m_diffuseReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("diffuseReflectance", Spectrum(0.5f)));
		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(0.2f)));
		
		std::string type = 
			boost::to_lower_copy(props.getString("type", "balanced"));
		if (type == "ward")
			m_modelType = EWard;
		else if (type == "ward-duer")
			m_modelType = EWardDuer;
		else if (type == "balanced")
			m_modelType = EBalanced;
		else
			Log(EError, "Specified an invalid model type \"%s\", must be "
				"\"ward\", \"ward-duer\", or \"balanced\"!", type.c_str());

		Float alpha = props.getFloat("alpha", 0.1f),
			  alphaU = props.getFloat("alphaU", alpha),
			  alphaV = props.getFloat("alphaV", alpha);

		m_alphaU = new ConstantFloatTexture(alphaU);
		if (alphaU == alphaV)
			m_alphaV = m_alphaU;
		else
			m_alphaV = new ConstantFloatTexture(alphaV);
	}

	Ward(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager) {
		m_modelType = (EModelType) stream->readUInt();
		m_diffuseReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaV = static_cast<Texture *>(manager->getInstance(stream));

		configure();
	}

	virtual ~Ward() { }

	void configure() {
		unsigned int extraFlags = 0;
		if (m_alphaU != m_alphaV)
			extraFlags |= EAnisotropic;

		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide | extraFlags);
		m_components.push_back(EDiffuseReflection | EFrontSide | extraFlags);

		/* Verify the input parameters and fix them if necessary */
		std::pair<Texture *, Texture *> result = ensureEnergyConservation(
			m_specularReflectance, m_diffuseReflectance,
			"specularReflectance", "diffuseReflectance", 1.0f);
		m_specularReflectance = result.first;
		m_diffuseReflectance = result.second;

		/* Compute weights that steer samples towards
		   the specular or diffuse components */
		Float dAvg = m_diffuseReflectance->getAverage().getLuminance(),
			  sAvg = m_specularReflectance->getAverage().getLuminance();
		m_specularSamplingWeight = sAvg / (dAvg + sAvg);

		m_usesRayDifferentials = 
			m_diffuseReflectance->usesRayDifferentials() ||
			m_specularReflectance->usesRayDifferentials() ||
			m_alphaU->usesRayDifferentials() ||
			m_alphaV->usesRayDifferentials();

		BSDF::configure();
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		if (Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 || measure != ESolidAngle)
			return Spectrum(0.0f);

		bool hasSpecular = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasDiffuse  = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		Spectrum result(0.0f);
		if (hasSpecular) {
			Vector H = bRec.wi+bRec.wo;
			Float alphaU = m_alphaU->getValue(bRec.its).average();
			Float alphaV = m_alphaV->getValue(bRec.its).average();
			
			Float factor1 = 0.0f;
			switch (m_modelType) {
				case EWard:
					factor1 = 1.0f / (4.0f * M_PI * alphaU * alphaV * 
						std::sqrt(Frame::cosTheta(bRec.wi)*Frame::cosTheta(bRec.wo)));
					break;
				case EWardDuer:
					factor1 = 1.0f / (4.0f * M_PI * alphaU * alphaV * 
						Frame::cosTheta(bRec.wi)*Frame::cosTheta(bRec.wo));
					break;
				case EBalanced:
					factor1 = dot(H,H) / (M_PI * alphaU * alphaV 
						* std::pow(Frame::cosTheta(H),4));
					break;
				default:
					Log(EError, "Unknown model type!");
			}

			Float factor2 = H.x / alphaU, factor3 = H.y / alphaV;
			Float exponent = -(factor2*factor2+factor3*factor3)/(H.z*H.z);
			Float specRef = factor1 * std::exp(exponent);
			/* Important to prevent numeric issues when evaluating the
			   sampling density of the Ward model in places where it takes
			   on miniscule values (Veach-MLT does this for instance) */
			if (specRef > 1e-10f)
				result += m_specularReflectance->getValue(bRec.its) * specRef;
		}

		if (hasDiffuse) 
			result += m_diffuseReflectance->getValue(bRec.its) * INV_PI;

		return result * Frame::cosTheta(bRec.wo);
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		if (Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 || measure != ESolidAngle)
			return 0.0f;

		bool hasSpecular = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasDiffuse  = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		Float diffuseProb = 0.0f, specProb = 0.0f;

		if (hasSpecular) {
			Float alphaU = m_alphaU->getValue(bRec.its).average();
			Float alphaV = m_alphaV->getValue(bRec.its).average();
			Vector H = normalize(bRec.wi+bRec.wo);
			Float factor1 = 1.0f / (4.0f * M_PI * alphaU * alphaV * 
				dot(H, bRec.wi) * std::pow(Frame::cosTheta(H), 3));
			Float factor2 = H.x / alphaU, factor3 = H.y / alphaV;

			Float exponent = -(factor2*factor2+factor3*factor3)/(H.z*H.z);
			specProb = factor1 * std::exp(exponent);
		}

		if (hasDiffuse) 
			diffuseProb = Frame::cosTheta(bRec.wo) * INV_PI;

		if (hasDiffuse && hasSpecular)
			return m_specularSamplingWeight * specProb + 
				   (1-m_specularSamplingWeight) * diffuseProb;
		else if (hasDiffuse)
			return diffuseProb;
		else if (hasSpecular)
			return specProb;
		else
			return 0.0f;
	}

	inline Spectrum sample(BSDFQueryRecord &bRec, Float &_pdf, const Point2 &_sample) const {
		Point2 sample(_sample);

		bool hasSpecular = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasDiffuse  = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		if (!hasSpecular && !hasDiffuse)
			return Spectrum(0.0f);

		bool choseSpecular = hasSpecular;

		if (hasDiffuse && hasSpecular) {
			if (sample.x <= m_specularSamplingWeight) {
				sample.x /= m_specularSamplingWeight;
			} else {
				sample.x = (sample.x - m_specularSamplingWeight)
					/ (1-m_specularSamplingWeight);
				choseSpecular = false;
			}
		}

		if (choseSpecular) {
			Float alphaU = m_alphaU->getValue(bRec.its).average();
			Float alphaV = m_alphaV->getValue(bRec.its).average();

			Float phiH = std::atan(alphaV/alphaU 
				* std::tan(2.0f * M_PI * sample.y));
			if (sample.y > 0.5f)
				phiH += M_PI;
			Float cosPhiH = std::cos(phiH);
			Float sinPhiH = std::sqrt(std::max((Float) 0.0f, 
				1.0f-cosPhiH*cosPhiH));

			Float thetaH = std::atan(std::sqrt(std::max((Float) 0.0f, 
				-std::log(sample.x) / (
					(cosPhiH*cosPhiH) / (alphaU*alphaU) +
					(sinPhiH*sinPhiH) / (alphaV*alphaV)
			))));
			Vector H = sphericalDirection(thetaH, phiH);
			bRec.wo = H * (2.0f * dot(bRec.wi, H)) - bRec.wi;

			bRec.sampledComponent = 1;
			bRec.sampledType = EGlossyReflection;

			if (Frame::cosTheta(bRec.wo) <= 0.0f)
				return Spectrum(0.0f);
		} else {
			bRec.wo = squareToHemispherePSA(sample);
			bRec.sampledComponent = 0;
			bRec.sampledType = EDiffuseReflection;
		}

		_pdf = pdf(bRec, ESolidAngle);

		if (_pdf == 0) 
			return Spectrum(0.0f);
		else
			return eval(bRec, ESolidAngle);
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		Float pdf=0;
		Spectrum result = Ward::sample(bRec, pdf, sample);

		if (result.isZero())
			return Spectrum(0.0f);
		else
			return result / pdf;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "alphaU")
				m_alphaU = static_cast<Texture *>(child);
			else if (name == "alphaV")
				m_alphaV = static_cast<Texture *>(child);
			else if (name == "diffuseReflectance")
				m_diffuseReflectance = static_cast<Texture *>(child);
			else if (name == "specularReflectance")
				m_specularReflectance = static_cast<Texture *>(child);
			else
				BSDF::addChild(name, child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeUInt(m_modelType);
		manager->serialize(stream, m_diffuseReflectance.get());
		manager->serialize(stream, m_specularReflectance.get());
		manager->serialize(stream, m_alphaU.get());
		manager->serialize(stream, m_alphaV.get());
	}

	//Shader *createShader(Renderer *renderer) const; 

	std::string toString() const {
		std::ostringstream oss;
		oss << "Ward[" << endl
			<< "  type = ";

		switch (m_modelType) {
			case EWard: oss << "ward," << endl; break;
			case EWardDuer: oss << "wardDuer," << endl; break;
			case EBalanced: oss << "balanced," << endl; break;
			default: Log(EError, "Unknown model type!");
		}

		oss << "  diffuseReflectance = " << indent(m_diffuseReflectance->toString()) << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  specularSamplingWeight = " << m_specularSamplingWeight << "," << endl
			<< "  alphaU = " << indent(m_alphaU->toString()) << "," << endl
			<< "  alphaV = " << indent(m_alphaV->toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	EModelType m_modelType;
	ref<Texture> m_diffuseReflectance;
	ref<Texture> m_specularReflectance;
	ref<Texture> m_alphaU;
	ref<Texture> m_alphaV;
	Float m_specularSamplingWeight;
};
#if 0
// ================ Hardware shader implementation ================ 

class WardShader : public Shader {
public:
	WardShader(Renderer *renderer, Ward::EModelType type, 
			const Texture *diffuseColor,
			const Texture *specularColor,
			Float alphaU, Float alphaV) : Shader(renderer, EBSDFShader), 
			m_modelType(type), m_diffuseReflectance(diffuseColor),
			m_specularReflectance(specularColor),
			m_alphaU(alphaU), m_alphaV(alphaV) {
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
		oss << "uniform int " << evalName << "_type;" << endl
			<< "uniform float " << evalName << "_alphaU;" << endl
			<< "uniform float " << evalName << "_alphaV;" << endl
			<< "uniform float " << evalName << "_ks;" << endl
			<< "uniform float " << evalName << "_kd;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (wi.z <= 0.0 || wo.z <= 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    vec3 H = normalize(wi + wo);" << endl
			<< "    float factor1;" << endl
			<< "    if (" << evalName << "_type == 1)" << endl
			<< "        factor1 = 1/(12.566 * " << evalName << "_alphaU * "  << endl
			<< "                    " << evalName << "_alphaV * sqrt(wi.z * wo.z));" << endl
			<< "    else if (" << evalName << "_type == 2)" << endl
			<< "        factor1 = 1/(12.566 * " << evalName << "_alphaU * " << endl
			<< "                    " << evalName << "_alphaV * wi.z * wo.z);" << endl
			<< "    else" << endl
			<< "        factor1 = dot(H, H)/(3.1415 * " << evalName << "_alphaU * "  << endl
			<< "                    " << evalName << "_alphaV * (H.z * H.z) * (H.z * H.z));" << endl
			<< "    float factor2 = H.x / " << evalName << "_alphaU;" << endl
			<< "    float factor3 = H.y / " << evalName << "_alphaV;" << endl
			<< "    float exponent = -(factor2*factor2 + factor3*factor3)/(H.z*H.z);" << endl
			<< "    float specRef = factor1 * exp(exponent) * " << evalName << "_ks;" << endl 
			<< "    return " << depNames[0] << "(uv) * (0.31831 * " << evalName << "_kd)" << endl
			<< "           + " << depNames[1] << "(uv) * specRef;" << endl
			<< "}" << endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (wi.z < 0.0 || wo.z < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return " << depNames[0] << "(uv) * (0.31831 * " << evalName << "_kd);" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_type"));
		parameterIDs.push_back(program->getParameterID(evalName + "_alphaU"));
		parameterIDs.push_back(program->getParameterID(evalName + "_alphaV"));
		parameterIDs.push_back(program->getParameterID(evalName + "_ks"));
		parameterIDs.push_back(program->getParameterID(evalName + "_kd"));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], (int) m_modelType);
		program->setParameter(parameterIDs[1], m_alphaU);
		program->setParameter(parameterIDs[2], m_alphaV);
	}

	MTS_DECLARE_CLASS()
private:
	Ward::EModelType m_modelType;
	ref<const Texture> m_diffuseReflectance;
	ref<const Texture> m_specularReflectance;
	ref<Shader> m_diffuseReflectanceShader;
	ref<Shader> m_specularReflectanceShader;
	Float m_alphaU, m_alphaV;
};

Shader *Ward::createShader(Renderer *renderer) const { 
	return new WardShader(renderer, m_modelType, m_diffuseReflectance.get(),
		m_specularReflectance.get(), m_alphaU, m_alphaV);
}

MTS_IMPLEMENT_CLASS(WardShader, false, Shader)
#endif
MTS_IMPLEMENT_CLASS_S(Ward, false, BSDF);
MTS_EXPORT_PLUGIN(Ward, "Anisotropic Ward BRDF");
MTS_NAMESPACE_END

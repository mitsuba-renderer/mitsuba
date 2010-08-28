#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/**
 * Modified Phong model based on the technical report
 * "Using the Modified Phong Reflectance Model for Physically Based Rendering"
 * by Eric P. Lafortune and Yves D. Willems
 */
class Phong : public BSDF {
public:
	Phong(const Properties &props) 
		: BSDF(props) {
		m_diffuseColor = new ConstantTexture(
			props.getSpectrum("diffuseColor", Spectrum(1.0f)));
		m_specularColor = new ConstantTexture(
			props.getSpectrum("specularColor", Spectrum(1.0f)));

		m_kd = props.getFloat("diffuseReflectance", 0.5f);
		m_ks = props.getFloat("specularReflectance", 0.2f);
		m_exponent = props.getFloat("exponent", 10.0f);

		if (m_kd * m_diffuseColor->getMaximum().max() + m_ks * m_specularColor->getMaximum().max() > 1.0f) {
			Log(EWarn, "%s: Energy conservation is violated!", props.getID().c_str());
			Log(EWarn, "Max. diffuse reflectance = %f * %f", m_kd, m_diffuseColor->getMaximum().max());
			Log(EWarn, "Max. specular reflectance = %f * %f", m_ks, m_specularColor->getMaximum().max());
			Float normalization = 1/(m_kd * m_diffuseColor->getMaximum().max() + m_ks * m_specularColor->getMaximum().max());
			Log(EWarn, "Reducing the albedo to %.1f%% of the original value", normalization * 100);
			m_kd *= normalization; m_ks *= normalization;
		}

		Float avgDiffReflectance = m_diffuseColor->getAverage().average() * m_kd;
		Float avgSpecularReflectance = m_specularColor->getAverage().average() * m_ks;
		
		m_specularSamplingWeight = props.getFloat("specularSamplingWeight", 
			avgSpecularReflectance / (avgDiffReflectance + avgSpecularReflectance));
		m_diffuseSamplingWeight = 1.0f - m_specularSamplingWeight;

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDiffuseReflection;
		m_type[1] = EGlossyReflection;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = false;
	}

	Phong(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager) {
		m_diffuseColor = static_cast<Texture *>(manager->getInstance(stream));
		m_specularColor = static_cast<Texture *>(manager->getInstance(stream));
		m_exponent = stream->readFloat();
		m_kd = stream->readFloat();
		m_ks = stream->readFloat();
		m_specularSamplingWeight = stream->readFloat();
		m_diffuseSamplingWeight = stream->readFloat();

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDiffuseReflection;
		m_type[1] = EGlossyReflection;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = 
			m_diffuseColor->usesRayDifferentials() ||
			m_specularColor->usesRayDifferentials();
	}

	virtual ~Phong() {
		delete[] m_type;
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return m_diffuseColor->getValue(its) * m_kd;
	}

	Spectrum f(const BSDFQueryRecord &bRec) const {
		Spectrum result(0.0f);

		if (bRec.wi.z <= 0 || bRec.wo.z <= 0)
			return result;

		bool hasDiffuse = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasGlossy   = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		Vector R = Vector(-bRec.wi.x, -bRec.wi.y, bRec.wi.z);
		Float alpha = dot(R, bRec.wo);

		if (hasGlossy) {
			Float specRef;
			if (alpha <= 0.0f)
				specRef = 0.0f;
			else
				specRef = (m_exponent + 2) * INV_TWOPI
					* std::pow(alpha, m_exponent) * m_ks;
			result += m_specularColor->getValue(bRec.its) * specRef;
		}

		if (hasDiffuse) 
			result += m_diffuseColor->getValue(bRec.its) * (INV_PI * m_kd);
		return result;
	}

	inline Float pdfSpec(const BSDFQueryRecord &bRec) const {
		Vector R = Vector(-bRec.wi.x, -bRec.wi.y, bRec.wi.z);
		Float alpha = dot(R, bRec.wo);
		Float specPdf = std::pow(alpha, m_exponent) * 
			(m_exponent + 1.0f) / (2.0f * M_PI);
		if (alpha <= 0) 
			specPdf = 0;
		return specPdf;
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		bool hasDiffuse = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasGlossy   = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		if (bRec.wo.z <= 0 || bRec.wi.z <= 0) 
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
		Vector R = Vector(-bRec.wi.x, -bRec.wi.y, bRec.wi.z);

		/* Sample a cosine lobe centered around the normal */
		Float sinAlpha = std::sqrt(1-std::pow(bRec.sample.y, 2/(m_exponent + 1)));
		Float cosAlpha = std::pow(bRec.sample.y, 1/(m_exponent + 1));
		Float phi = (2.0f * M_PI) * bRec.sample.x;
		Vector localDir = Vector(
			sinAlpha * std::cos(phi),
			sinAlpha * std::sin(phi),
			cosAlpha
		);

		/* Rotate into the correct coordinate system */
		bRec.wo = Frame(R).toWorld(localDir);
		bRec.sampledComponent = 1;
		bRec.sampledType = EGlossyReflection;

		if (bRec.wo.z <= 0) 
			return Spectrum(0.0f);

		return f(bRec) / pdf(bRec);
	}

	inline Float pdfDiffuse(const BSDFQueryRecord &bRec) const {
		return Frame::cosTheta(bRec.wo) * INV_PI;
	}

	inline Spectrum sampleDiffuse(BSDFQueryRecord &bRec) const {
		bRec.wo = squareToHemispherePSA(bRec.sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		return f(bRec) / pdf(bRec);
	}

	Spectrum sample(BSDFQueryRecord &bRec) const {
		if (bRec.wi.z <= 0)
			return Spectrum(0.0f);

		bool enableDiffuse = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool enableGlossy   = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		if (enableDiffuse && enableGlossy) {
			Point2 newSample = bRec.sample;

			if (bRec.sample.x <= m_specularSamplingWeight) {
				bRec.sample.x /= m_specularSamplingWeight;
				return sampleSpecular(bRec);
			} else {
				bRec.sample.x = (bRec.sample.x - m_specularSamplingWeight)
					/ m_diffuseSamplingWeight;
				return sampleDiffuse(bRec);
			}
		} else if (enableDiffuse) {
			return sampleDiffuse(bRec);
		} else if (enableGlossy) {
			return sampleSpecular(bRec);
		}

		return Spectrum(0.0f);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(Texture::m_theClass) && name == "diffuseColor") {
			m_diffuseColor = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_diffuseColor->usesRayDifferentials();
		} else if (child->getClass()->derivesFrom(Texture::m_theClass) && name == "specularColor") {
			m_specularColor = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_specularColor->usesRayDifferentials();
		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_diffuseColor.get());
		manager->serialize(stream, m_specularColor.get());
		stream->writeFloat(m_exponent);
		stream->writeFloat(m_kd);
		stream->writeFloat(m_ks);
		stream->writeFloat(m_specularSamplingWeight);
		stream->writeFloat(m_diffuseSamplingWeight);
	}

	Shader *createShader(Renderer *renderer) const; 

	std::string toString() const {
		std::ostringstream oss;
		oss << "Phong[" << endl
			<< "  diffuseColor = " << indent(m_diffuseColor->toString()) << "," << endl
			<< "  specularColor = " << indent(m_specularColor->toString()) << "," << endl
			<< "  exponent = " << m_exponent << endl
			<< "]";
		return oss.str();
	}


	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_diffuseColor;
	ref<Texture> m_specularColor;
	Float m_exponent;
	Float m_kd, m_ks;
	Float m_specularSamplingWeight;
	Float m_diffuseSamplingWeight;
};

// ================ Hardware shader implementation ================ 

class PhongShader : public Shader {
public:
	PhongShader(Renderer *renderer, 
			const Texture *diffuseColor,
			const Texture *specularColor,
			Float ks, Float kd,
			Float exponent) : Shader(renderer, EBSDFShader), 
			m_diffuseColor(diffuseColor),
			m_specularColor(specularColor),
			m_ks(ks), m_kd(kd),
			m_exponent(exponent) {
		m_diffuseColorShader = renderer->registerShaderForResource(m_diffuseColor.get());
		m_specularColorShader = renderer->registerShaderForResource(m_specularColor.get());
	}

	bool isComplete() const {
		return m_diffuseColorShader.get() != NULL &&
			   m_specularColorShader.get() != NULL;
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_diffuseColorShader.get());
		deps.push_back(m_specularColorShader.get());
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_diffuseColor.get());
		renderer->unregisterShaderForResource(m_specularColor.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform float " << evalName << "_exponent;" << endl
			<< "uniform float " << evalName << "_ks;" << endl
			<< "uniform float " << evalName << "_kd;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (wi.z < 0.0 || wo.z < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    vec3 R = vec3(-wi.x, -wi.y, wi.z);" << endl
			<< "    float alpha = dot(R, wo);" << endl 
			<< "    if (alpha < 0.0)" << endl 
			<< "       return vec3(0.0);" << endl 
			<< "    float specRef = pow(alpha, " << evalName << "_exponent) * " << endl
			<< "      (" << evalName << "_exponent + 2) * 0.15915 * " << evalName << "_ks;" << endl
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
		parameterIDs.push_back(program->getParameterID(evalName + "_exponent"));
		parameterIDs.push_back(program->getParameterID(evalName + "_ks"));
		parameterIDs.push_back(program->getParameterID(evalName + "_kd"));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_exponent);
		program->setParameter(parameterIDs[1], m_ks);
		program->setParameter(parameterIDs[2], m_kd);
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_diffuseColor;
	ref<const Texture> m_specularColor;
	ref<Shader> m_diffuseColorShader;
	ref<Shader> m_specularColorShader;
	Float m_ks, m_kd;
	Float m_exponent;
};

Shader *Phong::createShader(Renderer *renderer) const { 
	return new PhongShader(renderer, m_diffuseColor.get(),
		m_specularColor.get(), m_ks, m_kd, m_exponent);
}

MTS_IMPLEMENT_CLASS(PhongShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Phong, false, BSDF)
MTS_EXPORT_PLUGIN(Phong, "Modified Phong BRDF");
MTS_NAMESPACE_END

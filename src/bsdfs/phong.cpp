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
		m_diffuseReflectance = new ConstantTexture(
			props.getSpectrum("diffuseReflectance", Spectrum(.5f)));
		m_specularReflectance = new ConstantTexture(
			props.getSpectrum("specularReflectance", Spectrum(.2f)));
		m_kd = m_diffuseReflectance->getAverage().average();
		m_ks = m_specularReflectance->getAverage().average();
		m_specularSamplingWeight = props.getFloat("specularSamplingWeight", 
			m_ks / (m_kd+m_ks));
		m_diffuseSamplingWeight = 1.0f - m_specularSamplingWeight;
		m_exponent = props.getFloat("exponent", 10.0f);

		if (m_kd + m_ks > 1.0f)
			Log(EWarn, "Energy conservation is violated!");

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDiffuseReflection;
		m_type[1] = EGlossyReflection;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = false;
	}

	Phong(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager) {
		m_diffuseReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
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
			m_diffuseReflectance->usesRayDifferentials() ||
			m_specularReflectance->usesRayDifferentials();
	}

	virtual ~Phong() {
		delete[] m_type;
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return m_diffuseReflectance->getValue(its);
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
					* std::pow(alpha, m_exponent);
			result += m_specularReflectance->getValue(bRec.its) * specRef;
		}

		if (hasDiffuse) 
			result += m_diffuseReflectance->getValue(bRec.its) * INV_PI;
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

		if (m_diffuseSamplingWeight == 0) {
			return m_specularReflectance->getValue(bRec.its) * (m_exponent+2)/(m_exponent+1);
		} else {
			return f(bRec) / pdf(bRec);
		}
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
				bRec.sample.x = bRec.sample.x / m_specularSamplingWeight;
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
			<< "  diffuseReflectance = " << indent(m_diffuseReflectance->toString()) << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  exponent = " << m_exponent << endl
			<< "]";
		return oss.str();
	}


	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_diffuseReflectance;
	ref<Texture> m_specularReflectance;
	Float m_exponent;
	Float m_kd, m_ks;
	Float m_specularSamplingWeight;
	Float m_diffuseSamplingWeight;
};

// ================ Hardware shader implementation ================ 

class PhongShader : public Shader {
public:
	PhongShader(Renderer *renderer, 
			const Texture *diffuseReflectance,
			const Texture *specularReflectance,
			Float exponent) : Shader(renderer, EBSDFShader), 
			m_diffuseReflectance(diffuseReflectance),
			m_specularReflectance(specularReflectance),
			m_exponent(exponent) {
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
		oss << "uniform float " << evalName << "_exponent;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (wi.z < 0.0 || wo.z < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    vec3 R = vec3(-wi.x, -wi.y, wi.z);" << endl
			<< "    float alpha = dot(R, wo);" << endl 
			<< "    if (alpha < 0.0)" << endl 
			<< "       return vec3(0.0);" << endl 
			<< "    float specRef = pow(alpha, " << evalName << "_exponent) * " << endl
			<< "      (" << evalName << "_exponent + 2) * 0.15915;" << endl
			<< "    return " << depNames[0] << "(uv) * 0.31831" << endl
			<< "           + " << depNames[1] << "(uv) * specRef;" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_exponent"));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_exponent);
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_diffuseReflectance;
	ref<const Texture> m_specularReflectance;
	ref<Shader> m_diffuseReflectanceShader;
	ref<Shader> m_specularReflectanceShader;
	Float m_exponent;
};

Shader *Phong::createShader(Renderer *renderer) const { 
	return new PhongShader(renderer, m_diffuseReflectance.get(),
		m_specularReflectance.get(), m_exponent);
}

MTS_IMPLEMENT_CLASS(PhongShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Phong, false, BSDF)
MTS_EXPORT_PLUGIN(Phong, "Modified Phong BRDF");
MTS_NAMESPACE_END

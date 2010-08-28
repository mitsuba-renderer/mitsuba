#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/renderer.h>

MTS_NAMESPACE_BEGIN


/**
 * Simple diffuse transmitter
 */
class DiffuseTransmitter : public BSDF {
public:
	DiffuseTransmitter(const Properties &props) 
		: BSDF(props) {
		m_transmittance = new ConstantTexture(
			props.getSpectrum("transmittance", Spectrum(.5f)));
		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_combinedType = m_type[0] = EDiffuseTransmission;
		m_usesRayDifferentials = false;
	}

	DiffuseTransmitter(Stream *stream, InstanceManager *manager) 
		: BSDF(stream, manager) {
		m_transmittance = static_cast<Texture *>(manager->getInstance(stream));
		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_combinedType = m_type[0] = EDiffuseTransmission;
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

	Spectrum sample(BSDFQueryRecord &bRec) const {
		if (!(bRec.typeMask & m_combinedType))
			return Spectrum(0.0f);
		bRec.wo = squareToHemispherePSA(bRec.sample);
		if (bRec.wi.z > 0)
			bRec.wo.z *= -1;
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseTransmission;
		if (Frame::cosTheta(bRec.wo) == 0)
			return Spectrum(0.0f);
		return m_transmittance->getValue(bRec.its) / std::abs(Frame::cosTheta(bRec.wo));
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf) const {
		if (!(bRec.typeMask & m_combinedType)) 
			return Spectrum(0.0f);
		bRec.wo = squareToHemispherePSA(bRec.sample);
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
		if (child->getClass()->derivesFrom(Texture::m_theClass) && name == "transmittance") {
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
		oss << "DiffuseTransmitter[transmittance=" << m_transmittance->toString() << "]";
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

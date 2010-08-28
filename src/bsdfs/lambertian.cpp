#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/renderer.h>

MTS_NAMESPACE_BEGIN


/**
 * Simple one-sided Lambertian (i.e. perfectly diffuse) material
 */
class Lambertian : public BSDF {
public:
	Lambertian(const Properties &props) 
		: BSDF(props) {
		m_reflectance = new ConstantTexture(
			props.getSpectrum("reflectance", Spectrum(.5f)));
		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_combinedType = m_type[0] = EDiffuseReflection;
		m_usesRayDifferentials = false;
	}

	Lambertian(Stream *stream, InstanceManager *manager) 
		: BSDF(stream, manager) {
		m_reflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_combinedType = m_type[0] = EDiffuseReflection;
		m_usesRayDifferentials = m_reflectance->usesRayDifferentials();
	}

	virtual ~Lambertian() {
		delete[] m_type;
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return m_reflectance->getValue(its);
	}

	Spectrum f(const BSDFQueryRecord &bRec) const {
		if (!(bRec.typeMask & m_combinedType)
			|| bRec.wi.z <= 0 || bRec.wo.z <= 0)
			return Spectrum(0.0f);

		return m_reflectance->getValue(bRec.its) * INV_PI;
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		if (bRec.wi.z <= 0 || bRec.wo.z <= 0)
			return 0.0f;
		return Frame::cosTheta(bRec.wo) * INV_PI;
	}

	Spectrum sample(BSDFQueryRecord &bRec) const {
		if (!(bRec.typeMask & m_combinedType) || bRec.wi.z <= 0)
			return Spectrum(0.0f);
		bRec.wo = squareToHemispherePSA(bRec.sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		if (Frame::cosTheta(bRec.wo) == 0)
			return Spectrum(0.0f);
		return m_reflectance->getValue(bRec.its) / Frame::cosTheta(bRec.wo);
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf) const {
		if (!(bRec.typeMask & m_combinedType) || bRec.wi.z <= 0)
			return Spectrum(0.0f);
		bRec.wo = squareToHemispherePSA(bRec.sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		pdf = Frame::cosTheta(bRec.wo) * INV_PI;
		if (Frame::cosTheta(bRec.wo) == 0)
			return Spectrum(0.0f);
		return m_reflectance->getValue(bRec.its) * INV_PI;
	}
		
	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(Texture::m_theClass) && name == "reflectance") {
			m_reflectance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_reflectance->usesRayDifferentials();
		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_reflectance.get());
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Lambertian[reflectance=" << m_reflectance->toString() << "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_reflectance;
};

// ================ Hardware shader implementation ================ 

class LambertianShader : public Shader {
public:
	LambertianShader(Renderer *renderer, const Texture *reflectance) 
		: Shader(renderer, EBSDFShader), m_reflectance(reflectance) {
		m_reflectanceShader = renderer->registerShaderForResource(m_reflectance.get());
	}

	bool isComplete() const {
		return m_reflectanceShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_reflectance.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_reflectanceShader.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (wi.z < 0.0 || wo.z < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return " << depNames[0] << "(uv) * 0.31831;" << endl
			<< "}" << endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << evalName << "(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_reflectance;
	ref<Shader> m_reflectanceShader;
};

Shader *Lambertian::createShader(Renderer *renderer) const { 
	return new LambertianShader(renderer, m_reflectance.get());
}

MTS_IMPLEMENT_CLASS(LambertianShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Lambertian, false, BSDF)
MTS_EXPORT_PLUGIN(Lambertian, "Lambertian BRDF")
MTS_NAMESPACE_END

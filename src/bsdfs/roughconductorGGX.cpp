
#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>
#include "ior.h"

MTS_NAMESPACE_BEGIN

#include "roughGGX.h"

vec3 samplePhaseFunction_conductor(const vec3& wi, const float alpha_x, const float alpha_y, const Spectrum& m_eta, const Spectrum& m_k, Spectrum& weight)
{
	const float U1 = generateRandomNumber();
	const float U2 = generateRandomNumber();

	// sample D_wi
	// stretch to match configuration with alpha=1.0	
	const vec3 wi_11 = normalize(vec3(alpha_x * wi.x, alpha_y * wi.y, wi.z));

	// sample visible slope with alpha=1.0
	vec2 slope_11 = sampleP22_11(acosf(wi_11.z), U1, U2, alpha_x, alpha_y);

	// align with view direction
	const float phi = atan2(wi_11.y, wi_11.x);
	vec2 slope(cosf(phi)*slope_11.x - sinf(phi)*slope_11.y, sinf(phi)*slope_11.x + cos(phi)*slope_11.y, 0);

	// stretch back
	slope.x *= alpha_x;
	slope.y *= alpha_y;

	// compute normal
	vec3 wm;
	// if numerical instability
	if( (slope.x != slope.x) || !IsFiniteNumber(slope.x) ) 
	{
		if(wi.z > 0) wm = vec3(0.0f,0.0f,1.0f);
		else wm = normalize(vec3(wi.x, wi.y, 0.0f));
	}
	else
		wm = normalize(vec3(-slope.x, -slope.y, 1.0f));

	// reflect
	const vec3 wo = -wi + 2.0f * wm * dot(wi, wm);
	weight = fresnelConductorExact(dot(wi, wm), m_eta, m_k);

	return wo;
}

Spectrum evalPhaseFunction_conductor(const RayInfo& ray, const vec3& wo, const float alpha_x, const float alpha_y, const Spectrum& m_eta, const Spectrum& m_k)
{
	if(ray.w.z > 0.9999f)
		return Spectrum(0.0f);

	// half vector 
	const vec3 wh = normalize(-ray.w+wo);
	if(wh.z < 0.0f)
		return Spectrum(0.0f);
	
	// projected area
	float projectedArea;
	if(ray.w.z < -0.9999f)
		projectedArea = 1.0f;
	else 
		projectedArea = ray.Lambda * ray.w.z;

	// value
	const Spectrum value = fresnelConductorExact(dot(-ray.w, wh), m_eta, m_k) * std::max(0.0f, dot(-ray.w, wh)) * D_ggx(wh, alpha_x, alpha_y) / 4.0f / projectedArea / dot(-ray.w, wh);
	return value;
}

// MIS weights for bidirectional path tracing on the microsurface
float MISweight_conductor(const vec3& wi, const vec3& wo, const float alpha_x, const float alpha_y) 
{
	if(wi.x == -wo.x && wi.y == -wo.y && wi.z == -wo.z)
		return 1.0f;
	const vec3 wh = normalize(wi+wo);
	const float value = D_ggx( (wh.z>0) ? wh : -wh , alpha_x, alpha_y);
	return value;
}

Spectrum eval_conductor(const vec3& wi, const vec3& wo, const float alpha_x, const float alpha_y, const Spectrum& m_eta, const Spectrum& m_k, const int scatteringOrderMax)
{
	if(wi.z <= 0 || wo.z <= 0)
		return Spectrum(0.0f);

	// init
	RayInfo ray;
	ray.updateDirection(-wi, alpha_x, alpha_y);	
	ray.updateHeight(1.0f);
	Spectrum energy(1.0f);

	RayInfo ray_shadowing;
	ray_shadowing.updateDirection(wo, alpha_x, alpha_y);

	// eval single scattering	
	// half-vector
	const vec3 wh = normalize(wi+wo);
	const float D = D_ggx(wh, alpha_x, alpha_y);
	const float G2 = 1.0f / (1.0f + (-ray.Lambda-1.0f) + ray_shadowing.Lambda);
	Spectrum singleScattering = fresnelConductorExact(dot(-ray.w, wh), m_eta, m_k)  *  D * G2 / (4.0f * wi.z);
	
	// MIS weight 
	float wi_MISweight;

	// multiple scattering
	Spectrum multipleScattering(0.0f);
	
	// random walk
	int current_scatteringOrder = 0;	
	while(current_scatteringOrder < scatteringOrderMax)
	{
		// next height
		float U = generateRandomNumber();
		ray.updateHeight( sampleHeight(ray, U) );		
				
		// leave the microsurface?
		if( ray.h == std::numeric_limits<Float>::max() )
			break;
		else
			current_scatteringOrder++;

		// next event estimation 
		if( current_scatteringOrder > 1) // single scattering is already computed
		{
			Spectrum phasefunction = evalPhaseFunction_conductor(ray, wo, alpha_x, alpha_y, m_eta, m_k); 
			ray_shadowing.updateHeight(ray.h);
			float shadowing = ray_shadowing.G1;
			Spectrum I = energy * phasefunction * shadowing;

			// MIS
			const float MIS = wi_MISweight / ( wi_MISweight + MISweight_conductor(-ray.w, wo, alpha_x, alpha_y) );


			if ( IsFiniteNumber(I[0]) )
				multipleScattering += I * MIS;
		}

		// next direction
		Spectrum weight;
		ray.updateDirection(samplePhaseFunction_conductor(-ray.w, alpha_x, alpha_y, m_eta, m_k, weight), alpha_x, alpha_y);
		energy = energy * weight;
		ray.updateHeight(ray.h);

		if(current_scatteringOrder == 1)
			wi_MISweight = MISweight_conductor(wi, ray.w, alpha_x, alpha_y);

		// if NaN (should not happen, just in case)
		if( (ray.h != ray.h) || (ray.w.x != ray.w.x)) 
			return Spectrum(0.0f);
	}

	// 0.5f = MIS weight of singleScattering
	// multipleScattering already weighted by MIS
	return 0.5f*singleScattering + multipleScattering;
}

vec3 sample_conductor(const vec3& wi, const float alpha_x, const float alpha_y, const Spectrum& m_eta, const Spectrum& m_k, const int scatteringOrderMax, Spectrum& energy)
{
	energy = Spectrum(1.0f);

	// init
	RayInfo ray;
	ray.updateDirection(-wi, alpha_x, alpha_y);
	ray.updateHeight(1.0f);
		
	// random walk
	int current_scatteringOrder = 0;
	while(true)
	{
		// next height
		float U = generateRandomNumber();
		ray.updateHeight( sampleHeight(ray, U) );		

		// leave the microsurface?
		if( ray.h == std::numeric_limits<Float>::max() )
			break;
		else
			current_scatteringOrder++;

		// next direction
		Spectrum weight;
		ray.updateDirection(samplePhaseFunction_conductor(-ray.w, alpha_x, alpha_y, m_eta, m_k, weight), alpha_x, alpha_y);
		energy = energy * weight;
		ray.updateHeight(ray.h);

		// if NaN (should not happen, just in case)
		if( (ray.h != ray.h) || (ray.w.x != ray.w.x)) 
		{
			energy = Spectrum(0.0f);
			return vec3(0,0,1);
		}

		if( current_scatteringOrder > scatteringOrderMax )
		{
			energy = Spectrum(0.0f);
			return vec3(0,0,1);
		}
	}

	return ray.w;
}



class RoughConductor : public BSDF {
public:
	RoughConductor(const Properties &props) : BSDF(props) {
		ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));

		std::string materialName = props.getString("material", "Cu");

		Spectrum intEta, intK;
		if (boost::to_lower_copy(materialName) == "none") {
			intEta = Spectrum(0.0f);
			intK = Spectrum(1.0f);
		} else {
			intEta.fromContinuousSpectrum(InterpolatedSpectrum(
				fResolver->resolve("data/ior/" + materialName + ".eta.spd")));
			intK.fromContinuousSpectrum(InterpolatedSpectrum(
				fResolver->resolve("data/ior/" + materialName + ".k.spd")));
		}

		Float extEta = lookupIOR(props, "extEta", "air");

		m_eta = props.getSpectrum("eta", intEta) / extEta;
		m_k   = props.getSpectrum("k", intK) / extEta;

		// roughness
		float alphaU, alphaV;
		if (props.hasProperty("alpha")) {
			alphaU = alphaV = props.getFloat("alpha");
			if (props.hasProperty("alphaU") || props.hasProperty("alphaV"))
				SLog(EError, "Microfacet model: please specify either 'alpha' or 'alphaU'/'alphaV'.");
		} else if (props.hasProperty("alphaU") || props.hasProperty("alphaV")) {
			if (!props.hasProperty("alphaU") || !props.hasProperty("alphaV"))
				SLog(EError, "Microfacet model: both 'alphaU' and 'alphaV' must be specified.");
			if (props.hasProperty("alpha"))
				SLog(EError, "Microfacet model: please specify either 'alpha' or 'alphaU'/'alphaV'.");
			alphaU = props.getFloat("alphaU");
			alphaV = props.getFloat("alphaV");
		}
		m_alphaU = new ConstantFloatTexture(alphaU);
		if (alphaU == alphaV)
			m_alphaV = m_alphaU;
		else
			m_alphaV = new ConstantFloatTexture(alphaV);

		// scattering order
		m_scatteringOrderMax = props.getInteger("scatteringOrderMax", 10);
	}

	RoughConductor(Stream *stream, InstanceManager *manager)
	 : BSDF(stream, manager) {
		m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaV = static_cast<Texture *>(manager->getInstance(stream));
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_eta = Spectrum(stream);
		m_k = Spectrum(stream);
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_alphaU.get());
		manager->serialize(stream, m_alphaV.get());
		manager->serialize(stream, m_specularReflectance.get());
		m_eta.serialize(stream);
		m_k.serialize(stream);
	}

	void configure() {
		unsigned int extraFlags = 0;
		if (m_alphaU != m_alphaV)
			extraFlags |= EAnisotropic;

		if (!m_alphaU->isConstant() || !m_alphaV->isConstant() ||
			!m_specularReflectance->isConstant())
			extraFlags |= ESpatiallyVarying;

		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide | extraFlags);

		/* Verify the input parameters and fix them if necessary */
		m_specularReflectance = ensureEnergyConservation(
			m_specularReflectance, "specularReflectance", 1.0f);

		m_usesRayDifferentials =
			m_alphaU->usesRayDifferentials() ||
			m_alphaV->usesRayDifferentials() ||
			m_specularReflectance->usesRayDifferentials();

		BSDF::configure();
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		/* Stop if this component was not requested */
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);
		
		vec3 wi(bRec.wi.x, bRec.wi.y, bRec.wi.z);
		vec3 wo(bRec.wo.x, bRec.wo.y, bRec.wo.z);

		const float alpha_x = std::max(m_alphaU->eval(bRec.its).average(), (Float) 1e-4f);
		const float alpha_y = std::max(m_alphaV->eval(bRec.its).average(), (Float) 1e-4f);

		// start random walk from either wi or wo randomly (bidirectional path tracing)
		Spectrum res = (generateRandomNumber() > 0.5f) ? 
			2.0f * eval_conductor(wi, wo, alpha_x, alpha_y, m_eta, m_k, m_scatteringOrderMax) :
			2.0f * eval_conductor(wo, wi, alpha_x, alpha_y, m_eta, m_k, m_scatteringOrderMax)/Frame::cosTheta(bRec.wi)*Frame::cosTheta(bRec.wo);

		return res;
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return 0.0f;
		
		// Calculate the reflection half-vector 
		Vector wh = normalize(bRec.wo+bRec.wi);

		RayInfo ray;
		vec3 wi(bRec.wi.x, bRec.wi.y, bRec.wi.z);
		const float alpha_x = std::max(m_alphaU->eval(bRec.its).average(), (Float) 1e-4f);
		const float alpha_y = std::max(m_alphaV->eval(bRec.its).average(), (Float) 1e-4f);
		ray.updateDirection(wi, alpha_x, alpha_y);

		// single-scattering PDF + diffuse 
		// otherwise too many fireflies due to lack of multiple-scattering PDF
		// (MIS works even if the PDF is wrong and not normalized)
		return D_ggx(wh, alpha_x, alpha_y) / (1.0f + ray.Lambda) / (4.0f * Frame::cosTheta(bRec.wi)) + Frame::cosTheta(bRec.wo);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		Float pdf;
		return this->sample(bRec, pdf, sample);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;	

		vec3 wi(bRec.wi.x, bRec.wi.y, bRec.wi.z);

		const float alpha_x = std::max(m_alphaU->eval(bRec.its).average(), (Float) 1e-4f);
		const float alpha_y = std::max(m_alphaV->eval(bRec.its).average(), (Float) 1e-4f);

		Spectrum energy;
		vec3 wo = sample_conductor(wi, alpha_x, alpha_y, m_eta, m_k, m_scatteringOrderMax, energy);
		bRec.wo.x = wo.x;
		bRec.wo.y = wo.y;
		bRec.wo.z = wo.z;

		pdf = this->pdf(bRec, ESolidAngle);
		
		return energy;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "alpha")
				m_alphaU = m_alphaV = static_cast<Texture *>(child);
			else if (name == "alphaU")
				m_alphaU = static_cast<Texture *>(child);
			else if (name == "alphaV")
				m_alphaV = static_cast<Texture *>(child);
			else if (name == "specularReflectance")
				m_specularReflectance = static_cast<Texture *>(child);
			else
				BSDF::addChild(name, child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	Float getRoughness(const Intersection &its, int component) const {
		return 0.5f * (m_alphaU->eval(its).average()
			+ m_alphaV->eval(its).average());
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "RoughConductor[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  alphaU = " << indent(m_alphaU->toString()) << "," << endl
			<< "  alphaV = " << indent(m_alphaV->toString()) << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  eta = " << m_eta.toString() << "," << endl
			<< "  k = " << m_k.toString() << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_specularReflectance;
	ref<Texture> m_alphaU, m_alphaV;
	Spectrum m_eta, m_k;
	int m_scatteringOrderMax;
};

/**
 * GLSL port of the rough conductor shader. This version is much more
 * approximate -- it only supports the Ashikhmin-Shirley distribution,
 * does everything in RGB, and it uses the Schlick approximation to the
 * Fresnel reflectance of conductors. When the roughness is lower than
 * \alpha < 0.2, the shader clamps it to 0.2 so that it will still perform
 * reasonably well in a VPL-based preview.
 */
class RoughConductorShader : public Shader {
public:
	RoughConductorShader(Renderer *renderer, const Texture *specularReflectance,
			const Texture *alphaU, const Texture *alphaV, const Spectrum &eta,
			const Spectrum &k) : Shader(renderer, EBSDFShader),
			m_specularReflectance(specularReflectance), m_alphaU(alphaU), m_alphaV(alphaV){
		m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());
		m_alphaUShader = renderer->registerShaderForResource(m_alphaU.get());
		m_alphaVShader = renderer->registerShaderForResource(m_alphaV.get());

		/* Compute the reflectance at perpendicular incidence */
		m_R0 = fresnelConductorExact(1.0f, eta, k);
	}

	bool isComplete() const {
		return m_specularReflectanceShader.get() != NULL &&
			   m_alphaUShader.get() != NULL &&
			   m_alphaVShader.get() != NULL;
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_specularReflectanceShader.get());
		deps.push_back(m_alphaUShader.get());
		deps.push_back(m_alphaVShader.get());
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_specularReflectance.get());
		renderer->unregisterShaderForResource(m_alphaU.get());
		renderer->unregisterShaderForResource(m_alphaV.get());
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_R0", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_R0);
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform vec3 " << evalName << "_R0;" << endl
			<< endl
			<< "float " << evalName << "_D(vec3 m, float alphaU, float alphaV) {" << endl
			<< "    float ct = cosTheta(m), ds = 1-ct*ct;" << endl
			<< "    if (ds <= 0.0)" << endl
			<< "        return 0.0f;" << endl
			<< "    alphaU = 2 / (alphaU * alphaU) - 2;" << endl
			<< "    alphaV = 2 / (alphaV * alphaV) - 2;" << endl
			<< "    float exponent = (alphaU*m.x*m.x + alphaV*m.y*m.y)/ds;" << endl
			<< "    return sqrt((alphaU+2) * (alphaV+2)) * 0.15915 * pow(ct, exponent);" << endl
			<< "}" << endl
			<< endl
			<< "float " << evalName << "_G(vec3 m, vec3 wi, vec3 wo) {" << endl
			<< "    if ((dot(wi, m) * cosTheta(wi)) <= 0 || " << endl
			<< "        (dot(wo, m) * cosTheta(wo)) <= 0)" << endl
			<< "        return 0.0;" << endl
			<< "    float nDotM = cosTheta(m);" << endl
			<< "    return min(1.0, min(" << endl
			<< "        abs(2 * nDotM * cosTheta(wo) / dot(wo, m))," << endl
			<< "        abs(2 * nDotM * cosTheta(wi) / dot(wi, m))));" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_schlick(float ct) {" << endl
			<< "    float ctSqr = ct*ct, ct5 = ctSqr*ctSqr*ct;" << endl
			<< "    return " << evalName << "_R0 + (vec3(1.0) - " << evalName << "_R0) * ct5;" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "   if (cosTheta(wi) <= 0 || cosTheta(wo) <= 0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "   vec3 H = normalize(wi + wo);" << endl
			<< "   vec3 reflectance = " << depNames[0] << "(uv);" << endl
			<< "   float alphaU = max(0.2, " << depNames[1] << "(uv).r);" << endl
			<< "   float alphaV = max(0.2, " << depNames[2] << "(uv).r);" << endl
			<< "   float D = " << evalName << "_D(H, alphaU, alphaV)" << ";" << endl
			<< "   float G = " << evalName << "_G(H, wi, wo);" << endl
			<< "   vec3 F = " << evalName << "_schlick(1-dot(wi, H));" << endl
			<< "   return reflectance * F * (D * G / (4*cosTheta(wi)));" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return " << evalName << "_R0 * inv_pi * inv_pi * cosTheta(wo);"<< endl
			<< "}" << endl;
	}
	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_specularReflectance;
	ref<const Texture> m_alphaU;
	ref<const Texture> m_alphaV;
	ref<Shader> m_specularReflectanceShader;
	ref<Shader> m_alphaUShader;
	ref<Shader> m_alphaVShader;
	Spectrum m_R0;
};

Shader *RoughConductor::createShader(Renderer *renderer) const {
	return new RoughConductorShader(renderer,
		m_specularReflectance.get(), m_alphaU.get(), m_alphaV.get(), m_eta, m_k);
}

MTS_IMPLEMENT_CLASS(RoughConductorShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(RoughConductor, false, BSDF)
MTS_EXPORT_PLUGIN(RoughConductor, "Rough conductor BRDF");
MTS_NAMESPACE_END

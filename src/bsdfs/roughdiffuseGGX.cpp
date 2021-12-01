
#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

#include "roughGGX.h"







 // build orthonormal basis (Building an Orthonormal Basis from a 3D Unit Vector Without Normalization, [Frisvad2012])
 void buildOrthonormalBasis(vec3& omega_1, vec3& omega_2, const vec3& omega_3)
{
	if(omega_3.z < -0.9999999f) 
	{
	   omega_1 = vec3 ( 0.0f , -1.0f , 0.0f );
	   omega_2 = vec3 ( -1.0f , 0.0f , 0.0f );
	} else {
	   const float a = 1.0f /(1.0f + omega_3.z );
	   const float b = -omega_3.x*omega_3 .y*a ;
	   omega_1 = vec3 (1.0f - omega_3.x*omega_3. x*a , b , -omega_3.x );
	   omega_2 = vec3 (b , 1.0f - omega_3.y*omega_3.y*a , -omega_3.y );
	}
}

vec3 samplePhaseFunction_diffuse(const vec3& wm)
{
	const float U1 = generateRandomNumber();
	const float U2 = generateRandomNumber();

	// sample diffuse reflection
	vec3 w1, w2;
	buildOrthonormalBasis(w1, w2, wm);

	float r1 = 2.0f*U1 - 1.0f;
	float r2 = 2.0f*U2 - 1.0f;

	// concentric map code from
	// http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html
	float phi, r;
	if (r1 == 0 && r2 == 0) {
		r = phi = 0;
	} else if (r1*r1 > r2*r2) {
		r = r1;
		phi = (M_PI/4.0f) * (r2/r1);
	} else {
		r = r2;
		phi = (M_PI/2.0f) - (r1/r2) * (M_PI/4.0f);
	}
	float x = r*cosf(phi);
	float y = r*sinf(phi);
	float z = sqrtf(std::max(0.0f, 1.0f - x*x - y*y));
	vec3 wo = x*w1 + y*w2 + z*wm;

	return wo;
}

Spectrum eval_diffuse(const vec3& wi, const vec3& wo, const float alpha_x, const float alpha_y, const Spectrum& albedo, const int scatteringOrderMax)
{
	if(wi.z <= 0 || wo.z <= 0)
		return Spectrum(0.0f);

	// init
	RayInfo ray;
	ray.updateDirection(-wi, alpha_x, alpha_y);	

	RayInfo ray_shadowing;
	ray_shadowing.updateDirection(wo, alpha_x, alpha_y);

	Spectrum res(0.0f);

	ray.updateHeight(1.0f);
	Spectrum energy(1.0f);

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

		// sample VNDF
		vec3 wm = sampleVNDF(-ray.w, alpha_x, alpha_y);

		// next event estimation
		Spectrum phasefunction = albedo * std::max(0.0f, dot(wm, wo)) / M_PI;
		if(current_scatteringOrder == 1) // closed masking and shadowing (we compute G2 / G1 because G1 is already in the phase function)
		{			
			float G2_G1 = (1.0f + (-ray.Lambda-1.0f)) / (1.0f + (-ray.Lambda-1.0f) + ray_shadowing.Lambda);
			Spectrum I = energy * phasefunction * G2_G1;
			if ( IsFiniteNumber(I[0]) )
				res += I;
		}
		else
		{
			Spectrum phasefunction = albedo * std::max(0.0f, dot(wm, wo)) / M_PI;
			ray_shadowing.updateHeight(ray.h);
			float shadowing = ray_shadowing.G1;
			Spectrum I = energy * phasefunction * shadowing;
			if ( IsFiniteNumber(I[0]) )
				res += I;
		}

		// next direction
		ray.updateDirection(samplePhaseFunction_diffuse(wm), alpha_x, alpha_y);
		energy = energy * albedo;
		ray.updateHeight(ray.h);

		// if NaN (should not happen, just in case)
		if( (ray.h != ray.h) || (ray.w.x != ray.w.x)) 
			return Spectrum(0.0f);
	}

	return res;
}

vec3 sample_diffuse(const vec3& wi, const float alpha_x, const float alpha_y, const Spectrum& albedo, const int scatteringOrderMax, Spectrum& energy)
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
		if( ray.h == std::numeric_limits<Float>::max())
			break;
		else
			current_scatteringOrder++;

		// sample VNDF
		vec3 wm = sampleVNDF(-ray.w, alpha_x, alpha_y);

		// next direction
		ray.updateDirection(samplePhaseFunction_diffuse(wm), alpha_x, alpha_y);
		energy = energy * albedo;
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



class RoughDiffuse : public BSDF {
public:
	RoughDiffuse(const Properties &props) : BSDF(props) {
		ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

		/* For better compatibility with other models, support both
		   'reflectance' and 'diffuseReflectance' as parameter names */
		m_reflectance = new ConstantSpectrumTexture(props.getSpectrum(
			props.hasProperty("reflectance") ? "reflectance"
				: "diffuseReflectance", Spectrum(1.0f)));

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

	RoughDiffuse(Stream *stream, InstanceManager *manager)
	 : BSDF(stream, manager) {
		m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaV = static_cast<Texture *>(manager->getInstance(stream));
		m_reflectance = static_cast<Texture *>(manager->getInstance(stream));

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_alphaU.get());
		manager->serialize(stream, m_alphaV.get());
		manager->serialize(stream, m_reflectance.get());
	}

	void configure() {
		unsigned int extraFlags = 0;
		if (m_alphaU != m_alphaV)
			extraFlags |= EAnisotropic;

		if (!m_alphaU->isConstant() || !m_alphaV->isConstant() ||
			!m_reflectance->isConstant())
			extraFlags |= ESpatiallyVarying;

		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide | extraFlags);

		/* Verify the input parameters and fix them if necessary */
		m_reflectance = ensureEnergyConservation(
			m_reflectance, "reflectance", 1.0f);

		m_usesRayDifferentials =
			m_alphaU->usesRayDifferentials() ||
			m_alphaV->usesRayDifferentials() ||
			m_reflectance->usesRayDifferentials();

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
		
		Spectrum albedo = m_reflectance->eval(bRec.its);
		vec3 wi(bRec.wi.x, bRec.wi.y, bRec.wi.z);
		vec3 wo(bRec.wo.x, bRec.wo.y, bRec.wo.z);

		const float alpha_x = std::max(m_alphaU->eval(bRec.its).average(), (Float) 1e-4f);
		const float alpha_y = std::max(m_alphaV->eval(bRec.its).average(), (Float) 1e-4f);

		// start random walks from lower direction for eval
		Spectrum res = (wi.z < wo.z) ?
			eval_diffuse(wi, wo, alpha_x, alpha_y, albedo, m_scatteringOrderMax) :
			eval_diffuse(wo, wi, alpha_x, alpha_y, albedo, m_scatteringOrderMax)/Frame::cosTheta(bRec.wi)*Frame::cosTheta(bRec.wo);

		return res;
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return 0.0f;
		
		return 1.0f / M_PI * Frame::cosTheta(bRec.wo);
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

		Spectrum albedo = m_reflectance->eval(bRec.its);
		Spectrum energy;
		vec3 wo = sample_diffuse(wi, alpha_x, alpha_y, albedo, m_scatteringOrderMax, energy);
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
			else if (name == "reflectance")
				m_reflectance = static_cast<Texture *>(child);
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
			<< "  specularReflectance = " << indent(m_reflectance->toString()) << "," << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_reflectance;
	ref<Texture> m_alphaU, m_alphaV;
	int m_scatteringOrderMax;
};
// ================ Hardware shader implementation ================

class SmoothDiffuseShader : public Shader {
public:
	SmoothDiffuseShader(Renderer *renderer, const Texture *reflectance)
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
			<< "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return " << depNames[0] << "(uv) * inv_pi * cosTheta(wo);" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << evalName << "(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_reflectance;
	ref<Shader> m_reflectanceShader;
};

Shader *RoughDiffuse::createShader(Renderer *renderer) const {
	return new SmoothDiffuseShader(renderer, m_reflectance.get());
}

MTS_IMPLEMENT_CLASS(SmoothDiffuseShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(RoughDiffuse, false, BSDF)
MTS_EXPORT_PLUGIN(RoughDiffuse, "Smooth diffuse BRDF")
MTS_NAMESPACE_END
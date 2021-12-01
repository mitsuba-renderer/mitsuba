
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/hw/basicshader.h>
#include "ior.h"

MTS_NAMESPACE_BEGIN

#include "roughGGX.h"



vec3 refract(const vec3 &wi, const vec3 &wm, const float eta)  
{
	const float cos_theta_i = dot(wi, wm);
	const float cos_theta_t2 = 1.0f - (1.0f-cos_theta_i*cos_theta_i) / (eta*eta);
	const float cos_theta_t = -sqrtf(std::max(0.0f,cos_theta_t2));

	return wm * (dot(wi, wm) / eta + cos_theta_t) - wi / eta;
}

float Fresnel(const vec3& wi, const vec3& wm, const float eta) 
{	
	const float cos_theta_i = dot(wi, wm);
	const float cos_theta_t2 = 1.0f - (1.0f-cos_theta_i*cos_theta_i) / (eta*eta);

	// total internal reflection 
	if (cos_theta_t2 <= 0.0f) return 1.0f;

	const float cos_theta_t = sqrtf(cos_theta_t2);

	const float Rs = (cos_theta_i - eta * cos_theta_t) / (cos_theta_i + eta * cos_theta_t);
	const float Rp = (eta * cos_theta_i - cos_theta_t) / (eta * cos_theta_i + cos_theta_t);

	const float F = 0.5f * (Rs * Rs + Rp * Rp);
	return F;
}

// by convention, ray is always outside
float evalPhaseFunction_dielectric(const RayInfo& ray, const vec3& wo, const bool wo_outside, const float eta, const float alpha_x, const float alpha_y) 
{
	if(ray.w.z > 0.9999f)
		return 0.0f;

	// projected area
	float projectedArea;
	if(ray.w.z < -0.9999f)
		projectedArea = 1.0f;
	else 
		projectedArea = ray.Lambda * ray.w.z;

	if( wo_outside ) // reflection
	{
		// half vector 
		const vec3 wh = normalize(-ray.w+wo);
		if(wh.z < 0.0f)
			return 0.0f;

		// value
		const float value = Fresnel(-ray.w, wh, eta) * std::max(0.0f, dot(-ray.w, wh)) * D_ggx(wh, alpha_x, alpha_y) / 4.0f / projectedArea / dot(-ray.w, wh);

		return value;
	}
	else // transmission
	{
		vec3 wh = normalize(-ray.w+wo*eta);
		wh *= (wh.z>0)?1.0f:-1.0f;

		if(dot(wh, -ray.w) < 0)
			return 0;

		const float value = eta*eta * (1.0f-Fresnel(-ray.w, wh, eta)) *
				std::max(0.0f, dot(-ray.w, wh)) * D_ggx(wh, alpha_x, alpha_y) / projectedArea * 
				std::max(0.0f, -dot(wo, wh)) *
				1.0f / powf(dot(-ray.w, wh)+eta*dot(wo,wh), 2.0f);

		return value;	
	}
}

// by convention, wi is always outside
vec3 samplePhaseFunction_dielectric(const vec3& wi, const float alpha_x, const float alpha_y, const float eta, bool& wo_outside) 
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
	vec2 slope(cosf(phi)*slope_11.x - sinf(phi)*slope_11.y, sinf(phi)*slope_11.x + cos(phi)*slope_11.y, 0.0f);

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

	const float F = Fresnel(wi, wm, eta);

	if( generateRandomNumber() < F )
	{
		wo_outside = true;
		const vec3 wo = -wi + 2.0f * wm * dot(wi, wm); // reflect
		return wo;
	}
	else
	{
		wo_outside = false;
		const vec3 wo = refract(wi, wm, eta);
		return normalize(wo);
	}
}

// MIS weights for bidirectional path tracing on the microsurface
float MISweight_dielectric(const vec3& wi, const vec3& wo, const bool wo_outside, const float eta, const float alpha_x, const float alpha_y) 
{
	if( wo_outside ) // reflection
	{
		if(wi.x == -wo.x && wi.y == -wo.y && wi.z == -wo.z)
			return 1.0f;
		const vec3 wh = normalize(wi+wo);
		const float value = D_ggx( (wh.z>0) ? wh : -wh , alpha_x, alpha_y);
		return value;
	}
	else // transmission
	{
		const vec3 wh = normalize(wi+wo*eta);
		const float value = D_ggx( (wh.z>0) ? wh : -wh , alpha_x, alpha_y);
		return value;	
	}
}

float eval_dielectric(const vec3& wi, const vec3& wo, const bool wo_outside, const float alpha_x, const float alpha_y, const float eta, const int scatteringOrderMax)
{
	if( (wi.z <= 0) || (wo.z <= 0 && wo_outside) || (wo.z >= 0 && !wo_outside))
		return 0.0f;

	// init
	RayInfo ray;
	ray.updateDirection(-wi, alpha_x, alpha_y);	
	ray.updateHeight(1.0f);
	bool outside = true;

	RayInfo ray_shadowing;
	if(wo_outside)
		ray_shadowing.updateDirection(wo, alpha_x, alpha_y);
	else
		ray_shadowing.updateDirection(-wo, alpha_x, alpha_y);

	float singleScattering = 0;
	float multipleScattering = 0;

	float wi_MISweight;

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
		if(current_scatteringOrder == 1) // single scattering
		{
				float phasefunction = evalPhaseFunction_dielectric(ray, wo, wo_outside, eta, alpha_x, alpha_y); 
	
				// closed masking and shadowing (we compute G2 / G1 because G1 is already in the phase function)
				float G2_G1;
				if( wo_outside )
					G2_G1 = (1.0f + (-ray.Lambda-1.0f)) / (1.0f + (-ray.Lambda-1.0f) + ray_shadowing.Lambda);
				else
					G2_G1 = (1.0f + (-ray.Lambda-1.0f)) * (float) beta(1.0f + (-ray.Lambda-1.0f), 1.0f + ray_shadowing.Lambda);

				float I = phasefunction * G2_G1;
				if ( IsFiniteNumber(I) )
					singleScattering = I;
		}
		if(current_scatteringOrder > 1) // multiple scattering
		{
				float phasefunction;
				float MIS;
				if( outside )
				{
					phasefunction = evalPhaseFunction_dielectric(ray, wo, wo_outside, eta, alpha_x, alpha_y); 
					MIS = wi_MISweight / (wi_MISweight + MISweight_dielectric(-ray.w, wo, wo_outside, eta, alpha_x, alpha_y));
				}
				else
				{
					phasefunction = evalPhaseFunction_dielectric(ray, -wo, !wo_outside, 1.0f/eta, alpha_x, alpha_y); 
					MIS = wi_MISweight / (wi_MISweight + MISweight_dielectric(-ray.w, -wo, !wo_outside, 1.0f/eta, alpha_x, alpha_y));
				}

				if( outside == wo_outside )
					ray_shadowing.updateHeight(ray.h);
				else
					ray_shadowing.updateHeight(-ray.h);

				const float shadowing = ray_shadowing.G1;
				float I = phasefunction * shadowing;
				if ( IsFiniteNumber(I) )
					multipleScattering += I * MIS;
		}

		// next direction
		bool next_outside;
		vec3 w = samplePhaseFunction_dielectric(-ray.w, alpha_x, alpha_y, (outside ? eta:1.0f/eta), next_outside);
		if (next_outside)
		{
			ray.updateDirection(w, alpha_x, alpha_y);
			ray.updateHeight(ray.h);
		}
		else
		{
			outside = !outside;
			ray.updateDirection(-w, alpha_x, alpha_y);
			ray.updateHeight(-ray.h);
		}
		
		if(current_scatteringOrder == 1)
			wi_MISweight = MISweight_dielectric(wi, ray.w, outside, eta, alpha_x, alpha_y);

		// if NaN (should not happen, just in case)
		if( (ray.h != ray.h) || (ray.w.x != ray.w.x)) 
			return 0.0f;
	}

	// 0.5f = MIS weight of singleScattering
	// multipleScattering already weighted by MIS
	return 0.5f * singleScattering + multipleScattering;
}

vec3 sample_dielectric(const vec3& wi, const float alpha_x, const float alpha_y, const float eta, const int scatteringOrderMax, float& weight)
{
	// init
	RayInfo ray;
	ray.updateDirection(-wi, alpha_x, alpha_y);
	ray.updateHeight(1.0f);
	bool outside = true;
	
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
		bool next_outside;
		vec3 w = samplePhaseFunction_dielectric(-ray.w, alpha_x, alpha_y, (outside ? eta:1.0f/eta), next_outside);
		if (next_outside)
		{
			ray.updateDirection(w, alpha_x, alpha_y);
			ray.updateHeight(ray.h);
		}
		else
		{
			outside = !outside;
			ray.updateDirection(-w, alpha_x, alpha_y);
			ray.updateHeight(-ray.h);
		}
					
		// if NaN (should not happen, just in case)
		if( (ray.h != ray.h) || (ray.w.x != ray.w.x)) 
		{
			weight = 0.0f;
			return vec3(0,0,1);
		}

		if( current_scatteringOrder > scatteringOrderMax )
		{
			weight = 0.0f;
			return vec3(0,0,1);
		}

	}

	weight = 1.0f;
	return (outside) ? ray.w : -ray.w;
}




class RoughDielectric : public BSDF {
public:
	RoughDielectric(const Properties &props) : BSDF(props) {
		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));
		m_specularTransmittance = new ConstantSpectrumTexture(
			props.getSpectrum("specularTransmittance", Spectrum(1.0f)));

		/* Specifies the internal index of refraction at the interface */
		Float intIOR = lookupIOR(props, "intIOR", "bk7");

		/* Specifies the external index of refraction at the interface */
		Float extIOR = lookupIOR(props, "extIOR", "air");

		if (intIOR < 0 || extIOR < 0 || intIOR == extIOR)
			Log(EError, "The interior and exterior indices of "
				"refraction must be positive and differ!");

		m_eta = intIOR / extIOR;
		m_invEta = 1 / m_eta;
		
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

	RoughDielectric(Stream *stream, InstanceManager *manager)
	 : BSDF(stream, manager) {
		m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaV = static_cast<Texture *>(manager->getInstance(stream));
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_specularTransmittance = static_cast<Texture *>(manager->getInstance(stream));
		m_eta = stream->readFloat();
		m_invEta = 1 / m_eta;

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_alphaU.get());
		manager->serialize(stream, m_alphaV.get());
		manager->serialize(stream, m_specularReflectance.get());
		manager->serialize(stream, m_specularTransmittance.get());
		stream->writeFloat(m_eta);
	}

	void configure() {
		unsigned int extraFlags = 0;
		if (m_alphaU != m_alphaV)
			extraFlags |= EAnisotropic;

		if (!m_alphaU->isConstant() || !m_alphaV->isConstant())
			extraFlags |= ESpatiallyVarying;

		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide
			| EBackSide | EUsesSampler | extraFlags
			| (m_specularReflectance->isConstant() ? 0 : ESpatiallyVarying));
		m_components.push_back(EGlossyTransmission | EFrontSide
			| EBackSide | EUsesSampler | ENonSymmetric | extraFlags
			| (m_specularTransmittance->isConstant() ? 0 : ESpatiallyVarying));

		/* Verify the input parameters and fix them if necessary */
		m_specularReflectance = ensureEnergyConservation(
			m_specularReflectance, "specularReflectance", 1.0f);
		m_specularTransmittance = ensureEnergyConservation(
			m_specularTransmittance, "specularTransmittance", 1.0f);

		m_usesRayDifferentials =
			m_alphaU->usesRayDifferentials() ||
			m_alphaV->usesRayDifferentials() ||
			m_specularReflectance->usesRayDifferentials() ||
			m_specularTransmittance->usesRayDifferentials();

		BSDF::configure();
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (measure != ESolidAngle || Frame::cosTheta(bRec.wi) == 0)
			return Spectrum(0.0f);
		
		vec3 wi(bRec.wi.x, bRec.wi.y, bRec.wi.z);
		vec3 wo(bRec.wo.x, bRec.wo.y, bRec.wo.z);

		const float alpha_x = std::max(m_alphaU->eval(bRec.its).average(), (Float) 1e-4f);
		const float alpha_y = std::max(m_alphaV->eval(bRec.its).average(), (Float) 1e-4f);

		Float factor = (bRec.mode == ERadiance)
				? (Frame::cosTheta(bRec.wi) > 0 ? m_invEta : m_eta) : 1.0f;

		if ( Frame::cosTheta(bRec.wi) > 0 ) // outside
		{
			if( Frame::cosTheta(bRec.wo) >= 0 )
			{
				float value = (generateRandomNumber() > 0.5f) ? 
							2.0f * eval_dielectric(wi, wo, true , alpha_x, alpha_y, m_eta, m_scatteringOrderMax) :
							2.0f * eval_dielectric(wo, wi, true, alpha_x, alpha_y, m_eta, m_scatteringOrderMax) / Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo);
				return Spectrum(value);
			}
			else
			{
				float value = (generateRandomNumber() > 0.5f) ? 
							2.0f * eval_dielectric(wi, wo, false , alpha_x, alpha_y, m_eta, m_scatteringOrderMax) :
							2.0f * eval_dielectric(-wo, -wi, false, alpha_x, alpha_y, m_invEta, m_scatteringOrderMax) / Frame::cosTheta(bRec.wi) * Frame::cosTheta(-bRec.wo) / (factor*factor);
				return Spectrum(value) * ((Frame::cosTheta(bRec.wo) > 0) ? 1.0f : factor * factor);
			}			
		}
		else
		{
			if(Frame::cosTheta(bRec.wo) <= 0)
			{
				float value = (generateRandomNumber() > 0.5f) ? 
							2.0f * eval_dielectric(-wi, -wo, true , alpha_x, alpha_y, m_invEta, m_scatteringOrderMax) :
							2.0f * eval_dielectric(-wo, -wi, true, alpha_x, alpha_y, m_invEta, m_scatteringOrderMax) / Frame::cosTheta(-bRec.wi) * Frame::cosTheta(-bRec.wo);
				return Spectrum(value);
			}
			else
			{
				float value = (generateRandomNumber() > 0.5f) ? 
							2.0f * eval_dielectric(-wi, -wo, false , alpha_x, alpha_y, m_invEta, m_scatteringOrderMax) :
							2.0f * eval_dielectric(wo, wi, false , alpha_x, alpha_y, m_eta, m_scatteringOrderMax)/ Frame::cosTheta(-bRec.wi) * Frame::cosTheta(bRec.wo) / (factor*factor);
				return Spectrum(value) * ((Frame::cosTheta(bRec.wo) < 0) ? 1.0f : factor * factor);
			}
		}
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (measure != ESolidAngle)
			return 0.0f;

		/* Determine the type of interaction */
		bool hasReflection   = ((bRec.component == -1 || bRec.component == 0)
							  && (bRec.typeMask & EGlossyReflection)),
		     hasTransmission = ((bRec.component == -1 || bRec.component == 1)
							  && (bRec.typeMask & EGlossyTransmission)),
		     reflect         = Frame::cosTheta(bRec.wi)
				             * Frame::cosTheta(bRec.wo) > 0;

		Vector wh;
		Float dwh_dwo;

		if (reflect) {
			/* Zero probability if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 0)
				|| !(bRec.typeMask & EGlossyReflection))
				return 0.0f;

			/* Calculate the reflection half-vector */
			wh = normalize(bRec.wo+bRec.wi);

			/* Jacobian of the half-direction mapping */
			dwh_dwo = 1.0f / (4.0f * dot(bRec.wo, wh));
		} else {
			/* Zero probability if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 1)
				|| !(bRec.typeMask & EGlossyTransmission))
				return 0.0f;

			/* Calculate the transmission half-vector */
			Float eta = Frame::cosTheta(bRec.wi) > 0
				? m_eta : m_invEta;

			wh = normalize(bRec.wi + bRec.wo*eta);

			/* Jacobian of the half-direction mapping */
			Float sqrtDenom = dot(bRec.wi, wh) + eta * dot(bRec.wo, wh);
			dwh_dwo = (eta*eta * dot(bRec.wo, wh)) / (sqrtDenom*sqrtDenom);
		}

		/* Ensure that the half-vector points into the
		   same hemisphere as the macrosurface normal */
		wh *= math::signum(Frame::cosTheta(wh));

		RayInfo ray;
		Float s = math::signum(Frame::cosTheta(bRec.wi));
		vec3 wi(s*bRec.wi.x, s*bRec.wi.y, s*bRec.wi.z);
		const float alpha_x = std::max(m_alphaU->eval(bRec.its).average(), (Float) 1e-4f) ;
		const float alpha_y = std::max(m_alphaV->eval(bRec.its).average(), (Float) 1e-4f) ;
		ray.updateDirection(wi, alpha_x, alpha_y);
		
		Float prob = std::max(0.0f, dot(wh, wi)) * D_ggx(wh, alpha_x, alpha_y) / (1.0f + ray.Lambda) / Frame::cosTheta(wi);
		
		if (hasTransmission && hasReflection) {
			Float F = fresnelDielectricExt(dot(bRec.wi, wh), m_eta);
			prob *= reflect ? F : (1-F);
		}

		// single-scattering PDF + diffuse 
		// otherwise too many fireflies due to lack of multiple-scattering PDF
		// (MIS works even if the PDF is wrong and not normalized)
		return std::abs(prob * s) + Frame::cosTheta(bRec.wo); 
	}



	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &_sample) const {
		Float pdf;
		return this->sample(bRec, pdf, _sample);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample) const {
		Point2 sample(_sample);

		bool hasReflection = ((bRec.component == -1 || bRec.component == 0)
							  && (bRec.typeMask & EGlossyReflection)),
		     hasTransmission = ((bRec.component == -1 || bRec.component == 1)
							  && (bRec.typeMask & EGlossyTransmission)),
		     sampleReflection = hasReflection;

		if (!hasReflection && !hasTransmission)
			return Spectrum(0.0f);

		vec3 wi(bRec.wi.x, bRec.wi.y, bRec.wi.z);

		const float alpha_x = std::max(m_alphaU->eval(bRec.its).average(), (Float) 1e-4f);
		const float alpha_y = std::max(m_alphaV->eval(bRec.its).average(), (Float) 1e-4f);

		if ( Frame::cosTheta(bRec.wi) > 0 ) // outside
		{
			float weight;
			vec3 wo = sample_dielectric(wi, alpha_x, alpha_y, m_eta, m_scatteringOrderMax, weight);
							
			bRec.wo.x = wo.x;
			bRec.wo.y = wo.y;
			bRec.wo.z = wo.z;
			pdf = this->pdf(bRec, ESolidAngle);

			if( Frame::cosTheta(bRec.wo) > 0 ) // reflection
			{
				bRec.eta = 1.0f;
				bRec.sampledComponent = 0;
				bRec.sampledType = EGlossyReflection;
				const Spectrum R = m_specularReflectance->eval(bRec.its);
				return R * weight;
			}
			else // refraction
			{
				bRec.eta = m_eta;
				bRec.sampledComponent = 1;
				bRec.sampledType = EGlossyTransmission;
				Float factor = (bRec.mode == ERadiance) ? m_invEta : 1.0f; 
				const Spectrum T = m_specularTransmittance->eval(bRec.its);
				return T * factor * factor * weight;
			}
		}
		else // inside
		{
			float weight;
			vec3 wo = -sample_dielectric(-wi, alpha_x, alpha_y, m_invEta, m_scatteringOrderMax, weight);

			bRec.wo.x = wo.x;
			bRec.wo.y = wo.y;
			bRec.wo.z = wo.z;
			pdf = this->pdf(bRec, ESolidAngle);

			if( Frame::cosTheta(bRec.wo) > 0 ) // refraction
			{
				bRec.eta = m_invEta;
				bRec.sampledComponent = 1;
				bRec.sampledType = EGlossyTransmission; 
				Float factor = (bRec.mode == ERadiance) ? m_eta : 1.0f; 
				const Spectrum T = m_specularTransmittance->eval(bRec.its);
				return T * factor * factor * weight;
			}
			else // reflection
			{
				bRec.eta = 1.0f;
				bRec.sampledComponent = 0;
				bRec.sampledType = EGlossyReflection;
				const Spectrum R = m_specularReflectance->eval(bRec.its);
				return R * weight;
			}
		}
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
			else if (name == "specularTransmittance")
				m_specularTransmittance = static_cast<Texture *>(child);
			else
				BSDF::addChild(name, child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	Float getEta() const {
		return m_eta;
	}

	Float getRoughness(const Intersection &its, int component) const {
		return 0.5f * (m_alphaU->eval(its).average()
			+ m_alphaV->eval(its).average());
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "RoughDielectric[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  eta = " << m_eta << "," << endl
			<< "  alphaU = " << indent(m_alphaU->toString()) << "," << endl
			<< "  alphaV = " << indent(m_alphaV->toString()) << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  specularTransmittance = " << indent(m_specularTransmittance->toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_specularTransmittance;
	ref<Texture> m_specularReflectance;
	ref<Texture> m_alphaU, m_alphaV;
	Float m_eta, m_invEta;
	int m_scatteringOrderMax;
};

/* Fake glass shader -- it is really hopeless to visualize
   this material in the VPL renderer, so let's try to do at least
   something that suggests the presence of a transparent boundary */
class RoughDielectricShader : public Shader {
public:
	RoughDielectricShader(Renderer *renderer, Float eta) :
		Shader(renderer, EBSDFShader) {
		m_flags = ETransparent;
	}

	Float getAlpha() const {
		return 0.3f;
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return vec3(inv_pi * cosTheta(wo));" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << evalName << "(uv, wi, wo);" << endl
			<< "}" << endl;
	}


	MTS_DECLARE_CLASS()
};

Shader *RoughDielectric::createShader(Renderer *renderer) const {
	return new RoughDielectricShader(renderer, m_eta);
}

MTS_IMPLEMENT_CLASS(RoughDielectricShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(RoughDielectric, false, BSDF)
MTS_EXPORT_PLUGIN(RoughDielectric, "Rough dielectric BSDF");
MTS_NAMESPACE_END

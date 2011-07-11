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

#if !defined(__MICROFACET_H)
#define __MICROFACET_H

#include <mitsuba/core/quad.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/core/spline.h>
#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>

MTS_NAMESPACE_BEGIN

/**
 * Implements the microfacet distributions discussed in
 * "Microfacet Models for Refraction through Rough Surfaces"
 * by Bruce Walter, Stephen R. Marschner, Hongsong Li, and Kenneth E. Torrance
 */
class MicrofacetDistribution {
public:
	/// Supported distribution types
	enum EType {
		/// Beckmann distribution derived from Gaussian random surfaces
		EBeckmann         = 0,
		/// Long-tailed distribution proposed by Walter et al.
		EGGX              = 1,
		/// Classical Phong distribution
		EPhong            = 2,
		/// Anisotropic distribution by Ashikhmin and Shirley
		EAshikhminShirley = 3
	};

	/// Create a microfacet distribution of the specified type
	MicrofacetDistribution(EType type = EBeckmann)
		: m_type(type) { }

	/**
	 * \brief Create a microfacet distribution of the specified name
	 * (ggx/phong/beckmann/as)
	 */
	MicrofacetDistribution(const std::string &name) : m_type(EBeckmann) {
		std::string distr = boost::to_lower_copy(name);

		if (distr == "beckmann")
			m_type = EBeckmann;
		else if (distr == "phong")
			m_type = EPhong;
		else if (distr == "ggx")
			m_type = EGGX;
		else if (distr == "as")
			m_type = EAshikhminShirley;
		else 
			SLog(EError, "Specified an invalid distribution \"%s\", must be "
				"\"beckmann\", \"phong\", \"ggx\", or \"as\"!", distr.c_str());
	}

	/// Return the distribution type
	inline EType getType() const { return m_type; }

	/// Is this an anisotropic microfacet distribution?
	bool isAnisotropic() const {
		return m_type == EAshikhminShirley;
	}

	/**
	 * \brief Convert the roughness values so that they behave similarly to the
	 * Beckmann distribution.
	 *
	 * Also clamps to the minimal roughness 1e-4 to avoid numerical issues 
	 * (For lower roughness values, please switch to the smooth BSDF variants)
	 */
	Float transformRoughness(Float value) const {
		if (m_type == EPhong || m_type == EAshikhminShirley)
			value = 2 / (value * value) - 2;
		return std::max(value, (Float) 1e-4f);
	}
	
	/**
	 * \brief Implements the microfacet distribution function D
	 *
	 * \param m The microsurface normal
	 * \param alpha  The surface roughness
	 */
	inline Float eval(const Vector &m, Float alpha) const {
		return eval(m, alpha, alpha);
	}

	/**
	 * \brief Implements the microfacet distribution function D
	 *
	 * \param m The microsurface normal
	 * \param alphaU  The surface roughness in the tangent direction
	 * \param alphaV  The surface roughness in the bitangent direction
	 */
	Float eval(const Vector &m, Float alphaU, Float alphaV) const {
		if (Frame::cosTheta(m) <= 0)
			return 0.0f;
	
		Float result;
		switch (m_type) {
			case EBeckmann: {
					/* Beckmann distribution function for Gaussian random surfaces */
					const Float ex = Frame::tanTheta(m) / alphaU;
					result = std::exp(-(ex*ex)) / (M_PI * alphaU*alphaU * 
							std::pow(Frame::cosTheta(m), (Float) 4.0f));
				}
				break;
	
			case EGGX: {
					/* Empirical GGX distribution function for rough surfaces */
					const Float tanTheta = Frame::tanTheta(m),
						        cosTheta = Frame::cosTheta(m);
	
					const Float root = alphaU / (cosTheta*cosTheta * 
								(alphaU*alphaU + tanTheta*tanTheta));
	
					result = INV_PI * (root * root);
				}
				break;

			case EPhong: {
					/* Phong distribution function */
					result = (alphaU + 2) * INV_TWOPI 
							* std::pow(Frame::cosTheta(m), alphaU);
				}
				break;

			case EAshikhminShirley: {
					const Float cosTheta = Frame::cosTheta(m);
					const Float ds = 1 - cosTheta * cosTheta;
					if (ds < 0)
						return 0.0f;
					const Float exponent = (alphaU * m.x * m.x 
							+ alphaV * m.y * m.y) / ds;
					result = std::sqrt((alphaU + 2) * (alphaV + 2))
						* INV_TWOPI * std::pow(cosTheta, exponent);
				}
				break;

			default:
				SLog(EError, "Invalid distribution function!");
				return 0.0f;
		}

		/* Prevent potential numerical issues in other stages of the model */
		if (result < 1e-20f)
			result = 0;
	
		return result;
	}

	/**
	 * \brief Returns the density function associated with
	 * the \ref sample() function.
	 * \param m The microsurface normal
	 * \param alpha  The surface roughness 
	 */
	inline Float pdf(const Vector &m, Float alpha) const {
		return pdf(m, alpha, alpha);
	}

	/**
	 * \brief Returns the density function associated with
	 * the \ref sample() function.
	 * \param m The microsurface normal
	 * \param alphaU  The surface roughness in the tangent direction
	 * \param alphaV  The surface roughness in the bitangent direction
	 */
	Float pdf(const Vector &m, Float alphaU, Float alphaV) const {
		/* Usually, this is just D(m) * cos(theta_M) */
		if (m_type != EAshikhminShirley)
			return eval(m, alphaU, alphaV) * Frame::cosTheta(m);

		/* For the Ashikhmin-Shirley model, the sampling density
		   does not include the cos(theta_M) factor, and the
		   normalization is slightly different than in eval(). */
		const Float cosTheta = Frame::cosTheta(m);
		const Float ds = 1 - cosTheta * cosTheta;
		if (ds < 0)
			return 0.0f;
		const Float exponent = (alphaU * m.x * m.x 
				+ alphaV * m.y * m.y) / ds;
		Float result = std::sqrt((alphaU + 1) * (alphaV + 1))
			* INV_TWOPI * std::pow(cosTheta, exponent);

		/* Prevent potential numerical issues in other stages of the model */
		if (result < 1e-20f)
			result = 0;

		return result;
	}

	/// Helper routine: sample the first quadrant of the A&S distribution
	void sampleFirstQuadrant(Float alphaU, Float alphaV, Float u1, Float u2,
			Float &phi, Float &cosTheta) const {
		if (alphaU == alphaV)
			phi = M_PI * u1 * 0.5f;
		else
			phi = std::atan(
				std::sqrt((alphaU + 1.0f) / (alphaV + 1.0f)) *
				std::tan(M_PI * u1 * 0.5f));
		const Float cosPhi = std::cos(phi), sinPhi = std::sin(phi);
		cosTheta = std::pow(u2, 1.0f / 
			(alphaU * cosPhi * cosPhi + alphaV * sinPhi * sinPhi + 1.0f));
	}

	/**
	 * \brief Draw a sample from the microsurface normal distribution
	 *
	 * \param sample  A uniformly distributed 2D sample
	 * \param alpha  The surface roughness 
	 */
	inline Normal sample(const Point2 &sample, Float alpha) const {
		return MicrofacetDistribution::sample(sample, alpha, alpha);
	}

	/**
	 * \brief Draw a sample from the microsurface normal distribution
	 *
	 * \param sample  A uniformly distributed 2D sample
	 * \param alphaU  The surface roughness in the tangent direction
	 * \param alphaV  The surface roughness in the bitangent direction
	 */
	Normal sample(const Point2 &sample, Float alphaU, Float alphaV) const {
		/* The azimuthal component is always selected 
		   uniformly regardless of the distribution */
		Float phiM = (2.0f * M_PI) * sample.y,
			  thetaM = 0.0f;
	
		switch (m_type) {
			case EBeckmann: 
				thetaM = std::atan(std::sqrt(-alphaU*alphaU *
						 std::log(1.0f - sample.x)));
				break;
	
			case EGGX: 
				thetaM = std::atan(alphaU * std::sqrt(sample.x) /
						 std::sqrt(1.0f - sample.x));
				break;

			case EPhong:
				thetaM = std::acos(std::pow(sample.x, (Float) 1 / 
						 (alphaU + 2)));
				break;
	
			case EAshikhminShirley: {
					/* Sampling method based on code from PBRT */
					Float phi, cosTheta;
					if (sample.x < 0.25f) {
						sampleFirstQuadrant(alphaU, alphaV,
							4 * sample.x, sample.y, phi, cosTheta);
					} else if (sample.x < 0.5f) {
						sampleFirstQuadrant(alphaU, alphaV,
							4 * (0.5f - sample.x), sample.y, phi, cosTheta);
						phi = M_PI - phi;
					} else if (sample.x < 0.75f) {
						sampleFirstQuadrant(alphaU, alphaV,
							4 * (sample.x - 0.5f), sample.y, phi, cosTheta);
						phi += M_PI;
					} else {
						sampleFirstQuadrant(alphaU, alphaV,
							4 * (1 - sample.x), sample.y, phi, cosTheta);
						phi = 2 * M_PI - phi;
					}
					const Float sinTheta = std::sqrt(
						std::max((Float) 0, 1 - cosTheta*cosTheta));
					return Vector(
						sinTheta * std::cos(phi),
						sinTheta * std::sin(phi),
						cosTheta
					);
				}
				break;
			default: 
				SLog(EError, "Invalid distribution function!");
		}
	
		return Normal(sphericalDirection(thetaM, phiM));
	}

	/**
	 * \brief Draw a sample from an isotropic microsurface normal
	 * distribution and return the magnitude of its 'z' component.
	 *
	 * \param sample  A uniformly distributed number on [0,1]
	 * \param alphaU  The surface roughness 
	 */
	Float sampleIsotropic(Float sample, Float alpha) const {
		switch (m_type) {
			case EBeckmann: 
				return 1.0f / std::sqrt(1 + 
					std::abs(-alpha*alpha * std::log(1.0f - sample)));
	
			case EGGX: 
				return 1.0f / std::sqrt(1 + 
					alpha * alpha * sample / (1.0f - sample));

			case EPhong:
				return std::pow(sample, (Float) 1 / (alpha + 2));

			default: 
				SLog(EError, "Invalid distribution function!");
				return 0.0f;
		}
	}
	
	/**
	 * \brief Smith's shadow-masking function G1 for each
	 * of the supported microfacet distributions
	 *
	 * \param v An arbitrary direction
	 * \param m The microsurface normal
	 * \param alpha The surface roughness
	 */
	Float smithG1(const Vector &v, const Vector &m, Float alpha) const {
		const Float tanTheta = std::abs(Frame::tanTheta(v)); 

		/* perpendicular incidence -- no shadowing/masking */
		if (tanTheta == 0.0f)
			return 1.0f;
	
		/* Can't see the back side from the front and vice versa */
		if (dot(v, m) * Frame::cosTheta(v) <= 0)
			return 0.0f;
	
		switch (m_type) {
			case EPhong:
			case EBeckmann: {
					Float a;
					/* Approximation recommended by Bruce Walter: Use
					   the Beckmann shadowing-masking function with
					   specially chosen roughness value */
					if (m_type != EBeckmann)
						a = std::sqrt(0.5f * alpha + 1) / tanTheta;
					else
						a = 1.0f / (alpha * tanTheta);

					if (a >= 1.6f)
						return 1.0f;
	
					/* Use a fast and accurate (<0.35% rel. error) rational
					   approximation to the shadowing-masking function */
					const Float aSqr = a * a;
					return (3.535f * a + 2.181f * aSqr) 
						 / (1.0f + 2.276f * a + 2.577f * aSqr);
				}
				break;
	
			case EGGX: {
					const Float root = alpha * tanTheta;
					return 2.0f / (1.0f + std::sqrt(1.0f + root*root));
				}
				break;
	
			default:
				SLog(EError, "Invalid distribution function!");
				return 0.0f;
		}
	}

	/**
	 * \brief Shadow-masking function for each of the supported 
	 * microfacet distributions
	 *
	 * \param wi The incident direction
	 * \param wo The exitant direction
	 * \param m The microsurface normal
	 * \param alpha The surface roughness
	 */
	inline Float G(const Vector &wi, const Vector &wo, const Vector &m, Float alpha) const {
		return G(wi, wo, m, alpha, alpha);
	}

	/**
	 * \brief Shadow-masking function for each of the supported 
	 * microfacet distributions
	 *
	 * \param wi The incident direction
	 * \param wo The exitant direction
	 * \param m The microsurface normal
	 * \param alphaU  The surface roughness in the tangent direction
	 * \param alphaV  The surface roughness in the bitangent direction
	 */
	Float G(const Vector &wi, const Vector &wo, const Vector &m, Float alphaU, Float alphaV) const {
		if (m_type != EAshikhminShirley) {
			return smithG1(wi, m, alphaU)
				 * smithG1(wo, m, alphaU);
		} else {
			/* Can't see the back side from the front and vice versa */
			if (dot(wi, m) * Frame::cosTheta(wi) <= 0 ||
				dot(wo, m) * Frame::cosTheta(wo) <= 0)
				return 0.0f;

			/* Infinite groove shadowing/masking */
			const Float nDotM  = Frame::cosTheta(m),
						nDotWo = Frame::cosTheta(wo),
						nDotWi = Frame::cosTheta(wi),
						woDotM = dot(wo, m),
						wiDotM = dot(wi, m);

			return std::min((Float) 1, 
				std::min(std::abs(2 * nDotM * nDotWo / woDotM),
						 std::abs(2 * nDotM * nDotWi / wiDotM)));
		}
	}

	/**
	 * \brief Compute a spline representation for the overall Fresnel
	 * transmittance through a rough interface
	 *
	 * This function essentially computes the integral of 
	 *      1 - \int_{S^2} f(w_i, w_o) *  dw_o
	 * for incident directions 'wi' with a range of different inclinations
	 * (where f denotes a Cook-Torrance style reflectance model). It returns 
	 * a cubic spline interpolation parameterized by the cosine of the angle
	 * between 'wi' and the (macro-) surface normal.
	 *
	 * \remark This only works for isotropic microfacet distributions
	 */
	CubicSpline *computeRoughTransmittance(Float extIOR, Float intIOR, Float alpha, size_t resolution) const {
		if (isAnisotropic())
			SLog(EError, "MicrofacetDistribution::computeRoughTransmission(): only "
				"supports isotropic distributions!");

		NDIntegrator integrator(1, 2, 5000, 0, 1e-5f);
		CubicSpline *spline = new CubicSpline(resolution);
		size_t nEvals, nEvalsTotal = 0;
		ref<Timer> timer = new Timer();

		Float stepSize = (1.0f-2*Epsilon)/(resolution-1);
		for (size_t i=0; i<resolution; ++i) {
			Float z = stepSize * i + Epsilon;
			Vector wi(std::sqrt(std::max((Float) 0, 1-z*z)), 0, z);
			Float min[2] = {0, 0}, max[2] = {1, 1},
				  integral = 0, error = 0;

			integrator.integrateVectorized(
				boost::bind(&MicrofacetDistribution::integrand1, this,
					wi, extIOR, intIOR, alpha, _1, _2, _3),
				min, max, &integral, &error, &nEvals
			);

			spline->append(z, 1-integral);

			nEvalsTotal += nEvals;
		}
		SLog(EInfo, "Created a " SIZE_T_FMT "-node cubic spline approximation to the rough Frensel "
				"transmittance (integration took %i ms and " SIZE_T_FMT " function evaluations)",
				resolution, timer->getMilliseconds(), nEvalsTotal);
		spline->build();
		return spline;
	}

	/**
	 * \brief Compute a spline representation that gives the probability
	 * of choosing a reflection event when importance sampling wrt. the
	 * Fresnel coefficient between a sampled microsurface normal and the
	 * incident direction.
	 *
	 * This function is currently used by the plugin 'roughplastic'.
	 *
	 * Like \ref computeRoughTransmittance, the spline is parameterized by the
	 * cosine of the angle between the indident direction and the (macro-) 
	 * surface normal.
	 * 
	 * \remark This function only works for isotropic microfacet distributions
	 */
	CubicSpline *computeTransmissionProbability(Float extIOR, Float intIOR, 
			Float alpha, Float specularSamplingWeight, size_t resolution) const {
		if (isAnisotropic())
			SLog(EError, "MicrofacetDistribution::computeTransmissionProbability(): only "
				"supports isotropic distributions!");

		NDIntegrator integrator(1, 2, 5000, 0, 1e-5f);
		CubicSpline *spline = new CubicSpline(resolution);
		size_t nEvals, nEvalsTotal = 0;
		ref<Timer> timer = new Timer();

		Float stepSize = (1.0f-2*Epsilon)/(resolution-1);
		for (size_t i=0; i<resolution; ++i) {
			Float z = stepSize * i + Epsilon;
			Vector wi(std::sqrt(std::max((Float) 0, 1-z*z)), 0, z);
			Float min[2] = {0, 0}, max[2] = {1, 1},
				  integral = 0, error = 0;

			integrator.integrateVectorized(
				boost::bind(&MicrofacetDistribution::integrand2, this,
					wi, extIOR, intIOR, alpha, specularSamplingWeight, _1, _2, _3),
				min, max, &integral, &error, &nEvals
			);

			spline->append(z, integral);

			nEvalsTotal += nEvals;
		}
		SLog(EInfo, "Created a " SIZE_T_FMT "-node cubic spline approximation to the "
				"transmission probability (integration took %i ms and " SIZE_T_FMT 
				" function evaluations)", resolution, timer->getMilliseconds(), 
				nEvalsTotal);

		spline->build();
		return spline;
	}

	std::string toString() const {
		switch (m_type) {
			case EBeckmann: return "beckmann"; break;
			case EPhong: return "phong"; break;
			case EGGX: return "ggx"; break;
			case EAshikhminShirley: return "as"; break;
			default:
				SLog(EError, "Invalid distribution function");
				return "";
		}
	}
protected:
	/// Integrand helper function called by \ref computeRoughTransmission
	void integrand1(const Vector &wi, Float extIOR, Float intIOR, Float alpha,
			size_t nPts, const Float *in, Float *out) const {
		for (int i=0; i<(int) nPts; ++i) {
			Normal m = sample(Point2(in[2*i], in[2*i+1]), alpha);
			Vector wo = 2 * dot(wi, m) * Vector(m) - wi;
			if (Frame::cosTheta(wo) <= 0) {
				out[i] = 0;
				continue;
			}

			/* Calculate the specular reflection component */
			out[i] = std::abs(fresnel(dot(wi, m), extIOR, intIOR)
				* G(wi, wo, m, alpha) * dot(wi, m) /
				  (Frame::cosTheta(wi) * Frame::cosTheta(m)));
		}
	}

	/// Integrand helper function called by \ref computeTransmissionProbability
	void integrand2(const Vector &wi, Float extIOR, Float intIOR, Float alpha,
			Float specularSamplingWeight, size_t nPts, const Float *in, Float *out) const {
		for (int i=0; i<(int) nPts; ++i) {
			Normal m = sample(Point2(in[2*i], in[2*i+1]), alpha);
			Float probSpecular = fresnel(dot(wi, m), extIOR, intIOR);
			probSpecular = (probSpecular*specularSamplingWeight) /
				(probSpecular*specularSamplingWeight + 
				(1-probSpecular) * (1-specularSamplingWeight));
			out[i] = 1-probSpecular;
		}
	}

protected:
	EType m_type;
};

MTS_NAMESPACE_END

#endif /* __MICROFACET_H */

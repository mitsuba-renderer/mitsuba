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

#include <mitsuba/mitsuba.h>
#include <boost/algorithm/string.hpp>

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
		EBeckmann = 0,
		/// Classical Phong distribution
		EPhong    = 1,
		/// Long-tailed distribution proposed by Walter et al.
		EGGX      = 2
	};

	/// Create a microfacet distribution of the specified type
	MicrofacetDistribution(EType type = EBeckmann)
		: m_type(type) { }

	/**
	 * \brief Create a microfacet distribution of the specified name
	 * (ggx/phong/beckmann)
	 */
	MicrofacetDistribution(const std::string &name) {
		std::string distr = 
			boost::to_lower_copy(props.getString("distribution", "beckmann"));

		if (distr == "beckmann")
			m_type = EBeckmann;
		else if (distr == "phong")
			m_type = EPhong;
		else if (distr == "ggx")
			m_type = EGGX;
		else 
			SLog(EError, "Specified an invalid distribution \"%s\", must be "
				"\"beckmann\", \"phong\", or \"ggx\"!", distr.c_str());
	}

	/// Return the distribution type
	inline EType getType() const { return m_type; }

	/**
	 * \brief Convert the roughness values so that they behave similarly to the
	 * Beckmann distribution.
	 *
	 * Also clamps to the minimal roughness 1e-4 to avoid numerical issues 
	 * (For lower roughness values, please switch to the smooth BSDF variants)
	 */
	Float transformRoughness(Float value) const {
		if (m_type == EPhong)
			value = 2 / (value * value) - 2;
		return std::max(value, (Float) 1e-4f);
	}

	/**
	 * \brief Implements the microfacet distribution function D
	 *
	 * \param m The microsurface normal
	 * \param v An arbitrary direction
	 */
	Float eval(const Vector &m, Float alpha) const {
		if (Frame::cosTheta(m) <= 0)
			return 0.0f;
	
		Float result;
		switch (m_type) {
			case EBeckmann: {
					/* Beckmann distribution function for Gaussian random surfaces */
					const Float ex = Frame::tanTheta(m) / alpha;
					result = std::exp(-(ex*ex)) / (M_PI * alpha*alpha * 
							std::pow(Frame::cosTheta(m), (Float) 4.0f));
				}
				break;
	
			case EPhong: {
					/* Phong distribution function */
					result = (alpha + 2) * INV_TWOPI 
							* std::pow(Frame::cosTheta(m), alpha);
				}
				break;
	
			case EGGX: {
					/* Empirical GGX distribution function for rough surfaces */
					const Float tanTheta = Frame::tanTheta(m),
						        cosTheta = Frame::cosTheta(m);
	
					const Float root = alpha / (cosTheta*cosTheta * 
								(alpha*alpha + tanTheta*tanTheta));
	
					result = INV_PI * (root * root);
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
	 * \brief Sample microsurface normals according to 
	 * the selected distribution
	 *
	 * \param sample  A uniformly distributed 2D sample
	 * \param alpha   Surface roughness
	 */
	Normal sample(const Point2 &sample, Float alpha) const {
		/* The azimuthal component is always selected 
		   uniformly regardless of the distribution */
		Float phiM = (2.0f * M_PI) * sample.y,
			  thetaM = 0.0f;
	
		switch (m_type) {
			case EBeckmann: 
				thetaM = std::atan(std::sqrt(-alpha*alpha *
						 std::log(1.0f - sample.x)));
				break;
	
			case EPhong:
				thetaM = std::acos(std::pow(sample.x, (Float) 1 / 
						 (alpha + 2)));
				break;
	
			case EGGX: 
				thetaM = std::atan(alpha * std::sqrt(sample.x) /
						 std::sqrt(1.0f - sample.x));
				break;
	
			default: 
				SLog(EError, "Invalid distribution function!");
		}
	
		return Normal(sphericalDirection(thetaM, phiM));
	}
	
	/**
	 * \brief Smith's shadow-masking function G1 for each
	 * of the supported microfacet distributions
	 *
	 * \param m The microsurface normal
	 * \param v An arbitrary direction
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
				/* Approximation recommended by Bruce Walter: Use
				   the Beckmann shadowing-masking function with
				   specially chosen roughness value */
				alpha = std::sqrt(0.5f * alpha + 1) / tanTheta;
	
			case EBeckmann: {
					/* Use a fast and accurate (<0.35% rel. error) rational
					   approximation to the shadowing-masking function */
					const Float a = 1.0f / (alpha * tanTheta);
					const Float aSqr = a * a;
	
					if (a >= 1.6f)
						return 1.0f;
	
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
	
	std::string toString() const {
		switch (m_distribution) {
			case EBeckmann: return "beckmann"; break;
			case EPhong: return "phong"; break;
			case EGGX: return "ggx"; break;
			default:
				SLog(EError, "Invalid distribution function");
				return "";
		}
	}
private:
	EType type;
};

MTS_NAMESPACE_END

#endif /* __MICROFACET_H */

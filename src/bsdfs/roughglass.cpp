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
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

/*! \newpage\plugin{roughglass}{Rough dielectric/glass material}
 * \parameters{
 *     \parameter{alpha}{\Float}{Roughness value (Default: 0.1)}
 *     \parameter{intIOR}{\Float}{Interior index of refraction (Default: 1.5046)}
 *     \parameter{extIOR}{\Float}{Exterior index of refraction (Default: 1)}
 *     \parameter{specular\showbreak Reflectance}{\Spectrum}{Modulation)}
 *     \parameter{specular\showbreak Transmittance}{\Spectrum}{Modulation)}
 *     \lastparameter{distribution}{\String}{
 *       Specifies the microfacet distribution
 *       \begin{enumerate}[(i)]
 *           \item \code{beckmann}: Beckmann distribution derived from
 *               Gaussian random surfaces.
 *           \item \code{phong}: Phong distribution
 *           \item \code{ggx}: New distribution proposed by
 *               Walter et al., meant to better handle the long
 *               tails observed in measurements of rough glass.
 *       \end{enumerate}
 *     }
 * }
 *
 * Rough glass BSDF model based on
 * "Microfacet Models for Refraction through Rough Surfaces"
 * by Bruce Walter, Stephen R. Marschner, Hongsong Li
 * and Kenneth E. Torrance
 * The default settings are set to a borosilicate glass BK7/air interface.
 */


class RoughGlass : public BSDF {
public:
	//// Microfacet distribution types supported by the model
	enum EDistribution  {
		EBeckmann = 0x0000,
		EPhong    = 0x0001,
		EGGX      = 0x0002
	};

	RoughGlass(const Properties &props) 
		: BSDF(props) {
		m_specularReflectance = props.getSpectrum("specularReflectance", 
			Spectrum(1.0f));
		m_specularTransmittance = props.getSpectrum("specularTransmittance", 
			Spectrum(1.0f));

		if (props.hasProperty("alphaB")) {
			Log(EWarn, "Deprecation warning: the 'alphaB' parameter "
				"has been renamed to 'alpha'");

			m_alpha = props.getFloat("alphaB");
		} else {
			m_alpha = props.getFloat("alpha", .1f);
		}

		m_intIOR = props.getFloat("intIOR", 1.5046f);
		m_extIOR = props.getFloat("extIOR", 1.0f);

		if (m_intIOR == m_extIOR)
			Log(EError, "Indices of refraction must differ!");

		std::string distr = 
			boost::to_lower_copy(props.getString("distribution", "beckmann"));

		if (distr == "beckmann")
			m_distribution = EBeckmann;
		else if (distr == "phong")
			m_distribution = EPhong;
		else if (distr == "ggx")
			m_distribution = EGGX;
		else 
			Log(EError, "Specified an invalid distribution \"%s\", must be "
				"\"beckmann\", \"phong\", or \"ggx\"!", distr.c_str());

		if (m_distribution == EPhong) {
			/* Transform the Phong exponent to make it behave
			   similarly to the Beckmann microfacet distribution */
			m_alpha = 2 / (m_alpha * m_alpha) - 2;
			Assert(m_alpha > 0);
		}

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EGlossyReflection | EFrontSide | EBackSide;
		m_type[1] = EGlossyTransmission | EFrontSide | EBackSide;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = false;
	}

	RoughGlass(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager) {
		m_distribution = (EDistribution) stream->readInt();
		m_specularReflectance = Spectrum(stream);
		m_specularTransmittance = Spectrum(stream);
		m_alpha = stream->readFloat();
		m_intIOR = stream->readFloat();
		m_extIOR = stream->readFloat();

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EGlossyReflection | EFrontSide | EBackSide;
		m_type[1] = EGlossyTransmission | EFrontSide | EBackSide;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = false;
	}

	virtual ~RoughGlass() {
		delete[] m_type;
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return Spectrum(0.0f);
	}

	inline Vector reflect(const Vector &wi, const Normal &m) const {
		return 2 * dot(wi, m) * Vector(m) - wi;
	}

	inline bool refract(const Vector &wi, Vector &wo, const Normal &m) const {
		/* Determine the appropriate indices of refraction */
		Float etaI = m_extIOR, etaT = m_intIOR;
		if (Frame::cosTheta(wi) < 0)
			std::swap(etaI, etaT);

		Float eta = etaI / etaT, c = dot(wi, m);

		/* Using Snell's law, calculate the squared cosine of the
		   angle between the normal and the transmitted ray */
		Float cosThetaTSqr = 1 + eta * eta * (c*c-1);

		if (cosThetaTSqr < 0) 
			return false; // Total internal reflection

		/* Compute the transmitted direction */
		wo = m * (eta*c - signum(wi.z)
			   * std::sqrt(cosThetaTSqr)) - wi * eta;

		return true;
	}

	inline Float signum(Float value) const {
		return (value < 0) ? -1.0f : 1.0f;
	}


	/**
	 * \brief Implements the microfacet distribution function D
	 *
	 * \param m The microsurface normal
	 * \param v An arbitrary direction
	 */
	Float evalD(const Vector &m, Float alpha) const {
		if (Frame::cosTheta(m) < 0)
			return 0.0f;
	
		Float result;
		switch (m_distribution) {
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
					const Float tanTheta = Frame::tanTheta(m);
					const Float cosTheta = Frame::cosTheta(m);

					const Float root = alpha / (cosTheta*cosTheta * 
								(alpha*alpha + tanTheta*tanTheta));

					result = INV_PI * (root*root);
				}
				break;

			default :
				Log(EError, "Invalid distribution function!");
				return 0.0f;
		}

		/* Prevent potential numerical issues in other stages of the model */
		if (result < 1e-40)
			result = 0;

		return result;
	}

	/**
	 * \brief Smith's shadow-masking function G1 for each
	 * of the supported microfacet distributions
	 *
	 * \param m The microsurface normal
	 * \param v An arbitrary direction
	 */
	Float smithG1(const Vector &v, const Vector &m) const {
		const Float tanTheta = std::abs(Frame::tanTheta(v)); 
		Float alpha = m_alpha;

		/* perpendicular incidence -- no shadowing/masking */
		if (tanTheta == 0.0f)
			return 1.0f;
			
		/* Can't see the back side from the front and vice versa */
		if (dot(v, m) * Frame::cosTheta(v) <= 0)
			return 0.0f;

		switch (m_distribution) {
			case EPhong:
				/* Approximation recommended by Bruce Walter: Use
				   the Beckmann shadowing-masking function with
				   a modified alpha value */
				alpha = std::sqrt(0.5f * alpha + 1) / tanTheta;

			case EBeckmann: {
					const Float a = 1.0f / (alpha * tanTheta);
					const Float aSqr = a * a;

					if (a >= 1.6f)
						return 1.0f;

					return (3.535f * a + 2.181f * aSqr) 
						 / (1.0f + 2.276f * a + 2.577f * aSqr);
				}
				break;

			case EGGX: {
					const Float root = m_alpha * tanTheta;
					return 2.0f / (1.0f + std::sqrt(1.0f + root*root));
				}
				break;

			default:
				Log(EError, "Invalid distribution function!");
				return 0.0f;
		}
	}

	/**
	 * \brief Sample microsurface normals according to 
	 * the selected distribution
	 *
	 * \param sample  A uniformly distributed 2D sample
	 * \param alpha   Surface roughness
	 */
	Normal sampleD(const Point2 &sample, Float alpha) const {
		/* The azimuthal component is always selected 
		   uniformly regardless of the distribution */
		Float phiM = (2.0f * M_PI) * sample.y,
			  thetaM = 0.0f;

		switch (m_distribution) {
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
				Log(EError, "Invalid distribution function!");
		}

		return Normal(sphericalDirection(thetaM, phiM));
	}

	inline Spectrum fReflection(const BSDFQueryRecord &bRec) const {
		Float intIOR = m_intIOR, extIOR = m_extIOR;

		if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		if (Frame::cosTheta(bRec.wi) < 0)
			std::swap(intIOR, extIOR);

		/* Calculate the reflection half-vector (and possibly flip it
		   so that it lies inside the hemisphere around the normal) */
		Vector Hr = normalize(bRec.wo+bRec.wi) 
			* signum(Frame::cosTheta(bRec.wo));

		/* Fresnel factor */
		Float F = fresnel(dot(bRec.wi, Hr), m_extIOR, m_intIOR);

		/* Microsurface normal distribution */
		Float D = evalD(Hr, m_alpha);

		/* Smith's shadow-masking function */
		Float G = smithG1(bRec.wi, Hr) * smithG1(bRec.wo, Hr);

		/* Calculate the total amount of reflection */
		Float value = F * D * G / 
			(4.0f * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo));
		
		return m_specularReflectance * value; 
	}

	Spectrum fTransmission(const BSDFQueryRecord &bRec) const {
		if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
			return Spectrum(0.0f);

		Float etaI = m_extIOR, etaT = m_intIOR;
		if (Frame::cosTheta(bRec.wi) < 0)
			std::swap(etaI, etaT);

		/* Calculate the transmission half-vector */
		Vector Ht = -normalize(bRec.wi*etaI + bRec.wo*etaT);
		if (m_extIOR > m_intIOR)
			Ht *= -1;

		/* Fresnel factor */
		Float F = 1.0f - fresnel(dot(bRec.wi, Ht), m_extIOR, m_intIOR);

		/* Microsurface normal distribution */
		Float D = evalD(Ht, m_alpha);

		if (D == 0)
			return Spectrum(0.0f);

		/* Smith's shadow-masking function */
		Float G = smithG1(bRec.wi,  Ht) * smithG1(bRec.wo,  Ht);

		/* Calculate the total amount of transmission */
		Float value = F * D * G * std::abs((dot(bRec.wi, Ht)*dot(bRec.wo, Ht)) / 
			(Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo)));

		Float sqrtDenom = etaI * dot(bRec.wi, Ht) + etaT * dot(bRec.wo, Ht);
		value *= (etaT * etaT) / (sqrtDenom*sqrtDenom);

		if (bRec.quantity == ERadiance)
			value *= (etaI*etaI) / (etaT*etaT);

		return m_specularTransmittance * value;
	}


	inline Spectrum f(const BSDFQueryRecord &bRec) const {
		Spectrum result(0.0f);

		if (!bRec.typeMask & m_combinedType)
			return result;

		bool hasReflection = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasTransmission = (bRec.typeMask & EGlossyTransmission)
				&& (bRec.component == -1 || bRec.component == 1);
		
		if (hasReflection)
			result += fReflection(bRec);
		if (hasTransmission)
			result += fTransmission(bRec);

		return result;
	}
	
	inline Float pdfReflection(const BSDFQueryRecord &bRec, Float alpha) const {
		if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) < 0)
			return 0.0f;

		Vector Hr = normalize(bRec.wo + bRec.wi) 
			* signum(Frame::cosTheta(bRec.wi));

		/* Jacobian of the half-direction transform */
		Float dwhr_dwo = 1.0f / (4.0f * absDot(bRec.wo, Hr));

		return evalD(Hr, alpha) * std::abs(Frame::cosTheta(Hr)) * dwhr_dwo;
	}

	inline Float pdfTransmission(const BSDFQueryRecord &bRec, Float alpha) const {
		Float etaI = m_extIOR, etaT = m_intIOR;

		if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
			return 0.0f;

		if (Frame::cosTheta(bRec.wi) < 0)
			std::swap(etaI, etaT);

		Vector Ht = -normalize(bRec.wi*etaI + bRec.wo*etaT);
		if (m_extIOR > m_intIOR)
			Ht *= -1;

		/* Jacobian of the half-direction transform. */
		Float sqrtDenom = etaI * dot(bRec.wi, Ht) + etaT * dot(bRec.wo, Ht);
		Float dwht_dwo = (etaT*etaT * absDot(bRec.wo, Ht)) / (sqrtDenom*sqrtDenom);

		return evalD(Ht, alpha) * std::abs(Frame::cosTheta(Ht)) * dwht_dwo;
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		bool hasReflection = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasTransmission = (bRec.typeMask & EGlossyTransmission)
				&& (bRec.component == -1 || bRec.component == 1);
		
		/* Suggestion by Bruce Walter: sample using a slightly different 
		   value of alpha. This in practice limits the weights to 
		   values <= 4. See also \ref sample() */
		Float alpha = m_alpha * (1.2f - 0.2f * std::sqrt(
				std::abs(Frame::cosTheta(bRec.wi))));

		if (hasReflection && hasTransmission) {
			/* PDF for importance sampling according to approximate 
			   Fresnel coefficients (approximate, because we don't know 
			   the microfacet normal at this point) */
			Float fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);
			fr = std::min(std::max(fr, (Float) 0.05f), (Float) 0.95f);
			return fr * pdfReflection(bRec, alpha) +
				   (1-fr) * pdfTransmission(bRec, alpha);
		} else if (hasReflection) {
			return pdfReflection(bRec, alpha);
		} else if (hasTransmission) {
			return pdfTransmission(bRec, alpha);
		}

		return 0.0f;
	}

	inline Spectrum sampleReflection(BSDFQueryRecord &bRec, Float alpha, const Point2 &sample) const {
		/* Sample M, the microsurface normal */
		Normal m = sampleD(sample, alpha);

		/* Perfect specular reflection along the microsurface normal */
		bRec.wo = reflect(bRec.wi, m);

		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;

		/* Side check */
		if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		Float pdfValue = pdf(bRec);

		/* Guard against numerical imprecisions */
		if (pdfValue == 0)
			return Spectrum(0.0f);
		else
			return f(bRec) / pdfValue;
	}

	inline Spectrum sampleTransmission(BSDFQueryRecord &bRec, Float alpha, const Point2 &sample) const {
		/* Sample the microfacet normal */
		Vector m = sampleD(sample, alpha);

		/* Refract based on 'm' */
		if (!refract(bRec.wi, bRec.wo, m))
			return Spectrum(0.0f);

		bRec.sampledComponent = 1;
		bRec.sampledType = EGlossyTransmission;

		/* Side check */
		if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
			return Spectrum(0.0f);
			
		Float pdfValue = pdf(bRec);
		
		/* Guard against numerical imprecisions */
		if (pdfValue == 0)
			return Spectrum(0.0f);
		else
			return f(bRec) / pdfValue;
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &_sample) const {
		bool hasReflection = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasTransmission = (bRec.typeMask & EGlossyTransmission)
				&& (bRec.component == -1 || bRec.component == 1);
		Point2 sample(_sample);

		/* Suggestion by Bruce Walter: sample using a slightly different 
		   value of alpha. This in practice limits the weights to 
		   values <= 4. The change is of course also accounted for 
		   in \ref pdf(), hence no error is introduced. */
		Float alpha = m_alpha * (1.2f - 0.2f * std::sqrt(
				std::abs(Frame::cosTheta(bRec.wi))));

		if (hasReflection && hasTransmission) {
			/* PDF for importance sampling according to approximate 
			   Fresnel coefficients (approximate, because we don't know 
			   the microfacet normal at this point) */
			Float fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);
			fr = std::min(std::max(fr, (Float) 0.05f), (Float) 0.95f);
			if (sample.x < fr) {
				sample.x /= fr;
				return sampleReflection(bRec, alpha, sample);
			} else {
				sample.x = (sample.x - fr) / (1 - fr);
				return sampleTransmission(bRec, alpha, sample);
			}
		} else if (hasReflection) {
			return sampleReflection(bRec, alpha, sample);
		} else if (hasTransmission) {
			return sampleTransmission(bRec, alpha, sample);
		}

		return Spectrum(0.0f);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeInt(m_distribution);
		m_specularReflectance.serialize(stream);
		m_specularTransmittance.serialize(stream);
		stream->writeFloat(m_alpha);
		stream->writeFloat(m_intIOR);
		stream->writeFloat(m_extIOR);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "RoughGlass[" << endl
			<< "  distribution = ";
		switch (m_distribution) {
			case EBeckmann: oss << "beckmann," << endl; break;
			case EGGX: oss << "ggx," << endl; break;
			case EPhong: oss << "phong," << endl; break;
			default:
				Log(EError, "Invalid distribution function");
		}
		oss << "  specularReflectance = " << m_specularReflectance.toString() << "," << endl
			<< "  specularTransmittance = " << m_specularTransmittance.toString() << "," << endl
			<< "  intIOR = " << m_intIOR << "," << endl
			<< "  extIOR = " << m_extIOR << "," << endl
			<< "  alpha = " << m_alpha << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	EDistribution m_distribution;
	Spectrum m_specularReflectance;
	Spectrum m_specularTransmittance;
	Float m_alpha, m_intIOR, m_extIOR;
};

MTS_IMPLEMENT_CLASS_S(RoughGlass, false, BSDF)
MTS_EXPORT_PLUGIN(RoughGlass, "Rough glass BSDF");
MTS_NAMESPACE_END

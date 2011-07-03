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
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/consttexture.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

/*! \plugin{roughdielectric}{Rough dielectric material}
 * \parameters{
 *     \parameter{distribution}{\String}{
 *          Specifies the type of microfacet normal distribution 
 *          used to model the surface roughness.
 *       \begin{enumerate}[(i)]
 *           \item \code{beckmann}: Physically-based distribution derived from
 *               Gaussian random surfaces. This is the default.
 *           \item \code{phong}: Classical $\cos^p\theta$ distribution.
 *              The Phong exponent $p$ is obtained using a transformation that
 *              produces roughness similar to a Beckmann distribution of the same 
 *              parameter. Note that due to the underlying microfacet theory, 
 *              the use of this distribution here leads to more realistic 
 *              behavior than the separately available \pluginref{phong} plugin.
 *           \item \code{ggx}: New distribution proposed by
 *               Walter et al. meant to better handle the long
 *               tails observed in transmission measurements through
 *               ground glass. Renderings with this distribution may
 *               converge slowly.
 *       \end{enumerate}
 *       Default: \code{beckmann}
 *     }
 *     \parameter{alpha}{\Float\Or\Texture}{Roughness value of the
 *         unresolved surface microgeometry. When the Beckmann
 *         distribution is used, this parameter specifies the 
 *         \emph{root mean square} (RMS) slope of the microfacets.
 *         \default{0.1}
 *     }
 *     \parameter{intIOR}{\Float}{Interior index of refraction \default{1.5046}}
 *     \parameter{extIOR}{\Float}{Exterior index of refraction \default{1.0}}
 *     \parameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor used to modulate the reflectance component\default{1.0}}
 *     \lastparameter{specular\showbreak Transmittance}{\Spectrum\Or\Texture}{Optional
 *         factor used to modulate the transmittance component\default{1.0}}
 * }
 *
 *
 * This plugin implements a realistic microfacet scattering model for rendering
 * rough interfaces between dielectric materials, such as a transition from air to 
 * ground glass. Microfacet theory describes rough surfaces as an arrangement of 
 * unresolved and ideally specular facets, whose normal directions are given by 
 * a specially chosen \emph{microfacet distribution}. By accounting for shadowing
 * and masking effects between these facets, it is possible to reproduce the 
 * off-specular reflections peaks observed in real-world measurements of such 
 * materials.
 * \renderings{
 *     \rendering{Rough glass (Beckmann, $\alpha$=0.1)}{bsdf_roughdielectric_beckmann_0_1.jpg}
 *     \rendering{Ground glass (GGX, $\alpha$=0.304, \lstref{roughdielectric-roughglass})}{bsdf_roughdielectric_ggx_0_304.jpg}
 * }
 *
 * This plugin is essentially the ``roughened'' equivalent of the plugin
 * \pluginref{dielectric}. As the roughness value is decreased, it increasingly
 * approximates that model. Its implementation is based on the paper 
 * ``Microfacet Models for Refraction through Rough Surfaces'' by Walter et
 * al. \cite{Walter07Microfacet}. The model supports several types of microfacet
 * distributions and has a texturable roughness parameter. 
 * The default settings are set 
 * to a borosilicate glass BK7/air interface with a light amount of rougness 
 * modeled using a Beckmann distribution.
 *
 * When using this plugin, it is crucial that the scene contains
 * meaningful and mutally compatible index of refraction change---see
 * \figref{glass-explanation} for an example. Also, please note that
 * the importance sampling implementation of this model is close, but 
 * not perfect a perfect match to the underlying scattering distribution,
 * particularly for high roughness values and when the \texttt{GGX} 
 * model is used. Hence, such renderings may converge slowly.
 *
 * \begin{xml}[caption=Ground glass, label=lst:roughdielectric-roughglass]
 * <bsdf type="roughdielectric">
 *     <string name="distribution" value="ggx"/>
 *     <float name="alpha" value="0.304"/>
 *     <float name="intIOR" value="1.5046"/>
 *     <float name="extIOR" value="1.0"/>
 * </bsdf>
 * \end{xml}
 *
 * \begin{xml}[caption=Textured rougness, label=lst:roughdielectric-textured]
 * <bsdf type="roughdielectric">
 *     <string name="distribution" value="beckmann"/>
 *     <float name="intIOR" value="1.5046"/>
 *     <float name="extIOR" value="1.0"/>
 *
 *     <texture type="bitmap" name="alpha">
 *         <string name="filename" value="roughness.exr"/>
 *     </texture>
 * </bsdf>
 * \end{xml}
 */
class RoughDielectric : public BSDF {
public:
	//// Microfacet distribution types supported by the model
	enum EDistribution  {
		/// Beckmann distribution derived from Gaussian random surfaces
		EBeckmann = 0x0000,
		/// Classical Phong distribution
		EPhong    = 0x0001,
		/// Long-tailed distribution proposed by Walter et al.
		EGGX      = 0x0002
	};

	RoughDielectric(const Properties &props) 
		: BSDF(props) {
		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));
		m_specularTransmittance = new ConstantSpectrumTexture(
			props.getSpectrum("specularTransmittance", Spectrum(1.0f)));

		Float alpha;
		if (props.hasProperty("alphaB")) {
			Log(EWarn, "Deprecation warning: the 'alphaB' parameter "
				"has been renamed to 'alpha'");

			alpha = props.getFloat("alphaB");
		} else {
			alpha = props.getFloat("alpha", 0.1f);
		}

		m_intIOR = props.getFloat("intIOR", 1.5046f);
		m_extIOR = props.getFloat("extIOR", 1.0f);

		if (m_intIOR < 0 || m_extIOR < 0 || m_intIOR == m_extIOR)
			Log(EError, "The interior and exterior indices of "
				"refraction must be positive and differ!");

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
			alpha = 2 / (alpha * alpha) - 2;
			AssertEx(alpha > 0, "Oops -- unable to map to a "
				"valid Phong exponent.");
		}

		m_alpha = new ConstantFloatTexture(alpha);

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EGlossyReflection | EFrontSide | EBackSide | ECanUseSampler;
		m_type[1] = EGlossyTransmission | EFrontSide | EBackSide | ECanUseSampler;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = false;
	}

	RoughDielectric(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager) {
		m_distribution = (EDistribution) stream->readInt();
		m_alpha = static_cast<Texture *>(manager->getInstance(stream));
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_specularTransmittance = static_cast<Texture *>(manager->getInstance(stream));
		m_intIOR = stream->readFloat();
		m_extIOR = stream->readFloat();

		m_componentCount = 2;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EGlossyReflection | EFrontSide | EBackSide | ECanUseSampler;
		m_type[1] = EGlossyTransmission | EFrontSide | EBackSide | ECanUseSampler;
		m_combinedType = m_type[0] | m_type[1];
		m_usesRayDifferentials = 
			m_alpha->usesRayDifferentials() ||
			m_specularReflectance->usesRayDifferentials() ||
			m_specularTransmittance->usesRayDifferentials();
	}

	virtual ~RoughDielectric() {
		delete[] m_type;
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return Spectrum(0.0f);
	}

	/// Helper function: reflect \c wi with respect to a given surface normal
	inline Vector reflect(const Vector &wi, const Normal &m) const {
		return 2 * dot(wi, m) * Vector(m) - wi;
	}

	/// Helper function: refract \c wi with respect to a given surface normal
	inline bool refract(const Vector &wi, Vector &wo, const Normal &m, Float etaI, Float etaT) const {
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
		if (Frame::cosTheta(m) <= 0)
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
					const Float tanTheta = Frame::tanTheta(m),
						        cosTheta = Frame::cosTheta(m);

					const Float root = alpha / (cosTheta*cosTheta * 
								(alpha*alpha + tanTheta*tanTheta));

					result = INV_PI * (root * root);
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

		switch (m_distribution) {
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
				Log(EError, "Invalid distribution function!");
				return 0.0f;
		}
	}

	/**
	 * \brief Evaluate the BSDF f(wi, wo) or its adjoint version f^{*}(wi, wo)
	 */
	Spectrum f(const BSDFQueryRecord &bRec) const {
		/* Determine the type of interaction */
		bool reflect = Frame::cosTheta(bRec.wi) 
			* Frame::cosTheta(bRec.wo) > 0;

		/* Determine the appropriate indices of refraction */
		Float etaI = m_extIOR, etaT = m_intIOR;
		if (Frame::cosTheta(bRec.wi) < 0)
			std::swap(etaI, etaT);

		Vector H;
		if (reflect) {
			/* Stop if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 0)
				|| !(bRec.typeMask & EGlossyReflection))
				return Spectrum(0.0f);

			/* Calculate the reflection half-vector (and possibly flip it
			   so that it lies inside the hemisphere around the normal) */
			H = normalize(bRec.wo+bRec.wi) 
				* signum(Frame::cosTheta(bRec.wo));
		} else {
			/* Stop if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 1)
				|| !(bRec.typeMask & EGlossyTransmission))
				return Spectrum(0.0f);

			/* Calculate the transmission half-vector (and possibly flip it
			   when the surface normal points into the denser medium -- this
			   removes an assumption in the original paper) */
			H = (m_extIOR > m_intIOR ? (Float) 1 : (Float) -1)
				* normalize(bRec.wi*etaI + bRec.wo*etaT);
		}

		/* Evaluate the roughness */
		const Float alpha = 
			std::max(m_alpha->getValue(bRec.its).average(), (Float) 1e-4f);

		/* Microsurface normal distribution */
		const Float D = evalD(H, alpha);
		if (D == 0)
			return Spectrum(0.0f);

		/* Fresnel factor */
		const Float F = fresnel(dot(bRec.wi, H), m_extIOR, m_intIOR);

		/* Smith's shadow-masking function */
		const Float G = smithG1(bRec.wi, H, alpha) * smithG1(bRec.wo, H, alpha);

		if (reflect) {
			/* Calculate the total amount of reflection */
			Float value = F * D * G / 
				(4.0f * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo));

			return m_specularReflectance->getValue(bRec.its) * value; 
		} else {
			/* Calculate the total amount of transmission */
			Float sqrtDenom = etaI * dot(bRec.wi, H) + etaT * dot(bRec.wo, H);
			Float value = ((1 - F) * D * G * etaT * etaT * dot(bRec.wi, H)*dot(bRec.wo, H)) / 
				(Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) * sqrtDenom * sqrtDenom);

			/* Missing term in the original paper: account for the solid angle 
			   compression when tracing radiance -- this is necessary for
			   bidirectional method */
			if (bRec.quantity == ERadiance)
				value *= (etaI*etaI) / (etaT*etaT);

			return m_specularTransmittance->getValue(bRec.its) * std::abs(value);
		}
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		/* Determine the type of interaction */
		bool sampleReflection   = ((bRec.component == -1 || bRec.component == 0)
							  && (bRec.typeMask & EGlossyReflection)),
		     sampleTransmission = ((bRec.component == -1 || bRec.component == 1)
							  && (bRec.typeMask & EGlossyTransmission)),
		     reflect            = Frame::cosTheta(bRec.wi) 
				                * Frame::cosTheta(bRec.wo) > 0;

		/* Determine the appropriate indices of refraction */
		Float etaI = m_extIOR, etaT = m_intIOR;
		if (Frame::cosTheta(bRec.wi) < 0)
			std::swap(etaI, etaT);

		Vector H;
		Float dwh_dwo;

		if (reflect) {
			/* Zero probability if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 0)
				|| !(bRec.typeMask & EGlossyReflection))
				return 0.0f;

			/* Calculate the reflection half-vector (and possibly flip it
			   so that it lies inside the hemisphere around the normal) */
			H = normalize(bRec.wo+bRec.wi) 
				* signum(Frame::cosTheta(bRec.wo));
	
			/* Jacobian of the half-direction transform */
			dwh_dwo = 1.0f / (4.0f * dot(bRec.wo, H));
		} else {
			/* Zero probability if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 1)
				|| !(bRec.typeMask & EGlossyTransmission))
				return 0.0f;

			/* Calculate the transmission half-vector (and possibly flip it
			   when the surface normal points into the denser medium -- this
			   removes an assumption in the original paper) */
			H = (m_extIOR > m_intIOR ? (Float) 1 : (Float) -1)
				* normalize(bRec.wi*etaI + bRec.wo*etaT);

			/* Jacobian of the half-direction transform. */
			Float sqrtDenom = etaI * dot(bRec.wi, H) + etaT * dot(bRec.wo, H);
			dwh_dwo = (etaT*etaT * dot(bRec.wo, H)) / (sqrtDenom*sqrtDenom);
		}

		/* Evaluate the roughness */
		Float alpha = 
			std::max(m_alpha->getValue(bRec.its).average(), (Float) 1e-4f);

		/* Suggestion by Bruce Walter: sample using a slightly different 
		   value of alpha. This in practice limits the weights to 
		   values <= 4. See also \ref sample() */
		alpha = alpha * (1.2f - 0.2f * std::sqrt(
			std::abs(Frame::cosTheta(bRec.wi))));

		/* Microsurface normal distribution */
		Float prob = evalD(H, alpha);

		if (sampleTransmission && sampleReflection) {
			/* Please see the sample() methods if the 
			   following seems confusing */
			Float F;
			if (bRec.sampler) {
				/* We have access to extra random numbers, hence
				   the exact Fresnel term with respect to the
				   microfacet normal is sampled */
				F = fresnel(dot(bRec.wi, H), m_extIOR, m_intIOR);
			} else {
				/* Only a simple 2D sample is given, and hence the
				   best we can do is to sample a clamped Fresnel term 
				   that is taken with respect to the surface normal */
				F = std::min((Float) 0.9f, std::max((Float) 0.1f,
					fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR)));
			}
			prob *= reflect ? F : (1-F);
		}

		return std::abs(prob * Frame::cosTheta(H) * dwh_dwo);
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &_sample) const {
		Point2 sample(_sample);

		bool sampleReflection = ((bRec.component == -1 || bRec.component == 0)
							  && (bRec.typeMask & EGlossyReflection)),
		     sampleTransmission = ((bRec.component == -1 || bRec.component == 1)
							  && (bRec.typeMask & EGlossyTransmission)),
		     choseReflection = sampleReflection,
			 sampleExactFresnelTerm = false;

		Float sampleF = 1.0f;
		if (sampleReflection && sampleTransmission) {
			if (!bRec.sampler) {
				/* By default, the sample() method is given exactly two
				   uniformly chosen random numbers, which is a problem
				   when wanting to choose between the reflection and the
				   transmission component: by the time the microfacet normal
				   has been sampled, we have already "used up" both of them,
				   and no random bits are left. Therefore, the following
				   somewhat crude approach is taken here when no extra sampler
				   instance is provided:
				      1. Take the Fresnel term with respect to the surface
					     normal to be a good approximation to the microsurface
						 Fresnel term -- this will be less true for higher 
						 rougness values. To be safe, clamp it to some 
						 reasonable range.
					  2. Use this approximate term and a random number to
					     choose between reflection and refraction component.
					  3. Recycle the used random number! This is possible,
					     since we so far only used the knowledge that it is
						 smaller or larget than the purely deterministic value
						 from step 2.
				*/
				sampleF = std::min((Float) 0.9f, std::max((Float) 0.1f,
					fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR)));
				if (sample.x < sampleF) {
					sample.x /= sampleF;
				} else {
					sample.x = (sample.x - sampleF) / (1 - sampleF);
					choseReflection = false;
				}
			} else {
				sampleExactFresnelTerm = true;
			}
		} else if (!sampleReflection && !sampleTransmission) {
			return Spectrum(0.0f);
		}

		/* Evaluate the roughness */
		Float alpha = 
			std::max(m_alpha->getValue(bRec.its).average(), (Float) 1e-4f);

		/* Suggestion by Bruce Walter: sample using a slightly different 
		   value of alpha. This in practice limits the weights to 
		   values <= 4. See also \ref sample() */
		Float sampleAlpha = alpha * (1.2f - 0.2f * std::sqrt(
			std::abs(Frame::cosTheta(bRec.wi))));

		/* Sample M, the microsurface normal */
		const Normal m = sampleD(sample, sampleAlpha);
	
		if (sampleExactFresnelTerm) {
			sampleF = fresnel(dot(bRec.wi, m), m_extIOR, m_intIOR);
			if (bRec.sampler->next1D() > sampleF)
				choseReflection = false;
		}

		Spectrum result;
		if (choseReflection) {
			/* Perfect specular reflection based on the microsurface normal */
			bRec.wo = reflect(bRec.wi, m);
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum(0.0f);
		
			result = m_specularReflectance->getValue(bRec.its);
		} else {
			/* Determine the appropriate indices of refraction */
			Float etaI = m_extIOR, etaT = m_intIOR;
			if (Frame::cosTheta(bRec.wi) < 0)
				std::swap(etaI, etaT);

			/* Perfect specular transmission based on the microsurface normal */
			if (!refract(bRec.wi, bRec.wo, m, etaI, etaT))
				return Spectrum(0.0f);

			bRec.sampledComponent = 1;
			bRec.sampledType = EGlossyTransmission;
		
			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
				return Spectrum(0.0f);

			result = m_specularTransmittance->getValue(bRec.its)
				* ((bRec.quantity == ERadiance) ? ((etaI*etaI) / (etaT*etaT)) : (Float) 1);
		}

		Float numerator = evalD(m, alpha)
			* smithG1(bRec.wi, m, alpha)
			* smithG1(bRec.wo, m, alpha)
			* dot(bRec.wi, m);

		Float denominator = evalD(m, sampleAlpha)
			* Frame::cosTheta(m) 
			* Frame::cosTheta(bRec.wi) 
			* Frame::cosTheta(bRec.wo);

		if (!sampleExactFresnelTerm) {
			Float F = fresnel(dot(bRec.wi, m), m_extIOR, m_intIOR);
			if (!choseReflection) {
				sampleF = 1-sampleF;
				F = 1-F;
			}
			numerator *= F;
			if (sampleReflection && sampleTransmission) 
				denominator *= sampleF;
		}

		return  result * std::abs(numerator / denominator);
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float _pdf, const Point2 &_sample) const {
		Point2 sample(_sample);

		bool sampleReflection = ((bRec.component == -1 || bRec.component == 0)
							  && (bRec.typeMask & EGlossyReflection)),
		     sampleTransmission = ((bRec.component == -1 || bRec.component == 1)
							  && (bRec.typeMask & EGlossyTransmission)),
		     choseReflection = sampleReflection,
			 sampleExactFresnelTerm = false;

		if (sampleReflection && sampleTransmission) {
			if (!bRec.sampler) {
				/* By default, the sample() method is given exactly two
				   uniformly chosen random numbers, which is a problem
				   when wanting to choose between the reflection and the
				   transmission component: by the time the microfacet normal
				   has been sampled, we have already "used up" both of them,
				   and no random bits are left. Therefore, the following
				   somewhat crude approach is taken here when no extra sampler
				   instance is provided:
				      1. Take the Fresnel term with respect to the surface
					     normal to be a good approximation to the microsurface
						 Fresnel term -- this will be less true for higher 
						 rougness values. To be safe, clamp it to some 
						 reasonable range.
					  2. Use this approximate term and a random number to
					     choose between reflection and refraction component.
					  3. Recycle the used random number! This is possible,
					     since we so far only used the knowledge that it is
						 smaller or larget than the purely deterministic value
						 from step 2.
				*/
				Float sampleF = std::min((Float) 0.9f, std::max((Float) 0.1f,
					fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR)));
				if (sample.x < sampleF) {
					sample.x /= sampleF;
				} else {
					sample.x = (sample.x - sampleF) / (1 - sampleF);
					choseReflection = false;
				}
			} else {
				sampleExactFresnelTerm = true;
			}
		} else if (!sampleReflection && !sampleTransmission) {
			return Spectrum(0.0f);
		}

		/* Evaluate the roughness */
		Float alpha = 
			std::max(m_alpha->getValue(bRec.its).average(), (Float) 1e-4f);

		/* Suggestion by Bruce Walter: sample using a slightly different 
		   value of alpha. This in practice limits the weights to 
		   values <= 4. See also \ref sample() */
		Float sampleAlpha = alpha * (1.2f - 0.2f * std::sqrt(
			std::abs(Frame::cosTheta(bRec.wi))));

		/* Sample M, the microsurface normal */
		const Normal m = sampleD(sample, sampleAlpha);
	
		if (sampleExactFresnelTerm) {
			Float sampleF = fresnel(dot(bRec.wi, m), m_extIOR, m_intIOR);
			if (bRec.sampler->next1D() > sampleF)
				choseReflection = false;
		}

		if (choseReflection) {
			/* Perfect specular reflection based on the microsurface normal */
			bRec.wo = reflect(bRec.wi, m);
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum(0.0f);
		} else {
			/* Determine the appropriate indices of refraction */
			Float etaI = m_extIOR, etaT = m_intIOR;
			if (Frame::cosTheta(bRec.wi) < 0)
				std::swap(etaI, etaT);

			/* Perfect specular transmission based on the microsurface normal */
			if (!refract(bRec.wi, bRec.wo, m, etaI, etaT))
				return Spectrum(0.0f);

			bRec.sampledComponent = 1;
			bRec.sampledType = EGlossyTransmission;
		
			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
				return Spectrum(0.0f);
		}

		/* Guard against numerical imprecisions */
		_pdf = pdf(bRec);

		if (_pdf == 0) 
			return Spectrum(0.0f);
		else
			return f(bRec);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "alpha") {
			m_alpha = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_alpha->usesRayDifferentials();
		} else if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "specularReflectance") {
			m_specularReflectance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_specularReflectance->usesRayDifferentials();
		} else if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "specularTransmittance") {
			m_specularTransmittance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_specularTransmittance->usesRayDifferentials();
		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeInt(m_distribution);
		manager->serialize(stream, m_alpha.get());
		manager->serialize(stream, m_specularReflectance.get());
		manager->serialize(stream, m_specularTransmittance.get());
		stream->writeFloat(m_intIOR);
		stream->writeFloat(m_extIOR);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "RoughDielectric[" << endl
			<< "  distribution = ";
		switch (m_distribution) {
			case EBeckmann: oss << "beckmann," << endl; break;
			case EGGX: oss << "ggx," << endl; break;
			case EPhong: oss << "phong," << endl; break;
			default:
				Log(EError, "Invalid distribution function");
		}
		oss << "  alpha = " << indent(m_alpha->toString()) << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  specularTransmittance = " << indent(m_specularTransmittance->toString()) << "," << endl
			<< "  intIOR = " << m_intIOR << "," << endl
			<< "  extIOR = " << m_extIOR << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	EDistribution m_distribution;
	ref<Texture> m_specularTransmittance;
	ref<Texture> m_specularReflectance;
	ref<Texture> m_alpha;
	Float m_intIOR, m_extIOR;
};

/* Fake glass shader -- it is really hopeless to visualize
   this material in the VPL renderer, so let's try to do at least 
   something that suggests the presence of a transparent boundary */
class RoughDielectricShader : public Shader {
public:
	RoughDielectricShader(Renderer *renderer) :
		Shader(renderer, EBSDFShader) {
		m_flags = ETransparent;
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return vec3(0.08);" << endl
			<< "}" << endl;
		oss << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return vec3(0.08);" << endl
			<< "}" << endl;
	}
	MTS_DECLARE_CLASS()
};

Shader *RoughDielectric::createShader(Renderer *renderer) const { 
	return new RoughDielectricShader(renderer);
}

MTS_IMPLEMENT_CLASS(RoughDielectricShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(RoughDielectric, false, BSDF)
MTS_EXPORT_PLUGIN(RoughDielectric, "Rough glass BSDF");
MTS_NAMESPACE_END

/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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
#include <mitsuba/hw/basicshader.h>
#include "microfacet.h"
#include "ior.h"

MTS_NAMESPACE_BEGIN

/* Suggestion by Bruce Walter: sample the model using a slightly
   wider density function. This in practice limits the importance
   weights to values <= 4.
*/
#define ENLARGE_LOBE_TRICK 1

/*!\plugin{roughdielectric}{Rough dielectric material}
 * \order{5}
 * \icon{bsdf_roughdielectric}
 * \parameters{
 *     \parameter{distribution}{\String}{
 *          Specifies the type of microfacet normal distribution
 *          used to model the surface roughness.
 *       \begin{enumerate}[(i)]
 *           \item \code{beckmann}: Physically-based distribution derived from
 *               Gaussian random surfaces. This is the default.
 *           \item \code{ggx}: New distribution proposed by
 *              Walter et al. \cite{Walter07Microfacet}, which is meant to better handle
 *              the long tails observed in measurements of ground surfaces.
 *              Renderings with this distribution may converge slowly.
 *           \item \code{phong}: Classical $\cos^p\theta$ distribution.
 *              Due to the underlying microfacet theory,
 *              the use of this distribution here leads to more realistic
 *              behavior than the separately available \pluginref{phong} plugin.
 *           \item \code{as}: Anisotropic Phong-style microfacet distribution proposed by
 *              Ashikhmin and Shirley \cite{Ashikhmin2005Anisotropic}.\vspace{-3mm}
 *       \end{enumerate}
 *     }
 *     \parameter{alpha}{\Float\Or\Texture}{
 *         Specifies the roughness of the unresolved surface micro-geometry.
 *         When the Beckmann distribution is used, this parameter is equal to the
 *         \emph{root mean square} (RMS) slope of the microfacets. This
 *         parameter is only valid when \texttt{distribution=beckmann/phong/ggx}.
 *         \default{0.1}.
 *     }
 *     \parameter{alphaU, alphaV}{\Float\Or\Texture}{
 *         Specifies the anisotropic roughness values along the tangent and
 *         bitangent directions. These parameter are only valid when
 *         \texttt{distribution=as}. \default{0.1}.
 *     }
 *     \parameter{intIOR}{\Float\Or\String}{Interior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{bk7} / 1.5046}}
 *     \parameter{extIOR}{\Float\Or\String}{Exterior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{air} / 1.000277}}
 *     \parameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor that can be used to modulate the specular reflection component. Note
 *         that for physical realism, this parameter should never be touched. \default{1.0}}
 *     \parameter{specular\showbreak Transmittance}{\Spectrum\Or\Texture}{Optional
 *         factor that can be used to modulate the specular transmission component. Note
 *         that for physical realism, this parameter should never be touched. \default{1.0}}
 * }\vspace{4mm}
 *
 * This plugin implements a realistic microfacet scattering model for rendering
 * rough interfaces between dielectric materials, such as a transition from air to
 * ground glass. Microfacet theory describes rough surfaces as an arrangement of
 * unresolved and ideally specular facets, whose normal directions are given by
 * a specially chosen \emph{microfacet distribution}. By accounting for shadowing
 * and masking effects between these facets, it is possible to reproduce the important
 * off-specular reflections peaks observed in real-world measurements of such
 * materials.
 * \renderings{
 *     \rendering{Anti-glare glass (Beckmann, $\alpha=0.02$)}
 *     	   {bsdf_roughdielectric_beckmann_0_0_2.jpg}
 *     \rendering{Rough glass (Beckmann, $\alpha=0.1$)}
 *     	   {bsdf_roughdielectric_beckmann_0_1.jpg}
 * }
 *
 * This plugin is essentially the ``roughened'' equivalent of the (smooth) plugin
 * \pluginref{dielectric}. For very low values of $\alpha$, the two will
 * be identical, though scenes using this plugin will take longer to render
 * due to the additional computational burden of tracking surface roughness.
 *
 * The implementation is based on the paper ``Microfacet Models
 * for Refraction through Rough Surfaces'' by Walter et al.
 * \cite{Walter07Microfacet}. It supports several different types of microfacet
 * distributions and has a texturable roughness parameter. Exterior and
 * interior IOR values can be specified independently, where ``exterior''
 * refers to the side that contains the surface normal. Similar to the
 * \pluginref{dielectric} plugin, IOR values can either be specified
 * numerically, or based on a list of known materials (see
 * \tblref{dielectric-iors} for an overview). When no parameters are given,
 * the plugin activates the default settings, which describe a borosilicate
 * glass BK7/air interface with a light amount of roughness modeled using a
 * Beckmann distribution.
 *
 * To get an intuition about the effect of the surface roughness
 * parameter $\alpha$, consider the following approximate classification:
 * a value of $\alpha=0.001-0.01$ corresponds to a material
 * with slight imperfections on an
 * otherwise smooth surface finish, $\alpha=0.1$ is relatively rough,
 * and $\alpha=0.3-0.7$ is \emph{extremely} rough (e.g. an etched or ground
 * finish).
 *
 * Please note that when using this plugin, it is crucial that the scene contains
 * meaningful and mutually compatible index of refraction changes---see
 * \figref{glass-explanation} for an example of what this entails. Also, note that
 * the importance sampling implementation of this model is close, but
 * not always a perfect a perfect match to the underlying scattering distribution,
 * particularly for high roughness values and when the \texttt{ggx}
 * microfacet distribution is used. Hence, such renderings may
 * converge slowly.
 *
 * \subsubsection*{Technical details}
 * When rendering with the Ashikhmin-Shirley or Phong microfacet
 * distributions, a conversion is used to turn the specified
 * $\alpha$ roughness value into the exponents of these distributions.
 * This is done in a way, such that the different
 * distributions all produce a similar appearance for the same value of
 * $\alpha$.
 *
 * The Ashikhmin-Shirley microfacet distribution allows the specification
 * of two distinct roughness values along the tangent and bitangent
 * directions. This can be used to provide a material with a ``brushed''
 * appearance. The alignment of the anisotropy will follow the UV
 * parameterization of the underlying mesh in this case. This also means that
 * such an anisotropic material cannot be applied to triangle meshes that
 * are missing texture coordinates.\newpage
 *
 * \renderings{
 *     \rendering{Ground glass (GGX, $\alpha$=0.304,
 *     	   \lstref{roughdielectric-roughglass})}{bsdf_roughdielectric_ggx_0_304.jpg}
 *     \rendering{Textured roughness (\lstref{roughdielectric-textured})}
 *         {bsdf_roughdielectric_textured.jpg}
 * }
 *
 * \begin{xml}[caption=A material definition for ground glass, label=lst:roughdielectric-roughglass]
 * <bsdf type="roughdielectric">
 *     <string name="distribution" value="ggx"/>
 *     <float name="alpha" value="0.304"/>
 *     <string name="intIOR" value="bk7"/>
 *     <string name="extIOR" value="air"/>
 * </bsdf>
 * \end{xml}
 *
 * \begin{xml}[caption=A texture can be attached to the roughness parameter, label=lst:roughdielectric-textured]
 * <bsdf type="roughdielectric">
 *     <string name="distribution" value="beckmann"/>
 *     <float name="intIOR" value="1.5046"/>
 *     <float name="extIOR" value="1.0"/>
 *
 *     <texture name="alpha" type="bitmap">
 *         <string name="filename" value="roughness.exr"/>
 *     </texture>
 * </bsdf>
 * \end{xml}
 */
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

		m_distribution = MicrofacetDistribution(
			props.getString("distribution", "beckmann")
		);

		Float alpha = props.getFloat("alpha", 0.1f),
			  alphaU = props.getFloat("alphaU", alpha),
			  alphaV = props.getFloat("alphaV", alpha);

		m_alphaU = new ConstantFloatTexture(alphaU);
		if (alphaU == alphaV)
			m_alphaV = m_alphaU;
		else
			m_alphaV = new ConstantFloatTexture(alphaV);
	}

	RoughDielectric(Stream *stream, InstanceManager *manager)
	 : BSDF(stream, manager) {
		m_distribution = MicrofacetDistribution(
			(MicrofacetDistribution::EType) stream->readUInt()
		);
		m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaV = static_cast<Texture *>(manager->getInstance(stream));
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_specularTransmittance = static_cast<Texture *>(manager->getInstance(stream));
		m_eta = stream->readFloat();
		m_invEta = 1 / m_eta;

		configure();
	}

	void configure() {
		unsigned int extraFlags = 0;
		if (m_alphaU != m_alphaV) {
			extraFlags |= EAnisotropic;
			if (m_distribution.getType() !=
				MicrofacetDistribution::EAshikhminShirley)
				Log(EError, "Different roughness values along the tangent and "
						"bitangent directions are only supported when using the "
						"anisotropic Ashikhmin-Shirley microfacet distribution "
						"(named \"as\")");
		}

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
		if (measure != ESolidAngle)
			return Spectrum(0.0f);

		/* Determine the type of interaction */
		bool reflect = Frame::cosTheta(bRec.wi)
			* Frame::cosTheta(bRec.wo) > 0;

		Vector H;
		if (reflect) {
			/* Stop if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 0)
				|| !(bRec.typeMask & EGlossyReflection))
				return Spectrum(0.0f);

			/* Calculate the reflection half-vector */
			H = normalize(bRec.wo+bRec.wi);
		} else {
			/* Stop if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 1)
				|| !(bRec.typeMask & EGlossyTransmission))
				return Spectrum(0.0f);

			/* Calculate the transmission half-vector */
			Float eta = Frame::cosTheta(bRec.wi) > 0
				? m_eta : m_invEta;

			H = normalize(bRec.wi + bRec.wo*eta);
		}

		/* Ensure that the half-vector points into the
		   same hemisphere as the macrosurface normal */
		H *= math::signum(Frame::cosTheta(H));

		/* Evaluate the roughness */
		Float alphaU = m_distribution.transformRoughness(
					m_alphaU->eval(bRec.its).average()),
			  alphaV = m_distribution.transformRoughness(
					m_alphaV->eval(bRec.its).average());

		/* Evaluate the microsurface normal distribution */
		const Float D = m_distribution.eval(H, alphaU, alphaV);
		if (D == 0)
			return Spectrum(0.0f);

		/* Fresnel factor */
		const Float F = fresnelDielectricExt(dot(bRec.wi, H), m_eta);

		/* Smith's shadow-masking function */
		const Float G = m_distribution.G(bRec.wi, bRec.wo, H, alphaU, alphaV);

		if (reflect) {
			/* Calculate the total amount of reflection */
			Float value = F * D * G /
				(4.0f * std::abs(Frame::cosTheta(bRec.wi)));

			return m_specularReflectance->eval(bRec.its) * value;
		} else {
			Float eta = Frame::cosTheta(bRec.wi) > 0.0f ? m_eta : m_invEta;

			/* Calculate the total amount of transmission */
			Float sqrtDenom = dot(bRec.wi, H) + eta * dot(bRec.wo, H);
			Float value = ((1 - F) * D * G * eta * eta
				* dot(bRec.wi, H) * dot(bRec.wo, H)) /
				(Frame::cosTheta(bRec.wi) * sqrtDenom * sqrtDenom);

			/* Missing term in the original paper: account for the solid angle
			   compression when tracing radiance -- this is necessary for
			   bidirectional methods */
			Float factor = (bRec.mode == ERadiance)
				? (Frame::cosTheta(bRec.wi) > 0 ? m_invEta : m_eta) : 1.0f;

			return m_specularTransmittance->eval(bRec.its)
				* std::abs(value * factor * factor);
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

		Vector H;
		Float dwh_dwo;

		if (reflect) {
			/* Zero probability if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 0)
				|| !(bRec.typeMask & EGlossyReflection))
				return 0.0f;

			/* Calculate the reflection half-vector */
			H = normalize(bRec.wo+bRec.wi);

			/* Jacobian of the half-direction mapping */
			dwh_dwo = 1.0f / (4.0f * dot(bRec.wo, H));
		} else {
			/* Zero probability if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 1)
				|| !(bRec.typeMask & EGlossyTransmission))
				return 0.0f;

			/* Calculate the transmission half-vector */
			Float eta = Frame::cosTheta(bRec.wi) > 0
				? m_eta : m_invEta;

			H = normalize(bRec.wi + bRec.wo*eta);

			/* Jacobian of the half-direction mapping */
			Float sqrtDenom = dot(bRec.wi, H) + eta * dot(bRec.wo, H);
			dwh_dwo = (eta*eta * dot(bRec.wo, H)) / (sqrtDenom*sqrtDenom);
		}

		/* Ensure that the half-vector points into the
		   same hemisphere as the macrosurface normal */
		H *= math::signum(Frame::cosTheta(H));

		/* Evaluate the roughness */
		Float alphaU = m_alphaU->eval(bRec.its).average(),
			  alphaV = m_alphaV->eval(bRec.its).average();

#if ENLARGE_LOBE_TRICK == 1
		Float factor = (1.2f - 0.2f * std::sqrt(
			std::abs(Frame::cosTheta(bRec.wi))));
		alphaU *= factor; alphaV *= factor;
#endif

		alphaU = m_distribution.transformRoughness(alphaU);
		alphaV = m_distribution.transformRoughness(alphaV);

		/* Evaluate the microsurface normal sampling density */
		Float prob = m_distribution.pdf(H, alphaU, alphaV);

		if (hasTransmission && hasReflection) {
			Float F = fresnelDielectricExt(dot(bRec.wi, H), m_eta);
			prob *= reflect ? F : (1-F);
		}

		return std::abs(prob * dwh_dwo);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &_sample) const {
		Point2 sample(_sample);

		bool hasReflection = ((bRec.component == -1 || bRec.component == 0)
							  && (bRec.typeMask & EGlossyReflection)),
		     hasTransmission = ((bRec.component == -1 || bRec.component == 1)
							  && (bRec.typeMask & EGlossyTransmission)),
		     sampleReflection = hasReflection;

		if (!hasReflection && !hasTransmission)
			return Spectrum(0.0f);

		/* Evaluate the roughness */
		Float alphaU = m_alphaU->eval(bRec.its).average(),
		      alphaV = m_alphaV->eval(bRec.its).average(),
		      sampleAlphaU = alphaU,
		      sampleAlphaV = alphaV;

#if ENLARGE_LOBE_TRICK == 1
		Float factor = (1.2f - 0.2f * std::sqrt(
			std::abs(Frame::cosTheta(bRec.wi))));
		sampleAlphaU *= factor; sampleAlphaV *= factor;
#endif

		alphaU = m_distribution.transformRoughness(alphaU);
		alphaV = m_distribution.transformRoughness(alphaV);
		sampleAlphaU = m_distribution.transformRoughness(sampleAlphaU);
		sampleAlphaV = m_distribution.transformRoughness(sampleAlphaV);

		/* Sample M, the microsurface normal */
		Float microfacetPDF;
		const Normal m = m_distribution.sample(sample,
				sampleAlphaU, sampleAlphaV, microfacetPDF);

		if (microfacetPDF == 0)
			return Spectrum(0.0f);

		Float cosThetaT, numerator = 1.0f;
		Float F = fresnelDielectricExt(dot(bRec.wi, m), cosThetaT, m_eta);

		if (hasReflection && hasTransmission) {
			if (bRec.sampler->next1D() > F)
				sampleReflection = false;
		} else {
			numerator = hasReflection ? F : (1-F);
		}

		Spectrum result;
		if (sampleReflection) {
			/* Perfect specular reflection based on the microsurface normal */
			bRec.wo = reflect(bRec.wi, m);
			bRec.eta = 1.0f;
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum(0.0f);

			result = m_specularReflectance->eval(bRec.its);
		} else {
			if (cosThetaT == 0)
				return Spectrum(0.0f);

			/* Perfect specular transmission based on the microsurface normal */
			bRec.wo = refract(bRec.wi, m, m_eta, cosThetaT);
			bRec.eta = cosThetaT < 0 ? m_eta : m_invEta;
			bRec.sampledComponent = 1;
			bRec.sampledType = EGlossyTransmission;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
				return Spectrum(0.0f);

			/* Radiance must be scaled to account for the solid angle compression
			   that occurs when crossing the interface. */
			Float factor = (bRec.mode == ERadiance)
				? (cosThetaT < 0 ? m_invEta : m_eta) : 1.0f;

			result = m_specularTransmittance->eval(bRec.its) * (factor * factor);
		}

		numerator *= m_distribution.eval(m, alphaU, alphaV)
			* m_distribution.G(bRec.wi, bRec.wo, m, alphaU, alphaV)
			* dot(bRec.wi, m);

		Float denominator = microfacetPDF
			* Frame::cosTheta(bRec.wi);

		return result * std::abs(numerator / denominator);
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

		/* Evaluate the roughness */
		Float alphaU = m_alphaU->eval(bRec.its).average(),
		      alphaV = m_alphaV->eval(bRec.its).average(),
		      sampleAlphaU = alphaU,
		      sampleAlphaV = alphaV;

#if ENLARGE_LOBE_TRICK == 1
		Float factor = (1.2f - 0.2f * std::sqrt(
			std::abs(Frame::cosTheta(bRec.wi))));
		sampleAlphaU *= factor; sampleAlphaV *= factor;
#endif

		alphaU = m_distribution.transformRoughness(alphaU);
		alphaV = m_distribution.transformRoughness(alphaV);
		sampleAlphaU = m_distribution.transformRoughness(sampleAlphaU);
		sampleAlphaV = m_distribution.transformRoughness(sampleAlphaV);

		/* Sample M, the microsurface normal */
		Float microfacetPDF;
		const Normal m = m_distribution.sample(sample,
				sampleAlphaU, sampleAlphaV, microfacetPDF);

		if (microfacetPDF == 0)
			return Spectrum(0.0f);

		pdf = microfacetPDF;

		Float cosThetaT, numerator = 1.0f;
		Float F = fresnelDielectricExt(dot(bRec.wi, m), cosThetaT, m_eta);

		if (hasReflection && hasTransmission) {
			if (bRec.sampler->next1D() > F) {
				sampleReflection = false;
				pdf *= 1-F;
			} else {
				pdf *= F;
			}
		} else {
			numerator = hasReflection ? F : (1-F);
		}

		Spectrum result;
		Float dwh_dwo;

		if (sampleReflection) {
			/* Perfect specular reflection based on the microsurface normal */
			bRec.wo = reflect(bRec.wi, m);
			bRec.eta = 1.0f;
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum(0.0f);

			result = m_specularReflectance->eval(bRec.its);

			/* Jacobian of the half-direction mapping */
			dwh_dwo = 1.0f / (4.0f * dot(bRec.wo, m));
		} else {
			if (cosThetaT == 0)
				return Spectrum(0.0f);

			/* Perfect specular transmission based on the microsurface normal */
			bRec.wo = refract(bRec.wi, m, m_eta, cosThetaT);
			bRec.eta = cosThetaT < 0 ? m_eta : m_invEta;
			bRec.sampledComponent = 1;
			bRec.sampledType = EGlossyTransmission;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
				return Spectrum(0.0f);

			/* Radiance must be scaled to account for the solid angle compression
			   that occurs when crossing the interface. */
			Float factor = (bRec.mode == ERadiance)
				? (cosThetaT < 0 ? m_invEta : m_eta) : 1.0f;

			result = m_specularTransmittance->eval(bRec.its) * (factor * factor);

			/* Jacobian of the half-direction mapping */
			Float sqrtDenom = dot(bRec.wi, m) + bRec.eta * dot(bRec.wo, m);
			dwh_dwo = (bRec.eta*bRec.eta * dot(bRec.wo, m)) / (sqrtDenom*sqrtDenom);
		}

		numerator *= m_distribution.eval(m, alphaU, alphaV)
			* m_distribution.G(bRec.wi, bRec.wo, m, alphaU, alphaV)
			* dot(bRec.wi, m);

		Float denominator = microfacetPDF * Frame::cosTheta(bRec.wi);

		pdf *= std::abs(dwh_dwo);

		return result * std::abs(numerator / denominator);
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

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeUInt((uint32_t) m_distribution.getType());
		manager->serialize(stream, m_alphaU.get());
		manager->serialize(stream, m_alphaV.get());
		manager->serialize(stream, m_specularReflectance.get());
		manager->serialize(stream, m_specularTransmittance.get());
		stream->writeFloat(m_eta);
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
			<< "  distribution = " << m_distribution.toString() << "," << endl
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
	MicrofacetDistribution m_distribution;
	ref<Texture> m_specularTransmittance;
	ref<Texture> m_specularReflectance;
	ref<Texture> m_alphaU, m_alphaV;
	Float m_eta, m_invEta;
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

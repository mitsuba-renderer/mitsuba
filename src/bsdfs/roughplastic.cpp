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
#include <mitsuba/hw/basicshader.h>
#include "microfacet.h"
#include "rtrans.h"
#include "ior.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{roughplastic}{Rough plastic material}
 * \order{8}
 * \icon{bsdf_roughplastic}
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
 *              \vspace{-4mm}
 *       \end{enumerate}
 *     }
 *     \parameter{alpha}{\Float}{
 *         Specifies the roughness of the unresolved surface micro-geometry. 
 *         When the Beckmann distribution is used, this parameter is equal to the 
 *         \emph{root mean square} (RMS) slope of the microfacets. 
 *         \default{0.1}. 
 *     }
 *     \parameter{intIOR}{\Float\Or\String}{Interior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{polypropylene} / 1.49}}
 *     \parameter{extIOR}{\Float\Or\String}{Exterior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{air} / 1.000277}}
 *     \parameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor that can be used to modulate the specular reflection component. Note 
 *         that for physical realism, this parameter should never be touched. \default{1.0}}
 *     \parameter{diffuse\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor used to modulate the diffuse reflection component\default{0.5}}
 *     \parameter{preserveColors}{\Boolean}{
 *         Account for color shifts due to internal scattering? See the main text
 *         for a detailed description.\default{Don't account for them and
 *         preserve the colors, i.e. \code{true}}
 * }\vspace{-1mm}
 * \renderings{
 *     \rendering{Beckmann, $\alpha=0.1$}{bsdf_roughplastic_beckmann}
 *     \rendering{GGX, $\alpha=0.3$}{bsdf_roughplastic_ggx}
 * }\vspace{-1mm}
 *
 * This plugin implements a realistic microfacet scattering model for rendering
 * rough dielectric materials with internal scattering, such as plastic. It can 
 * be interpreted as a fancy version of the Cook-Torrance model and should be 
 * preferred over empirical models like \pluginref{phong} and \pluginref{ward} 
 * when possible.
 * \renderings{
 *     \setcounter{subfigure}{2}
 *     \rendering{Beckmann, $\alpha=0.05$, diffuseReflectance=0}
 *         {bsdf_roughplastic_beckmann_lacquer}
 * }
 *
 * Microfacet theory describes rough surfaces as an arrangement of unresolved and 
 * ideally specular facets, whose normal directions are given by a specially
 * chosen \emph{microfacet distribution}. 
 * By accounting for shadowing and masking effects between these facets, it is 
 * possible to reproduce the important off-specular reflections peaks observed 
 * in real-world measurements of such materials.
 *
 * This plugin is essentially the ``roughened'' equivalent of the (smooth) plugin
 * \pluginref{plastic}. For very low values of $\alpha$, the two will
 * be very similar, though scenes using this plugin will take longer to render 
 * due to the additional computational burden of tracking surface roughness.
 *
 * The model uses the integrated specular reflectance to interpolate between the 
 * specular and diffuse components (i.e. any light that is not scattered
 * specularly is assumed to contribute to the diffuse component).
 * Similar to the \pluginref{dielectric} plugin, IOR values 
 * can either be specified numerically, or based on a list of known materials 
 * (see \tblref{dielectric-iors} for an overview). 
 * 
 * The implementation is based on the paper ``Microfacet Models
 * for Refraction through Rough Surfaces'' by Walter et al. 
 * \cite{Walter07Microfacet}. It supports several different types of microfacet
 * distributions. Note that the choices are a bit more restricted here---in 
 * comparison to other rough scattering models in Mitsuba,
 * the roughness cannot be textured, and anisotropic microfacet 
 * distributions are not allowed.
 *
 * When no parameters are given, the plugin activates the defaults, 
 * which describe a white polypropylene plastic material with a light amount
 * of roughness modeled using the Beckmann distribution.
 *
 * To get an intuition about the effect of the surface roughness
 * parameter $\alpha$, consider the following approximate differentiation: 
 * a value of $\alpha=0.001-0.01$ corresponds to a material 
 * with slight imperfections on an
 * otherwise smooth surface finish, $\alpha=0.1$ is relatively rough,
 * and $\alpha=0.3-0.7$ is \emph{extremely} rough (e.g. an etched or ground
 * finish). Values significantly above that are probably not too realistic.
*
 * When rendering with the Phong microfacet 
 * distributions, a conversion is used to turn the specified 
 * $\alpha$ roughness value into the Phong exponent.
 * This is done in a way, such that the different 
 * distributions all produce a similar appearance for 
 * the same value of $\alpha$.
 * \begin{xml}[caption={A material definition for rough, black laquer.}, label=lst:roughplastic-lacquer]
 * <bsdf type="roughplastic">
 *     <string name="distribution" value="beckmann"/>
 *     <float name="alpha" value="0.05"/>
 *     <float name="intIOR" value="1.61"/>
 *     <spectrum name="diffuseReflectance" value="0"/>
 * </bsdf>
 * \end{xml}
 *
 */
class RoughPlastic : public BSDF {
public:
	RoughPlastic(const Properties &props) : BSDF(props) {
		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));
		m_diffuseReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("diffuseReflectance", Spectrum(0.5f)));

		/* Specifies the internal index of refraction at the interface */
		m_intIOR = lookupIOR(props, "intIOR", "polypropylene");

		/* Specifies the external index of refraction at the interface */
		m_extIOR = lookupIOR(props, "extIOR", "air");

		if (m_intIOR < 0 || m_extIOR < 0 || m_intIOR == m_extIOR)
			Log(EError, "The interior and exterior indices of "
				"refraction must be positive and differ!");

		m_distribution = MicrofacetDistribution(
			props.getString("distribution", "beckmann")
		);

		if (m_distribution.isAnisotropic())
			Log(EError, "The 'roughplastic' plugin currently does not support "
				"anisotropic microfacet distributions!");

		m_preserveColors = props.getBoolean("preserveColors", true);

		m_alpha = new ConstantFloatTexture(
			props.getFloat("alpha", 0.1f));

		m_specularSamplingWeight = 0.0f;
	}

	RoughPlastic(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager) {
		m_distribution = MicrofacetDistribution(
			(MicrofacetDistribution::EType) stream->readUInt()
		);
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_diffuseReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_alpha = static_cast<Texture *>(manager->getInstance(stream));
		m_intIOR = stream->readFloat();
		m_extIOR = stream->readFloat();
		m_preserveColors = stream->readBool();

		configure();
	}

	void configure() {
		m_components.clear();

		bool constAlpha = m_alpha->isConstant();

		m_components.push_back(EGlossyReflection | EFrontSide 
			| ((constAlpha && m_specularReflectance->isConstant())
				? 0 : ESpatiallyVarying));
		m_components.push_back(EDiffuseReflection | EFrontSide 
			| ((constAlpha && m_diffuseReflectance->isConstant()) 
				? 0 : ESpatiallyVarying));

		/* Verify the input parameters and fix them if necessary */
		m_specularReflectance = ensureEnergyConservation(
			m_specularReflectance, "specularReflectance", 1.0f);
		m_diffuseReflectance = ensureEnergyConservation(
			m_diffuseReflectance, "diffuseReflectance", 1.0f);

		/* Compute weights that further steer samples towards
		   the specular or diffuse components */
		Float dAvg = m_diffuseReflectance->getAverage().getLuminance(),
			  sAvg = m_specularReflectance->getAverage().getLuminance();
		m_specularSamplingWeight = sAvg / (dAvg + sAvg);
			
		Float eta = m_intIOR / m_extIOR;

		if (!m_roughTransmittance.get()) {
			/* Load precomputed data used to compute the rough
			   transmittance through the dielectric interface */
			m_roughTransmittance = new RoughTransmittance(
				m_distribution.getType());

			m_roughTransmittance->checkEta(eta);
			m_roughTransmittance->checkAlpha(m_alpha->getMinimum().average());
			m_roughTransmittance->checkAlpha(m_alpha->getMaximum().average());

			/* Reduce the rough transmittance data to a 2D slice */
			m_roughTransmittance->setEta(eta);

			/* If possible, even reduce it to a 1D slice */
			if (constAlpha) 
				m_roughTransmittance->setAlpha(
					m_alpha->getValue(Intersection()).average());
		}

		m_usesRayDifferentials = 
			m_specularReflectance->usesRayDifferentials() ||
			m_diffuseReflectance->usesRayDifferentials();

		BSDF::configure();
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return m_diffuseReflectance->getValue(its);
	}

	/// Helper function: reflect \c wi with respect to a given surface normal
	inline Vector reflect(const Vector &wi, const Normal &m) const {
		return 2 * dot(wi, m) * Vector(m) - wi;
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		bool hasSpecular = (bRec.typeMask & EGlossyReflection) &&
			(bRec.component == -1 || bRec.component == 0);
		bool hasDiffuse = (bRec.typeMask & EDiffuseReflection) &&
			(bRec.component == -1 || bRec.component == 1);
			
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			(!hasSpecular && !hasDiffuse))
			return Spectrum(0.0f);

		/* Evaluate the roughness texture */
		Float alpha = m_alpha->getValue(bRec.its).average();
		Float alphaT = m_distribution.transformRoughness(alpha);
	
		Spectrum result(0.0f);
		if (hasSpecular) {
			/* Calculate the reflection half-vector */
			const Vector H = normalize(bRec.wo+bRec.wi);

			/* Evaluate the microsurface normal distribution */
			const Float D = m_distribution.eval(H, alphaT);

			/* Fresnel term */
			const Float F = fresnel(dot(bRec.wi, H), m_extIOR, m_intIOR);

			/* Smith's shadow-masking function */
			const Float G = m_distribution.G(bRec.wi, bRec.wo, H, alphaT);

			/* Calculate the specular reflection component */
			Float value = F * D * G / 
				(4.0f * Frame::cosTheta(bRec.wi));

			result += m_specularReflectance->getValue(bRec.its) * value; 
		}

		if (hasDiffuse) { 
			Spectrum diff = m_diffuseReflectance->getValue(bRec.its);
			Float T12 = m_roughTransmittance->eval(Frame::cosTheta(bRec.wi), alpha);
			Float T21 = m_roughTransmittance->eval(Frame::cosTheta(bRec.wo), alpha);
			Float Fdr = 1-m_roughTransmittance->evalDiffuse(alpha);

			if (m_preserveColors)
				diff /= 1-Fdr;
			else
				diff /= Spectrum(1.0f) - diff * Fdr;

			result += diff * (INV_PI * Frame::cosTheta(bRec.wo) * T12 * T21);
		}

		return result;
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		bool hasSpecular = (bRec.typeMask & EGlossyReflection) &&
			(bRec.component == -1 || bRec.component == 0);
		bool hasDiffuse = (bRec.typeMask & EDiffuseReflection) &&
			(bRec.component == -1 || bRec.component == 1);

		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			(!hasSpecular && !hasDiffuse))
			return 0.0f;

		/* Evaluate the roughness texture */
		Float alpha = m_alpha->getValue(bRec.its).average();
		Float alphaT = m_distribution.transformRoughness(alpha);

		/* Calculate the reflection half-vector */
		const Vector H = normalize(bRec.wo+bRec.wi);

		Float probDiffuse, probSpecular;
		if (hasSpecular && hasDiffuse) {
			/* Find the probability of sampling the specular component */
			probSpecular = 1-m_roughTransmittance->eval(Frame::cosTheta(bRec.wi), alpha);

			/* Reallocate samples */
			probSpecular = (probSpecular*m_specularSamplingWeight) /
				(probSpecular*m_specularSamplingWeight + 
				(1-probSpecular) * (1-m_specularSamplingWeight));

			probDiffuse = 1 - probSpecular;
		} else {
			probDiffuse = probSpecular = 1.0f;
		}

		Float result = 0.0f;
		if (hasSpecular) {
			/* Jacobian of the half-direction transform */
			const Float dwh_dwo = 1.0f / (4.0f * dot(bRec.wo, H));

			/* Evaluate the microsurface normal distribution */
			const Float prob = m_distribution.pdf(H, alphaT);

			result = prob * dwh_dwo * probSpecular;
		}

		if (hasDiffuse) 
			result += Frame::cosTheta(bRec.wo) * INV_PI * probDiffuse;

		return result;
	}

	inline Spectrum sample(BSDFQueryRecord &bRec, Float &_pdf, const Point2 &_sample) const {
		bool hasSpecular = (bRec.typeMask & EGlossyReflection) &&
			(bRec.component == -1 || bRec.component == 0);
		bool hasDiffuse = (bRec.typeMask & EDiffuseReflection) &&
			(bRec.component == -1 || bRec.component == 1);
		
		if (Frame::cosTheta(bRec.wi) <= 0 || (!hasSpecular && !hasDiffuse))
			return Spectrum(0.0f);

		bool choseSpecular = hasSpecular;
		Point2 sample(_sample);

		/* Evaluate the roughness texture */
		Float alpha = m_alpha->getValue(bRec.its).average();
		Float alphaT = m_distribution.transformRoughness(alpha);

		Float probSpecular;
		if (hasSpecular && hasDiffuse) {
			/* Find the probability of sampling the specular component */
			probSpecular = 1 - m_roughTransmittance->eval(Frame::cosTheta(bRec.wi), alpha);

			/* Reallocate samples */
			probSpecular = (probSpecular*m_specularSamplingWeight) /
				(probSpecular*m_specularSamplingWeight + 
				(1-probSpecular) * (1-m_specularSamplingWeight));

			if (sample.x <= probSpecular) {
				sample.x /= probSpecular;
			} else {
				sample.x = (sample.x - probSpecular) / (1 - probSpecular);
				choseSpecular = false;
			}
		}

		if (choseSpecular) {
			/* Perfect specular reflection based on the microsurface normal */
			Normal m = m_distribution.sample(sample, alphaT);
			bRec.wo = reflect(bRec.wi, m);
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum(0.0f);
		} else {
			bRec.sampledComponent = 1;
			bRec.sampledType = EDiffuseReflection;
			bRec.wo = squareToHemispherePSA(sample);
		}

		/* Guard against numerical imprecisions */
		_pdf = pdf(bRec, ESolidAngle);

		if (_pdf == 0) 
			return Spectrum(0.0f);
		else
			return eval(bRec, ESolidAngle) / _pdf;
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		Float pdf;
		return RoughPlastic::sample(bRec, pdf, sample);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeUInt((uint32_t) m_distribution.getType());
		manager->serialize(stream, m_specularReflectance.get());
		manager->serialize(stream, m_diffuseReflectance.get());
		manager->serialize(stream, m_alpha.get());
		stream->writeFloat(m_intIOR);
		stream->writeFloat(m_extIOR);
		stream->writeBool(m_preserveColors);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "alpha")
				m_alpha = static_cast<Texture *>(child);
			else if (name == "specularReflectance") 
				m_specularReflectance = static_cast<Texture *>(child);
			else if (name == "diffuseReflectance")
				m_diffuseReflectance = static_cast<Texture *>(child);
			else 
				BSDF::addChild(name, child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "RoughPlastic[" << endl
			<< "  name = \"" << getName() << "\"," << endl
			<< "  distribution = " << m_distribution.toString() << "," << endl
			<< "  alpha = " << indent(m_alpha->toString()) << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  diffuseReflectance = " << indent(m_diffuseReflectance->toString()) << "," << endl
			<< "  specularSamplingWeight = " << m_specularSamplingWeight << "," << endl
			<< "  diffuseSamplingWeight = " << (1-m_specularSamplingWeight) << "," << endl
			<< "  preserveColors = " << m_preserveColors << "," << endl
			<< "  intIOR = " << m_intIOR << "," << endl
			<< "  extIOR = " << m_extIOR << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	MicrofacetDistribution m_distribution;
	ref<RoughTransmittance> m_roughTransmittance;
	ref<Texture> m_diffuseReflectance;
	ref<Texture> m_specularReflectance;
	ref<Texture> m_alpha;
	Float m_intIOR, m_extIOR;
	Float m_specularSamplingWeight;
	bool m_preserveColors;
};

/**
 * GLSL port of the rough plastic shader. This version is much more
 * approximate -- it only supports the Beckmann distribution, 
 * does everything in RGB, uses a cheaper shadowing-masking term, and 
 * it also makes use of the Schlick approximation to the Fresnel 
 * reflectance of dielectrics. When the roughness is lower than 
 * \alpha < 0.2, the shader clamps it to 0.2 so that it will still perform
 * reasonably well in a VPL-based preview. There is no support for
 * non-linear effects due to internal scattering.
 */
class RoughPlasticShader : public Shader {
public:
	RoughPlasticShader(Renderer *renderer, const Texture *specularReflectance,
			const Texture *diffuseReflectance, const Texture *alpha, Float extIOR, 
			Float intIOR) : Shader(renderer, EBSDFShader), 
			m_specularReflectance(specularReflectance), 
			m_diffuseReflectance(diffuseReflectance), 
			m_alpha(alpha), m_extIOR(extIOR), m_intIOR(intIOR) {
		m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());
		m_diffuseReflectanceShader = renderer->registerShaderForResource(m_diffuseReflectance.get());
		m_alphaShader = renderer->registerShaderForResource(m_alpha.get());
		m_R0 = fresnel(1.0f, m_extIOR, m_intIOR);
	}

	bool isComplete() const {
		return m_specularReflectanceShader.get() != NULL &&
			m_diffuseReflectanceShader.get() != NULL &&
			m_alphaShader.get() != NULL;
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_specularReflectanceShader.get());
		deps.push_back(m_diffuseReflectanceShader.get());
		deps.push_back(m_alphaShader.get());
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_specularReflectance.get());
		renderer->unregisterShaderForResource(m_diffuseReflectance.get());
		renderer->unregisterShaderForResource(m_alpha.get());
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_R0", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[1], m_R0);
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform float " << evalName << "_R0;" << endl
			<< endl
			<< "float " << evalName << "_D(vec3 m, float alpha) {" << endl
			<< "    float ct = cosTheta(m);" << endl
			<< "    if (cosTheta(m) <= 0.0)" << endl
			<< "        return 0.0;" << endl
			<< "    float ex = tanTheta(m) / alpha;" << endl
			<< "    return exp(-(ex*ex)) / (pi * alpha * alpha *" << endl
			<< "               pow(cosTheta(m), 4.0));" << endl
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
			<< endl
			<< "float " << evalName << "_schlick(float ct) {" << endl
			<< "    float ctSqr = ct*ct, ct5 = ctSqr*ctSqr*ct;" << endl
			<< "    return " << evalName << "_R0 + (1.0 - " << evalName << "_R0) * ct5;" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) <= 0 || cosTheta(wo) <= 0)" << endl
			<< "        return vec3(0.0);" << endl
			<< "    vec3 H = normalize(wi + wo);" << endl
			<< "    vec3 specRef = " << depNames[0] << "(uv);" << endl
			<< "    vec3 diffuseRef = " << depNames[1] << "(uv);" << endl
			<< "    float alpha = " << depNames[2] << "(uv)[0];" << endl
			<< "    float D = " << evalName << "_D(H, alpha)" << ";" << endl
			<< "    float G = " << evalName << "_G(H, wi, wo);" << endl
			<< "    float F = " << evalName << "_schlick(1-dot(wi, H));" << endl
			<< "    return specRef    * (F * D * G / (4*cosTheta(wi))) + " << endl
			<< "           diffuseRef * ((1-F) * cosTheta(wo) * 0.31831);" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    vec3 diffuseRef = " << depNames[1] << "(uv);" << endl
			<< "    return diffuseRef * 0.31831 * cosTheta(wo);"<< endl
			<< "}" << endl;
	}
	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_specularReflectance;
	ref<const Texture> m_diffuseReflectance;
	ref<const Texture> m_alpha;
	ref<Shader> m_specularReflectanceShader;
	ref<Shader> m_diffuseReflectanceShader;
	ref<Shader> m_alphaShader;
	Float m_extIOR, m_intIOR, m_R0;
};

Shader *RoughPlastic::createShader(Renderer *renderer) const { 
	return new RoughPlasticShader(renderer,
		m_specularReflectance.get(), m_diffuseReflectance.get(),
		m_alpha.get(), m_extIOR, m_intIOR);
}

MTS_IMPLEMENT_CLASS(RoughPlasticShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(RoughPlastic, false, BSDF)
MTS_EXPORT_PLUGIN(RoughPlastic, "Rough plastic BRDF");
MTS_NAMESPACE_END

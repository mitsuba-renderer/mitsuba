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

#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>
#include "microfacet.h"
#include "ior.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{roughconductor}{Rough conductor material}
 * \order{7}
 * \icon{bsdf_roughconductor}
 * \parameters{
 *     \parameter{distribution}{\String}{
 *          Specifies the type of microfacet normal distribution
 *          used to model the surface roughness.
 *       \begin{enumerate}[(i)]
 *           \item \code{beckmann}: Physically-based distribution derived from
 *               Gaussian random surfaces. This is the default.\vspace{-1mm}
 *           \item \code{ggx}: New distribution proposed by
 *              Walter et al. \cite{Walter07Microfacet}, which is meant to better handle
 *              the long tails observed in measurements of ground surfaces.
 *              Renderings with this distribution may converge slowly.\vspace{-1mm}
 *           \item \code{phong}: Classical $\cos^p\theta$ distribution.
 *              Due to the underlying microfacet theory,
 *              the use of this distribution here leads to more realistic
 *              behavior than the separately available \pluginref{phong} plugin.\vspace{-1mm}
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
 *     \parameter{material}{\String}{Name of a material preset, see
 *           \tblref{conductor-iors}.\!\default{\texttt{Cu} / copper}}
 *     \parameter{eta, k}{\Spectrum}{Real and imaginary components of the material's index of
 *             refraction \default{based on the value of \texttt{material}}}
 *     \parameter{extEta}{\Float\Or\String}{
 *           Real-valued index of refraction of the surrounding dielectric,
 *           or a material name of a dielectric \default{\code{air}}
 *     }
 *     \parameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor that can be used to modulate the specular reflection component. Note
 *         that for physical realism, this parameter should never be touched. \default{1.0}}
 * }
 * \vspace{4mm}
 * This plugin implements a realistic microfacet scattering model for rendering
 * rough conducting materials, such as metals. It can be interpreted as a fancy
 * version of the Cook-Torrance model and should be preferred over
 * heuristic models like \pluginref{phong} and \pluginref{ward} when possible.
 * \renderings{
 *     \rendering{Rough copper (Beckmann, $\alpha=0.1$)}
 *     	   {bsdf_roughconductor_copper.jpg}
 *     \rendering{Vertically brushed aluminium (Ashikhmin-Shirley,
 *         $\alpha_u=0.05,\ \alpha_v=0.3$), see
 *         \lstref{roughconductor-aluminium}}
 *         {bsdf_roughconductor_anisotropic_aluminium.jpg}
 * }
 *
 * Microfacet theory describes rough
 * surfaces as an arrangement of unresolved and ideally specular facets, whose
 * normal directions are given by a specially chosen \emph{microfacet distribution}.
 * By accounting for shadowing and masking effects between these facets, it is
 * possible to reproduce the important off-specular reflections peaks observed
 * in real-world measurements of such materials.
 *
 * This plugin is essentially the ``roughened'' equivalent of the (smooth) plugin
 * \pluginref{conductor}. For very low values of $\alpha$, the two will
 * be identical, though scenes using this plugin will take longer to render
 * due to the additional computational burden of tracking surface roughness.
 *
 * The implementation is based on the paper ``Microfacet Models
 * for Refraction through Rough Surfaces'' by Walter et al.
 * \cite{Walter07Microfacet}. It supports several different types of microfacet
 * distributions and has a texturable roughness parameter.
 * To facilitate the tedious task of specifying spectrally-varying index of
 * refraction information, this plugin can access a set of measured materials
 * for which visible-spectrum information was publicly available
 * (see \tblref{conductor-iors} for the full list).
 * There is also a special material profile named \code{none}, which disables
 * the computation of Fresnel reflectances and produces an idealized
 * 100% reflecting mirror.
 *
 * When no parameters are given, the plugin activates the default settings,
 * which describe copper with a light amount of roughness modeled using a
 * Beckmann distribution.
 *
 * To get an intuition about the effect of the surface roughness
 * parameter $\alpha$, consider the following approximate classification:
 * a value of $\alpha=0.001-0.01$ corresponds to a material
 * with slight imperfections on an
 * otherwise smooth surface finish, $\alpha=0.1$ is relatively rough,
 * and $\alpha=0.3-0.7$ is \emph{extremely} rough (e.g. an etched or ground
 * finish). Values significantly above that are probably not too realistic.
 * \vspace{4mm}
 * \begin{xml}[caption={A material definition for brushed aluminium}, label=lst:roughconductor-aluminium]
 * <bsdf type="roughconductor">
 *     <string name="material" value="Al"/>
 *     <string name="distribution" value="as"/>
 *     <float name="alphaU" value="0.05"/>
 *     <float name="alphaV" value="0.3"/>
 * </bsdf>
 * \end{xml}
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
 * are missing texture coordinates.
 *
 * When using this plugin, you should ideally compile Mitsuba with support for
 * spectral rendering to get the most accurate results. While it also works
 * in RGB mode, the computations will be more approximate in nature.
 * Also note that this material is one-sided---that is, observed from the
 * back side, it will be completely black. If this is undesirable,
 * consider using the \pluginref{twosided} BRDF adapter.
 */
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

	RoughConductor(Stream *stream, InstanceManager *manager)
	 : BSDF(stream, manager) {
		m_distribution = MicrofacetDistribution(
			(MicrofacetDistribution::EType) stream->readUInt()
		);
		m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaV = static_cast<Texture *>(manager->getInstance(stream));
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_eta = Spectrum(stream);
		m_k = Spectrum(stream);

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

	/// Helper function: reflect \c wi with respect to a given surface normal
	inline Vector reflect(const Vector &wi, const Normal &m) const {
		return 2 * dot(wi, m) * Vector(m) - wi;
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		/* Stop if this component was not requested */
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) < 0 ||
			Frame::cosTheta(bRec.wo) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		/* Calculate the reflection half-vector */
		Vector H = normalize(bRec.wo+bRec.wi);

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
		const Spectrum F = fresnelConductorExact(dot(bRec.wi, H), m_eta, m_k);

		/* Smith's shadow-masking function */
		const Float G = m_distribution.G(bRec.wi, bRec.wo, H, alphaU, alphaV);

		/* Calculate the total amount of reflection */
		Float value = D * G / (4.0f * Frame::cosTheta(bRec.wi));

		return m_specularReflectance->eval(bRec.its) * F * value;
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) < 0 ||
			Frame::cosTheta(bRec.wo) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return 0.0f;

		/* Calculate the reflection half-vector */
		Vector H = normalize(bRec.wo+bRec.wi);

		/* Evaluate the roughness */
		Float alphaU = m_distribution.transformRoughness(
					m_alphaU->eval(bRec.its).average()),
			  alphaV = m_distribution.transformRoughness(
					m_alphaV->eval(bRec.its).average());

		return m_distribution.pdf(H, alphaU, alphaV)
			/ (4 * absDot(bRec.wo, H));
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		/* Evaluate the roughness */
		Float alphaU = m_distribution.transformRoughness(
					m_alphaU->eval(bRec.its).average()),
			  alphaV = m_distribution.transformRoughness(
					m_alphaV->eval(bRec.its).average());

		/* Sample M, the microsurface normal */
		Float microfacetPDF;
		const Normal m = m_distribution.sample(sample,
			alphaU, alphaV, microfacetPDF);

		if (microfacetPDF == 0)
			return Spectrum(0.0f);

		/* Perfect specular reflection based on the microsurface normal */
		bRec.wo = reflect(bRec.wi, m);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;

		/* Side check */
		if (Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		const Spectrum F = fresnelConductorExact(dot(bRec.wi, m),
				m_eta, m_k);

		Float numerator = m_distribution.eval(m, alphaU, alphaV)
			* m_distribution.G(bRec.wi, bRec.wo, m, alphaU, alphaV)
			* dot(bRec.wi, m);

		Float denominator = microfacetPDF
			* Frame::cosTheta(bRec.wi);

		return m_specularReflectance->eval(bRec.its) * F
				* (numerator / denominator);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		/* Evaluate the roughness */
		Float alphaU = m_distribution.transformRoughness(
					m_alphaU->eval(bRec.its).average()),
			  alphaV = m_distribution.transformRoughness(
					m_alphaV->eval(bRec.its).average());

		/* Sample M, the microsurface normal */
		const Normal m = m_distribution.sample(sample,
			alphaU, alphaV, pdf);

		if (pdf == 0)
			return Spectrum(0.0f);

		/* Perfect specular reflection based on the microsurface normal */
		bRec.wo = reflect(bRec.wi, m);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;

		/* Side check */
		if (Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		const Spectrum F = fresnelConductorExact(dot(bRec.wi, m),
				m_eta, m_k);

		Float numerator = m_distribution.eval(m, alphaU, alphaV)
			* m_distribution.G(bRec.wi, bRec.wo, m, alphaU, alphaV)
			* dot(bRec.wi, m);

		Float denominator = pdf * Frame::cosTheta(bRec.wi);

		/* Jacobian of the half-direction mapping */
		pdf /= 4.0f * dot(bRec.wo, m);

		return m_specularReflectance->eval(bRec.its) * F
				* (numerator / denominator);
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

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeUInt((uint32_t) m_distribution.getType());
		manager->serialize(stream, m_alphaU.get());
		manager->serialize(stream, m_alphaV.get());
		manager->serialize(stream, m_specularReflectance.get());
		m_eta.serialize(stream);
		m_k.serialize(stream);
	}

	Float getRoughness(const Intersection &its, int component) const {
		return 0.5f * (m_alphaU->eval(its).average()
			+ m_alphaV->eval(its).average());
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "RoughConductor[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  distribution = " << m_distribution.toString() << "," << endl
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
	MicrofacetDistribution m_distribution;
	ref<Texture> m_specularReflectance;
	ref<Texture> m_alphaU, m_alphaV;
	Spectrum m_eta, m_k;
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

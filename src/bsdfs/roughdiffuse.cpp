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
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/renderer.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{diffuse}{Rough diffuse material}
 * \order{2}
 * \parameters{
 *     \parameter{reflectance}{\Spectrum\Or\Texture}{
 *         Specifies the reflectance / albedo of the material \linebreak(Default: 0.5)
 *     }
 *     \lastparameter{alpha}{\Spectrum\Or\Texture}{
 *         Specifies the roughness of the unresolved surface microgeometry. 
 *         This parameter is approximately equal to the \emph{root mean square} 
 *         (RMS) slope of the microfacets.\default{0.1}
 *     }
 * }
 *
 * \renderings{
 *     \rendering{Homogeneous reflectance, see \lstref{diffuse-uniform}}{bsdf_diffuse_plain}
 *     \rendering{Textured reflectance, see \lstref{diffuse-textured}}{bsdf_diffuse_textured}
 * }
 * 
 * This reflectance model describes scattering from a rough diffuse material,
 * such as plaster, sand, clay, or concrete.
 * The underlying theory was developed by Oren and Nayar 
 * \cite{Oren1994Generalization}, who model the microscopic surface structure as 
 * an arrangement of unresolved planar facets with different slopes, where each facet
 * is an ideal diffuse reflector. The model takes into account shadowing,
 * masking, as well as interreflections between the facets.
 *
 * Since the original publication in 1994, this approach has been shown to 
 * be a very good match for many real-world materials, in particular 
 * compared to Lambertian scattering, which does not take surface 
 * roughness into account.
 *
 * To get an intuition about the effect of the surface roughness
 * parameter $\alpha$, consider the following approximate differentiation: 
 * a value of $\alpha=0.001-0.01$ corresponds to a material 
 * with slight imperfections on an otherwise smooth surface (for such small
 * values, the model will behave almost identically to \pluginref{diffuse}), $\alpha=0.1$ 
 * is relatively rough, and $\alpha=0.3-0.5$ is \emph{extremely} rough 
 * (e.g. an etched or ground finish). 
 *
 * Note that this material is one-sided---that is, observed from the 
 * back side, it will be completely black. If this is undesirable, 
 * consider using the \pluginref{twosided} BRDF adapter plugin.
 * \vspace{4mm}
 *
 * \begin{xml}[caption={A diffuse material, whose reflectance is specified as an sRGB color}, label=lst:diffuse-uniform]
 * <bsdf type="diffuse">
 *     <srgb name="reflectance" value="#6d7185"/>
 * </bsdf>
 * \end{xml}
 */
class RoughDiffuse : public BSDF {
public:
	RoughDiffuse(const Properties &props) 
		: BSDF(props) {
		/* For better compatibility with other models, support both
		   'reflectance' and 'diffuseReflectance' as parameter names */
		m_reflectance = new ConstantSpectrumTexture(props.getSpectrum(
			props.hasProperty("reflectance") ? "reflectance" 
				: "diffuseReflectance", Spectrum(.5f)));

		m_alpha = new ConstantFloatTexture(props.getFloat("alpha", 0.1f));
		m_components.push_back(EGlossyReflection | EFrontSide);
		m_usesRayDifferentials = false;
	}

	RoughDiffuse(Stream *stream, InstanceManager *manager) 
		: BSDF(stream, manager) {
		m_reflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_alpha = static_cast<Texture *>(manager->getInstance(stream));
		m_components.push_back(EGlossyReflection | EFrontSide);
		m_usesRayDifferentials = m_reflectance->usesRayDifferentials();
	}

	virtual ~RoughDiffuse() { }

	void configure() {
		BSDF::configure();

		/* Verify the input parameter and fix them if necessary */
		m_reflectance = ensureEnergyConservation(m_reflectance, "reflectance", 1.0f);
	}

	Spectrum eval(const BSDFQueryRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EGlossyReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0 
			|| Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		/* Conversion from Beckmann-style RMS roughness to
		   Oren-Nayar-style slope-area variance. The factor
		   of 1/sqrt(2) was found to be a perfect fit up
		   to extreme roughness values (>.5), after which 
		   the match is not as good anymore */
		const Float conversionFactor = 1 / std::sqrt((Float) 2);

		Float sigma = m_alpha->getValue(bRec.its).average() 
			* conversionFactor;

		Float sigma2 = sigma*sigma;
		Float A = 10.f - (sigma2 / (2.0f * (sigma2 + 0.33f)));
		Float B = 0.45f * sigma2 / (sigma2 + 0.09f);

		return m_reflectance->getValue(bRec.its)
			* (INV_PI * Frame::cosTheta(bRec.wo));
	}

	Float pdf(const BSDFQueryRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EGlossyReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0 
			|| Frame::cosTheta(bRec.wo) <= 0)
			return 0.0f;

		return Frame::cosTheta(bRec.wo) * INV_PI;
	}

	Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const {
		if (!(bRec.typeMask & EGlossyReflection) || Frame::cosTheta(bRec.wi) <= 0) 
			return Spectrum(0.0f);

		bRec.wo = squareToHemispherePSA(sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;
		return m_reflectance->getValue(bRec.its);
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (!(bRec.typeMask & EGlossyReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);
		
		bRec.wo = squareToHemispherePSA(sample);
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;
		pdf = Frame::cosTheta(bRec.wo) * INV_PI;
		return m_reflectance->getValue(bRec.its) 
			* (INV_PI * Frame::cosTheta(bRec.wo));
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) 
				&& (name == "reflectance" || name == "diffuseReflectance")) {
			m_reflectance = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_reflectance->usesRayDifferentials();
		} else if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) 
				&& name == "alpha") {
			m_alpha = static_cast<Texture *>(child);
			m_usesRayDifferentials |= m_reflectance->usesRayDifferentials();
		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_reflectance.get());
		manager->serialize(stream, m_alpha.get());
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "RoughDiffuse[" << endl
			<< "  name = \"" << getName() << "\"," << endl
			<< "  reflectance = " << indent(m_reflectance->toString()) << "," << endl
			<< "  alpha = " << indent(m_alpha->toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_reflectance;
	ref<Texture> m_alpha;
};

// ================ Hardware shader implementation ================ 

class RoughDiffuseShader : public Shader {
public:
	RoughDiffuseShader(Renderer *renderer, const Texture *reflectance, const Texture *alpha) 
		: Shader(renderer, EBSDFShader), m_reflectance(reflectance), m_alpha(alpha) {
		m_reflectanceShader = renderer->registerShaderForResource(m_reflectance.get());
		m_alphaShader = renderer->registerShaderForResource(m_alpha.get());
	}

	bool isComplete() const {
		return m_reflectanceShader.get() != NULL &&
			m_alphaShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_reflectance.get());
		renderer->unregisterShaderForResource(m_alpha.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_reflectanceShader.get());
		deps.push_back(m_alphaShader.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return " << depNames[0] << "(uv) * 0.31831 * cosTheta(wo);" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << evalName << "(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_reflectance;
	ref<const Texture> m_alpha;
	ref<Shader> m_reflectanceShader;
	ref<Shader> m_alphaShader;
};

Shader *RoughDiffuse::createShader(Renderer *renderer) const { 
	return new RoughDiffuseShader(renderer, m_reflectance.get(), m_alpha.get());
}

MTS_IMPLEMENT_CLASS(RoughDiffuseShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(RoughDiffuse, false, BSDF)
MTS_EXPORT_PLUGIN(RoughDiffuse, "Rough diffuse BRDF")
MTS_NAMESPACE_END

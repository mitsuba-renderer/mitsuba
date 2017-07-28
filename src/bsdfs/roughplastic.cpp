/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

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
#include <mitsuba/core/warp.h>
#include "microfacet.h"
#include "rtrans.h"
#include "ior.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{roughplastic}{Rough plastic material}
 * \order{9}
 * \icon{bsdf_roughplastic}
 * \parameters{
 *     \parameter{distribution}{\String}{
 *          Specifies the type of microfacet normal distribution
 *          used to model the surface roughness.
 *          \vspace{-1mm}
 *       \begin{enumerate}[(i)]
 *           \item \code{beckmann}: Physically-based distribution derived from
 *               Gaussian random surfaces. This is the default.\vspace{-1.5mm}
 *           \item \code{ggx}: The GGX \cite{Walter07Microfacet} distribution (also known as
 *               Trowbridge-Reitz \cite{Trowbridge19975Average} distribution)
 *               was designed to better approximate the long tails observed in measurements
 *               of ground surfaces, which are not modeled by the Beckmann distribution.
 *           \vspace{-1.5mm}
 *           \item \code{phong}: Classical Phong distribution.
 *              In most cases, the \code{ggx} and \code{beckmann} distributions
 *              should be preferred, since they provide better importance sampling
 *              and accurate shadowing/masking computations.
 *              \vspace{-4mm}
 *       \end{enumerate}
 *     }
 *     \parameter{alpha}{\Float\Or\Texture}{
 *         Specifies the roughness of the unresolved surface micro-geometry.
 *         When the Beckmann distribution is used, this parameter is equal to the
 *         \emph{root mean square} (RMS) slope of the microfacets.
 *         \default{0.1}.
 *     }
 *
 *     \parameter{intIOR}{\Float\Or\String}{Interior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{polypropylene} / 1.49}}
 *     \parameter{extIOR}{\Float\Or\String}{Exterior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{air} / 1.000277}}
 *     \parameter{sampleVisible}{\Boolean}{
 *         Enables an improved importance sampling technique. Refer to
 *         pages \pageref{plg:roughconductor} and \pageref{sec:visiblenormal-sampling}
 *         for details. \default{\code{true}}
 *     }
 *     \parameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor that can be used to modulate the specular reflection component. Note
 *         that for physical realism, this parameter should never be touched. \default{1.0}}
 *     \parameter{diffuse\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor used to modulate the diffuse reflection component\default{0.5}}
 *     \parameter{nonlinear}{\Boolean}{
 *         Account for nonlinear color shifts due to internal scattering? See the
 *         \pluginref{plastic} plugin for details.\default{Don't account for them and
 *         preserve the texture colors, i.e. \code{false}}
 *     }
 * }
 *
 * \vspace{3mm}
 * This plugin implements a realistic microfacet scattering model for rendering
 * rough dielectric materials with internal scattering, such as plastic. It can
 * be interpreted as a fancy version of the Cook-Torrance model and should be
 * preferred over heuristic models like \pluginref{phong} and \pluginref{ward}
 * when possible.
 *
 * Microfacet theory describes rough surfaces as an arrangement of
 * unresolved and ideally specular facets, whose normal directions are given by
 * a specially chosen \emph{microfacet distribution}. By accounting for shadowing
 * and masking effects between these facets, it is possible to reproduce the important
 * off-specular reflections peaks observed in real-world measurements of such
 * materials.
 *
 * \renderings{
 *     \rendering{Beckmann, $\alpha=0.1$}{bsdf_roughplastic_beckmann}
 *     \rendering{GGX, $\alpha=0.3$}{bsdf_roughplastic_ggx}
 * }
 *
 * This plugin is essentially the ``roughened'' equivalent of the (smooth) plugin
 * \pluginref{plastic}. For very low values of $\alpha$, the two will
 * be identical, though scenes using this plugin will take longer to render
 * due to the additional computational burden of tracking surface roughness.
 *
 * For convenience, this model allows to specify IOR values either numerically,
 * or based on a list of known materials (see \tblref{dielectric-iors} on
 * \tblpage{dielectric-iors} for an overview).
 * When no parameters are given, the plugin activates the defaults,
 * which describe a white polypropylene plastic material with a light amount
 * of roughness modeled using the Beckmann distribution.
 *
 * Like the \pluginref{plastic} material, this model internally simulates the
 * interaction of light with a diffuse base surface coated by a thin dielectric
 * layer (where the coating layer is now \emph{rough}). This is a convenient
 * abstraction rather than a restriction. In other words, there are many
 * materials that can be rendered with this model, even if they might not
 * fit this description perfectly well.
 *
 * The simplicity of this setup makes it possible to account for interesting
 * nonlinear effects due to internal scattering, which is controlled by
 * the \texttt{nonlinear} parameter. For more details, please refer to the description
 * of this parameter given in the \pluginref{plastic} plugin section
 * on \pluginpage{plastic}.
 *
 * To get an intuition about the effect of the surface roughness parameter
 * $\alpha$, consider the following approximate classification: a value of
 * $\alpha=0.001-0.01$ corresponds to a material with slight imperfections
 * on an otherwise smooth surface finish, $\alpha=0.1$ is relatively rough,
 * and $\alpha=0.3-0.7$ is \emph{extremely} rough (e.g. an etched or ground
 * finish). Values significantly above that are probably not too realistic.
 *
 * \renderings{
 *     \medrendering{Diffuse textured rendering}{bsdf_plastic_diffuse}
 *     \medrendering{Textured rough plastic model and \code{nonlinear=false}}{bsdf_roughplastic_preserve}
 *     \medrendering{Textured rough plastic model and \code{nonlinear=true}}{bsdf_roughplastic_nopreserve}
 *     \caption{
 *        When asked to do so, this model can account for subtle nonlinear color shifts due
 *        to internal scattering processes. The above images show a textured
 *        object first rendered using \pluginref{diffuse}, then
 *        \pluginref{roughplastic} with the default parameters, and finally using
 *        \pluginref{roughplastic} and support for nonlinear color shifts.
 *     }
 * }
 * \renderings{
 *     \rendering{Wood material with smooth horizontal stripes}{bsdf_roughplastic_roughtex1}
 *     \rendering{A material with imperfections at a much smaller scale than what
 *       is modeled e.g. using a bump map.}{bsdf_roughplastic_roughtex2}\vspace{-3mm}
 *     \caption{
 *         The ability to texture the roughness parameter makes it possible
 *         to render materials with a structured finish, as well as
 *         ``smudgy'' objects.
 *     }
 * }
 * \vspace{2mm}
 * \begin{xml}[caption={A material definition for black plastic material with
 *    a spatially varying roughness.},
 *    label=lst:roughplastic-varyingalpha]
 * <bsdf type="roughplastic">
 *     <string name="distribution" value="beckmann"/>
 *     <float name="intIOR" value="1.61"/>
 *     <spectrum name="diffuseReflectance" value="0"/>
 *     <!-- Fetch roughness values from a texture and slightly reduce them -->
 *     <texture type="scale" name="alpha">
 *         <texture name="alpha" type="bitmap">
 *             <string name="filename" value="roughness.png"/>
 *         </texture>
 *         <float name="scale" value="0.6"/>
 *     </texture>
 * </bsdf>
 * \end{xml}
 *
 * \subsubsection*{Technical details}
 * The implementation of this model is partly based on the paper ``Microfacet
 * Models for Refraction through Rough Surfaces'' by Walter et al.
 * \cite{Walter07Microfacet}. Several different types of microfacet
 * distributions are supported. Note that the choices are slightly more
 * restricted here---in comparison to other rough scattering models in
 * Mitsuba, anisotropic distributions are not allowed.
 *
 * The implementation of this model makes heavy use of a \emph{rough
 * Fresnel transmittance} function, which is a generalization of the
 * usual Fresnel transmittion coefficient to microfacet surfaces. Unfortunately,
 * this function is normally prohibitively expensive, since each
 * evaluation involves a numerical integration over the sphere.
 *
 * To avoid this performance issue, Mitsuba ships with data files
 * (contained in the \code{data/microfacet} directory) containing precomputed
 * values of this function over a large range of parameter values. At runtime,
 * the relevant parts are extracted using tricubic interpolation.
 *
 * When rendering with the Phong microfacet distribution, a conversion is
 * used to turn the specified Beckmann-equivalent $\alpha$ roughness value
 * into the exponent parameter of this distribution. This is done in a way,
 * such that the same value $\alpha$ will produce a similar appearance across
 * different microfacet distributions.
 */
class RoughPlastic : public BSDF {
public:
    RoughPlastic(const Properties &props) : BSDF(props) {
        m_specularReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("specularReflectance", Spectrum(1.0f)));
        m_diffuseReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("diffuseReflectance", Spectrum(0.5f)));

        /* Specifies the internal index of refraction at the interface */
        Float intIOR = lookupIOR(props, "intIOR", "polypropylene");

        /* Specifies the external index of refraction at the interface */
        Float extIOR = lookupIOR(props, "extIOR", "air");

        if (intIOR < 0 || extIOR < 0 || intIOR == extIOR)
            Log(EError, "The interior and exterior indices of "
                "refraction must be positive and differ!");

        m_eta = intIOR / extIOR;

        m_nonlinear = props.getBoolean("nonlinear", false);

        MicrofacetDistribution distr(props);
        m_type = distr.getType();
        m_sampleVisible = distr.getSampleVisible();

        if (distr.isAnisotropic())
            Log(EError, "The 'roughplastic' plugin currently does not support "
                "anisotropic microfacet distributions!");

        m_alpha = new ConstantFloatTexture(distr.getAlpha());

        m_specularSamplingWeight = 0.0f;
    }

    RoughPlastic(Stream *stream, InstanceManager *manager)
     : BSDF(stream, manager) {
        m_type = (MicrofacetDistribution::EType) stream->readUInt();
        m_sampleVisible = stream->readBool();
        m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_diffuseReflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_alpha = static_cast<Texture *>(manager->getInstance(stream));
        m_eta = stream->readFloat();
        m_nonlinear = stream->readBool();

        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        stream->writeUInt((uint32_t) m_type);
        stream->writeBool(m_sampleVisible);
        manager->serialize(stream, m_specularReflectance.get());
        manager->serialize(stream, m_diffuseReflectance.get());
        manager->serialize(stream, m_alpha.get());
        stream->writeFloat(m_eta);
        stream->writeBool(m_nonlinear);
    }

    void configure() {
        bool constAlpha = m_alpha->isConstant();

        m_components.clear();

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

        m_invEta2 = 1.0f / (m_eta*m_eta);

        if (!m_externalRoughTransmittance.get()) {
            /* Load precomputed data used to compute the rough
               transmittance through the dielectric interface */
            m_externalRoughTransmittance = new RoughTransmittance(m_type);

            m_externalRoughTransmittance->checkEta(m_eta);
            m_externalRoughTransmittance->checkAlpha(m_alpha->getMinimum().average());
            m_externalRoughTransmittance->checkAlpha(m_alpha->getMaximum().average());

            /* Reduce the rough transmittance data to a 2D slice */
            m_internalRoughTransmittance = m_externalRoughTransmittance->clone();
            m_externalRoughTransmittance->setEta(m_eta);
            m_internalRoughTransmittance->setEta(1/m_eta);

            /* If possible, even reduce it to a 1D slice */
            if (constAlpha)
                m_externalRoughTransmittance->setAlpha(
                    m_alpha->eval(Intersection()).average());
        }

        m_usesRayDifferentials =
            m_specularReflectance->usesRayDifferentials() ||
            m_diffuseReflectance->usesRayDifferentials() ||
            m_alpha->usesRayDifferentials();

        BSDF::configure();
    }

    Spectrum getDiffuseReflectance(const Intersection &its) const {
        /* Evaluate the roughness texture */
        Float alpha = m_alpha->eval(its).average();
        Float Ftr = m_externalRoughTransmittance->evalDiffuse(alpha);

        return m_diffuseReflectance->eval(its) * Ftr;
    }

    Spectrum getSpecularReflectance(const Intersection &its) const {
        return m_specularReflectance->eval(its);
    }

    /// Helper function: reflect \c wi with respect to a given surface normal
    inline Vector reflect(const Vector &wi, const Normal &m) const {
        return 2 * dot(wi, m) * Vector(m) - wi;
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        bool hasSpecular = (bRec.typeMask & EGlossyReflection) &&
            (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse = (bRec.typeMask & EDiffuseReflection) &&
            (bRec.component == -1 || bRec.component == 1);

        if (measure != ESolidAngle ||
            Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 ||
            (!hasSpecular && !hasDiffuse))
            return Spectrum(0.0f);

        /* Construct the microfacet distribution matching the
           roughness values at the current surface position. */
        MicrofacetDistribution distr(
            m_type,
            m_alpha->eval(bRec.its).average(),
            m_sampleVisible
        );

        Spectrum result(0.0f);
        if (hasSpecular) {
            /* Calculate the reflection half-vector */
            const Vector H = normalize(bRec.wo+bRec.wi);

            /* Evaluate the microfacet normal distribution */
            const Float D = distr.eval(H);

            /* Fresnel term */
            const Float F = fresnelDielectricExt(dot(bRec.wi, H), m_eta);

            /* Smith's shadow-masking function */
            const Float G = distr.G(bRec.wi, bRec.wo, H);

            /* Calculate the specular reflection component */
            Float value = F * D * G /
                (4.0f * Frame::cosTheta(bRec.wi));

            result += m_specularReflectance->eval(bRec.its) * value;
        }

        if (hasDiffuse) {
            Spectrum diff = m_diffuseReflectance->eval(bRec.its);
            Float T12 = m_externalRoughTransmittance->eval(Frame::cosTheta(bRec.wi), distr.getAlpha());
            Float T21 = m_externalRoughTransmittance->eval(Frame::cosTheta(bRec.wo), distr.getAlpha());
            Float Fdr = 1-m_internalRoughTransmittance->evalDiffuse(distr.getAlpha());

            if (m_nonlinear)
                diff /= Spectrum(1.0f) - diff * Fdr;
            else
                diff /= 1-Fdr;

            result += diff * (INV_PI * Frame::cosTheta(bRec.wo) * T12 * T21 * m_invEta2);
        }

        return result;
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        bool hasSpecular = (bRec.typeMask & EGlossyReflection) &&
            (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse = (bRec.typeMask & EDiffuseReflection) &&
            (bRec.component == -1 || bRec.component == 1);

        if (measure != ESolidAngle ||
            Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 ||
            (!hasSpecular && !hasDiffuse))
            return 0.0f;

        /* Construct the microfacet distribution matching the
           roughness values at the current surface position. */
        MicrofacetDistribution distr(
            m_type,
            m_alpha->eval(bRec.its).average(),
            m_sampleVisible
        );

        /* Calculate the reflection half-vector */
        const Vector H = normalize(bRec.wo+bRec.wi);

        Float probDiffuse, probSpecular;
        if (hasSpecular && hasDiffuse) {
            /* Find the probability of sampling the specular component */
            probSpecular = 1-m_externalRoughTransmittance->eval(Frame::cosTheta(bRec.wi), distr.getAlpha());

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
            /* Jacobian of the half-direction mapping */
            const Float dwh_dwo = 1.0f / (4.0f * dot(bRec.wo, H));

            /* Evaluate the microfacet model sampling density function */
            const Float prob = distr.pdf(bRec.wi, H);

            result = prob * dwh_dwo * probSpecular;
        }

        if (hasDiffuse)
            result += probDiffuse * warp::squareToCosineHemispherePdf(bRec.wo);

        return result;
    }

    inline Spectrum sample(BSDFSamplingRecord &bRec, Float &_pdf, const Point2 &_sample) const {
        bool hasSpecular = (bRec.typeMask & EGlossyReflection) &&
            (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse = (bRec.typeMask & EDiffuseReflection) &&
            (bRec.component == -1 || bRec.component == 1);

        if (Frame::cosTheta(bRec.wi) <= 0 || (!hasSpecular && !hasDiffuse))
            return Spectrum(0.0f);

        bool choseSpecular = hasSpecular;
        Point2 sample(_sample);

        /* Construct the microfacet distribution matching the
           roughness values at the current surface position. */
        MicrofacetDistribution distr(
            m_type,
            m_alpha->eval(bRec.its).average(),
            m_sampleVisible
        );

        Float probSpecular;
        if (hasSpecular && hasDiffuse) {
            /* Find the probability of sampling the specular component */
            probSpecular = 1 - m_externalRoughTransmittance->eval(Frame::cosTheta(bRec.wi), distr.getAlpha());

            /* Reallocate samples */
            probSpecular = (probSpecular*m_specularSamplingWeight) /
                (probSpecular*m_specularSamplingWeight +
                (1-probSpecular) * (1-m_specularSamplingWeight));

            if (sample.y < probSpecular) {
                sample.y /= probSpecular;
            } else {
                sample.y = (sample.y - probSpecular) / (1 - probSpecular);
                choseSpecular = false;
            }
        }

        if (choseSpecular) {
            /* Perfect specular reflection based on the microfacet normal */
            Normal m = distr.sample(bRec.wi, sample);
            bRec.wo = reflect(bRec.wi, m);
            bRec.sampledComponent = 0;
            bRec.sampledType = EGlossyReflection;

            /* Side check */
            if (Frame::cosTheta(bRec.wo) <= 0)
                return Spectrum(0.0f);
        } else {
            bRec.sampledComponent = 1;
            bRec.sampledType = EDiffuseReflection;
            bRec.wo = warp::squareToCosineHemisphere(sample);
        }
        bRec.eta = 1.0f;

        /* Guard against numerical imprecisions */
        _pdf = pdf(bRec, ESolidAngle);

        if (_pdf == 0)
            return Spectrum(0.0f);
        else
            return eval(bRec, ESolidAngle) / _pdf;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        Float pdf;
        return RoughPlastic::sample(bRec, pdf, sample);
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

    Float getRoughness(const Intersection &its, int component) const {
        Assert(component == 0 || component == 1);

        if (component == 0)
            return m_alpha->eval(its).average();
        else
            return std::numeric_limits<Float>::infinity();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "RoughPlastic[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  distribution = " << MicrofacetDistribution::distributionName(m_type) << "," << endl
            << "  sampleVisible = " << m_sampleVisible << "," << endl
            << "  alpha = " << indent(m_alpha->toString()) << "," << endl
            << "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
            << "  diffuseReflectance = " << indent(m_diffuseReflectance->toString()) << "," << endl
            << "  specularSamplingWeight = " << m_specularSamplingWeight << "," << endl
            << "  diffuseSamplingWeight = " << (1-m_specularSamplingWeight) << "," << endl
            << "  eta = " << m_eta << "," << endl
            << "  nonlinear = " << m_nonlinear << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    MicrofacetDistribution::EType m_type;
    ref<RoughTransmittance> m_externalRoughTransmittance;
    ref<RoughTransmittance> m_internalRoughTransmittance;
    ref<Texture> m_diffuseReflectance;
    ref<Texture> m_specularReflectance;
    ref<Texture> m_alpha;
    Float m_eta, m_invEta2;
    Float m_specularSamplingWeight;
    bool m_nonlinear;
    bool m_sampleVisible;
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
            const Texture *diffuseReflectance, const Texture *alpha, Float eta)
        : Shader(renderer, EBSDFShader),
            m_specularReflectance(specularReflectance),
            m_diffuseReflectance(diffuseReflectance),
            m_alpha(alpha) {
        m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());
        m_diffuseReflectanceShader = renderer->registerShaderForResource(m_diffuseReflectance.get());
        m_alphaShader = renderer->registerShaderForResource(m_alpha.get());
        m_R0 = fresnelDielectricExt(1.0f, eta);
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
        program->setParameter(parameterIDs[0], m_R0);
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
            << "    float alpha = max(0.2, " << depNames[2] << "(uv)[0]);" << endl
            << "    float D = " << evalName << "_D(H, alpha)" << ";" << endl
            << "    float G = " << evalName << "_G(H, wi, wo);" << endl
            << "    float F = " << evalName << "_schlick(1-dot(wi, H));" << endl
            << "    return specRef    * (F * D * G / (4*cosTheta(wi))) + " << endl
            << "           diffuseRef * ((1-F) * cosTheta(wo) * inv_pi);" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
            << "        return vec3(0.0);" << endl
            << "    vec3 diffuseRef = " << depNames[1] << "(uv);" << endl
            << "    return diffuseRef * inv_pi * cosTheta(wo);"<< endl
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
    Float m_R0;
};

Shader *RoughPlastic::createShader(Renderer *renderer) const {
    return new RoughPlasticShader(renderer,
        m_specularReflectance.get(), m_diffuseReflectance.get(),
        m_alpha.get(), m_eta);
}

MTS_IMPLEMENT_CLASS(RoughPlasticShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(RoughPlastic, false, BSDF)
MTS_EXPORT_PLUGIN(RoughPlastic, "Rough plastic BRDF");
MTS_NAMESPACE_END

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
#include "microfacet.h"
#include "rtrans.h"
#include "ior.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{roughcoating}{Rough dielectric coating}
 * \order{11}
 * \icon{bsdf_roughcoating}
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
 *     \parameter{sampleVisible}{\Boolean}{
 *         Enables an improved importance sampling technique. Refer to
 *         pages \pageref{plg:roughconductor} and \pageref{sec:visiblenormal-sampling}
 *         for details. \default{\code{true}}
 *     }
 *     \parameter{intIOR}{\Float\Or\String}{Interior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{bk7} / 1.5046}}
 *     \parameter{extIOR}{\Float\Or\String}{Exterior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{air} / 1.000277}}
 *     \parameter{thickness}{\Float}{Denotes the thickness of the layer (to
 *      model absorption --- should be specified in inverse units of \code{sigmaA})\default{1}}
 *     \parameter{sigmaA}{\Spectrum\Or\Texture}{The absorption coefficient of the
 *      coating layer. \default{0, i.e. there is no absorption}}
 *     \parameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor that can be used to modulate the specular reflection component. Note
 *         that for physical realism, this parameter should never be touched. \default{1.0}}
 *     \parameter{\Unnamed}{\BSDF}{A nested BSDF model that should be coated.}
 * }
 * \renderings{
 *     \rendering{Rough gold coated with a \emph{smooth} varnish layer}
 *         {bsdf_roughcoating_gold_smooth}
 *     \rendering{Rough gold coated with a \emph{rough} ($\alpha\!=\!0.03$) varnish layer}
 *         {bsdf_roughcoating_gold_rough}
 * }
 *
 * This plugin implements a \emph{very} approximate\footnote{
 * The model only accounts for roughness
 * in the specular reflection and Fresnel transmittance through the interface.
 * The interior model receives incident illumination
 * that is transformed \emph{as if} the coating was smooth. While
 * that's not quite correct, it is a convenient workaround when the
 * \pluginref{coating} plugin produces specular highlights that are too sharp.}
 * model that simulates a rough dielectric coating. It is essentially the
 * roughened version of \pluginref{coating}.
 * Any BSDF in Mitsuba can be coated using this plugin and multiple coating
 * layers can even be applied in sequence, which allows designing interesting
 * custom materials. The coating layer can optionally be tinted (i.e. filled
 * with an absorbing medium), in which case this model also accounts for the
 * directionally dependent absorption within the layer.
 *
 * Note that the plugin discards illumination that undergoes internal
 * reflection within the coating. This can lead to a noticeable energy
 * loss for materials that reflect much of their energy near or below the critical
 * angle (i.e. diffuse or very rough materials).
 *
 * The implementation here is motivated by the paper
 * ``Arbitrarily Layered Micro-Facet Surfaces'' by Weidlich and
 * Wilkie \cite{Weidlich2007Arbitrarily}, though the implementation
 * works differently.
 */
class RoughCoating : public BSDF {
public:
    /// \sa refractTo()
    enum EDestination {
        EInterior = 0,
        EExterior = 1
    };

    RoughCoating(const Properties &props) : BSDF(props) {
        /* Specifies the internal index of refraction at the interface */
        Float intIOR = lookupIOR(props, "intIOR", "bk7");

        /* Specifies the external index of refraction at the interface */
        Float extIOR = lookupIOR(props, "extIOR", "air");

        if (intIOR < 0 || extIOR < 0 || intIOR == extIOR)
            Log(EError, "The interior and exterior indices of "
                "refraction must be positive and differ!");

        m_eta = intIOR / extIOR;
        m_invEta = 1 / m_eta;

        /* Specifies the absorption within the layer */
        m_sigmaA = new ConstantSpectrumTexture(
            props.getSpectrum("sigmaA", Spectrum(0.0f)));

        /* Specifies the layer's thickness using the inverse units of sigmaA */
        m_thickness = props.getFloat("thickness", 1);

        /* Specifies a multiplier for the specular reflectance component */
        m_specularReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("specularReflectance", Spectrum(1.0f)));

        MicrofacetDistribution distr(props);
        m_type = distr.getType();
        m_sampleVisible = distr.getSampleVisible();

        if (distr.isAnisotropic())
            Log(EError, "The 'roughplastic' plugin currently does not support "
                "anisotropic microfacet distributions!");

        m_alpha = new ConstantFloatTexture(distr.getAlpha());

        m_specularSamplingWeight = 0.0f;
    }

    RoughCoating(Stream *stream, InstanceManager *manager)
     : BSDF(stream, manager) {
        m_type = (MicrofacetDistribution::EType) stream->readUInt();
        m_sampleVisible = stream->readBool();
        m_nested = static_cast<BSDF *>(manager->getInstance(stream));
        m_sigmaA = static_cast<Texture *>(manager->getInstance(stream));
        m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_alpha = static_cast<Texture *>(manager->getInstance(stream));
        m_eta = stream->readFloat();
        m_thickness = stream->readFloat();
        m_invEta = 1 / m_eta;

        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        stream->writeUInt((uint32_t) m_type);
        stream->writeBool(m_sampleVisible);
        manager->serialize(stream, m_nested.get());
        manager->serialize(stream, m_sigmaA.get());
        manager->serialize(stream, m_specularReflectance.get());
        manager->serialize(stream, m_alpha.get());
        stream->writeFloat(m_eta);
        stream->writeFloat(m_thickness);
    }

    void configure() {
        unsigned int extraFlags = 0;
        if (!m_sigmaA->isConstant() || !m_alpha->isConstant())
            extraFlags |= ESpatiallyVarying;

        m_components.clear();
        for (int i=0; i<m_nested->getComponentCount(); ++i)
            m_components.push_back(m_nested->getType(i) | extraFlags);

        m_components.push_back(EGlossyReflection | EFrontSide | EBackSide
            | (m_specularReflectance->isConstant() ? 0 : ESpatiallyVarying));

        m_usesRayDifferentials = m_nested->usesRayDifferentials()
            || m_sigmaA->usesRayDifferentials()
            || m_alpha->usesRayDifferentials()
            || m_specularReflectance->usesRayDifferentials();

        /* Compute weights that further steer samples towards
           the specular or nested components */
        Float avgAbsorption = (m_sigmaA->getAverage()
             *(-2*m_thickness)).exp().average();

        m_specularSamplingWeight = 1.0f / (avgAbsorption + 1.0f);

        /* Verify the input parameters and fix them if necessary */
        m_specularReflectance = ensureEnergyConservation(
            m_specularReflectance, "specularReflectance", 1.0f);

        if (!m_roughTransmittance.get()) {
            /* Load precomputed data used to compute the rough
               transmittance through the dielectric interface */
            m_roughTransmittance = new RoughTransmittance(m_type);

            m_roughTransmittance->checkEta(m_eta);
            m_roughTransmittance->checkAlpha(m_alpha->getMinimum().average());
            m_roughTransmittance->checkAlpha(m_alpha->getMaximum().average());

            /* Reduce the rough transmittance data to a 2D slice */
            m_roughTransmittance->setEta(m_eta);

            /* If possible, even reduce it to a 1D slice */
            if (m_alpha->isConstant())
                m_roughTransmittance->setAlpha(
                    m_alpha->eval(Intersection()).average());
        }

        BSDF::configure();
    }

    /// Helper function: reflect \c wi with respect to a given surface normal
    inline Vector reflect(const Vector &wi, const Normal &m) const {
        return 2 * dot(wi, m) * Vector(m) - wi;
    }

    /// Refraction in local coordinates
    Vector refractTo(EDestination dest, const Vector &wi) const {
        Float cosThetaI = Frame::cosTheta(wi);
        Float invEta = (dest == EInterior) ? m_invEta : m_eta;

        bool entering = cosThetaI > 0.0f;

        /* Using Snell's law, calculate the squared sine of the
           angle between the normal and the transmitted ray */
        Float sinThetaTSqr = invEta*invEta * Frame::sinTheta2(wi);

        if (sinThetaTSqr >= 1.0f) {
            /* Total internal reflection */
            return Vector(0.0f);
        } else {
            Float cosThetaT = std::sqrt(1.0f - sinThetaTSqr);

            /* Retain the directionality of the vector */
            return Vector(invEta*wi.x, invEta*wi.y,
                entering ? cosThetaT : -cosThetaT);
        }
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        bool hasNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
            && (bRec.component == -1 || bRec.component < (int) m_components.size()-1);
        bool hasSpecular = (bRec.typeMask & EGlossyReflection)
            && (bRec.component == -1 || bRec.component == (int) m_components.size()-1)
            && measure == ESolidAngle;

        /* Construct the microfacet distribution matching the
           roughness values at the current surface position. */
        MicrofacetDistribution distr(
            m_type,
            m_alpha->eval(bRec.its).average(),
            m_sampleVisible
        );

        Spectrum result(0.0f);
        if (hasSpecular && Frame::cosTheta(bRec.wo) * Frame::cosTheta(bRec.wi) > 0) {
            /* Calculate the reflection half-vector */
            const Vector H = normalize(bRec.wo+bRec.wi)
                * math::signum(Frame::cosTheta(bRec.wo));

            /* Evaluate the microfacet normal distribution */
            const Float D = distr.eval(H);

            /* Fresnel term */
            const Float F = fresnelDielectricExt(absDot(bRec.wi, H), m_eta);

            /* Smith's shadow-masking function */
            const Float G = distr.G(bRec.wi, bRec.wo, H);

            /* Calculate the specular reflection component */
            Float value = F * D * G /
                (4.0f * std::abs(Frame::cosTheta(bRec.wi)));

            result += m_specularReflectance->eval(bRec.its) * value;
        }

        if (hasNested) {
            BSDFSamplingRecord bRecInt(bRec);
            bRecInt.wi = refractTo(EInterior, bRec.wi);
            bRecInt.wo = refractTo(EInterior, bRec.wo);

            Spectrum nestedResult = m_nested->eval(bRecInt, measure) *
                m_roughTransmittance->eval(std::abs(Frame::cosTheta(bRec.wi)), distr.getAlpha()) *
                m_roughTransmittance->eval(std::abs(Frame::cosTheta(bRec.wo)), distr.getAlpha());

            Spectrum sigmaA = m_sigmaA->eval(bRec.its) * m_thickness;
            if (!sigmaA.isZero())
                nestedResult *= (-sigmaA *
                    (1/std::abs(Frame::cosTheta(bRecInt.wi)) +
                     1/std::abs(Frame::cosTheta(bRecInt.wo)))).exp();

            /* Solid angle compression & irradiance conversion factors */
            if (measure == ESolidAngle)
                nestedResult *= m_invEta * m_invEta *
                    Frame::cosTheta(bRec.wo) / Frame::cosTheta(bRecInt.wo);

            result += nestedResult;
        }

        return result;
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        bool hasNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
            && (bRec.component == -1 || bRec.component < (int) m_components.size()-1);
        bool hasSpecular = (bRec.typeMask & EGlossyReflection)
            && (bRec.component == -1 || bRec.component == (int) m_components.size()-1)
            && measure == ESolidAngle;

        /* Calculate the reflection half-vector */
        const Vector H = normalize(bRec.wo+bRec.wi)
                * math::signum(Frame::cosTheta(bRec.wo));

        /* Construct the microfacet distribution matching the
           roughness values at the current surface position. */
        MicrofacetDistribution distr(
            m_type,
            m_alpha->eval(bRec.its).average(),
            m_sampleVisible
        );

        Float probNested, probSpecular;
        if (hasSpecular && hasNested) {
            /* Find the probability of sampling the specular component */
            probSpecular = 1-m_roughTransmittance->eval(
                std::abs(Frame::cosTheta(bRec.wi)), distr.getAlpha());

            /* Reallocate samples */
            probSpecular = (probSpecular*m_specularSamplingWeight) /
                (probSpecular*m_specularSamplingWeight +
                (1-probSpecular) * (1-m_specularSamplingWeight));

            probNested = 1 - probSpecular;
        } else {
            probNested = probSpecular = 1.0f;
        }

        Float result = 0.0f;
        if (hasSpecular && Frame::cosTheta(bRec.wo) * Frame::cosTheta(bRec.wi) > 0) {
            /* Jacobian of the half-direction mapping */
            const Float dwh_dwo = 1.0f / (4.0f * absDot(bRec.wo, H));

            /* Evaluate the microfacet model sampling density function */
            const Float prob = distr.pdf(bRec.wi, H);

            result = prob * dwh_dwo * probSpecular;
        }

        if (hasNested) {
            BSDFSamplingRecord bRecInt(bRec);
            bRecInt.wi = refractTo(EInterior, bRec.wi);
            bRecInt.wo = refractTo(EInterior, bRec.wo);

            Float prob = m_nested->pdf(bRecInt, measure);

            if (measure == ESolidAngle) {
                prob *= m_invEta * m_invEta * Frame::cosTheta(bRec.wo)
                      / Frame::cosTheta(bRecInt.wo);
            }

            result += prob * probNested;
        }

        return result;
    }

    inline Spectrum sample(BSDFSamplingRecord &bRec, Float &_pdf, const Point2 &_sample) const {
        bool hasNested = (bRec.typeMask & m_nested->getType() & BSDF::EAll)
            && (bRec.component == -1 || bRec.component < (int) m_components.size()-1);
        bool hasSpecular = (bRec.typeMask & EGlossyReflection)
            && (bRec.component == -1 || bRec.component == (int) m_components.size()-1);

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
        if (hasSpecular && hasNested) {
            /* Find the probability of sampling the diffuse component */
            probSpecular = 1 - m_roughTransmittance->eval(
                std::abs(Frame::cosTheta(bRec.wi)), distr.getAlpha());

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
            bRec.sampledComponent = (int) m_components.size() - 1;
            bRec.sampledType = EGlossyReflection;
            bRec.eta = 1.0f;

            /* Side check */
            if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(bRec.wi) <= 0)
                return Spectrum(0.0f);
        } else {
            Vector wiBackup = bRec.wi;
            bRec.wi = refractTo(EInterior, bRec.wi);
            Spectrum result = m_nested->sample(bRec, _pdf, sample);
            bRec.wi = wiBackup;
            if (result.isZero())
                return Spectrum(0.0f);
            bRec.wo = refractTo(EExterior, bRec.wo);
            if (bRec.wo.isZero())
                return Spectrum(0.0f);
        }

        /* Guard against numerical imprecisions */
        EMeasure measure = getMeasure(bRec.sampledType);
        _pdf = pdf(bRec, measure);

        if (_pdf == 0)
            return Spectrum(0.0f);
        else
            return eval(bRec, measure) / _pdf;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        Float pdf;
        return RoughCoating::sample(bRec, pdf, sample);
    }

    Float getRoughness(const Intersection &its, int component) const {
        return component < (int) m_components.size() - 1
            ? m_nested->getRoughness(its, component)
            : m_alpha->eval(its).average();
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
            if (m_nested != NULL)
                Log(EError, "Only a single nested BRDF can be added!");
            m_nested = static_cast<BSDF *>(child);
        } else if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
            if (name == "sigmaA")
                m_sigmaA = static_cast<Texture *>(child);
            else if (name == "alpha")
                m_alpha = static_cast<Texture *>(child);
            else
                BSDF::addChild(name, child);
        } else {
            BSDF::addChild(name, child);
        }
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "RoughCoating[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  distribution = " << MicrofacetDistribution::distributionName(m_type) << "," << endl
            << "  sampleVisible = " << m_sampleVisible << "," << endl
            << "  alpha = " << indent(m_alpha->toString()) << "," << endl
            << "  sigmaA = " << indent(m_sigmaA->toString()) << "," << endl
            << "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
            << "  specularSamplingWeight = " << m_specularSamplingWeight << "," << endl
            << "  diffuseSamplingWeight = " << (1-m_specularSamplingWeight) << "," << endl
            << "  eta = " << m_eta << "," << endl
            << "  nested = " << indent(m_nested.toString()) << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    MicrofacetDistribution::EType m_type;
    ref<RoughTransmittance> m_roughTransmittance;
    ref<Texture> m_sigmaA;
    ref<Texture> m_alpha;
    ref<Texture> m_specularReflectance;
    ref<BSDF> m_nested;
    Float m_eta, m_invEta;
    Float m_specularSamplingWeight;
    Float m_thickness;
    bool m_sampleVisible;
};

/**
 * GLSL port of the rough coating shader. This version is much more
 * approximate -- it only supports the Beckmann distribution,
 * does everything in RGB, uses a cheaper shadowing-masking term, and
 * it also makes use of the Schlick approximation to the Fresnel
 * reflectance of dielectrics. When the roughness is lower than
 * \alpha < 0.2, the shader clamps it to 0.2 so that it will still perform
 * reasonably well in a VPL-based preview.
 */
class RoughCoatingShader : public Shader {
public:
    RoughCoatingShader(Renderer *renderer, const BSDF *nested,
                const Texture *sigmaA, const Texture *alpha,
                Float eta) : Shader(renderer, EBSDFShader),
            m_nested(nested), m_sigmaA(sigmaA), m_alpha(alpha), m_eta(eta) {
        m_nestedShader = renderer->registerShaderForResource(m_nested.get());
        m_sigmaAShader = renderer->registerShaderForResource(m_sigmaA.get());
        m_alphaShader = renderer->registerShaderForResource(m_alpha.get());
        m_R0 = fresnelDielectricExt(1.0f, eta);
    }

    bool isComplete() const {
        return m_nestedShader.get() != NULL
            && m_sigmaAShader.get() != NULL
            && m_alphaShader.get() != NULL;
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(m_nestedShader.get());
        deps.push_back(m_sigmaAShader.get());
        deps.push_back(m_alphaShader.get());
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(m_nested.get());
        renderer->unregisterShaderForResource(m_sigmaA.get());
        renderer->unregisterShaderForResource(m_alpha.get());
    }

    void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
        parameterIDs.push_back(program->getParameterID(evalName + "_R0", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_eta", false));
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
        program->setParameter(parameterIDs[0], m_R0);
        program->setParameter(parameterIDs[1], m_eta);
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "uniform float " << evalName << "_R0;" << endl
            << "uniform float " << evalName << "_eta;" << endl
            << endl
            << "float " << evalName << "_schlick(float ct) {" << endl
            << "    float ctSqr = ct*ct, ct5 = ctSqr*ctSqr*ct;" << endl
            << "    return " << evalName << "_R0 + (1.0 - " << evalName << "_R0) * ct5;" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_refract(vec3 wi, out float T) {" << endl
            << "    float cosThetaI = cosTheta(wi);" << endl
            << "    bool entering = cosThetaI > 0.0;" << endl
            << "    float invEta = 1.0 / " << evalName << "_eta;" << endl
            << "    float sinThetaTSqr =  invEta * invEta * sinTheta2(wi);" << endl
            << "    if (sinThetaTSqr >= 1.0) {" << endl
            << "        T = 0.0; /* Total internal reflection */" << endl
            << "        return vec3(0.0);" << endl
            << "    } else {" << endl
            << "        float cosThetaT = sqrt(1.0 - sinThetaTSqr);" << endl
            << "        T = 1.0 - " << evalName << "_schlick(1.0 - abs(cosThetaI));" << endl
            << "        return vec3(invEta*wi.x, invEta*wi.y, entering ? cosThetaT : -cosThetaT);" << endl
            << "    }" << endl
            << "}" << endl
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
            << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    float T12, T21;" << endl
            << "    vec3 wiPrime = " << evalName << "_refract(wi, T12);" << endl
            << "    vec3 woPrime = " << evalName << "_refract(wo, T21);" << endl
            << "    vec3 nested = " << depNames[0] << "(uv, wiPrime, woPrime);" << endl
            << "    vec3 sigmaA = " << depNames[1] << "(uv);" << endl
            << "    vec3 result = nested * " << evalName << "_eta * " << evalName << "_eta" << endl
            << "                  * T12 * T21 * (cosTheta(wi)*cosTheta(wo)) /" << endl
            << "                  (cosTheta(wiPrime)*cosTheta(woPrime));" << endl
            << "    if (sigmaA != vec3(0.0))" << endl
            << "        result *= exp(-sigmaA * (1/abs(cosTheta(wiPrime)) + " << endl
            << "                                 1/abs(cosTheta(woPrime))));" << endl
            << "    if (cosTheta(wi)*cosTheta(wo) > 0) {" << endl
            << "        vec3 H = normalize(wi + wo);" << endl
            << "        float alpha = max(0.2, " << depNames[2] << "(uv)[0]);" << endl
            << "        float D = " << evalName << "_D(H, alpha)" << ";" << endl
            << "        float G = " << evalName << "_G(H, wi, wo);" << endl
            << "        float F = " << evalName << "_schlick(1-dot(wi, H));" << endl
            << "        result += vec3(F * D * G / (4*cosTheta(wi)));" << endl
            << "    }" << endl
            << "    return result;" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    return " << depNames[0] << "_diffuse(uv, wi, wo);" << endl
            << "}" << endl;
    }

    MTS_DECLARE_CLASS()
private:
    ref<const BSDF> m_nested;
    ref<Shader> m_nestedShader;
    ref<const Texture> m_sigmaA;
    ref<Shader> m_sigmaAShader;
    ref<const Texture> m_alpha;
    ref<Shader> m_alphaShader;
    Float m_R0, m_eta;
};

Shader *RoughCoating::createShader(Renderer *renderer) const {
    return new RoughCoatingShader(renderer, m_nested.get(),
        m_sigmaA.get(), m_alpha.get(), m_eta);
}

MTS_IMPLEMENT_CLASS(RoughCoatingShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(RoughCoating, false, BSDF)
MTS_EXPORT_PLUGIN(RoughCoating, "Rough coating BSDF");
MTS_NAMESPACE_END

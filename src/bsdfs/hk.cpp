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
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/hw/basicshader.h>
#include "../medium/materials.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{hk}{Hanrahan-Krueger BSDF}
 * \icon{bsdf_hk}
 *
 * \parameters{
 *     \parameter{material}{\String}{Name of a material preset, see
 *           \tblref{medium-coefficients}. \default{\texttt{skin1}}}
 *     \parameter{sigmaS}{\Spectrum\Or\Texture}{Specifies the scattering coefficient
 *      of the internal layer. \default{based on \code{material}}}
 *     \parameter{sigmaA}{\Spectrum\Or\Texture}{Specifies the absorption coefficient
 *      of the internal layer. \default{based on \code{material}}}
 *     \parameter{sigmaT \& albedo}{\Spectrum\Or\Texture}{
 *      Optional: Alternatively, the scattering and absorption coefficients may also be
 *      specified using the extinction coefficient \code{sigmaT} and the
 *      single-scattering albedo. Note that only one of the parameter passing
 *      conventions can be used at a time (i.e. use either \code{sigmaS\&sigmaA}
 *      \emph{or} \code{sigmaT\&albedo})}
 *     \parameter{thickness}{\Float}{Denotes the thickness of the layer.
 *      (should be specified in inverse units of \code{sigmaA} and \code{sigmaS})\default{1}}
 *     \parameter{\Unnamed}{\Phase}{A nested phase function instance that represents
 *      the type of scattering interactions occurring within the layer}
 * }
 *
 * \renderings{
 *     \rendering{An index-matched scattering layer with parameters $\sigma_s=2$, $\sigma_a=0.1$, thickness$=0.1$}{bsdf_hk_1}
 *     \rendering{Example of the HK model with a dielectric coating (and the \code{ketchup} material preset, see \lstref{hk-coated})}{bsdf_hk_2}
 *     \caption{
 *     \label{fig:hk-example}
 *     Renderings using the uncoated and coated form of the Hanrahan-Krueger model.
 *     }
 *     \vspace{3mm}
 * }
 *
 * This plugin provides an implementation of the Hanrahan-Krueger BSDF
 * \cite{Hanrahan1993Reflection} for simulating single scattering in thin
 * index-matched layers filled with a random scattering medium.
 * In addition, the implementation also accounts for attenuated
 * light that passes through the medium without undergoing any scattering events.
 *
 * This BSDF requires a phase function to model scattering interactions within the
 * random medium. When no phase function is explicitly specified, it uses an
 * isotropic one ($g=0$) by default. A sample usage for instantiating the
 * plugin is given on the next page:\newpage
 * \begin{xml}
 * <bsdf type="hk">
 *     <spectrum name="sigmaS" value="2"/>
 *     <spectrum name="sigmaA" value="0.1"/>
 *     <float name="thickness" value="0.1"/>
 *
 *     <phase type="hg">
 *         <float name="g" value="0.8"/>
 *     </phase>
 * </bsdf>
 * \end{xml}
 *
 * When used in conjuction with the \pluginref{coating} plugin, it is possible
 * to model refraction and reflection at the layer boundaries when the indices
 * of refraction are mismatched. The combination of these two plugins then
 * reproduces the full model as it was originally proposed by Hanrahan and
 * Krueger \cite{Hanrahan1993Reflection}.
 *
 * Note that this model does not account for light that undergoes multiple
 * scattering events within the layer. This leads to energy loss,
 * particularly at grazing angles, which can be seen in the left-hand image of
 * \figref{hk-example}.
 *
 * \begin{xml}[caption=A thin dielectric layer with measured ketchup scattering parameters, label=lst:hk-coated]
 * <bsdf type="coating">
 *     <float name="extIOR" value="1.0"/>
 *     <float name="intIOR" value="1.5"/>
 *
 *     <bsdf type="hk">
 *         <string name="material" value="ketchup"/>
 *         <float name="thickness" value="0.01"/>
 *     </bsdf>
 * </bsdf>
 * \end{xml}
 *
 * Note that when \texttt{sigmaS} = \texttt{sigmaA}$\ = 0$, or when \texttt{thickness=0},
 * any geometry associated with this BSDF becomes invisible, as light will pass through
 * unchanged.
 *
 * The implementation in Mitsuba is based on code by Tom Kazimiers and Marios Papas.
 * Marios Papas has kindly verified the implementation of the coated and uncoated variants
 * against both a path tracer and a separate reference implementation.
 */
class HanrahanKrueger : public BSDF {
public:
    HanrahanKrueger(const Properties &props) : BSDF(props) {
        Spectrum sigmaS, sigmaA, g;
        lookupMaterial(props, sigmaS, sigmaA, g, NULL);
        sigmaS *= Spectrum(1.0f) - g;

        /* Scattering coefficient of the layer */
        m_sigmaS = new ConstantSpectrumTexture(
            props.getSpectrum("sigmaS", sigmaS));

        /* Absorption coefficient of the layer */
        m_sigmaA = new ConstantSpectrumTexture(
            props.getSpectrum("sigmaA", sigmaA));

        /* Slab thickness in inverse units of sigmaS and sigmaA */
        m_thickness = props.getFloat("thickness", 1);

        if (props.hasProperty("sigmaT"))
            m_sigmaT = new ConstantSpectrumTexture(
                props.getSpectrum("sigmaT"));
        if (props.hasProperty("albedo"))
            m_albedo = new ConstantSpectrumTexture(
                props.getSpectrum("albedo"));
    }

    HanrahanKrueger(Stream *stream, InstanceManager *manager)
     : BSDF(stream, manager) {
        m_phase = static_cast<PhaseFunction *>(manager->getInstance(stream));
        m_sigmaS = static_cast<Texture *>(manager->getInstance(stream));
        m_sigmaA = static_cast<Texture *>(manager->getInstance(stream));
        m_thickness = stream->readFloat();
        configure();
    }

    void configure() {
        if (m_phase == NULL)
            m_phase = static_cast<PhaseFunction *> (PluginManager::getInstance()->
                    createObject(MTS_CLASS(PhaseFunction), Properties("isotropic")));

        if (m_sigmaT != NULL || m_albedo != NULL) {
            /* Support for the alternative scattering/absorption
             * coefficient parameter passing convention */
            if (m_sigmaT == NULL || m_albedo == NULL)
                SLog(EError, "Please provide *both* sigmaT & albedo!");

            m_sigmaS = new SpectrumProductTexture(m_sigmaT, m_albedo);
            m_sigmaA = new SpectrumSubtractionTexture(m_sigmaT, m_sigmaS);
            m_sigmaT = NULL;
            m_albedo = NULL;
        }

        int extraFlags = m_sigmaS->isConstant() && m_sigmaA->isConstant() ? 0 : ESpatiallyVarying;
        m_components.clear();
        m_components.push_back(EGlossyReflection   | EFrontSide | EBackSide | EUsesSampler | extraFlags);

        if (m_thickness != std::numeric_limits<Float>::infinity()) {
            m_components.push_back(EGlossyTransmission | EFrontSide | EBackSide | EUsesSampler | extraFlags);
            m_components.push_back(EDeltaTransmission  | EFrontSide | EBackSide | EUsesSampler | extraFlags);
        }

        m_usesRayDifferentials = m_sigmaS->usesRayDifferentials()
            || m_sigmaA->usesRayDifferentials();

        BSDF::configure();
    }

    Spectrum getDiffuseReflectance(const Intersection &its) const {
        Spectrum sigmaA = m_sigmaA->eval(its),
                 sigmaS = m_sigmaS->eval(its),
                 sigmaT = sigmaA + sigmaS,
                 albedo;
        for (int i = 0; i < SPECTRUM_SAMPLES; i++)
            albedo[i] = sigmaT[i] > 0 ? (sigmaS[i]/sigmaT[i]) : (Float) 0;
        return albedo; /* Very approximate .. */
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        Spectrum sigmaA = m_sigmaA->eval(bRec.its),
                 sigmaS = m_sigmaS->eval(bRec.its),
                 sigmaT = sigmaA + sigmaS,
                 tauD = sigmaT * m_thickness,
                 result(0.0f);

        if (measure == EDiscrete) {
            /* Figure out if the specular transmission is specifically requested */
            bool hasSpecularTransmission = (bRec.typeMask & EDeltaTransmission)
                && (bRec.component == -1 || bRec.component == 2);

            /* Return the attenuated light if requested */
            if (hasSpecularTransmission &&
                std::abs(1+dot(bRec.wi, bRec.wo)) < DeltaEpsilon)
                result = (-tauD/std::abs(Frame::cosTheta(bRec.wi))).exp();
        } else if (measure == ESolidAngle) {
            /* Sample single scattering events */
            bool hasGlossyReflection = (bRec.typeMask & EGlossyReflection)
                && (bRec.component == -1 || bRec.component == 0);
            bool hasGlossyTransmission = (bRec.typeMask & EGlossyTransmission)
                && (bRec.component == -1 || bRec.component == 1);

            Spectrum albedo;
            for (int i = 0; i < SPECTRUM_SAMPLES; i++)
                albedo[i] = sigmaT[i] > 0 ? (sigmaS[i]/sigmaT[i]) : (Float) 0;

            const Float cosThetaI = Frame::cosTheta(bRec.wi),
                        cosThetaO = Frame::cosTheta(bRec.wo),
                        dp = cosThetaI*cosThetaO;

            bool reflection = dp > 0, transmission = dp < 0;

            /* ==================================================================== */
            /*                        Reflection component                          */
            /* ==================================================================== */

            if (hasGlossyReflection && reflection) {
                MediumSamplingRecord dummy;
                PhaseFunctionSamplingRecord pRec(dummy,bRec.wi,bRec.wo);
                const Float phaseVal = m_phase->eval(pRec);

                result = albedo * (phaseVal*cosThetaI/(cosThetaI+cosThetaO)) *
                    (Spectrum(1.0f)-((-1.0f/std::abs(cosThetaI)-1.0f/std::abs(cosThetaO)) * tauD).exp());
            }

            /* ==================================================================== */
            /*                       Transmission component                         */
            /* ==================================================================== */

            if (hasGlossyTransmission && transmission
                    && m_thickness < std::numeric_limits<Float>::infinity()) {
                MediumSamplingRecord dummy;
                PhaseFunctionSamplingRecord pRec(dummy,bRec.wi,bRec.wo);
                const Float phaseVal = m_phase->eval(pRec);

                /* Hanrahan etal 93 Single Scattering transmission term */
                if (std::abs(cosThetaI + cosThetaO) < Epsilon) {
                    /* avoid division by zero */
                    result += albedo * phaseVal*tauD/std::abs(cosThetaO) *
                                ((-tauD/std::abs(cosThetaO)).exp());
                } else {
                    /* Guaranteed to be positive even if |cosThetaO| > |cosThetaI| */
                    result += albedo * phaseVal*std::abs(cosThetaI)/(std::abs(cosThetaI)-std::abs(cosThetaO)) *
                        ((-tauD/std::abs(cosThetaI)).exp() - (-tauD/std::abs(cosThetaO)).exp());
                }
            }
            return result * std::abs(cosThetaO);
        }
        return result;
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        bool hasSingleScattering = (bRec.typeMask & EGlossy)
            && (bRec.component == -1 || bRec.component == 0 || bRec.component == 1);
        bool hasSpecularTransmission = (bRec.typeMask & EDeltaTransmission)
            && (bRec.component == -1 || bRec.component == 2);

        const Spectrum sigmaA = m_sigmaA->eval(bRec.its),
                 sigmaS = m_sigmaS->eval(bRec.its),
                 sigmaT = sigmaA + sigmaS,
                 tauD = sigmaT * m_thickness;

        Float probSpecularTransmission = (-tauD/std::abs(Frame::cosTheta(bRec.wi))).exp().average();

        if (measure == EDiscrete) {
            bool hasSpecularTransmission = (bRec.typeMask & EDeltaTransmission)
                && (bRec.component == -1 || bRec.component == 2);
            /* Return the attenuated light if requested */
            if (hasSpecularTransmission &&
                std::abs(1+dot(bRec.wi, bRec.wo)) < DeltaEpsilon)
                return hasSingleScattering ? probSpecularTransmission : 1.0f;
        } else if (hasSingleScattering && measure == ESolidAngle) {
            bool hasGlossyReflection = (bRec.typeMask & EGlossyReflection)
                && (bRec.component == -1 || bRec.component == 0);
            bool hasGlossyTransmission = (bRec.typeMask & EGlossyTransmission)
                && (bRec.component == -1 || bRec.component == 1);
            bool reflection = Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0;

            if ((!hasGlossyReflection && reflection) ||
                (!hasGlossyTransmission && !reflection))
                return 0.0f;

            /* Sampled according to the phase function lobe(s) */
            MediumSamplingRecord dummy;
            PhaseFunctionSamplingRecord pRec(dummy, bRec.wi, bRec.wo);
            Float pdf = m_phase->pdf(pRec);
            if (hasSpecularTransmission)
                pdf *= 1-probSpecularTransmission;
            return pdf;
        }
        return 0.0f;
    }

    inline Spectrum sample(BSDFSamplingRecord &bRec, Float &_pdf, const Point2 &_sample) const {
        AssertEx(bRec.sampler != NULL, "The BSDFSamplingRecord needs to have a sampler!");

        bool hasSpecularTransmission = (bRec.typeMask & EDeltaTransmission)
            && (bRec.component == -1 || bRec.component == 2);
        bool hasSingleScattering = (bRec.typeMask & EGlossy)
            && (bRec.component == -1 || bRec.component == 0 || bRec.component == 1);

        const Spectrum sigmaA = m_sigmaA->eval(bRec.its),
                 sigmaS = m_sigmaS->eval(bRec.its),
                 sigmaT = sigmaA + sigmaS,
                 tauD = sigmaT * m_thickness;

        /* Probability for a specular transmission is approximated by the average (per wavelength)
         * probability of a photon exiting without a scattering event or an absorption event */
        Float probSpecularTransmission = (-tauD/std::abs(Frame::cosTheta(bRec.wi))).exp().average();

        bool choseSpecularTransmission = hasSpecularTransmission;

        Point2 sample(_sample);
        if (hasSpecularTransmission && hasSingleScattering) {
            if (sample.x > probSpecularTransmission) {
                sample.x = (sample.x - probSpecularTransmission) / (1 - probSpecularTransmission);
                choseSpecularTransmission = false;
            }
        }

        bRec.eta = 1.0f;
        if (choseSpecularTransmission) {
            /* The specular transmission component was sampled */
            bRec.sampledComponent = 2;
            bRec.sampledType = EDeltaTransmission;

            bRec.wo = -bRec.wi;

            _pdf = hasSingleScattering ? probSpecularTransmission : 1.0f;
            return eval(bRec, EDiscrete) / _pdf;
        } else {
            /* The glossy transmission/scattering component should be sampled */
            bool hasGlossyReflection = (bRec.typeMask & EGlossyReflection)
                && (bRec.component == -1 || bRec.component == 0);
            bool hasGlossyTransmission = (bRec.typeMask & EGlossyTransmission)
                && (bRec.component == -1 || bRec.component == 1);

            /* Sample According to the phase function lobes */
            PhaseFunctionSamplingRecord pRec(MediumSamplingRecord(), bRec.wi, bRec.wo);
            m_phase->sample(pRec, _pdf, bRec.sampler);

            /* Store the sampled direction */
            bRec.wo = pRec.wo;

            bool reflection = Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0;
            if ((!hasGlossyReflection && reflection) ||
                (!hasGlossyTransmission && !reflection))
                return Spectrum(0.0f);

            /* Notify that the scattering component was sampled */
            bRec.sampledComponent = reflection ? 0 : 1;
            bRec.sampledType = EGlossy;

            _pdf *= (hasSpecularTransmission ? (1 - probSpecularTransmission) : 1.0f);

            /* Guard against numerical imprecisions */
            if (_pdf == 0)
                return Spectrum(0.0f);
            else
                return eval(bRec, ESolidAngle) / _pdf;

        }
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        Float pdf;
        return HanrahanKrueger::sample(bRec, pdf, sample);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        manager->serialize(stream, m_phase.get());
        manager->serialize(stream, m_sigmaS.get());
        manager->serialize(stream, m_sigmaA.get());
        stream->writeFloat(m_thickness);
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        const Class *cClass = child->getClass();

        if (cClass->derivesFrom(MTS_CLASS(PhaseFunction))) {
            Assert(m_phase == NULL);
            m_phase = static_cast<PhaseFunction *>(child);
        } else if (cClass->derivesFrom(MTS_CLASS(Texture))) {
            if (name == "sigmaS")
                m_sigmaS = static_cast<Texture *>(child);
            else if (name == "sigmaA")
                m_sigmaA = static_cast<Texture *>(child);
            else if (name == "sigmaT")
                m_sigmaT = static_cast<Texture *>(child);
            else if (name == "albedo")
                m_albedo = static_cast<Texture *>(child);
            else
                BSDF::addChild(name, child);
        } else {
            BSDF::addChild(name, child);
        }
    }

    Float getRoughness(const Intersection &its, int component) const {
        /* For lack of a better value, treat this material as diffuse
           in Manifold Exploration */
        return std::numeric_limits<Float>::infinity();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "HanrahanKrueger[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  sigmaS = " << indent(m_sigmaS->toString()) << "," << endl
            << "  sigmaA = " << indent(m_sigmaA->toString()) << "," << endl
            << "  phase = " << indent(m_phase->toString()) << "," << endl
            << "  thickness = " << m_thickness << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    ref<PhaseFunction> m_phase;
    ref<Texture> m_sigmaS;
    ref<Texture> m_sigmaA;
    Float m_thickness;
    /* Temporary fields */
    ref<Texture> m_sigmaT;
    ref<Texture> m_albedo;
};


// ================ Hardware shader implementation ================

/**
 * This is a relatively approximate GLSL shader for the HK model.
 * It assumes that the layer is infinitely thick (i.e. there is no
 * transmission) and that the phase function is isotropic
 */
class HanrahanKruegerShader : public Shader {
public:
    HanrahanKruegerShader(Renderer *renderer, const Texture *sigmaS, const Texture *sigmaA)
        : Shader(renderer, EBSDFShader), m_sigmaS(sigmaS), m_sigmaA(sigmaA) {
        m_sigmaSShader = renderer->registerShaderForResource(m_sigmaS.get());
        m_sigmaAShader = renderer->registerShaderForResource(m_sigmaA.get());
    }

    bool isComplete() const {
        return m_sigmaSShader.get() != NULL
            && m_sigmaAShader.get() != NULL;
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(m_sigmaS.get());
        renderer->unregisterShaderForResource(m_sigmaA.get());
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(m_sigmaSShader.get());
        deps.push_back(m_sigmaAShader.get());
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    vec3 sigmaS = " << depNames[0] << "(uv);" << endl
            << "    vec3 sigmaA = " << depNames[1] << "(uv);" << endl
            << "    vec3 albedo = sigmaS/(sigmaS + sigmaA);" << endl
            << "    float cosThetaI = abs(cosTheta(wi));" << endl
            << "    float cosThetaO = abs(cosTheta(wo));" << endl
            << "    return albedo * (0.079577*cosThetaI*cosThetaO/(cosThetaI + cosThetaO));" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    vec3 sigmaS = " << depNames[0] << "(uv);" << endl
            << "    vec3 sigmaA = " << depNames[1] << "(uv);" << endl
            << "    vec3 albedo = sigmaS/(sigmaS + sigmaA);" << endl
            << "    float cosThetaO = abs(cosTheta(wo));" << endl
            << "    return albedo * 0.079577 * cosThetaO;" << endl
            << "}" << endl;
    }

    MTS_DECLARE_CLASS()
private:
    ref<const Texture> m_sigmaS;
    ref<const Texture> m_sigmaA;
    ref<Shader> m_sigmaSShader;
    ref<Shader> m_sigmaAShader;
};

Shader *HanrahanKrueger::createShader(Renderer *renderer) const {
    return new HanrahanKruegerShader(renderer, m_sigmaS.get(), m_sigmaA.get());
}

MTS_IMPLEMENT_CLASS(HanrahanKruegerShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(HanrahanKrueger, false, BSDF)
MTS_EXPORT_PLUGIN(HanrahanKrueger, "Hanrahan-Krueger BSDF");
MTS_NAMESPACE_END

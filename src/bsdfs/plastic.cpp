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
#include "ior.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{plastic}{Smooth plastic material}
 * \order{8}
 * \icon{bsdf_plastic}
 * \parameters{
 *     \parameter{intIOR}{\Float\Or\String}{Interior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{polypropylene} / 1.49}}
 *     \parameter{extIOR}{\Float\Or\String}{Exterior index of refraction specified
 *      numerically or using a known material name. \default{\texttt{air} / 1.000277}}
 *     \parameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor that can be used to modulate the specular reflection component. Note that
 *         for physical realism, this parameter should never be touched. \default{1.0}}
 *     \parameter{diffuse\showbreak Reflectance}{\Spectrum\Or\Texture}{Optional
 *         factor used to modulate the diffuse reflection component\default{0.5}}
 *     \parameter{nonlinear}{\Boolean}{
 *         Account for nonlinear color shifts due to internal scattering? See the
 *         main text for details.\default{Don't account for them and
 *         preserve the texture colors, i.e. \code{false}}
 *     }
 * }
 *
 * \renderings{
 *     \rendering{A rendering with the default parameters}{bsdf_plastic_default}
 *     \rendering{A rendering with custom parameters (\lstref{plastic-shiny})}
 *         {bsdf_plastic_shiny}
 * }
 *
 * \vspace{3mm}
 * This plugin describes a smooth plastic-like material with internal scattering.
 * It uses the Fresnel reflection and transmission coefficients to provide
 * direction-dependent specular and diffuse components.
 * Since it is simple, realistic, and fast, this model is often a better choice
 * than the \pluginref{phong}, \pluginref{ward}, and \pluginref{roughplastic}
 * plugins when rendering smooth plastic-like materials.
 *
 * For convenience, this model allows to specify IOR values either numerically,
 * or based on a list of known materials (see \tblref{dielectric-iors} for
 * an overview).
 *
 * Note that this plugin is quite similar to what one would get by applying the
 * \pluginref{coating} plugin to the \pluginref{diffuse} material. The main
 * difference is that this plugin is significantly faster, while at the same
 * time causing less variance. Furthermore, it accounts for multiple
 * interreflections inside the material (read on for details), which avoids
 * a serious energy loss problem of the aforementioned plugin
 * combination.
 * \newpage
 *
 * \begin{xml}[caption=A shiny material whose diffuse reflectance is
 *     specified using sRGB, label=lst:plastic-shiny]
 * <bsdf type="plastic">
 *     <srgb name="diffuseReflectance" value="#18455c"/>
 *     <float name="intIOR" value="1.9"/>
 * </bsdf>
 * \end{xml}
 *
 * \renderings{
 *     \medrendering{Diffuse textured rendering}{bsdf_plastic_diffuse}
 *     \medrendering{Plastic model, \code{nonlinear=false}}{bsdf_plastic_preserve}
 *     \medrendering{Plastic model, \code{nonlinear=true}}{bsdf_plastic_nopreserve}
 *     \caption{
 *        \label{fig:plastic-nonlinear}
 *        When asked to do so, this model can account for subtle nonlinear color shifts due
 *        to internal scattering processes. The above images show a textured
 *        object first rendered using \pluginref{diffuse}, then
 *        \pluginref{plastic} with the default parameters, and finally using
 *        \pluginref{plastic} and support for nonlinear color shifts.
 *     }
 * }
 *
 * \subsubsection*{Internal scattering}
 * Internally, this is model simulates the interaction of light with a diffuse
 * base surface coated by a thin dielectric layer. This is a convenient
 * abstraction rather than a restriction. In other words, there are many
 * materials that can be rendered with this model, even if they might not
 * fit this description perfectly well.
 *
 * \begin{figure}[h]
 * \setcounter{subfigure}{0}
 * \centering
 * \subfloat[At the boundary, incident illumination is partly \mbox{reflected} and refracted]
 *      {\includegraphics[width=4.9cm]{images/bsdf_plastic_intscat_1.pdf}}\hfill
 * \subfloat[The refracted portion scatters diffusely at the base layer]
 *     {\includegraphics[width=4.9cm]{images/bsdf_plastic_intscat_2.pdf}}\hfill
 * \subfloat[Some of the illumination undergoes further internal scattering events]
 *     {\includegraphics[width=4.9cm]{images/bsdf_plastic_intscat_3.pdf}}
 * \caption{
 *     \label{fig:plastic-intscat}
 *     An illustration of the scattering events that are internally
 *     handled by this plugin}
 * \end{figure}
 *
 * Given illumination that is incident upon such a material, a portion
 * of the illumination is specularly reflected at the material
 * boundary, which results in a sharp reflection in the mirror direction
 * (\subfigref{plastic-intscat}{a}).
 * The remaining illumination refracts into the material, where it
 * scatters from the diffuse base layer. (\subfigref{plastic-intscat}{b}).
 * While some of the diffusely scattered illumination is able to
 * directly refract outwards again, the remainder is reflected from the
 * interior side of the dielectric boundary and will in fact remain
 * trapped inside the material for some number of internal scattering
 * events until it is finally able to escape (\subfigref{plastic-intscat}{c}).
 *
 * Due to the mathematical simplicity of this setup, it is possible to work
 * out the correct form of the model without actually having to simulate
 * the potentially large number of internal scattering events.
 *
 * Note that due to the internal scattering, the diffuse color of the
 * material is in practice slightly different from the color of the
 * base layer on its own---in particular, the material color will tend to shift towards
 * darker colors with higher saturation. Since this can be counter-intuitive when
 * using bitmap textures, these color shifts are disabled by default. Specify
 * the parameter \code{nonlinear=true} to enable them. \figref{plastic-nonlinear}
 * illustrates the resulting change. This effect is also seen in real life,
 * for instance a piece of wood will look slightly darker after coating it
 * with a layer of varnish.
 */
class SmoothPlastic : public BSDF {
public:
    SmoothPlastic(const Properties &props) : BSDF(props) {
        /* Specifies the internal index of refraction at the interface */
        Float intIOR = lookupIOR(props, "intIOR", "polypropylene");

        /* Specifies the external index of refraction at the interface */
        Float extIOR = lookupIOR(props, "extIOR", "air");

        if (intIOR < 0 || extIOR < 0)
            Log(EError, "The interior and exterior indices of "
                "refraction must be positive!");

        m_eta = intIOR / extIOR;

        m_specularReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("specularReflectance", Spectrum(1.0f)));
        m_diffuseReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("diffuseReflectance", Spectrum(0.5f)));

        m_nonlinear = props.getBoolean("nonlinear", false);

        m_specularSamplingWeight = 0.0f;
    }

    SmoothPlastic(Stream *stream, InstanceManager *manager)
            : BSDF(stream, manager) {
        m_eta = stream->readFloat();
        m_nonlinear = stream->readBool();
        m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_diffuseReflectance = static_cast<Texture *>(manager->getInstance(stream));
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        stream->writeFloat(m_eta);
        stream->writeBool(m_nonlinear);
        manager->serialize(stream, m_specularReflectance.get());
        manager->serialize(stream, m_diffuseReflectance.get());
    }

    void configure() {
        /* Verify the input parameters and fix them if necessary */
        m_specularReflectance = ensureEnergyConservation(
            m_specularReflectance, "specularReflectance", 1.0f);
        m_diffuseReflectance = ensureEnergyConservation(
            m_diffuseReflectance, "diffuseReflectance", 1.0f);

        /* Numerically approximate the diffuse Fresnel reflectance */
        m_fdrInt = fresnelDiffuseReflectance(1/m_eta, false);
        m_fdrExt = fresnelDiffuseReflectance(m_eta, false);

        /* Compute weights that further steer samples towards
           the specular or diffuse components */
        Float dAvg = m_diffuseReflectance->getAverage().getLuminance(),
              sAvg = m_specularReflectance->getAverage().getLuminance();

        m_specularSamplingWeight = sAvg / (dAvg + sAvg);

        m_invEta2 = 1/(m_eta*m_eta);

        m_usesRayDifferentials =
            m_specularReflectance->usesRayDifferentials() ||
            m_diffuseReflectance->usesRayDifferentials();

        m_components.clear();
        m_components.push_back(EDeltaReflection | EFrontSide
            | (m_specularReflectance->isConstant() ? 0 : ESpatiallyVarying));
        m_components.push_back(EDiffuseReflection | EFrontSide
            | (m_diffuseReflectance->isConstant() ? 0 : ESpatiallyVarying));

        BSDF::configure();
    }

    Spectrum getDiffuseReflectance(const Intersection &its) const {
        return m_diffuseReflectance->eval(its) * (1-m_fdrExt);
    }

    Spectrum getSpecularReflectance(const Intersection &its) const {
        return m_specularReflectance->eval(its);
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
            if (name == "specularReflectance")
                m_specularReflectance = static_cast<Texture *>(child);
            else if (name == "diffuseReflectance")
                m_diffuseReflectance = static_cast<Texture *>(child);
            else
                BSDF::addChild(name, child);
        } else {
            BSDF::addChild(name, child);
        }
    }

    /// Reflection in local coordinates
    inline Vector reflect(const Vector &wi) const {
        return Vector(-wi.x, -wi.y, wi.z);
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        bool hasSpecular   = (bRec.typeMask & EDeltaReflection)
                && (bRec.component == -1 || bRec.component == 0)
                && measure == EDiscrete;
        bool hasDiffuse = (bRec.typeMask & EDiffuseReflection)
                && (bRec.component == -1 || bRec.component == 1)
                && measure == ESolidAngle;

        if (Frame::cosTheta(bRec.wo) <= 0 || Frame::cosTheta(bRec.wi) <= 0)
            return Spectrum(0.0f);

        Float Fi = fresnelDielectricExt(Frame::cosTheta(bRec.wi), m_eta);

        if (hasSpecular) {
            /* Check if the provided direction pair matches an ideal
               specular reflection; tolerate some roundoff errors */
            if (std::abs(dot(reflect(bRec.wi), bRec.wo)-1) < DeltaEpsilon)
                return m_specularReflectance->eval(bRec.its) * Fi;
        } else if (hasDiffuse) {
            Float Fo = fresnelDielectricExt(Frame::cosTheta(bRec.wo), m_eta);

            Spectrum diff = m_diffuseReflectance->eval(bRec.its);

            if (m_nonlinear)
                diff /= Spectrum(1.0f) - diff * m_fdrInt;
            else
                diff /= 1 - m_fdrInt;

            return diff * (warp::squareToCosineHemispherePdf(bRec.wo)
                * m_invEta2 * (1-Fi) * (1-Fo));
        }

        return Spectrum(0.0f);
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        bool hasSpecular   = (bRec.typeMask & EDeltaReflection)
                && (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse = (bRec.typeMask & EDiffuseReflection)
                && (bRec.component == -1 || bRec.component == 1);

        if (Frame::cosTheta(bRec.wo) <= 0 || Frame::cosTheta(bRec.wi) <= 0)
            return 0.0f;

        Float probSpecular = hasSpecular ? 1.0f : 0.0f;
        if (hasSpecular && hasDiffuse) {
            Float Fi = fresnelDielectricExt(Frame::cosTheta(bRec.wi), m_eta);
            probSpecular = (Fi*m_specularSamplingWeight) /
                (Fi*m_specularSamplingWeight +
                (1-Fi) * (1-m_specularSamplingWeight));
        }

        if (hasSpecular && measure == EDiscrete) {
            /* Check if the provided direction pair matches an ideal
               specular reflection; tolerate some roundoff errors */
            if (std::abs(dot(reflect(bRec.wi), bRec.wo)-1) < DeltaEpsilon)
                return probSpecular;
        } else if (hasDiffuse && measure == ESolidAngle) {
            return warp::squareToCosineHemispherePdf(bRec.wo) * (1-probSpecular);
        }

        return 0.0f;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        bool hasSpecular   = (bRec.typeMask & EDeltaReflection)
                && (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse = (bRec.typeMask & EDiffuseReflection)
                && (bRec.component == -1 || bRec.component == 1);

        if ((!hasDiffuse && !hasSpecular) || Frame::cosTheta(bRec.wi) <= 0)
            return Spectrum(0.0f);

        Float Fi = fresnelDielectricExt(Frame::cosTheta(bRec.wi), m_eta);

        bRec.eta = 1.0f;
        if (hasDiffuse && hasSpecular) {
            Float probSpecular = (Fi*m_specularSamplingWeight) /
                (Fi*m_specularSamplingWeight +
                (1-Fi) * (1-m_specularSamplingWeight));

            /* Importance sample wrt. the Fresnel reflectance */
            if (sample.x < probSpecular) {
                bRec.sampledComponent = 0;
                bRec.sampledType = EDeltaReflection;
                bRec.wo = reflect(bRec.wi);

                return m_specularReflectance->eval(bRec.its)
                    * Fi / probSpecular;
            } else {
                bRec.sampledComponent = 1;
                bRec.sampledType = EDiffuseReflection;
                bRec.wo = warp::squareToCosineHemisphere(Point2(
                    (sample.x - probSpecular) / (1 - probSpecular),
                    sample.y
                ));
                Float Fo = fresnelDielectricExt(Frame::cosTheta(bRec.wo), m_eta);

                Spectrum diff = m_diffuseReflectance->eval(bRec.its);
                if (m_nonlinear)
                    diff /= Spectrum(1.0f) - diff*m_fdrInt;
                else
                    diff /= 1 - m_fdrInt;

                return diff * (m_invEta2 * (1-Fi) * (1-Fo) / (1-probSpecular));
            }
        } else if (hasSpecular) {
            bRec.sampledComponent = 0;
            bRec.sampledType = EDeltaReflection;
            bRec.wo = reflect(bRec.wi);
            return m_specularReflectance->eval(bRec.its) * Fi;
        } else {
            bRec.sampledComponent = 1;
            bRec.sampledType = EDiffuseReflection;
            bRec.wo = warp::squareToCosineHemisphere(sample);
            Float Fo = fresnelDielectricExt(Frame::cosTheta(bRec.wo), m_eta);

            Spectrum diff = m_diffuseReflectance->eval(bRec.its);
            if (m_nonlinear)
                diff /= Spectrum(1.0f) - diff*m_fdrInt;
            else
                diff /= 1 - m_fdrInt;

            return diff * (m_invEta2 * (1-Fi) * (1-Fo));
        }
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
        bool hasSpecular   = (bRec.typeMask & EDeltaReflection)
                && (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse = (bRec.typeMask & EDiffuseReflection)
                && (bRec.component == -1 || bRec.component == 1);

        if ((!hasDiffuse && !hasSpecular) || Frame::cosTheta(bRec.wi) <= 0)
            return Spectrum(0.0f);

        Float Fi = fresnelDielectricExt(Frame::cosTheta(bRec.wi), m_eta);

        bRec.eta = 1.0f;
        if (hasDiffuse && hasSpecular) {
            Float probSpecular = (Fi*m_specularSamplingWeight) /
                (Fi*m_specularSamplingWeight +
                (1-Fi) * (1-m_specularSamplingWeight));

            /* Importance sample wrt. the Fresnel reflectance */
            if (sample.x < probSpecular) {
                bRec.sampledComponent = 0;
                bRec.sampledType = EDeltaReflection;
                bRec.wo = reflect(bRec.wi);

                pdf = probSpecular;
                return m_specularReflectance->eval(bRec.its)
                    * Fi / probSpecular;
            } else {
                bRec.sampledComponent = 1;
                bRec.sampledType = EDiffuseReflection;
                bRec.wo = warp::squareToCosineHemisphere(Point2(
                    (sample.x - probSpecular) / (1 - probSpecular),
                    sample.y
                ));
                Float Fo = fresnelDielectricExt(Frame::cosTheta(bRec.wo), m_eta);

                Spectrum diff = m_diffuseReflectance->eval(bRec.its);
                if (m_nonlinear)
                    diff /= Spectrum(1.0f) - diff*m_fdrInt;
                else
                    diff /= 1 - m_fdrInt;

                pdf = (1-probSpecular) *
                    warp::squareToCosineHemispherePdf(bRec.wo);

                return diff * (m_invEta2 * (1-Fi) * (1-Fo) / (1-probSpecular));
            }
        } else if (hasSpecular) {
            bRec.sampledComponent = 0;
            bRec.sampledType = EDeltaReflection;
            bRec.wo = reflect(bRec.wi);
            pdf = 1;
            return m_specularReflectance->eval(bRec.its) * Fi;
        } else {
            bRec.sampledComponent = 1;
            bRec.sampledType = EDiffuseReflection;
            bRec.wo = warp::squareToCosineHemisphere(sample);
            Float Fo = fresnelDielectricExt(Frame::cosTheta(bRec.wo), m_eta);

            Spectrum diff = m_diffuseReflectance->eval(bRec.its);
            if (m_nonlinear)
                diff /= Spectrum(1.0f) - diff*m_fdrInt;
            else
                diff /= 1 - m_fdrInt;

            pdf = warp::squareToCosineHemispherePdf(bRec.wo);

            return diff * (m_invEta2 * (1-Fi) * (1-Fo));
        }
    }

    Float getRoughness(const Intersection &its, int component) const {
        Assert(component == 0 || component == 1);

        if (component == 0)
            return 0.0f;
        else
            return std::numeric_limits<Float>::infinity();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "SmoothPlastic[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
            << "  diffuseReflectance = " << indent(m_diffuseReflectance->toString()) << "," << endl
            << "  specularSamplingWeight = " << m_specularSamplingWeight << "," << endl
            << "  diffuseSamplingWeight = " << (1-m_specularSamplingWeight) << "," << endl
            << "  nonlinear = " << m_nonlinear << "," << endl
            << "  eta = " << m_eta << "," << endl
            << "  fdrInt = " << m_fdrInt << "," << endl
            << "  fdrExt = " << m_fdrExt << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    Float m_fdrInt, m_fdrExt, m_eta, m_invEta2;
    ref<Texture> m_diffuseReflectance;
    ref<Texture> m_specularReflectance;
    Float m_specularSamplingWeight;
    bool m_nonlinear;
};

/**
 * Smooth plastic shader -- it is really hopeless to visualize
 * this material in the VPL renderer, so let's try to do at least
 * something that suggests the presence of a specularly-reflecting
 * dielectric material with something diffuse underneath.
 */
class SmoothPlasticShader : public Shader {
public:
    SmoothPlasticShader(Renderer *renderer, const Texture *specularReflectance,
            const Texture *diffuseReflectance, Float eta) : Shader(renderer, EBSDFShader),
            m_specularReflectance(specularReflectance),
            m_diffuseReflectance(diffuseReflectance) {
        m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());
        m_diffuseReflectanceShader = renderer->registerShaderForResource(m_diffuseReflectance.get());
        m_alpha = 0.4f;
        m_R0 = fresnelDielectricExt(1.0f, eta);
    }

    bool isComplete() const {
        return m_specularReflectanceShader.get() != NULL &&
            m_diffuseReflectanceShader.get() != NULL;
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(m_specularReflectanceShader.get());
        deps.push_back(m_diffuseReflectanceShader.get());
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(m_specularReflectance.get());
        renderer->unregisterShaderForResource(m_diffuseReflectance.get());
    }

    void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
        parameterIDs.push_back(program->getParameterID(evalName + "_alpha", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_R0", false));
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
        program->setParameter(parameterIDs[0], m_alpha);
        program->setParameter(parameterIDs[1], m_R0);
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "uniform float " << evalName << "_alpha;" << endl
            << "uniform float " << evalName << "_R0;" << endl
            << endl
            << "float " << evalName << "_D(vec3 m) {" << endl
            << "    float ct = cosTheta(m);" << endl
            << "    if (cosTheta(m) <= 0.0)" << endl
            << "        return 0.0;" << endl
            << "    float ex = tanTheta(m) / " << evalName << "_alpha;" << endl
            << "    return exp(-(ex*ex)) / (pi * " << evalName << "_alpha" << endl
            << "        * " << evalName << "_alpha * pow(cosTheta(m), 4.0));" << endl
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
            << "    float D = " << evalName << "_D(H)" << ";" << endl
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
    ref<Shader> m_specularReflectanceShader;
    ref<Shader> m_diffuseReflectanceShader;
    Float m_alpha, m_R0;
};

Shader *SmoothPlastic::createShader(Renderer *renderer) const {
    return new SmoothPlasticShader(renderer,
        m_specularReflectance.get(), m_diffuseReflectance.get(), m_eta);
}

MTS_IMPLEMENT_CLASS(SmoothPlasticShader, false, Shader)

MTS_IMPLEMENT_CLASS_S(SmoothPlastic, false, BSDF)
MTS_EXPORT_PLUGIN(SmoothPlastic, "Smooth plastic BRDF");
MTS_NAMESPACE_END

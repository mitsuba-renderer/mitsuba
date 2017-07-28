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
#include "ior.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{thindielectric}{Thin dielectric material}
 * \order{4}
 * \parameters{
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
 * }
 *
 * This plugin models a \emph{thin} dielectric material that is embedded inside another
 * dielectric---for instance, glass surrounded by air. The interior of the material
 * is assumed to be so thin that its effect on transmitted rays is negligible,
 * Hence, light exits such a material without any form of angular deflection
 * (though there is still specular reflection).
 *
 * This model should be used for things like glass windows that were modeled using only a
 * single sheet of triangles or quads. On the other hand, when the window consists of
 * proper closed geometry, \pluginref{dielectric} is the right choice. This is illustrated below:
 *
 * \begin{figure}[h]
 * \setcounter{subfigure}{0}
 * \centering
 * \hfill
 * \subfloat[The \pluginref{dielectric} plugin models a single transition from one index of refraction to another]
 *     {\includegraphics[width=4.5cm]{images/bsdf_dielectric_figure.pdf}}\hfill
 * \subfloat[The \pluginref{thindielectric} plugin models a pair of interfaces causing a transient index of refraction change]
 *      {\includegraphics[width=4.5cm]{images/bsdf_thindielectric_figure.pdf}}\hfill
 * \subfloat[Windows modeled using a single sheet of geometry are the most frequent application of this BSDF]
 *      {\fbox{\includegraphics[width=4.5cm]{images/bsdf_thindielectric_window.jpg}}}\hspace*\fill
 * \caption{
 *     \label{fig:thindielectric-diff}
 *     An illustration of the difference between the \pluginref{dielectric} and \pluginref{thindielectric} plugins}
 * \end{figure}
 *
 * The implementation correctly accounts for multiple internal reflections
 * inside the thin dielectric at \emph{no significant extra cost}, i.e. paths
 * of the type $R, TRT, TR^3T, ..$ for reflection and $TT, TR^2, TR^4T, ..$ for
 * refraction, where $T$ and $R$ denote individual reflection and refraction
 * events, respectively.
 */
class ThinDielectric : public BSDF {
public:
    ThinDielectric(const Properties &props) : BSDF(props) {
        /* Specifies the internal index of refraction at the interface */
        Float intIOR = lookupIOR(props, "intIOR", "bk7");

        /* Specifies the external index of refraction at the interface */
        Float extIOR = lookupIOR(props, "extIOR", "air");

        if (intIOR < 0 || extIOR < 0)
            Log(EError, "The interior and exterior indices of "
                "refraction must be positive!");

        m_eta = intIOR / extIOR;

        m_specularReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("specularReflectance", Spectrum(1.0f)));
        m_specularTransmittance = new ConstantSpectrumTexture(
            props.getSpectrum("specularTransmittance", Spectrum(1.0f)));
    }

    ThinDielectric(Stream *stream, InstanceManager *manager)
            : BSDF(stream, manager) {
        m_eta = stream->readFloat();
        m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_specularTransmittance = static_cast<Texture *>(manager->getInstance(stream));
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        stream->writeFloat(m_eta);
        manager->serialize(stream, m_specularReflectance.get());
        manager->serialize(stream, m_specularTransmittance.get());
    }

    void configure() {
        /* Verify the input parameters and fix them if necessary */
        m_specularReflectance = ensureEnergyConservation(
            m_specularReflectance, "specularReflectance", 1.0f);
        m_specularTransmittance = ensureEnergyConservation(
            m_specularTransmittance, "specularTransmittance", 1.0f);

        m_components.clear();
        m_components.push_back(EDeltaReflection | EFrontSide | EBackSide
            | (m_specularReflectance->isConstant() ? 0 : ESpatiallyVarying));
        m_components.push_back(ENull | EFrontSide | EBackSide
            | (m_specularTransmittance->isConstant() ? 0 : ESpatiallyVarying));

        m_usesRayDifferentials = false;

        m_usesRayDifferentials =
            m_specularReflectance->usesRayDifferentials() ||
            m_specularTransmittance->usesRayDifferentials();

        BSDF::configure();
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
            if (name == "specularReflectance")
                m_specularReflectance = static_cast<Texture *>(child);
            else if (name == "specularTransmittance")
                m_specularTransmittance = static_cast<Texture *>(child);
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

    /// Transmission in local coordinates
    inline Vector transmit(const Vector &wi) const {
        return -wi;
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
                && (bRec.component == -1 || bRec.component == 0) && measure == EDiscrete;
        bool sampleTransmission = (bRec.typeMask & ENull)
                && (bRec.component == -1 || bRec.component == 1) && measure == EDiscrete;

        Float R = fresnelDielectricExt(std::abs(Frame::cosTheta(bRec.wi)), m_eta), T = 1-R;

        // Account for internal reflections: R' = R + TRT + TR^3T + ..
        if (R < 1)
            R += T*T * R / (1-R*R);

        if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0) {
            if (!sampleReflection || std::abs(dot(reflect(bRec.wi), bRec.wo)-1) > DeltaEpsilon)
                return Spectrum(0.0f);

            return m_specularReflectance->eval(bRec.its) * R;
        } else {
            if (!sampleTransmission || std::abs(dot(transmit(bRec.wi), bRec.wo)-1) > DeltaEpsilon)
                return Spectrum(0.0f);

            return m_specularTransmittance->eval(bRec.its) * (1 - R);
        }
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
                && (bRec.component == -1 || bRec.component == 0) && measure == EDiscrete;
        bool sampleTransmission = (bRec.typeMask & ENull)
                && (bRec.component == -1 || bRec.component == 1) && measure == EDiscrete;

        Float R = fresnelDielectricExt(std::abs(Frame::cosTheta(bRec.wi)), m_eta), T = 1-R;

        // Account for internal reflections: R' = R + TRT + TR^3T + ..
        if (R < 1)
            R += T*T * R / (1-R*R);

        if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0) {
            if (!sampleReflection || std::abs(dot(reflect(bRec.wi), bRec.wo)-1) > DeltaEpsilon)
                return 0.0f;

            return sampleTransmission ? R : 1.0f;
        } else {
            if (!sampleTransmission || std::abs(dot(transmit(bRec.wi), bRec.wo)-1) > DeltaEpsilon)
                return 0.0f;

            return sampleReflection ? 1-R : 1.0f;
        }
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
        bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
                && (bRec.component == -1 || bRec.component == 0);
        bool sampleTransmission = (bRec.typeMask & ENull)
                && (bRec.component == -1 || bRec.component == 1);

        Float R = fresnelDielectricExt(std::abs(Frame::cosTheta(bRec.wi)), m_eta), T = 1-R;

        // Account for internal reflections: R' = R + TRT + TR^3T + ..
        if (R < 1)
            R += T*T * R / (1-R*R);

        if (sampleTransmission && sampleReflection) {
            if (sample.x <= R) {
                bRec.sampledComponent = 0;
                bRec.sampledType = EDeltaReflection;
                bRec.wo = reflect(bRec.wi);
                bRec.eta = 1.0f;
                pdf = R;

                return m_specularReflectance->eval(bRec.its);
            } else {
                bRec.sampledComponent = 1;
                bRec.sampledType = ENull;
                bRec.wo = transmit(bRec.wi);
                bRec.eta = 1.0f;
                pdf = 1-R;

                return m_specularTransmittance->eval(bRec.its);
            }
        } else if (sampleReflection) {
            bRec.sampledComponent = 0;
            bRec.sampledType = EDeltaReflection;
            bRec.wo = reflect(bRec.wi);
            bRec.eta = 1.0f;
            pdf = 1.0f;

            return m_specularReflectance->eval(bRec.its) * R;
        } else if (sampleTransmission) {
            bRec.sampledComponent = 1;
            bRec.sampledType = ENull;
            bRec.wo = transmit(bRec.wi);
            bRec.eta = 1.0f;
            pdf = 1.0f;

            return m_specularTransmittance->eval(bRec.its) * (1-R);
        }

        return Spectrum(0.0f);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
                && (bRec.component == -1 || bRec.component == 0);
        bool sampleTransmission = (bRec.typeMask & ENull)
                && (bRec.component == -1 || bRec.component == 1);

        Float R = fresnelDielectricExt(Frame::cosTheta(bRec.wi), m_eta), T = 1-R;

        // Account for internal reflections: R' = R + TRT + TR^3T + ..
        if (R < 1)
            R += T*T * R / (1-R*R);

        if (sampleTransmission && sampleReflection) {
            if (sample.x <= R) {
                bRec.sampledComponent = 0;
                bRec.sampledType = EDeltaReflection;
                bRec.wo = reflect(bRec.wi);
                bRec.eta = 1.0f;

                return m_specularReflectance->eval(bRec.its);
            } else {
                bRec.sampledComponent = 1;
                bRec.sampledType = ENull;
                bRec.wo = transmit(bRec.wi);
                bRec.eta = 1.0f;

                return m_specularTransmittance->eval(bRec.its);
            }
        } else if (sampleReflection) {
            bRec.sampledComponent = 0;
            bRec.sampledType = EDeltaReflection;
            bRec.wo = reflect(bRec.wi);
            bRec.eta = 1.0f;

            return m_specularReflectance->eval(bRec.its) * R;
        } else if (sampleTransmission) {
            bRec.sampledComponent = 1;
            bRec.sampledType = ENull;
            bRec.wo = transmit(bRec.wi);
            bRec.eta = 1.0f;

            return m_specularTransmittance->eval(bRec.its) * (1-R);
        }

        return Spectrum(0.0f);
    }

    Float getEta() const {
        /* The rrelative IOR across this interface is 1, since the internal
           material is thin: it begins and ends here. */
        return 1.0f;
    }

    Float getRoughness(const Intersection &its, int component) const {
        return 0.0f;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "ThinDielectric[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  eta = " << m_eta << "," << endl
            << "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
            << "  specularTransmittance = " << indent(m_specularTransmittance->toString()) << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    Float m_eta;
    ref<Texture> m_specularTransmittance;
    ref<Texture> m_specularReflectance;
};

/* Fake glass shader -- it is really hopeless to visualize
   this material in the VPL renderer, so let's try to do at least
   something that suggests the presence of a transparent boundary */
class ThinDielectricShader : public Shader {
public:
    ThinDielectricShader(Renderer *renderer) :
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
            << "        return vec3(0.0);" << endl
            << "    return vec3(inv_pi * cosTheta(wo));" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    return " << evalName << "(uv, wi, wo);" << endl
            << "}" << endl;
    }


    MTS_DECLARE_CLASS()
};

Shader *ThinDielectric::createShader(Renderer *renderer) const {
    return new ThinDielectricShader(renderer);
}

MTS_IMPLEMENT_CLASS(ThinDielectricShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(ThinDielectric, false, BSDF)
MTS_EXPORT_PLUGIN(ThinDielectric, "Thin dielectric BSDF");
MTS_NAMESPACE_END

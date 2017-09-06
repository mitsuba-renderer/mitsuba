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
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/*! \plugin{difftrans}{Diffuse transmitter}
 * \icon{bsdf_difftrans}
 *
 * \parameters{
 *     \parameter{transmittance}{\Spectrum\Or\Texture}{
 *       Specifies the diffuse transmittance of the material
 *       \default{0.5}
 *     }
 * }
 *
 * \renderings{
 *     \rendering{The model with default parameters}{bsdf_difftrans.jpg}
 * }
 *
 * This BSDF models a non-reflective material, where any entering light loses
 * its directionality and is diffusely scattered from the other side. This
 * model can be combined\footnote{For instance using the
 * \pluginref{mixturebsdf} plugin.} with a surface reflection model to
 * describe translucent substances that have internal multiple scattering
 * processes (e.g. plant leaves).
 */
class DiffuseTransmitter : public BSDF {
public:
    DiffuseTransmitter(const Properties &props)
        : BSDF(props) {
        /* For better compatibility with other models, support both
           'transmittance' and 'diffuseTransmittance' as parameter names */
        m_transmittance = new ConstantSpectrumTexture(props.getSpectrum(
            props.hasProperty("transmittance") ? "transmittance"
                : "diffuseTransmittance", Spectrum(.5f)));
        m_usesRayDifferentials = false;
    }

    DiffuseTransmitter(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager) {
        m_transmittance = static_cast<Texture *>(manager->getInstance(stream));
        m_usesRayDifferentials = m_transmittance->usesRayDifferentials();
        configure();
    }

    void configure() {
        /* Verify the input parameters and fix them if necessary */
        m_transmittance = ensureEnergyConservation(m_transmittance, "transmittance", 1.0f);

        m_components.clear();
        m_components.push_back(EDiffuseTransmission | EFrontSide | EBackSide
            | (m_transmittance->isConstant() ? 0 : ESpatiallyVarying));
        BSDF::configure();
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        if (!(bRec.typeMask & EDiffuseTransmission) || measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
            return Spectrum(0.0f);

        return m_transmittance->eval(bRec.its)
            * (INV_PI * std::abs(Frame::cosTheta(bRec.wo)));
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        if (!(bRec.typeMask & EDiffuseTransmission) || measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
            return 0.0f;

        return std::abs(Frame::cosTheta(bRec.wo)) * INV_PI;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        if (!(bRec.typeMask & EDiffuseTransmission))
            return Spectrum(0.0f);
        bRec.wo = warp::squareToCosineHemisphere(sample);
        if (Frame::cosTheta(bRec.wi) > 0)
            bRec.wo.z *= -1;
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EDiffuseTransmission;
        return m_transmittance->eval(bRec.its);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
        if (!(bRec.typeMask & m_combinedType))
            return Spectrum(0.0f);
        bRec.wo = warp::squareToCosineHemisphere(sample);
        if (Frame::cosTheta(bRec.wi) > 0)
            bRec.wo.z *= -1;
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EDiffuseTransmission;
        pdf = std::abs(Frame::cosTheta(bRec.wo)) * INV_PI;
        return m_transmittance->eval(bRec.its);
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) &&
                (name == "transmittance" || name == "diffuseTransmittance")) {
            m_transmittance = static_cast<Texture *>(child);
            m_usesRayDifferentials |= m_transmittance->usesRayDifferentials();
        } else {
            BSDF::addChild(name, child);
        }
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        manager->serialize(stream, m_transmittance.get());
    }

    Float getRoughness(const Intersection &its, int component) const {
        return std::numeric_limits<Float>::infinity();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "DiffuseTransmitter[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  transmittance = " << indent(m_transmittance->toString()) << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    ref<Texture> m_transmittance;
};

// ================ Hardware shader implementation ================

class DiffuseTransmitterShader : public Shader {
public:
    DiffuseTransmitterShader(Renderer *renderer, const Texture *reflectance)
        : Shader(renderer, EBSDFShader), m_transmittance(reflectance) {
        m_transmittanceShader = renderer->registerShaderForResource(m_transmittance.get());
    }

    bool isComplete() const {
        return m_transmittanceShader.get() != NULL;
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(m_transmittance.get());
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(m_transmittanceShader.get());
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (cosTheta(wi) * cosTheta(wo) >= 0.0)" << endl
            << "        return vec3(0.0);" << endl
            << "    return " << depNames[0] << "(uv) * inv_pi * abs(cosTheta(wo));" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    return " << evalName << "(uv, wi, wo);" << endl
            << "}" << endl;
    }

    MTS_DECLARE_CLASS()
private:
    ref<const Texture> m_transmittance;
    ref<Shader> m_transmittanceShader;
};

Shader *DiffuseTransmitter::createShader(Renderer *renderer) const {
    return new DiffuseTransmitterShader(renderer, m_transmittance.get());
}

MTS_IMPLEMENT_CLASS(DiffuseTransmitterShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(DiffuseTransmitter, false, BSDF)
MTS_EXPORT_PLUGIN(DiffuseTransmitter, "Diffuse transmitter")
MTS_NAMESPACE_END

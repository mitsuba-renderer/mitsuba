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
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

/*!\plugin{ward}{Anisotropic Ward BRDF}
 * \order{15}
 * \parameters{
 *     \parameter{variant}{\String}{
 *         Determines the variant of the Ward model to use:
 *         \begin{enumerate}[(i)]
 *             \item \code{ward}: The original model by Ward \cite{Ward1992Measuring}
 *             --- suffers from energy loss at grazing angles.
 *             \item \code{ward-duer}: Corrected Ward model with lower energy loss
 *             at grazing angles \cite{Dur2006Improved}.
 *             Does not always conserve energy.
 *             \item \code{balanced}: Improved version of the \code{ward-duer}
 *             model with energy balance at all angles \cite{Geisler2010New}.
 *         \end{enumerate}
 *         Default: \texttt{balanced}
 *     }
 *     \parameter{alphaU, alphaV}{\Float\Or\Texture}{
 *         Specifies the anisotropic roughness values along the tangent and
 *         bitangent directions.
 *         \default{0.1}.
 *     }
 *     \parameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{
 *         Specifies the weight of the specular reflectance component.\default{0.2}}
 *     \parameter{diffuse\showbreak Reflectance}{\Spectrum\Or\Texture}{
 *         Specifies the weight of the diffuse reflectance component\default{0.5}}
 * }
 * \renderings{
 *     \rendering{$\alpha_u=0.1,\ \alpha_v=0.3$}{bsdf_ward_01_03}
 *     \rendering{$\alpha_u=0.3,\ \alpha_v=0.1$}{bsdf_ward_03_01}
 * }

 * This plugin implements the anisotropic Ward reflectance model and
 * several extensions. They are described in the papers
 * \begin{enumerate}[(i)]
 *    \item ``Measuring and Modeling Anisotropic Reflection''
 *      by Greg Ward \cite{Ward1992Measuring}
 *    \item ``Notes on the Ward BRDF'' by Bruce Walter \cite{Walter2005Notes}
 *    \item ``An Improved Normalization for the Ward Reflectance Model''
 *      by Arne D\"ur \cite{Dur2006Improved}
 *    \item ``A New Ward BRDF Model with Bounded Albedo'' by
 *      Geisler-Moroder et al. \cite{Geisler2010New}
 * \end{enumerate}
 *
 * Like the Phong BRDF, the Ward model does not take the Fresnel reflectance
 * of the material into account. In an experimental study by Ngan et al.
 * \cite{Ngan2005Experimental}, the Ward model performed noticeably worse than
 * models based on microfacet theory.
 *
 * For this reason, it is usually preferable to switch to a microfacet model
 * that incorporates knowledge about the material's index of refraction. In Mitsuba,
 * two such alternatives to \pluginref{ward} are given by the plugins
 * \pluginref{roughconductor} and \pluginref{roughplastic} (depending on the
 * material type).
 *
 * When using this plugin, note that the diffuse and specular reflectance
 * components should add up to a value less than or equal to one (for each
 * color channel). Otherwise, they will automatically be scaled appropriately
 * to ensure energy conservation.
 */
class Ward : public BSDF {
public:
    /// Supported model types
    enum EModelVariant {
        /// The original Ward model
        EWard = 0,
        /// Ward model with correction by Arne Duer
        EWardDuer = 1,
        /// Energy-balanced Ward model
        EBalanced = 2
    };

    Ward(const Properties &props)
        : BSDF(props) {
        m_diffuseReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("diffuseReflectance", Spectrum(0.5f)));
        m_specularReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("specularReflectance", Spectrum(0.2f)));

        std::string type =
            boost::to_lower_copy(props.getString("variant", "balanced"));
        if (type == "ward")
            m_modelVariant = EWard;
        else if (type == "ward-duer")
            m_modelVariant = EWardDuer;
        else if (type == "balanced")
            m_modelVariant = EBalanced;
        else
            Log(EError, "Specified an invalid model type \"%s\", must be "
                "\"ward\", \"ward-duer\", or \"balanced\"!", type.c_str());

        Float alpha = props.getFloat("alpha", 0.1f),
              alphaU = props.getFloat("alphaU", alpha),
              alphaV = props.getFloat("alphaV", alpha);

        m_alphaU = new ConstantFloatTexture(alphaU);
        if (alphaU == alphaV)
            m_alphaV = m_alphaU;
        else
            m_alphaV = new ConstantFloatTexture(alphaV);
        m_specularSamplingWeight = 0.0f;
    }

    Ward(Stream *stream, InstanceManager *manager)
     : BSDF(stream, manager) {
        m_modelVariant = (EModelVariant) stream->readUInt();
        m_diffuseReflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
        m_alphaV = static_cast<Texture *>(manager->getInstance(stream));

        configure();
    }

    void configure() {
        unsigned int extraFlags = 0;
        if (m_alphaU != m_alphaV)
            extraFlags |= EAnisotropic;

        m_components.clear();
        m_components.push_back(EGlossyReflection | EFrontSide | extraFlags
            | ((!m_specularReflectance->isConstant() || !m_alphaU->isConstant()
              || !m_alphaV->isConstant()) ? ESpatiallyVarying : 0));
        m_components.push_back(EDiffuseReflection | EFrontSide | extraFlags
            | (m_diffuseReflectance->isConstant() ? 0 : ESpatiallyVarying));

        /* Verify the input parameters and fix them if necessary */
        std::pair<Texture *, Texture *> result = ensureEnergyConservation(
            m_specularReflectance, m_diffuseReflectance,
            "specularReflectance", "diffuseReflectance", 1.0f);
        m_specularReflectance = result.first;
        m_diffuseReflectance = result.second;

        /* Compute weights that steer samples towards
           the specular or diffuse components */
        Float dAvg = m_diffuseReflectance->getAverage().getLuminance(),
              sAvg = m_specularReflectance->getAverage().getLuminance();
        m_specularSamplingWeight = sAvg / (dAvg + sAvg);

        m_usesRayDifferentials =
            m_diffuseReflectance->usesRayDifferentials() ||
            m_specularReflectance->usesRayDifferentials() ||
            m_alphaU->usesRayDifferentials() ||
            m_alphaV->usesRayDifferentials();

        BSDF::configure();
    }

    Spectrum getDiffuseReflectance(const Intersection &its) const {
        return m_diffuseReflectance->eval(its);
    }


    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        if (Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 || measure != ESolidAngle)
            return Spectrum(0.0f);

        bool hasSpecular = (bRec.typeMask & EGlossyReflection)
                && (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse  = (bRec.typeMask & EDiffuseReflection)
                && (bRec.component == -1 || bRec.component == 1);

        Spectrum result(0.0f);
        if (hasSpecular) {
            Vector H = bRec.wi+bRec.wo;
            Float alphaU = m_alphaU->eval(bRec.its).average();
            Float alphaV = m_alphaV->eval(bRec.its).average();

            Float factor1 = 0.0f;
            switch (m_modelVariant) {
                case EWard:
                    factor1 = 1.0f / (4.0f * M_PI * alphaU * alphaV *
                        std::sqrt(Frame::cosTheta(bRec.wi)*Frame::cosTheta(bRec.wo)));
                    break;
                case EWardDuer:
                    factor1 = 1.0f / (4.0f * M_PI * alphaU * alphaV *
                        Frame::cosTheta(bRec.wi)*Frame::cosTheta(bRec.wo));
                    break;
                case EBalanced:
                    factor1 = dot(H,H) / (M_PI * alphaU * alphaV
                        * std::pow(Frame::cosTheta(H),4));
                    break;
                default:
                    Log(EError, "Unknown model type!");
            }

            Float factor2 = H.x / alphaU, factor3 = H.y / alphaV;
            Float exponent = -(factor2*factor2+factor3*factor3)/(H.z*H.z);
            Float specRef = factor1 * math::fastexp(exponent);
            /* Important to prevent numeric issues when evaluating the
               sampling density of the Ward model in places where it takes
               on miniscule values (Veach-MLT does this for instance) */
            if (specRef > 1e-10f)
                result += m_specularReflectance->eval(bRec.its) * specRef;
        }

        if (hasDiffuse)
            result += m_diffuseReflectance->eval(bRec.its) * INV_PI;

        return result * Frame::cosTheta(bRec.wo);
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        if (Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 || measure != ESolidAngle)
            return 0.0f;

        bool hasSpecular = (bRec.typeMask & EGlossyReflection)
                && (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse  = (bRec.typeMask & EDiffuseReflection)
                && (bRec.component == -1 || bRec.component == 1);

        Float diffuseProb = 0.0f, specProb = 0.0f;

        if (hasSpecular) {
            Float alphaU = m_alphaU->eval(bRec.its).average();
            Float alphaV = m_alphaV->eval(bRec.its).average();
            Vector H = normalize(bRec.wi+bRec.wo);
            Float factor1 = 1.0f / (4.0f * M_PI * alphaU * alphaV *
                dot(H, bRec.wi) * std::pow(Frame::cosTheta(H), 3));
            Float factor2 = H.x / alphaU, factor3 = H.y / alphaV;

            Float exponent = -(factor2*factor2+factor3*factor3)/(H.z*H.z);
            specProb = factor1 * math::fastexp(exponent);
        }

        if (hasDiffuse)
            diffuseProb = warp::squareToCosineHemispherePdf(bRec.wo);

        if (hasDiffuse && hasSpecular)
            return m_specularSamplingWeight * specProb +
                   (1-m_specularSamplingWeight) * diffuseProb;
        else if (hasDiffuse)
            return diffuseProb;
        else if (hasSpecular)
            return specProb;
        else
            return 0.0f;
    }

    inline Spectrum sample(BSDFSamplingRecord &bRec, Float &_pdf, const Point2 &_sample) const {
        Point2 sample(_sample);

        bool hasSpecular = (bRec.typeMask & EGlossyReflection)
                && (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse  = (bRec.typeMask & EDiffuseReflection)
                && (bRec.component == -1 || bRec.component == 1);

        if (!hasSpecular && !hasDiffuse)
            return Spectrum(0.0f);

        bool choseSpecular = hasSpecular;

        if (hasDiffuse && hasSpecular) {
            if (sample.x <= m_specularSamplingWeight) {
                sample.x /= m_specularSamplingWeight;
            } else {
                sample.x = (sample.x - m_specularSamplingWeight)
                    / (1-m_specularSamplingWeight);
                choseSpecular = false;
            }
        }

        if (choseSpecular) {
            Float alphaU = m_alphaU->eval(bRec.its).average();
            Float alphaV = m_alphaV->eval(bRec.its).average();

            Float phiH = std::atan(alphaV/alphaU
                * std::tan(2.0f * M_PI * sample.y));
            if (sample.y > 0.5f)
                phiH += M_PI;
            Float cosPhiH = std::cos(phiH);
            Float sinPhiH = math::safe_sqrt(1.0f-cosPhiH*cosPhiH);

            Float thetaH = std::atan(math::safe_sqrt(
                -math::fastlog(sample.x) / (
                    (cosPhiH*cosPhiH) / (alphaU*alphaU) +
                    (sinPhiH*sinPhiH) / (alphaV*alphaV)
            )));
            Vector H = sphericalDirection(thetaH, phiH);
            bRec.wo = H * (2.0f * dot(bRec.wi, H)) - bRec.wi;

            bRec.sampledComponent = 1;
            bRec.sampledType = EGlossyReflection;

            if (Frame::cosTheta(bRec.wo) <= 0.0f)
                return Spectrum(0.0f);
        } else {
            bRec.wo = warp::squareToCosineHemisphere(sample);
            bRec.sampledComponent = 0;
            bRec.sampledType = EDiffuseReflection;
        }
        bRec.eta = 1.0f;

        _pdf = pdf(bRec, ESolidAngle);

        if (_pdf == 0)
            return Spectrum(0.0f);
        else
            return eval(bRec, ESolidAngle) / _pdf;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        Float pdf;
        return Ward::sample(bRec, pdf, sample);
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
            if (name == "alphaU")
                m_alphaU = static_cast<Texture *>(child);
            else if (name == "alphaV")
                m_alphaV = static_cast<Texture *>(child);
            else if (name == "diffuseReflectance")
                m_diffuseReflectance = static_cast<Texture *>(child);
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

        stream->writeUInt(m_modelVariant);
        manager->serialize(stream, m_diffuseReflectance.get());
        manager->serialize(stream, m_specularReflectance.get());
        manager->serialize(stream, m_alphaU.get());
        manager->serialize(stream, m_alphaV.get());
    }

    Float getRoughness(const Intersection &its, int component) const {
        Assert(component == 0 || component == 1);

        if (component == 0)
            return 0.5f * (m_alphaU->eval(its).average()
                + m_alphaV->eval(its).average());
        else
            return std::numeric_limits<Float>::infinity();
    }

    Shader *createShader(Renderer *renderer) const;

    std::string toString() const {
        std::ostringstream oss;
        oss << "Ward[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  variant = ";

        switch (m_modelVariant) {
            case EWard: oss << "ward," << endl; break;
            case EWardDuer: oss << "wardDuer," << endl; break;
            case EBalanced: oss << "balanced," << endl; break;
            default: Log(EError, "Unknown model type!");
        }

        oss << "  diffuseReflectance = " << indent(m_diffuseReflectance->toString()) << "," << endl
            << "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
            << "  specularSamplingWeight = " << m_specularSamplingWeight << "," << endl
            << "  alphaU = " << indent(m_alphaU->toString()) << "," << endl
            << "  alphaV = " << indent(m_alphaV->toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    EModelVariant m_modelVariant;
    ref<Texture> m_diffuseReflectance;
    ref<Texture> m_specularReflectance;
    ref<Texture> m_alphaU;
    ref<Texture> m_alphaV;
    Float m_specularSamplingWeight;
};

// ================ Hardware shader implementation ================

/**
 * GLSL port of the Ward shader. This version only implements the variant
 * with energy balance. When the roughness is lower than
 * \alpha < 0.2, the shader clamps it to 0.2 so that it will still perform
 * reasonably well in a VPL-based preview.
 */
class WardShader : public Shader {
public:
    WardShader(Renderer *renderer,
            const Texture *diffuseColor,
            const Texture *specularColor,
            const Texture *alphaU,
            const Texture *alphaV) : Shader(renderer, EBSDFShader),
            m_diffuseReflectance(diffuseColor),
            m_specularReflectance(specularColor),
            m_alphaU(alphaU), m_alphaV(alphaV) {
        m_diffuseReflectanceShader = renderer->registerShaderForResource(m_diffuseReflectance.get());
        m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());
        m_alphaUShader = renderer->registerShaderForResource(m_alphaU.get());
        m_alphaVShader = renderer->registerShaderForResource(m_alphaV.get());
    }

    bool isComplete() const {
        return m_diffuseReflectanceShader.get() != NULL &&
               m_specularReflectanceShader.get() != NULL &&
               m_alphaU.get() != NULL &&
               m_alphaV.get() != NULL;
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(m_diffuseReflectanceShader.get());
        deps.push_back(m_specularReflectanceShader.get());
        deps.push_back(m_alphaUShader.get());
        deps.push_back(m_alphaVShader.get());
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(m_diffuseReflectance.get());
        renderer->unregisterShaderForResource(m_specularReflectance.get());
        renderer->unregisterShaderForResource(m_alphaU.get());
        renderer->unregisterShaderForResource(m_alphaV.get());
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (wi.z <= 0.0 || wo.z <= 0.0)" << endl
            << "        return vec3(0.0);" << endl
            << "    vec3 H = wi + wo;" << endl
            << "    float cosSqr = H.z * H.z;" << endl
            << "    float alphaU = max(0.3, " << depNames[2] << "(uv)[0]);" << endl
            << "    float alphaV = max(0.3, " << depNames[3] << "(uv)[0]);" << endl
            << "    float factor1 = dot(H, H)/(3.1415*alphaU*alphaV*cosSqr*cosSqr);"  << endl
            << "    float factor2 = H.x / alphaU, factor3 = H.y / alphaV;" << endl
            << "    float exponent = -(factor2*factor2 + factor3*factor3)/(H.z*H.z);" << endl
            << "    float specRef = factor1 * exp(exponent);" << endl
            << "    return (" << depNames[0] << "(uv) * inv_pi" << endl
            << "           + " << depNames[1] << "(uv) * specRef) * cosTheta(wo);" << endl
            << "}" << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (wi.z <= 0.0 || wo.z <= 0.0)" << endl
            << "        return vec3(0.0);" << endl
            << "    return " << depNames[0] << "(uv) * (inv_pi * cosTheta(wo));" << endl
            << "}" << endl;
    }

    MTS_DECLARE_CLASS()
private:
    ref<const Texture> m_diffuseReflectance;
    ref<const Texture> m_specularReflectance;
    ref<const Texture> m_alphaU;
    ref<const Texture> m_alphaV;
    ref<Shader> m_diffuseReflectanceShader;
    ref<Shader> m_specularReflectanceShader;
    ref<Shader> m_alphaUShader;
    ref<Shader> m_alphaVShader;
};

Shader *Ward::createShader(Renderer *renderer) const {
    return new WardShader(renderer, m_diffuseReflectance.get(),
        m_specularReflectance.get(), m_alphaU.get(), m_alphaV.get());
}

MTS_IMPLEMENT_CLASS(WardShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Ward, false, BSDF);
MTS_EXPORT_PLUGIN(Ward, "Anisotropic Ward BRDF");
MTS_NAMESPACE_END

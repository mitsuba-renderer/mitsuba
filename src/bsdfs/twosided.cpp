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
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{twosided}{Two-sided BRDF adapter}
 * \order{19}
 * \parameters{
 *     \parameter{\Unnamed}{\BSDF}{A nested BRDF that should
 *     be turned into a two-sided scattering model. If two BRDFs
 *     are specified, they will be placed on the front and back side, respectively.}
 * }
 *
 * \renderings{
 *     \unframedrendering{From this angle, the Cornell box scene shows visible back-facing geometry}
 *         {bsdf_twosided_before}
 *     \unframedrendering{Applying the \pluginref{twosided} plugin fixes the rendering}
 *         {bsdf_twosided_after}
 * }
 *
 * By default, all non-transmissive scattering models in Mitsuba
 * are \emph{one-sided} --- in other words, they absorb all light
 * that is received on the interior-facing side of any associated
 * surfaces. Holes and visible back-facing parts are thus exposed
 * as black regions.
 *
 * Usually, this is a good idea, since it will reveal modeling
 * issues early on. But sometimes one is forced to deal with
 * improperly closed geometry, where the one-sided behavior is
 * bothersome. In that case, this plugin can be used to turn
 * one-sided scattering models into proper two-sided versions of
 * themselves. The plugin has no parameters other than a required
 * nested BSDF specification. It is also possible to supply two
 * different BRDFs that should be placed on the front and back
 * side, respectively.
 * \vspace{4mm}
 *
 * \begin{xml}[caption=A two-sided diffuse material]
 * <bsdf type="twosided">
 *     <bsdf type="diffuse">
 *          <spectrum name="reflectance" value="0.4"/>
 *     </bsdf>
 * </bsdf>
 * \end{xml}
 */
class TwoSidedBRDF : public BSDF {
public:
    TwoSidedBRDF(const Properties &props)
        : BSDF(props) { }

    TwoSidedBRDF(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager) {
        m_nestedBRDF[0] = static_cast<BSDF *>(manager->getInstance(stream));
        m_nestedBRDF[1] = static_cast<BSDF *>(manager->getInstance(stream));
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        manager->serialize(stream, m_nestedBRDF[0].get());
        manager->serialize(stream, m_nestedBRDF[1].get());
    }

    void configure() {
        if (!m_nestedBRDF[0])
            Log(EError, "A nested one-sided material is required!");
        if (!m_nestedBRDF[1])
            m_nestedBRDF[1] = m_nestedBRDF[0];


        m_usesRayDifferentials = m_nestedBRDF[0]->usesRayDifferentials()
            || m_nestedBRDF[1]->usesRayDifferentials();

        m_components.clear();

        for (int i=0; i<m_nestedBRDF[0]->getComponentCount(); ++i)
            m_components.push_back((m_nestedBRDF[0]->getType(i) & ~EBackSide) | EFrontSide);

        for (int i=0; i<m_nestedBRDF[1]->getComponentCount(); ++i)
            m_components.push_back((m_nestedBRDF[1]->getType(i) & ~EFrontSide) | EBackSide);

        BSDF::configure();
        if (m_combinedType & BSDF::ETransmission)
            Log(EError, "Only materials without "
                "a transmission component can be nested!");
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        BSDFSamplingRecord b(bRec);

        if (Frame::cosTheta(b.wi) > 0) {
            return m_nestedBRDF[0]->eval(b, measure);
        } else {
            if (b.component != -1)
                b.component -= m_nestedBRDF[0]->getComponentCount();
            b.wi.z *= -1;
            b.wo.z *= -1;
            return m_nestedBRDF[1]->eval(b, measure);
        }
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        BSDFSamplingRecord b(bRec);

        if (b.wi.z > 0) {
            return m_nestedBRDF[0]->pdf(b, measure);
        } else {
            if (b.component != -1)
                b.component -= m_nestedBRDF[0]->getComponentCount();
            b.wi.z *= -1;
            b.wo.z *= -1;
            return m_nestedBRDF[1]->pdf(b, measure);
        }
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        bool flipped = false;

        if (Frame::cosTheta(bRec.wi) < 0) {
            bRec.wi.z *= -1;
            flipped = true;
            if (bRec.component != -1)
                bRec.component -= m_nestedBRDF[0]->getComponentCount();
        }

        Spectrum result = m_nestedBRDF[flipped ? 1 : 0]->sample(bRec, sample);

        if (flipped) {
            bRec.wi.z *= -1;
            if (bRec.component != -1)
                bRec.component += m_nestedBRDF[0]->getComponentCount();
            if (!result.isZero()) {
                bRec.wo.z *= -1;
                bRec.sampledComponent += m_nestedBRDF[0]->getComponentCount();
            }
        }

        return result;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
        bool flipped = false;
        if (Frame::cosTheta(bRec.wi) < 0) {
            bRec.wi.z *= -1;
            flipped = true;
            if (bRec.component != -1)
                bRec.component -= m_nestedBRDF[0]->getComponentCount();
        }

        Spectrum result = m_nestedBRDF[flipped ? 1 : 0]->sample(bRec, pdf, sample);

        if (flipped) {
            bRec.wi.z *= -1;

            if (bRec.component != -1)
                bRec.component += m_nestedBRDF[0]->getComponentCount();
            if (!result.isZero() && pdf != 0) {
                bRec.wo.z *= -1;
                bRec.sampledComponent += m_nestedBRDF[0]->getComponentCount();
            }
        }
        return result;
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(BSDF::m_theClass)) {
            if (m_nestedBRDF[0] == NULL)
                m_nestedBRDF[0] = static_cast<BSDF *>(child);
            else if (m_nestedBRDF[1] == NULL)
                m_nestedBRDF[1] = static_cast<BSDF *>(child);
            else
                Log(EError, "No more than two nested BRDFs can be added!");
        } else {
            BSDF::addChild(name, child);
        }
    }

    Spectrum getDiffuseReflectance(const Intersection &its) const {
        if (its.wi.z > 0)
            return m_nestedBRDF[0]->getDiffuseReflectance(its);
        else
            return m_nestedBRDF[1]->getDiffuseReflectance(its);
    }

    Float getRoughness(const Intersection &its, int component) const {
        if (component < m_nestedBRDF[0]->getComponentCount()) {
            return m_nestedBRDF[0]->getRoughness(its, component);
        } else {
            return m_nestedBRDF[1]->getRoughness(its, component
                  -m_nestedBRDF[0]->getComponentCount());
        }
    }

    Float getEta() const {
        return 1.0f;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "TwoSided[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  nestedBRDF[0] = " << indent(m_nestedBRDF[0].toString()) << "," << endl
            << "  nestedBRDF[1] = " << indent(m_nestedBRDF[1].toString()) << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
protected:
    ref<BSDF> m_nestedBRDF[2];
};


// ================ Hardware shader implementation ================

class TwoSidedShader : public Shader {
public:
    TwoSidedShader(Renderer *renderer,
            const ref<BSDF> *nestedBRDF) : Shader(renderer, EBSDFShader) {
        m_nestedBRDF[0] = nestedBRDF[0].get();
        m_nestedBRDF[1] = nestedBRDF[1].get();

        m_nestedBRDFShader[0] = renderer->registerShaderForResource(m_nestedBRDF[0]);
        if (m_nestedBRDF[0] != m_nestedBRDF[1])
            m_nestedBRDFShader[1] = renderer->registerShaderForResource(m_nestedBRDF[1]);
        else
            m_nestedBRDFShader[1] = NULL;
    }

    bool isComplete() const {
        return m_nestedBRDFShader[0].get() != NULL &&
               (m_nestedBRDF[0] == m_nestedBRDF[1] || m_nestedBRDFShader[1].get() != NULL);
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(m_nestedBRDFShader[0].get());
        if (m_nestedBRDF[0] != m_nestedBRDF[1])
            deps.push_back(m_nestedBRDFShader[1].get());
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(m_nestedBRDF[0]);
        if (m_nestedBRDF[0] != m_nestedBRDF[1])
            renderer->unregisterShaderForResource(m_nestedBRDF[1]);
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (cosTheta(wi) <= 0.0) {" << endl
            << "        wi.z *= -1; wo.z *= -1;" << endl
            << "        return " << (depNames.size() == 2 ? depNames[1] : depNames[0]) << "(uv, wi, wo);" << endl
            << "    } else {" << endl
            << "        return " << depNames[0] << "(uv, wi, wo);" << endl
            << "    }" << endl
            << "}" << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (cosTheta(wi) <= 0.0) {" << endl
            << "        wi.z *= -1; wo.z *= -1;" << endl
            << "        return " << (depNames.size() == 2 ? depNames[1] : depNames[0]) << "_diffuse(uv, wi, wo);" << endl
            << "    } else {" << endl
            << "        return " << depNames[0] << "_diffuse(uv, wi, wo);" << endl
            << "    }" << endl
            << "}" << endl;
    }

    MTS_DECLARE_CLASS()
private:
    const BSDF *m_nestedBRDF[2];
    ref<Shader> m_nestedBRDFShader[2];
};

Shader *TwoSidedBRDF::createShader(Renderer *renderer) const {
    return new TwoSidedShader(renderer, m_nestedBRDF);
}

MTS_IMPLEMENT_CLASS(TwoSidedShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(TwoSidedBRDF, false, BSDF)
MTS_EXPORT_PLUGIN(TwoSidedBRDF, "Two-sided BRDF adapter");
MTS_NAMESPACE_END

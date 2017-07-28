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
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

/*! \plugin{blendbsdf}{Blended material}
 * \order{17}
 * \parameters{
 *     \parameter{weight}{\Float\Or\Texture}{A floating point value or texture
 *      with values between zero and one. The extreme values zero and one activate the
 *      first and second nested BSDF respectively, and inbetween values
 *      interpolate accordingly. \default{0.5}}
 *     \parameter{\Unnamed}{\BSDF}{Two nested BSDF instances that should be
 *     mixed according to the specified blending weight}
 * }
 *
 * \renderings{
 *     \rendering{A material created by blending between dark rough plastic and
 *     smooth gold based on a binary bitmap texture (\lstref{blendbsdf})}{bsdf_blendbsdf}
 * }
 *
 * This plugin implements a ``blend'' material, which represents
 * linear combinations of two BSDF instances. It is conceptually very similar
 * to the \pluginref{mixturebsdf} plugin. The main difference is that
 * \pluginref{blendbsdf} can interpolate based on a texture rather than a set
 * of constants.
 *
 * Any surface scattering model in Mitsuba (be it smooth, rough, reflecting, or
 * transmitting) can be mixed with others in this manner to synthesize new models.
 *
 * \begin{xml}[caption=Description of the material shown above,
 *     label=lst:blendbsdf]
 * <bsdf type="blendbsdf">
 *     <texture name="weight" type="bitmap">
 *         <string name="wrapMode" value="repeat"/>
 *         <string name="filename" value="pattern.png"/>
 *     </texture>
 *
 *     <bsdf type="conductor">
 *         <string name="material" value="Au"/>
 *     </bsdf>
 *
 *     <bsdf type="roughplastic">
 *         <spectrum name="diffuseReflectance" value="0"/>
 *     </bsdf>
 * </bsdf>
 * \end{xml}
 */

class BlendBSDF : public BSDF {
public:
    BlendBSDF(const Properties &props)
        : BSDF(props) {
        m_weight = new ConstantFloatTexture(props.getFloat("weight", 0.5f));
    }

    BlendBSDF(Stream *stream, InstanceManager *manager)
     : BSDF(stream, manager) {
        m_weight = static_cast<Texture *>(manager->getInstance(stream));
        m_bsdfs.push_back(static_cast<BSDF *>(manager->getInstance(stream)));
        m_bsdfs.push_back(static_cast<BSDF *>(manager->getInstance(stream)));
        m_bsdfs[0]->incRef();
        m_bsdfs[1]->incRef();
        configure();
    }

    virtual ~BlendBSDF() {
        for (size_t i=0; i<m_bsdfs.size(); ++i)
            m_bsdfs[i]->decRef();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        Assert(m_bsdfs.size() == 2);
        manager->serialize(stream, m_weight.get());
        manager->serialize(stream, m_bsdfs[0]);
        manager->serialize(stream, m_bsdfs[1]);
    }

    void configure() {
        m_usesRayDifferentials = false;
        size_t componentCount = 0;

        if (m_bsdfs.size() != 2)
            Log(EError, "BSDF count mismatch: expected two nested BSDF instances!");

        for (size_t i=0; i<m_bsdfs.size(); ++i)
            componentCount += m_bsdfs[i]->getComponentCount();

        m_components.reserve(componentCount);
        m_components.clear();
        m_indices.reserve(componentCount);
        m_indices.clear();
        m_offsets.reserve(m_bsdfs.size());
        m_offsets.clear();

        int offset = 0;
        for (size_t i=0; i<m_bsdfs.size(); ++i) {
            const BSDF *bsdf = m_bsdfs[i];
            m_offsets.push_back(offset);

            for (int j=0; j<bsdf->getComponentCount(); ++j) {
                int componentType = bsdf->getType(j);
                m_components.push_back(componentType);
                m_indices.push_back(std::make_pair((int) i, j));
            }

            offset += bsdf->getComponentCount();
            m_usesRayDifferentials |= bsdf->usesRayDifferentials();
        }
        BSDF::configure();
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        Float weight = std::min((Float) 1.0f, std::max((Float) 0.0f,
            m_weight->eval(bRec.its).average()));

        if (bRec.component == -1) {
            return
                m_bsdfs[0]->eval(bRec, measure) * (1-weight) +
                m_bsdfs[1]->eval(bRec, measure) * weight;
        } else {
            /* Pick out an individual component */
            int idx = m_indices[bRec.component].first;
            if (idx == 0)
                weight = 1-weight;
            BSDFSamplingRecord bRec2(bRec);
            bRec2.component = m_indices[bRec.component].second;
            return m_bsdfs[idx]->eval(bRec2, measure) * weight;
        }
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        Spectrum result;

        Float weight = std::min((Float) 1.0f, std::max((Float) 0.0f,
            m_weight->eval(bRec.its).average()));

        if (bRec.component == -1) {
            return
                m_bsdfs[0]->pdf(bRec, measure) * (1-weight) +
                m_bsdfs[1]->pdf(bRec, measure) * weight;
        } else {
            /* Pick out an individual component */
            int idx = m_indices[bRec.component].first;
            if (idx == 0)
                weight = 1-weight;
            BSDFSamplingRecord bRec2(bRec);
            bRec2.component = m_indices[bRec.component].second;
            return m_bsdfs[idx]->pdf(bRec2, measure) * weight;
        }
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &_sample) const {
        Point2 sample(_sample);

        Float weights[2];
        weights[1] = std::min((Float) 1.0f, std::max((Float) 0.0f,
            m_weight->eval(bRec.its).average()));
        weights[0] = 1-weights[1];

        if (bRec.component == -1) {
            size_t entry;
            if (sample.x < weights[0]) {
                entry = 0; sample.x /= weights[0];
            } else {
                entry = 1; sample.x = (sample.x - weights[0]) / weights[1];
            }

            Float pdf;
            Spectrum result = m_bsdfs[entry]->sample(bRec, pdf, sample);
            if (result.isZero()) // sampling failed
                return result;

            result *= weights[entry] * pdf;
            pdf *= weights[entry];

            EMeasure measure = BSDF::getMeasure(bRec.sampledType);
            for (size_t i=0; i<m_bsdfs.size(); ++i) {
                if (entry == i)
                    continue;
                pdf += m_bsdfs[i]->pdf(bRec, measure) * weights[i];
                result += m_bsdfs[i]->eval(bRec, measure) * weights[i];
            }

            bRec.sampledComponent += m_offsets[entry];
            return result/pdf;
        } else {
            /* Pick out an individual component */
            int requestedComponent = bRec.component;
            int bsdfIndex = m_indices[requestedComponent].first;
            bRec.component = m_indices[requestedComponent].second;
            Spectrum result = m_bsdfs[bsdfIndex]->sample(bRec, sample)
                * weights[bsdfIndex];
            bRec.component = bRec.sampledComponent = requestedComponent;
            return result;
        }
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample) const {
        Point2 sample(_sample);

        Float weights[2];
        weights[1] = std::min((Float) 1.0f, std::max((Float) 0.0f,
            m_weight->eval(bRec.its).average()));
        weights[0] = 1-weights[1];

        if (bRec.component == -1) {
            size_t entry;
            if (sample.x < weights[0]) {
                entry = 0; sample.x /= weights[0];
            } else {
                entry = 1; sample.x = (sample.x - weights[0]) / weights[1];
            }

            Spectrum result = m_bsdfs[entry]->sample(bRec, pdf, sample);
            if (result.isZero()) // sampling failed
                return result;

            result *= weights[entry] * pdf;
            pdf *= weights[entry];

            EMeasure measure = BSDF::getMeasure(bRec.sampledType);
            for (size_t i=0; i<m_bsdfs.size(); ++i) {
                if (entry == i)
                    continue;
                pdf += m_bsdfs[i]->pdf(bRec, measure) * weights[i];
                result += m_bsdfs[i]->eval(bRec, measure) * weights[i];
            }

            bRec.sampledComponent += m_offsets[entry];
            return result/pdf;
        } else {
            /* Pick out an individual component */
            int requestedComponent = bRec.component;
            int bsdfIndex = m_indices[requestedComponent].first;
            bRec.component = m_indices[requestedComponent].second;
            Spectrum result = m_bsdfs[bsdfIndex]->sample(bRec, pdf, sample)
                * weights[bsdfIndex];
            bRec.component = bRec.sampledComponent = requestedComponent;
            return result;
        }
    }

    Float getRoughness(const Intersection &its, int component) const {
        int bsdfIndex = m_indices[component].first;
        component = m_indices[component].second;
        return m_bsdfs[bsdfIndex]->getRoughness(its, component);
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
            BSDF *bsdf = static_cast<BSDF *>(child);
            m_bsdfs.push_back(bsdf);
            bsdf->incRef();
        } else if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "weight") {
            m_weight = static_cast<Texture *>(child);
        } else {
            BSDF::addChild(name, child);
        }
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "BlendBSDF[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  weight = " << indent(m_weight->toString()) << endl
            << "  bsdfs = {" << endl;
        for (size_t i=0; i<m_bsdfs.size(); ++i)
            oss << "    " << indent(m_bsdfs[i]->toString(), 2) << "," << endl;
        oss << "  }"
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    std::vector<BSDF *> m_bsdfs;
    ref<Texture> m_weight;
    std::vector<std::pair<int, int> > m_indices;
    std::vector<int> m_offsets;
};


MTS_IMPLEMENT_CLASS_S(BlendBSDF, false, BSDF)
MTS_EXPORT_PLUGIN(BlendBSDF, "Blend BSDF")
MTS_NAMESPACE_END

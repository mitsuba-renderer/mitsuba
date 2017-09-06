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

/*! \plugin{mixturebsdf}{Mixture material}
 * \order{16}
 * \parameters{
 *     \parameter{weights}{\String}{A comma-separated list of BSDF weights}
 *     \parameter{\Unnamed}{\BSDF}{Multiple BSDF instances that should be
 *     mixed according to the specified weights}
 * }
 * \renderings{
 *     \medrendering{Smooth glass}{bsdf_mixturebsdf_smooth}
 *     \medrendering{Rough glass}{bsdf_mixturebsdf_rough}
 *     \medrendering{An mixture of 70% smooth glass and 30% rough glass
 *     results in a more realistic smooth material with imperfections
 *     (\lstref{mixture-example})}{bsdf_mixturebsdf_result}
 * }
 *
 * This plugin implements a ``mixture'' material, which represents
 * linear combinations of multiple BSDF instances. Any surface scattering
 * model in Mitsuba (be it smooth, rough, reflecting, or transmitting) can
 * be mixed with others in this manner to synthesize new models. There
 * is no limit on how many models can be mixed, but their combination
 * weights must be non-negative and sum to a value of one or less to ensure
 * energy balance. When they sum to less than one, the material will
 * absorb a proportional amount of the incident illlumination.
 *
 * \vspace{4mm}
 * \begin{xml}[caption={A material definition for a mixture of 70% smooth
 *     and 30% rough glass},
 *     label=lst:mixture-example]
 * <bsdf type="mixturebsdf">
 *     <string name="weights" value="0.7, 0.3"/>
 *
 *     <bsdf type="dielectric"/>
 *
 *     <bsdf type="roughdielectric">
 *         <float name="alpha" value="0.3"/>
 *     </bsdf>
 * </bsdf>
 * \end{xml}
 */

class MixtureBSDF : public BSDF {
public:
    MixtureBSDF(const Properties &props)
        : BSDF(props) {
        /* Parse the weight parameter */
        std::vector<std::string> weights =
            tokenize(props.getString("weights", ""), " ,;");
        if (weights.size() == 0)
            Log(EError, "No weights were supplied!");
        m_weights.resize(weights.size());

        char *end_ptr = NULL;
        for (size_t i=0; i<weights.size(); ++i) {
            Float weight = (Float) strtod(weights[i].c_str(), &end_ptr);
            if (*end_ptr != '\0')
                SLog(EError, "Could not parse the BSDF weights!");
            if (weight < 0)
                SLog(EError, "Invalid BSDF weight!");
            m_weights[i] = weight;
        }
    }

    MixtureBSDF(Stream *stream, InstanceManager *manager)
     : BSDF(stream, manager) {
        size_t bsdfCount = stream->readSize();
        m_weights.resize(bsdfCount);
        for (size_t i=0; i<bsdfCount; ++i) {
            m_weights[i] = stream->readFloat();
            BSDF *bsdf = static_cast<BSDF *>(manager->getInstance(stream));
            bsdf->incRef();
            m_bsdfs.push_back(bsdf);
        }
        configure();
    }

    virtual ~MixtureBSDF() {
        for (size_t i=0; i<m_bsdfs.size(); ++i)
            m_bsdfs[i]->decRef();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        stream->writeSize(m_bsdfs.size());
        for (size_t i=0; i<m_bsdfs.size(); ++i) {
            stream->writeFloat(m_weights[i]);
            manager->serialize(stream, m_bsdfs[i]);
        }
    }

    void configure() {
        m_usesRayDifferentials = false;
        size_t componentCount = 0;

        if (m_bsdfs.size() != m_weights.size())
            Log(EError, "BSDF count mismatch: " SIZE_T_FMT " bsdfs, but specified " SIZE_T_FMT " weights",
                m_bsdfs.size(), m_bsdfs.size());

        Float totalWeight = 0;
        for (size_t i=0; i<m_weights.size(); ++i)
            totalWeight += m_weights[i];

        if (totalWeight <= 0)
            Log(EError, "The weights must sum to a value greater than zero!");

        if (m_ensureEnergyConservation && totalWeight > 1) {
            std::ostringstream oss;
            Float scale = 1.0f / totalWeight;
            oss << "The BSDF" << endl << toString() << endl
                << "potentially violates energy conservation, since the weights "
                << "sum to " << totalWeight << ", which is greater than one! "
                << "They will be re-scaled to avoid potential issues. Specify "
                << "the parameter ensureEnergyConservation=false to prevent "
                << "this from happening.";
            Log(EWarn, "%s", oss.str().c_str());
            for (size_t i=0; i<m_weights.size(); ++i)
                m_weights[i] *= scale;
        }

        for (size_t i=0; i<m_bsdfs.size(); ++i)
            componentCount += m_bsdfs[i]->getComponentCount();

        m_pdf = DiscreteDistribution(m_bsdfs.size());
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
            m_pdf.append(m_weights[i]);
        }
        m_pdf.normalize();
        BSDF::configure();
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        Spectrum result(0.0f);

        if (bRec.component == -1) {
            for (size_t i=0; i<m_bsdfs.size(); ++i)
                result += m_bsdfs[i]->eval(bRec, measure) * m_weights[i];
        } else {
            /* Pick out an individual component */
            int idx = m_indices[bRec.component].first;
            BSDFSamplingRecord bRec2(bRec);
            bRec2.component = m_indices[bRec.component].second;
            return m_bsdfs[idx]->eval(bRec2, measure) * m_weights[idx];
        }

        return result;
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        Float result = 0.0f;

        if (bRec.component == -1) {
            for (size_t i=0; i<m_bsdfs.size(); ++i)
                result += m_bsdfs[i]->pdf(bRec, measure) * m_pdf[i];
        } else {
            /* Pick out an individual component */
            int idx = m_indices[bRec.component].first;
            BSDFSamplingRecord bRec2(bRec);
            bRec2.component = m_indices[bRec.component].second;
            return m_bsdfs[idx]->pdf(bRec2, measure);
        }

        return result;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &_sample) const {
        Point2 sample(_sample);
        if (bRec.component == -1) {
            /* Choose a component based on the normalized weights */
            size_t entry = m_pdf.sampleReuse(sample.x);

            Float pdf;
            Spectrum result = m_bsdfs[entry]->sample(bRec, pdf, sample);
            if (result.isZero()) // sampling failed
                return result;

            result *= m_weights[entry] * pdf;
            pdf *= m_pdf[entry];

            EMeasure measure = BSDF::getMeasure(bRec.sampledType);
            for (size_t i=0; i<m_bsdfs.size(); ++i) {
                if (entry == i)
                    continue;
                pdf += m_bsdfs[i]->pdf(bRec, measure) * m_pdf[i];
                result += m_bsdfs[i]->eval(bRec, measure) * m_weights[i];
            }

            bRec.sampledComponent += m_offsets[entry];
            return result / pdf;
        } else {
            /* Pick out an individual component */
            int requestedComponent = bRec.component;
            int bsdfIndex = m_indices[requestedComponent].first;
            bRec.component = m_indices[requestedComponent].second;
            Spectrum result = m_bsdfs[bsdfIndex]->sample(bRec, sample)
                * m_weights[bsdfIndex];
            bRec.component = bRec.sampledComponent = requestedComponent;
            return result;
        }
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample) const {
        Point2 sample(_sample);
        if (bRec.component == -1) {
            /* Choose a component based on the normalized weights */
            size_t entry = m_pdf.sampleReuse(sample.x);

            Spectrum result = m_bsdfs[entry]->sample(bRec, pdf, sample);
            if (result.isZero()) // sampling failed
                return result;

            result *= m_weights[entry] * pdf;
            pdf *= m_pdf[entry];

            EMeasure measure = BSDF::getMeasure(bRec.sampledType);
            for (size_t i=0; i<m_bsdfs.size(); ++i) {
                if (entry == i)
                    continue;
                pdf += m_bsdfs[i]->pdf(bRec, measure) * m_pdf[i];
                result += m_bsdfs[i]->eval(bRec, measure) * m_weights[i];
            }

            bRec.sampledComponent += m_offsets[entry];
            return result/pdf;
        } else {
            /* Pick out an individual component */
            int requestedComponent = bRec.component;
            int bsdfIndex = m_indices[requestedComponent].first;
            bRec.component = m_indices[requestedComponent].second;
            Spectrum result = m_bsdfs[bsdfIndex]->sample(bRec, pdf, sample)
                * m_weights[bsdfIndex];
            bRec.component = bRec.sampledComponent = requestedComponent;
            return result;
        }
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
            BSDF *bsdf = static_cast<BSDF *>(child);
            m_bsdfs.push_back(bsdf);
            bsdf->incRef();
        } else {
            BSDF::addChild(name, child);
        }
    }

    Float getRoughness(const Intersection &its, int component) const {
        int bsdfIndex = m_indices[component].first;
        component = m_indices[component].second;
        return m_bsdfs[bsdfIndex]->getRoughness(its, component);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "MixtureBSDF[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  weights = {";
        for (size_t i=0; i<m_bsdfs.size(); ++i) {
            oss << " " << m_weights[i];
            if (i + 1 < m_bsdfs.size())
                oss << ",";
        }
        oss << " }," << endl
            << "  bsdfs = {" << endl;
        for (size_t i=0; i<m_bsdfs.size(); ++i)
            oss << "    " << indent(m_bsdfs[i]->toString(), 2) << "," << endl;
        oss << "  }" << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    std::vector<Float> m_weights;
    std::vector<std::pair<int, int> > m_indices;
    std::vector<int> m_offsets;
    std::vector<BSDF *> m_bsdfs;
    DiscreteDistribution m_pdf;
};

// ================ Hardware shader implementation ================

class MixtureBSDFShader : public Shader {
public:
    MixtureBSDFShader(Renderer *renderer, const std::vector<BSDF *> &bsdfs, const std::vector<Float> &weights)
        : Shader(renderer, EBSDFShader), m_bsdfs(bsdfs), m_weights(weights), m_complete(false) {
        m_bsdfShader.resize(bsdfs.size());
        for (size_t i=0; i<bsdfs.size(); ++i) {
            ref<Shader> shader = renderer->registerShaderForResource(bsdfs[i]);
            if (shader) {
                shader->incRef();
                /* At least one shader has a hardware implementation */
                m_complete = true;
            }
            m_bsdfShader[i] = shader;
        }
    }

    bool isComplete() const {
        return m_complete;
    }

    void cleanup(Renderer *renderer) {
        for (size_t i=0; i<m_bsdfs.size(); ++i) {
            renderer->unregisterShaderForResource(m_bsdfs[i]);
            if (m_bsdfShader[i])
                m_bsdfShader[i]->decRef();
        }
        m_bsdfShader.clear();
    }

    void putDependencies(std::vector<Shader *> &deps) {
        for (size_t i=0; i<m_bsdfs.size(); ++i) {
            if (m_bsdfShader[i])
                deps.push_back(m_bsdfShader[i]);
        }
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        Assert(m_complete);
        int ctr = 0;
        for (size_t i=0; i<m_bsdfs.size(); ++i) {
            if (!m_bsdfShader[i])
                continue;
            oss << "uniform float " << evalName << "_weight_" << ctr++ << ";" << endl;
        }
        oss << endl;
        oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    return ";
        ctr = 0;
        for (size_t i=0; i<m_bsdfs.size(); ++i) {
            if (!m_bsdfShader[i])
                continue;
            oss << endl << "      ";
            if (ctr != 0)
                oss << "+ ";
            else
                oss << "  ";
            oss << depNames[ctr] << "(uv, wi, wo) * "
                              << evalName << "_weight_" << ctr;
            ctr++;
        }
        oss << ";" << endl << "}" << endl << endl;
        oss << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    return ";
        ctr = 0;
        for (size_t i=0; i<m_bsdfs.size(); ++i) {
            if (!m_bsdfShader[i])
                continue;
            oss << endl << "      ";
            if (ctr != 0)
                oss << "+ ";
            else
                oss << "  ";
            oss << depNames[ctr] << "_diffuse(uv, wi, wo) * "
                              << evalName << "_weight_" << ctr;
            ctr++;
        }
        oss << ";" << endl << "}" << endl;
    }

    void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
        int ctr = 0;
        for (size_t i=0; i<m_bsdfs.size(); ++i) {
            if (!m_bsdfShader[i])
                continue;
            parameterIDs.push_back(
                program->getParameterID(formatString("%s_weight_%i", evalName.c_str(), ctr++), false));
        }
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
        int ctr = 0;
        for (size_t i=0; i<m_bsdfs.size(); ++i) {
            if (!m_bsdfShader[i])
                continue;
            program->setParameter(parameterIDs[ctr++], m_weights[i]);
        }
    }

    MTS_DECLARE_CLASS()
private:
    std::vector<Shader *> m_bsdfShader;
    const std::vector<BSDF *> &m_bsdfs;
    const std::vector<Float> &m_weights;
    bool m_complete;
};

Shader *MixtureBSDF::createShader(Renderer *renderer) const {
    return new MixtureBSDFShader(renderer, m_bsdfs, m_weights);
}

MTS_IMPLEMENT_CLASS(MixtureBSDFShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(MixtureBSDF, false, BSDF)
MTS_EXPORT_PLUGIN(MixtureBSDF, "Mixture BSDF")
MTS_NAMESPACE_END

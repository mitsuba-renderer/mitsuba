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

/* 
 * File:   pssmlt_sampler_base.h
 * Author: swl
 *
 * Created on April 13, 2017, 11:17 PM
 */

#ifndef PSSMLT_SAMPLER_BASE_H
#define PSSMLT_SAMPLER_BASE_H

#include <mitsuba/render/sampler.h>
#include <mitsuba/core/random.h>

/// Use Kelemen-style mutations in random number space?
#define MUTATION_STYLE 2
// 0: KELEMEN  1: VEACH  2: TOSHIYA

MTS_NAMESPACE_BEGIN

class PSSMLTSamplerBase : public Sampler {
public:
    // Construct a new MLT sampler

    PSSMLTSamplerBase(Float mutationSize) : Sampler(Properties()) {
        m_random = new Random();
        m_s1 = 1.0f/1024.0f;
        m_s2 = 1.0f/64.0f;
        configure();
    }

    /// Set whether the current step should be large

    inline void setLargeStep(bool value) {
        m_largeStep = value;
    }

    /// Check if the current step is a large step

    inline bool isLargeStep() const {
        return m_largeStep;
    }

    /// 1D mutation routine
    inline void updateStepSize(const Float& s1, const Float& s2) {
        m_s1 = s1;
        m_s2 = s2;
        m_logRatio = -math::fastlog(m_s2 / m_s1);
    }
    
    inline void updateStrength(const Float& s) {
        m_s1 = m_s2 * s;
        m_logRatio = -math::fastlog(m_s2 / m_s1);
    }
    
    inline Float getStrength() const {
        return m_s1 / m_s2;
    }

    inline Float mutate(Float value) {
#if MUTATION_STYLE == 0
        Float sample = m_random->nextFloat();
        bool add;

        if (sample < 0.5f) {
            add = true;
            sample *= 2.0f;
        } else {
            add = false;
            sample = 2.0f * (sample - 0.5f);
        }

        Float dv = m_s2 * math::fastexp(sample * m_logRatio);
        if (add) {
            value += dv;
            if (value > 1)
                value -= 1;
        } else {
            value -= dv;
            if (value < 0)
                value += 1;
        }
#elif MUTATION_STYLE == 1
        Float tmp1 = std::sqrt(-2 * std::log(1 - m_random->nextFloat()));
        Float dv = tmp1 * std::cos(2 * M_PI * m_random->nextFloat());
        value = modulo(value + 1e-2f * dv, 1.0f);
#else
        Float sample = m_random->nextFloat();
        if (sample < 0.5f) {
            sample *= 2.0f;
            value += pow(sample, m_s2 / m_s1 + 1.f);
        } else {
            sample = 2.0f * (sample - 0.5f);
            value -= pow(sample, m_s2 / m_s1 + 1.f);
        }
        value -= floor(value);
#endif
        return value;
    }

    /// Replace the underlying random number generator

    inline void setRandom(Random *random) {
        m_random = random;
    }

    /// Return the underlying random number generator

    inline Random *getRandom() {
        return m_random;
    }

    /* The following functions do nothing in this implementation */
    virtual void advance() {
    }

    virtual void generate(const Point2i &pos) {
    }

    /* The following functions are unsupported by this implementation */
    void request1DArray(size_t size) {
        Log(EError, "request1DArray(): Unsupported!");
    }

    void request2DArray(size_t size) {
        Log(EError, "request2DArray(): Unsupported!");
    }

    void setSampleIndex(size_t sampleIndex) {
        Log(EError, "setSampleIndex(): Unsupported!");
    }

    PSSMLTSamplerBase(PSSMLTSamplerBase *sampler) : Sampler(Properties()),
    m_random(sampler->m_random) {
        m_s1 = sampler->m_s1;
        m_s2 = sampler->m_s2;
        configure();
    }

    PSSMLTSamplerBase(Stream *stream, InstanceManager *manager)
    : Sampler(stream, manager) {
        m_random = static_cast<Random *> (manager->getInstance(stream));
        m_s1 = stream->readFloat();
        m_s2 = stream->readFloat();
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Sampler::serialize(stream, manager);
        manager->serialize(stream, m_random.get());
        stream->writeFloat(m_s1);
        stream->writeFloat(m_s2);
    }

    void configure() {
        m_logRatio = -math::fastlog(m_s2 / m_s1);
        m_time = 0;
        m_largeStepTime = 0;
        m_largeStep = false;
        m_sampleIndex = 0;
        m_sampleCount = 0;
    }

    void reset() {
        m_time = m_sampleIndex = m_largeStepTime = 0;
        m_u.clear();
    }
    
    void accept() {
        if (m_largeStep)
            m_largeStepTime = m_time;
        m_time++;
        m_backup.clear();
        m_sampleIndex = 0;
    }

    void reject() {
        for (size_t i = 0; i < m_backup.size(); ++i)
            m_u[m_backup[i].first] = m_backup[i].second;
        m_backup.clear();
        m_sampleIndex = 0;
    }

    Float primarySample(size_t i) {
        while (i >= m_u.size())
            m_u.push_back(SampleStruct(m_random->nextFloat()));

        if (m_u[i].modify < m_time) {
            if (m_largeStep) {
                m_backup.push_back(std::pair<size_t, SampleStruct>(i, m_u[i]));
                m_u[i].modify = m_time;
                m_u[i].value = m_random->nextFloat();
            } else {
                if (m_u[i].modify < m_largeStepTime) {
                    m_u[i].modify = m_largeStepTime;
                    m_u[i].value = m_random->nextFloat();
                }

                while (m_u[i].modify + 1 < m_time) {
                    m_u[i].value = mutate(m_u[i].value);
                    m_u[i].modify++;
                }

                m_backup.push_back(std::pair<size_t, SampleStruct>(i, m_u[i]));

                m_u[i].value = mutate(m_u[i].value);
                m_u[i].modify++;
            }
        }

        return m_u[i].value;
    }

    ref<Sampler> clone() {
        ref<PSSMLTSamplerBase> sampler = new PSSMLTSamplerBase(this);
        sampler->m_sampleCount = m_sampleCount;
        sampler->m_sampleIndex = m_sampleIndex;
        sampler->m_random = new Random(m_random);
        return sampler.get();
    }

    Float next1D() {
        return primarySample(m_sampleIndex++);
    }

    Point2 next2D() {
        /// Enforce a specific order of evaluation
        Float value1 = primarySample(m_sampleIndex++);
        Float value2 = primarySample(m_sampleIndex++);
        return Point2(value1, value2);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "PSSMLTSampler[" << endl
                << "  sampleCount = " << m_sampleCount << endl
                << "]";
        return oss.str();
    }

protected:
    /// Virtual destructor
    virtual ~PSSMLTSamplerBase() {
    }
protected:

    struct SampleStruct {
        Float value;
        size_t modify;

        inline SampleStruct(Float value) : value(value), modify(0) {
        }
    };

    ref<Random> m_random;
    Float m_s1, m_s2, m_logRatio;
    bool m_largeStep;
    std::vector<std::pair<size_t, SampleStruct> > m_backup;
    std::vector<SampleStruct> m_u;
    size_t m_time, m_largeStepTime;
    Float m_probLargeStep;
};

MTS_NAMESPACE_END

#endif /* PSSMLT_SAMPLER_BASE_H */


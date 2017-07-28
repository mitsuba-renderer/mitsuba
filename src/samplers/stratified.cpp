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

#include <mitsuba/render/sampler.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{stratified}{Stratified sampler}
 * \order{2}
 * \parameters{
 *     \parameter{sampleCount}{\Integer}{
 *       Number of samples per pixel; should be a perfect square
 *       (e.g. 1, 4, 9, 16, 25, etc.), or it will be rounded up to the
 *       next one \default{4}
 *     }
 *     \parameter{dimension}{\Integer}{
 *       Effective dimension, up to which stratified samples are provided. The
 *       number here is to be interpreted as the number of subsequent 1D or 2D sample
 *       requests that can be satisfied using ``good'' samples. Higher high values
 *       increase both storage and computational costs.
 *       \default{4}
 *     }
 * }
 * \renderings{
 *     \unframedrendering{A projection of the first 1024 points
 *     onto the first two dimensions.}{sampler_stratified}
 *     \unframedrendering{The same samples shown together with the
 *     underlying strata for illustrative purposes}{sampler_stratified_strata}
 * }
 *
 * The stratified sample generator divides the domain into a discrete number
 * of strata and produces a sample within each one of them. This generally leads to less
 * sample clumping when compared to the independent sampler, as well as better
 * convergence. Due to internal storage costs, stratified samples are only provided up to a
 * certain dimension, after which independent sampling takes over.
 *
 * Like the \pluginref{independent} sampler, multicore and network renderings
 * will generally produce different images in subsequent runs due to the nondeterminism
 * introduced by the operating system scheduler.
 */
class StratifiedSampler : public Sampler {
public:
    StratifiedSampler() : Sampler(Properties()) { }

    StratifiedSampler(const Properties &props) : Sampler(props) {
        /* Sample count (will be rounded up to the next perfect square) */
        size_t desiredSampleCount = props.getSize("sampleCount", 4);

        size_t i = 1;
        while (i * i < desiredSampleCount)
            ++i;
        m_sampleCount = i*i;

        if (m_sampleCount != desiredSampleCount) {
            Log(EWarn, "Sample count should be a perfect square -- rounding to "
                    SIZE_T_FMT, m_sampleCount);
        }

        m_resolution = (int) i;

        /* Dimension, up to which which stratified samples are guaranteed to be available. */
        m_maxDimension = props.getInteger("dimension", 4);

        m_sampleCount = m_resolution*m_resolution;
        m_permutations1D = new uint32_t*[m_maxDimension];
        m_permutations2D = new uint32_t*[m_maxDimension];

        for (int i=0; i<m_maxDimension; i++) {
            m_permutations1D[i] = new uint32_t[m_sampleCount];
            m_permutations2D[i] = new uint32_t[m_sampleCount];
        }

        m_invResolution = 1 / (Float) m_resolution;
        m_invResolutionSquare = 1 / (Float) m_sampleCount;
        m_random = new Random();
    }

    StratifiedSampler(Stream *stream, InstanceManager *manager)
     : Sampler(stream, manager) {
        m_maxDimension = stream->readInt();
        m_resolution = stream->readInt();
        m_random = static_cast<Random *>(manager->getInstance(stream));
        m_permutations1D = new uint32_t*[m_maxDimension];
        m_permutations2D = new uint32_t*[m_maxDimension];
        for (int i=0; i<m_maxDimension; i++) {
            m_permutations1D[i] = new uint32_t[m_sampleCount];
            m_permutations2D[i] = new uint32_t[m_sampleCount];
        }
        m_invResolution = 1.0f / m_resolution;
        m_invResolutionSquare = 1.0f / m_sampleCount;
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Sampler::serialize(stream, manager);

        stream->writeInt(m_maxDimension);
        stream->writeInt(m_resolution);
        manager->serialize(stream, m_random.get());
    }

    virtual ~StratifiedSampler() {
        for (int i=0; i<m_maxDimension; i++) {
            delete[] m_permutations1D[i];
            delete[] m_permutations2D[i];
        }
        delete[] m_permutations1D;
        delete[] m_permutations2D;
    }

    ref<Sampler> clone() {
        ref<StratifiedSampler> sampler = new StratifiedSampler();
        sampler->m_sampleCount = m_sampleCount;
        sampler->m_maxDimension = m_maxDimension;
        sampler->m_resolution = m_resolution;
        sampler->m_invResolution = m_invResolution;
        sampler->m_invResolutionSquare = m_invResolutionSquare;
        sampler->m_random = new Random(m_random);

        sampler->m_permutations1D = new uint32_t*[m_maxDimension];
        sampler->m_permutations2D = new uint32_t*[m_maxDimension];
        for (int i=0; i<m_maxDimension; i++) {
            sampler->m_permutations1D[i] = new uint32_t[m_sampleCount];
            sampler->m_permutations2D[i] = new uint32_t[m_sampleCount];
        }
        for (size_t i=0; i<m_req1D.size(); ++i)
            sampler->request1DArray(m_req1D[i]);
        for (size_t i=0; i<m_req2D.size(); ++i)
            sampler->request2DArray(m_req2D[i]);
        return sampler.get();
    }

    void generate(const Point2i &) {
        for (int i=0; i<m_maxDimension; i++) {
            for (size_t j=0; j<m_sampleCount; j++)
                m_permutations1D[i][j] = (uint32_t) j;
            m_random->shuffle(&m_permutations1D[i][0], &m_permutations1D[i][m_sampleCount]);

            for (size_t j=0; j<m_sampleCount; j++)
                m_permutations2D[i][j] = (uint32_t) j;
            m_random->shuffle(&m_permutations2D[i][0], &m_permutations2D[i][m_sampleCount]);
        }

        for (size_t i=0; i<m_req1D.size(); i++)
            latinHypercube(m_random, m_sampleArrays1D[i], m_req1D[i] * m_sampleCount, 1);
        for (size_t i=0; i<m_req2D.size(); i++)
            latinHypercube(m_random, reinterpret_cast<Float *>(m_sampleArrays2D[i]),
                m_req2D[i] * m_sampleCount, 2);

        m_sampleIndex = 0;
        m_dimension1D = m_dimension2D = 0;
        m_dimension1DArray = m_dimension2DArray = 0;
    }

    void setSampleIndex(size_t sampleIndex) {
        m_sampleIndex = sampleIndex;
        m_dimension1D = m_dimension2D = 0;
        m_dimension1DArray = m_dimension2DArray = 0;
    }

    void advance() {
        m_sampleIndex++;
        m_dimension1D = m_dimension2D = 0;
        m_dimension1DArray = m_dimension2DArray = 0;
    }

    Float next1D() {
        Assert(m_sampleIndex < m_sampleCount);
        if (m_dimension1D < m_maxDimension) {
            int k = m_permutations1D[m_dimension1D++][m_sampleIndex];
            return (k + m_random->nextFloat()) * m_invResolutionSquare;
        } else {
            return m_random->nextFloat();
        }
    }

    Point2 next2D() {
        Assert(m_sampleIndex < m_sampleCount);
        if (m_dimension2D < m_maxDimension) {
            int k = m_permutations2D[m_dimension2D++][m_sampleIndex];
            int x = k % m_resolution;
            int y = k / m_resolution;
            return Point2(
                (x + m_random->nextFloat()) * m_invResolution,
                (y + m_random->nextFloat()) * m_invResolution
            );
        } else {
            return Point2(
                m_random->nextFloat(),
                m_random->nextFloat()
            );
        }
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "StratifiedSampler[" << endl
            << "  resolution = " << m_resolution << "," << endl
            << "  sampleCount = " << m_sampleCount << "," << endl
            << "  dimension = " << m_maxDimension << "," << endl
            << "  sampleIndex = " << m_sampleIndex << "," << endl
            << "  dimension = " << m_maxDimension << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    ref<Random> m_random;
    int m_resolution;
    int m_maxDimension;
    Float m_invResolution, m_invResolutionSquare;
    uint32_t **m_permutations1D, **m_permutations2D;
    int m_dimension1D, m_dimension2D;
};

MTS_IMPLEMENT_CLASS_S(StratifiedSampler, false, Sampler)
MTS_EXPORT_PLUGIN(StratifiedSampler, "Stratified sampler");
MTS_NAMESPACE_END

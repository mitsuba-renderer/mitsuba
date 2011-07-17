/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

/**
 * Stratified sample generator. Given a resolution $R$ and a depth $D$, it 
 * generates $R*R$ samples, each of which can be queried for up to $D$ 1- 
 * or 2-dimensional vectors by an integrator. The returned 1D/2D-vectors 
 * of a particular depth have the property of being stratified over all 
 * $R*R$ samples. When the maximum depth is exceeded, independent sampling 
 * takes over.
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

		/* Depth, up to which which stratified samples are guaranteed to be available. */
		m_depth = props.getInteger("depth", 3);

		m_sampleCount = m_resolution*m_resolution;
		m_permutations1D = new uint32_t*[m_depth];
		m_permutations2D = new uint32_t*[m_depth];

		for (int i=0; i<m_depth; i++) {
			m_permutations1D[i] = new uint32_t[m_sampleCount];
			m_permutations2D[i] = new uint32_t[m_sampleCount];
		}

		m_invResolution = 1 / (Float) m_resolution;
		m_invResolutionSquare = 1 / (Float) m_sampleCount;
		m_random = new Random();
	}

	StratifiedSampler(Stream *stream, InstanceManager *manager) 
	 : Sampler(stream, manager) {
		m_depth = stream->readInt();
		m_resolution = stream->readInt();
		m_random = static_cast<Random *>(manager->getInstance(stream));
		m_permutations1D = new uint32_t*[m_depth];
		m_permutations2D = new uint32_t*[m_depth];
		for (int i=0; i<m_depth; i++) {
			m_permutations1D[i] = new uint32_t[m_sampleCount];
			m_permutations2D[i] = new uint32_t[m_sampleCount];
		}
		m_invResolution = 1.0f / m_resolution;
		m_invResolutionSquare = 1.0f / m_sampleCount;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Sampler::serialize(stream, manager);

		stream->writeInt(m_depth);
		stream->writeInt(m_resolution);
		manager->serialize(stream, m_random.get());
	}

	virtual ~StratifiedSampler() {
		for (int i=0; i<m_depth; i++) {
			delete[] m_permutations1D[i];
			delete[] m_permutations2D[i];
		}
		delete[] m_permutations1D;
		delete[] m_permutations2D;
	}

	ref<Sampler> clone() {
		ref<StratifiedSampler> sampler = new StratifiedSampler();
		sampler->m_sampleCount = m_sampleCount;
		sampler->m_depth = m_depth;
		sampler->m_resolution = m_resolution;
		sampler->m_invResolution = m_invResolution;
		sampler->m_invResolutionSquare = m_invResolutionSquare;
		sampler->m_random = new Random(m_random);

		sampler->m_permutations1D = new uint32_t*[m_depth];
		sampler->m_permutations2D = new uint32_t*[m_depth];
		for (int i=0; i<m_depth; i++) {
			sampler->m_permutations1D[i] = new uint32_t[m_sampleCount];
			sampler->m_permutations2D[i] = new uint32_t[m_sampleCount];
		}
		for (size_t i=0; i<m_req1D.size(); ++i)
			sampler->request2DArray(m_req1D[i]);
		for (size_t i=0; i<m_req2D.size(); ++i)
			sampler->request2DArray(m_req2D[i]);
		return sampler.get();
	}

	void generate() {
		for (int i=0; i<m_depth; i++) {
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
		m_sampleDepth1D = m_sampleDepth2D = 0;
		m_sampleDepth1DArray = m_sampleDepth2DArray = 0;
	}

	void setSampleIndex(size_t sampleIndex) {
		m_sampleIndex = sampleIndex;
		m_sampleDepth1D = m_sampleDepth2D = 0;
		m_sampleDepth1DArray = m_sampleDepth2DArray = 0;
	}

	void advance() {
		m_sampleIndex++;
		m_sampleDepth1D = m_sampleDepth2D = 0;
		m_sampleDepth1DArray = m_sampleDepth2DArray = 0;
	}

	Float next1D() {
		Assert(m_sampleIndex < m_sampleCount);
		if (m_sampleDepth1D < m_depth) {
			int k = m_permutations1D[m_sampleDepth1D++][m_sampleIndex];
			return (k + m_random->nextFloat()) * m_invResolutionSquare;
		} else {
			return m_random->nextFloat();
		}
	}

	Point2 next2D() {
		Assert(m_sampleIndex < m_sampleCount);
		if (m_sampleDepth2D < m_depth) {
			int k = m_permutations2D[m_sampleDepth2D++][m_sampleIndex];
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

	Float independent1D() {
		return m_random->nextFloat();
	}

	Point2 independent2D() {
		return Point2(
			m_random->nextFloat(),
			m_random->nextFloat()
		);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "StratifiedSampler[" << endl
			<< "  resolution = " << m_resolution << "," << endl
			<< "  sampleCount = " << m_sampleCount << "," << endl
			<< "  depth = " << m_depth << "," << endl
			<< "  sampleIndex = " << m_sampleIndex << "," << endl
			<< "  sampleDepth = " << m_depth << endl
			<< "]"; 
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<Random> m_random;
	int m_resolution;
	int m_depth;
	Float m_invResolution, m_invResolutionSquare;
	uint32_t **m_permutations1D, **m_permutations2D;
	int m_sampleDepth1D, m_sampleDepth2D;
};

MTS_IMPLEMENT_CLASS_S(StratifiedSampler, false, Sampler)
MTS_EXPORT_PLUGIN(StratifiedSampler, "Stratified sampling");
MTS_NAMESPACE_END

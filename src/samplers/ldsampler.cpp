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
 * Adapted version of the low discrepancy sampler in PBRT.
 * Provides samples up to a specified depth, after which independent 
 * sampling takes over.
 */
class LowDiscrepancySampler : public Sampler {
public:
	LowDiscrepancySampler() : Sampler(Properties()) { }

	LowDiscrepancySampler(const Properties &props) : Sampler(props) {
		/* Sample count (will be rounded up to the next power of two) */
		m_sampleCount = props.getSize("sampleCount", 4);

		/* Depth, up to which which low discrepancy samples are guaranteed to be available. */
		m_depth = props.getInteger("depth", 3);

		if (!isPowerOfTwo(m_sampleCount)) {
			m_sampleCount = roundToPowerOfTwo(m_sampleCount);
			Log(EWarn, "Sample count should be a power of two -- rounding to "
					SIZE_T_FMT, m_sampleCount);
		}

		m_samples1D = new Float*[m_depth];
		m_samples2D = new Point2*[m_depth];

		for (int i=0; i<m_depth; i++) {
			m_samples1D[i] = new Float[m_sampleCount];
			m_samples2D[i] = new Point2[m_sampleCount];
		}

		m_random = new Random();
	}

	LowDiscrepancySampler(Stream *stream, InstanceManager *manager) 
	 : Sampler(stream, manager) {
		m_random = static_cast<Random *>(manager->getInstance(stream));
		m_depth = stream->readInt();

		m_samples1D = new Float*[m_depth];
		m_samples2D = new Point2*[m_depth];
		for (int i=0; i<m_depth; i++) {
			m_samples1D[i] = new Float[(size_t) m_sampleCount];
			m_samples2D[i] = new Point2[(size_t) m_sampleCount];
		}
	}

	virtual ~LowDiscrepancySampler() {
		for (int i=0; i<m_depth; i++) {
			delete[] m_samples1D[i];
			delete[] m_samples2D[i];
		}
		delete[] m_samples1D;
		delete[] m_samples2D;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Sampler::serialize(stream, manager);
		manager->serialize(stream, m_random.get());
		stream->writeInt(m_depth);
	}

	ref<Sampler> clone() {
		ref<LowDiscrepancySampler> sampler = new LowDiscrepancySampler();

		sampler->m_sampleCount = m_sampleCount;
		sampler->m_depth = m_depth;
		sampler->m_random = new Random(m_random);
		sampler->m_samples1D = new Float*[m_depth];
		sampler->m_samples2D = new Point2*[m_depth];
		for (int i=0; i<m_depth; i++) {
			sampler->m_samples1D[i] = new Float[m_sampleCount];
			sampler->m_samples2D[i] = new Point2[m_sampleCount];
		}
		for (size_t i=0; i<m_req1D.size(); ++i)
			sampler->request2DArray(m_req1D[i]);
		for (size_t i=0; i<m_req2D.size(); ++i)
			sampler->request2DArray(m_req2D[i]);

		return sampler.get();
	}

	inline Float vanDerCorput(uint32_t n, uint32_t scramble) {
		n = (n << 16) | (n >> 16);
		n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
		n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
		n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
		n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
		n ^= scramble;
		return (Float) n / (Float) 0x100000000LL;
	}

	inline Float sobol2(uint32_t n, uint32_t scramble) {
		for (uint32_t v = 1 << 31; n != 0; n >>= 1, v ^= v >> 1)
			if (n & 0x1) scramble ^= v;
		return (Float) scramble / (Float) 0x100000000LL;
	}

	inline void sample02(uint32_t n, uint32_t scramble[2], Point2 &sample) {
		sample.x = vanDerCorput(n, scramble[0]);
		sample.y = sobol2(n, scramble[1]);
	}

	inline void generate1D(Float *samples, size_t sampleCount) {
		uint32_t scramble = m_random->nextULong() & 0xFFFFFFFF;
		for (size_t i = 0; i < sampleCount; ++i)
			samples[i] = vanDerCorput((uint32_t) i, scramble);
		m_random->shuffle(samples, samples + sampleCount);
	}

	inline void generate2D(Point2 *samples, size_t sampleCount) {
		union {
			uint64_t qword;
			uint32_t dword[2];
		} scramble;
		scramble.qword = m_random->nextULong();
		for (size_t i = 0; i < sampleCount; ++i)
			sample02((uint32_t) i, scramble.dword, samples[i]);
		m_random->shuffle(samples, samples + sampleCount);
	}

	void generate() {
		for (int i=0; i<m_depth; ++i) {
			generate1D(m_samples1D[i], m_sampleCount);
			generate2D(m_samples2D[i], m_sampleCount);
		}
		
		for (size_t i=0; i<m_req1D.size(); i++)
			generate1D(m_sampleArrays1D[i], m_sampleCount * m_req1D[i]);

		for (size_t i=0; i<m_req2D.size(); i++)
			generate2D(m_sampleArrays2D[i], m_sampleCount * m_req2D[i]);

		m_sampleIndex = 0;
		m_sampleDepth1D = m_sampleDepth2D = 0;
		m_sampleDepth1DArray = m_sampleDepth2DArray = 0;
	}

	void advance() {
		m_sampleIndex++;
		m_sampleDepth1D = m_sampleDepth2D = 0;
		m_sampleDepth1DArray = m_sampleDepth2DArray = 0;
	}
	
	void setSampleIndex(size_t sampleIndex) {
		m_sampleIndex = sampleIndex;
		m_sampleDepth1D = m_sampleDepth2D = 0;
		m_sampleDepth1DArray = m_sampleDepth2DArray = 0;
	}

	Float next1D() {
		Assert(m_sampleIndex < m_sampleCount);
		if (m_sampleDepth1D < m_depth)
			return m_samples1D[m_sampleDepth1D++][m_sampleIndex];
		else
			return m_random->nextFloat();
	}

	Point2 next2D() {
		Assert(m_sampleIndex < m_sampleCount);
		if (m_sampleDepth2D < m_depth)
			return m_samples2D[m_sampleDepth2D++][m_sampleIndex];
		else 
			return Point2(m_random->nextFloat(), m_random->nextFloat());
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
		oss << "LowDiscrepancySampler[" << endl
			<< "  sampleCount = " << m_sampleCount << "," << endl
			<< "  depth = " << m_depth << endl
			<< "]"; 
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<Random> m_random;
	int m_depth;
	int m_sampleDepth1D;
	int m_sampleDepth2D;
	Float **m_samples1D;
	Point2 **m_samples2D;
};

MTS_IMPLEMENT_CLASS_S(LowDiscrepancySampler, false, Sampler)
MTS_EXPORT_PLUGIN(LowDiscrepancySampler, "Low discrepancy sampler");
MTS_NAMESPACE_END

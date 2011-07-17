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
 * Deterministic 2D sample generator based on the Hammersley sequence. Internally
 * uses a table of prime numbers to provide elements of the sequence up to a
 * depth of 1000.
 * Because of the high correlation amongst neighboring pixels, this
 * sampler, by itself, is not meant to be used as a source of random numbers 
 * for sample-based integrators such as <tt>direct</tt>, <tt>volpath</tt> etc. 
 */
class HammersleySequence : public Sampler {
public:
	HammersleySequence() : Sampler(Properties()) {
	}

	HammersleySequence(Stream *stream, InstanceManager *manager) 
	 : Sampler(stream, manager) {
		m_invSamplesPerPixel = 1.0f / m_sampleCount;
	}

	HammersleySequence(const Properties &props) : Sampler(props) {
		/* Number of samples per pixel when used with a sampling-based integrator */
		m_sampleCount = props.getSize("sampleCount", 1);
		m_invSamplesPerPixel = 1.0f / m_sampleCount;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Sampler::serialize(stream, manager);
	}

	ref<Sampler> clone() {
		ref<HammersleySequence> sampler = new HammersleySequence();
		sampler->m_sampleCount = m_sampleCount;
		sampler->m_invSamplesPerPixel = m_invSamplesPerPixel;
		sampler->m_sampleIndex = m_sampleIndex;
		sampler->m_sampleDepth = m_sampleDepth;
		return sampler.get();
	}

	void generate() {
		m_sampleIndex = 0;
		m_sampleDepth = 0;
	}

	void advance() {
		m_sampleIndex++;
		m_sampleDepth = 0;
	}

	void setSampleIndex(size_t sampleIndex) {
		m_sampleDepth = 0;
		m_sampleIndex = sampleIndex;
	}

	inline Float nextValue() {
		if (m_sampleDepth == 0) {
			m_sampleDepth++;
			return m_sampleIndex * m_invSamplesPerPixel;
		} else {
			return radicalInverse(primeTable[(m_sampleDepth++)-1], m_sampleIndex);
		}
	}

	Float next1D() {
		Assert(m_sampleIndex < m_sampleCount);
		if (m_sampleDepth >= primeTableSize)
			Log(EError, "Lookup depth exceeds the prime number table size!");
		return nextValue();
	}

	Point2 next2D() {
		Assert(m_sampleIndex < m_sampleCount);
		if (m_sampleDepth >= primeTableSize)
			Log(EError, "Lookup depth exceeds the prime number table size!");
		return Point2(nextValue(), nextValue());
	}

	Float independent1D() {
		Log(EError, "independent1D() - unsupported (this is a purely deterministic sampler)!");
		return 0.0f;
	}

	Point2 independent2D() {
		Log(EError, "independent2D() - unsupported (this is a purely deterministic sampler)!");
		return Point2();
	}

	void request1DArray(unsigned int size) {
		Log(EError, "request1DArray() is not supported by QMC samplers!");
	}
	
	void request2DArray(unsigned int size) {
		Log(EError, "request2DArray() is not supported by QMC samplers!");
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HammersleySequence[" << endl
			<< "  sampleCount = " << m_sampleCount << "," << endl
			<< "  sampleIndex = " << m_sampleIndex << "," << endl
			<< "  sampleDepth = " << m_sampleDepth << endl
			<< "]"; 
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	int m_sampleDepth;
	Float m_invSamplesPerPixel;
};

MTS_IMPLEMENT_CLASS_S(HammersleySequence, false, Sampler)
MTS_EXPORT_PLUGIN(HammersleySequence, "Hammersley sequence");
MTS_NAMESPACE_END

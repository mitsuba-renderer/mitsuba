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
 * Independent sampling - returns independent uniformly distributed
 * random numbers on <tt>[0, 1)x[0, 1)</tt>.
 */
class IndependentSampler : public Sampler {
public:
	IndependentSampler() : Sampler(Properties()) {
	}

	IndependentSampler(const Properties &props) : Sampler(props) {
		/* Number of samples per pixel when used with a sampling-based integrator */
		m_sampleCount = props.getSize("sampleCount", 4);
		m_random = new Random();
	}

	IndependentSampler(Stream *stream, InstanceManager *manager) 
	 : Sampler(stream, manager) {
		m_random = static_cast<Random *>(manager->getInstance(stream));
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Sampler::serialize(stream, manager);
		manager->serialize(stream, m_random.get());
	}

	ref<Sampler> clone() {
		ref<IndependentSampler> sampler = new IndependentSampler();
		sampler->m_sampleCount = m_sampleCount;
		sampler->m_random = new Random(m_random);
		for (size_t i=0; i<m_req1D.size(); ++i)
			sampler->request2DArray(m_req1D[i]);
		for (size_t i=0; i<m_req2D.size(); ++i)
			sampler->request2DArray(m_req2D[i]);
		return sampler.get();
	}

	void generate() {
		for (size_t i=0; i<m_req1D.size(); i++)
			for (size_t j=0; j<m_sampleCount * m_req1D[i]; ++j)
				m_sampleArrays1D[i][j] = m_random->nextFloat();
		for (size_t i=0; i<m_req2D.size(); i++)
			for (size_t j=0; j<m_sampleCount * m_req2D[i]; ++j)
				m_sampleArrays2D[i][j] = Point2(
					m_random->nextFloat(),
					m_random->nextFloat());
		m_sampleIndex = 0;
		m_sampleDepth1DArray = m_sampleDepth2DArray = 0;
	}
	
	Float next1D() {
		return m_random->nextFloat();
	}

	Point2 next2D() {
		/// Enforce a specific order of evaluation
		Float value1 = m_random->nextFloat();
		Float value2 = m_random->nextFloat();
		return Point2(value1, value2);
	}
	
	Float independent1D() {
		return m_random->nextFloat();
	}

	Point2 independent2D() {
		/// Enforce a specific order of evaluation
		Float value1 = m_random->nextFloat();
		Float value2 = m_random->nextFloat();
		return Point2(value1, value2);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "IndependentSampler[" << endl
			<< "  sampleCount = " << m_sampleCount << endl
			<< "]"; 
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<Random> m_random;
};

MTS_IMPLEMENT_CLASS_S(IndependentSampler, false, Sampler)
MTS_EXPORT_PLUGIN(IndependentSampler, "Independent sampling");
MTS_NAMESPACE_END

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
#include <mitsuba/core/qmc.h>
#include "sobolseq.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{sobol}{Sobol QMC sampler}
 * \order{6}
 * \parameters{
 *     \parameter{sampleCount}{\Integer}{
 *       Number of samples per pixel \default{4}
 *     }
 *     \parameter{scramble}{\Integer}{
 *       This parameter can be used to set a scramble value to break up temporally coherent
 *       noise patterns. For stills, this is irrelevant. When rendering
 *       an animation, simply set it to the current frame index. \default{0}
 *     }
 * }
 *
 * \renderings{
 *     \unframedrendering{A projection of the first 1024 points
 *     onto the first two dimensions.}{sampler_sobol_nonscrambled_0}
 *     \unframedrendering{A projection of the first 1024 points
 *     onto the 32 and 33th dimension.}{sampler_sobol_nonscrambled_32}
 * }
 *
 * This plugin implements a Quasi-Monte Carlo (QMC) sample generator based on the
 * Sobol sequence. QMC number sequences are designed to reduce sample clumping
 * across integration dimensions, which can lead to a higher order of
 * convergence in renderings. Because of the deterministic character of the samples,
 * errors will manifest as grid or moir\'e patterns rather than random noise, but
 * these diminish as the number of samples is increased.
 *
 * The Sobol sequence in particular provides a relatively good point set that can
 * be computed extremely efficiently. One downside is the susceptibility to pattern artifacts
 * in the generated image. To minimize these artifacts, it is advisable to use
 * a number of samples per pixel that is a power of two.
 *
 * Because everything that happens inside this sampler is completely
 * deterministic and independent of operating system scheduling behavior, subsequent
 * runs of Mitsuba will always compute the same image, and this even holds when rendering
 * with multiple threads and/or machines.

 * The plugin relies on a fast implementation of the Sobol sequence by Leonhard
 * Gr\"unschlo\ss\ using direction numbers provided by Joe and Kuo
 * \cite{Joe2008Constructing}.
 * These direction numbers are given up to a dimension of 1024. Because of this
 * upper bound, the maximum path depth of the integrator must be limited (e.g. to 100), or
 * rendering might fail with the following error message: \emph{Lookup dimension
 * exceeds the direction number table size! You may have to reduce the 'maxDepth'
 * parameter of your integrator}.
 *
 * Note that this sampler generates a $(0,2)$-sequence in the first two dimensions,
 * and therefore the point plot shown in (a) happens to match the corresponding plots of
 * \pluginref{ldsampler}. In higher dimensions, however, they behave rather differently.
 *
 * When this sampler is used to perform parallel block-based renderings,
 * the sequence is internally enumerated using a scheme proposed and implemented
 * by Gr\"unschlo\ss\ et al. \cite{Grunschloss2010Enumerating}.
 * \remarks{
 *   \item This sampler is incompatible with Metropolis Light Transport (all variants).
 * }
 */
class SobolSampler : public Sampler {
public:
	SobolSampler() : Sampler(Properties()) { }

	SobolSampler(const Properties &props) : Sampler(props) {
		/* Number of samples per pixel when used with a sampling-based integrator */
		m_sampleCount = props.getSize("sampleCount", 4);

		/* Scramble value, which can be used to break up temporally coherent
		   noise patterns when rendering the frames of an animation. */
		m_scramble = (uint64_t) props.getSize("scramble", 0);

		/* When scrambling was requested, run TEA to turn the frame
		   number into a pseudorandom uniformly distributed 64-bit integer */
		if (m_scramble) {
			union {
				uint64_t ui64;
				uint32_t v[2];
			} u = {m_scramble};
			m_scramble = sampleTEA(u.v[0], u.v[1]);
		}

		m_resolution = 1; m_logResolution = 0;
		m_arrayStartDim = m_arrayEndDim = 5;
		m_pixelPosition = Point2i(0);
	}

	SobolSampler(Stream *stream, InstanceManager *manager)
	 : Sampler(stream, manager) {
		m_scramble = stream->readULong();
		m_resolution = stream->readFloat();
		m_logResolution = stream->readUInt();
		m_arrayStartDim = stream->readUInt();
		m_arrayEndDim = stream->readUInt();
		m_pixelPosition = Point2i(0);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Sampler::serialize(stream, manager);
		stream->writeULong(m_scramble);
		stream->writeFloat(m_resolution);
		stream->writeUInt(m_logResolution);
		stream->writeUInt(m_arrayStartDim);
		stream->writeUInt(m_arrayEndDim);
	}

	ref<Sampler> clone() {
		ref<SobolSampler> sampler = new SobolSampler();
		sampler->m_sampleCount = m_sampleCount;
		sampler->m_sampleIndex = m_sampleIndex;
		sampler->m_sobolSampleIndex = m_sobolSampleIndex;
		sampler->m_dimension = m_dimension;
		sampler->m_scramble = m_scramble;
		sampler->m_resolution = m_resolution;
		sampler->m_logResolution = m_logResolution;
		sampler->m_pixelPosition = m_pixelPosition;
		sampler->m_arrayStartDim = m_arrayStartDim;
		sampler->m_arrayEndDim = m_arrayEndDim;
		for (size_t i=0; i<m_req1D.size(); ++i)
			sampler->request1DArray(m_req1D[i]);
		for (size_t i=0; i<m_req2D.size(); ++i)
			sampler->request2DArray(m_req2D[i]);
		return sampler.get();
	}

	void setFilmResolution(const Vector2i &res, bool bucketed) {
		if (!bucketed) {
			m_resolution = 1;
			m_logResolution = 0;
		} else {
			uint32_t resolution = math::roundToPowerOfTwo(
				(uint32_t) std::max(res.x, res.y));

			m_resolution = (Float) resolution;
			m_logResolution = math::log2i(resolution);
		}
	}

	template <typename Iterator> void shuffle(uint32_t seed, Iterator it1, Iterator it2) {
		uint32_t ctr = 0;
		for (Iterator it = it2 - 1; it > it1; --it) {
			size_t size = it - it1;
			size_t value = std::max((size_t) (size * (sampleTEA(seed, ctr++)
					/ (Float) std::numeric_limits<uint64_t>::max())), size-1);
			std::iter_swap(it, it1 + value);
		}
	}


	void generate(const Point2i &pos) {
		m_pixelPosition = pos;
		setSampleIndex(0);

		/* Dimensions reserved to sample array requests */
		m_arrayStartDim = 5;
		m_arrayEndDim = m_arrayStartDim +
				static_cast<uint32_t>(m_req1D.size() + 2 * m_req2D.size());

		uint32_t dim = m_arrayStartDim;
		for (size_t i=0; i<m_req1D.size(); i++) {
			for (size_t j=0; j<m_sampleCount * m_req1D[i]; ++j) {
				uint64_t idx = sobol::look_up(m_logResolution, (uint32_t) j, m_pixelPosition.x, m_pixelPosition.y, m_scramble);
				m_sampleArrays1D[i][j] = sobol::sample(idx, dim, m_scramble);
			}
			dim += 1;
		}

		for (size_t i=0; i<m_req2D.size(); i++) {
			for (size_t j=0; j<m_sampleCount * m_req2D[i]; ++j) {
				uint64_t idx = sobol::look_up(m_logResolution, (uint32_t) j, m_pixelPosition.x, m_pixelPosition.y, m_scramble);
				m_sampleArrays2D[i][j] = Point2(
					sobol::sample(idx, dim, m_scramble),
					sobol::sample(idx, dim+1, m_scramble));
			}
			dim += 2;
		}
	}

	void advance() {
		setSampleIndex(m_sampleIndex + 1);
	}

	void setSampleIndex(size_t sampleIndex) {
		m_dimension = 0;
		m_dimension1DArray = m_dimension2DArray = 0;
		m_sampleIndex = sampleIndex;

		if (m_logResolution > 1 && m_pixelPosition.x >= 0) {
			/* Find the next sample that is located in the current pixel */
			m_sobolSampleIndex = sobol::look_up(m_logResolution, (uint32_t) m_sampleIndex,
					m_pixelPosition.x, m_pixelPosition.y, m_scramble);
		} else {
			m_sobolSampleIndex = (uint64_t) m_sampleIndex;
		}
	}

	Float next1D() {
		/* Skip over dimensions that were reserved to arrays */
		if (m_dimension >= m_arrayStartDim && m_dimension < m_arrayEndDim)
			m_dimension = m_arrayEndDim;

		if (m_dimension >= sobol::Matrices::num_dimensions)
			Log(EError, "Lookup dimension exceeds the direction number table size! You "
				"may have to reduce the 'maxDepth' parameter of your integrator.");

		return sobol::sample(m_sobolSampleIndex, m_dimension++, m_scramble);
	}

	Point2 next2D() {
		Float value1, value2;

		/* Skip over dimensions that were reserved to arrays */
		if (m_dimension + 1 >= m_arrayStartDim && m_dimension < m_arrayEndDim)
			m_dimension = m_arrayEndDim;

		if (m_dimension + 1 >= sobol::Matrices::num_dimensions)
			Log(EError, "Lookup dimension exceeds the direction number table size! You "
				"may have to reduce the 'maxDepth' parameter of your integrator.");

		if (m_dimension == 0 && m_sobolSampleIndex != (uint64_t) m_sampleIndex) {
			value1 = sobol::sample(m_sobolSampleIndex, m_dimension++, m_scramble) * m_resolution - m_pixelPosition.x;
			value2 = sobol::sample(m_sobolSampleIndex, m_dimension++, m_scramble) * m_resolution - m_pixelPosition.y;
		} else {
			value1 = sobol::sample(m_sobolSampleIndex, m_dimension++, m_scramble);
			value2 = sobol::sample(m_sobolSampleIndex, m_dimension++, m_scramble);
		}

		return Point2(value1, value2);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SobolSampler[" << endl
			<< "  sampleCount = " << m_sampleCount << "," << endl
			<< "  sampleIndex = " << m_sampleIndex << "," << endl
			<< "  resolution = " << m_resolution << "," << endl
			<< "  scramble = " << m_scramble << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	uint32_t m_dimension;
	uint64_t m_scramble;
	uint64_t m_sobolSampleIndex;
	Float m_resolution;
	uint32_t m_logResolution;
	uint32_t m_arrayStartDim;
	uint32_t m_arrayEndDim;
	Point2i m_pixelPosition;
};

MTS_IMPLEMENT_CLASS_S(SobolSampler, false, Sampler)
MTS_EXPORT_PLUGIN(SobolSampler, "Sobol QMC sampler");
MTS_NAMESPACE_END

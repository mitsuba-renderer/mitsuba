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
#include <mitsuba/core/lock.h>
#include <mitsuba/core/qmc.h>
#include "faure.h"

/* This implementation limits the maximum pixel resolution and instead
   tiles a block of identical sequences across the screen. This is
   important to avoid running out of single precision bits very quickly ..
   The tile size should not be too small, or visible patterns will emerge */

#define MAX_RESOLUTION 128

MTS_NAMESPACE_BEGIN

/*!\plugin{hammersley}{Hammersley QMC sampler}
 * \order{5}
 * \parameters{
 *     \parameter{sampleCount}{\Integer}{
 *        Number of samples per pixel \default{4}
 *     }
 *     \parameter{scramble}{\Integer}{
 *        This plugin can operate in one of three scrambling modes:
 *        \begin{enumerate}[(i)]
 *        \item When set to \code{0}, the implementation will provide the standard Hammersley sequence.
 *
 *        \item When set to \code{-1}, the implementation will compute
 *        a scrambled variant of the Hammersley sequence based on permutations by
 *        Faure \cite{Faure1992Good}, which has better equidistribution properties
 *        in high dimensions.
 *
 *        \item When set to a value greater than one, a random permutation is chosen based
 *        on this number. This is useful to break up temporally coherent noise when rendering
 *        the frames of an animation --- in this case, simply set the parameter to the current frame index.
 *        \end{enumerate}
 *        Default: \code{-1}, i.e. use the Faure permutations. Note that permutations rely on a
 *        precomputed table that consumes approximately 7 MiB of additional memory at run time.
 *     }
 * }
  * \renderings{
 *     \unframedrendering{Projection of the first 1024 points
 *     of the Faure-scrambled sequence onto the first two dimensions.}{sampler_hammersley_0}
 *     \unframedrendering{Projection of the first 1024 points
 *     of the Faure-scrambled sequence onto the 32th and 33th dim.}{sampler_hammersley_32}
 * }
 * This plugin implements a Quasi-Monte Carlo (QMC) sample generator based on the
 * Hammersley sequence. QMC number sequences are designed to reduce sample clumping
 * across integration dimensions, which can lead to a higher order of
 * convergence in renderings. Because of the deterministic character of the samples,
 * errors will manifest as grid or moir\'e patterns rather than random noise, but
 * these diminish as the number of samples is increased.
 *
 * The Hammerlsey sequence is closely related to the Halton sequence and yields a very
 * high quality point set that is slightly more regular (and has lower discrepancy),
 * especially in the first few dimensions. As is the case with the Halton sequence,
 * the points should be scrambled to reduce patterns that manifest due to correlations
 * in higher dimensions. Please refer to the \pluginref{halton} page for more information
 * on how this works.
 *
 * Note that this sampler will cause odd-looking intermediate results when combined with rendering
 * techniques that trace paths starting at light source (e.g. \pluginref{ptracer})---these vanish
 * by the time the rendering process finishes.
 *
 * \remarks{
 *   \item This sampler is incompatible with Metropolis Light Transport (all variants).
 *   It interoperates poorly with Bidirectional Path Tracing and Energy Redistribution
 *   Path Tracing, hence these should not be used together. The \pluginref{sobol} QMC
 *   sequence is an alternative for the latter two cases, and \pluginref{ldsampler}
 *   works as well.
 * }
 */
class HammersleySampler : public Sampler {
public:
	HammersleySampler() : Sampler(Properties()) { }

	HammersleySampler(const Properties &props) : Sampler(props) {
		/* Number of samples per pixel */
		m_sampleCount = props.getSize("sampleCount", 4);

		/* Scramble value, which can be used to break up temporally coherent
		   noise patterns when rendering the frames of an animation. */
		m_scramble = props.getInteger("scramble", -1);

		setFilmResolution(Vector2i(1), false);

		m_arrayStartDim = m_arrayEndDim = 5;
	}

	HammersleySampler(Stream *stream, InstanceManager *manager)
	 : Sampler(stream, manager) {
		m_arrayStartDim = stream->readUInt();
		m_arrayEndDim = stream->readUInt();
		m_offset = stream->readULong();
		m_stride = stream->readULong();
		m_scramble = stream->readInt();
		m_logHeight = stream->readUInt();
		m_samplesPerBatch = stream->readSize();
		m_factor = stream->readFloat();
		m_resolution = Vector2i(stream);
		m_pixelPosition = Point2i(0);
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Sampler::serialize(stream, manager);
		stream->writeUInt(m_arrayStartDim);
		stream->writeUInt(m_arrayEndDim);
		stream->writeULong(m_offset);
		stream->writeULong(m_stride);
		stream->writeInt(m_scramble);
		stream->writeUInt(m_logHeight);
		stream->writeSize(m_samplesPerBatch);
		stream->writeFloat(m_factor);
		m_resolution.serialize(stream);
	}

	ref<Sampler> clone() {
		ref<HammersleySampler> sampler = new HammersleySampler();
		sampler->m_sampleCount = m_sampleCount;
		sampler->m_samplesPerBatch = m_samplesPerBatch;
		sampler->m_factor = m_factor;
		sampler->m_sampleIndex = m_sampleIndex;
		sampler->m_dimension = m_dimension;
		sampler->m_arrayStartDim = m_arrayStartDim;
		sampler->m_arrayEndDim = m_arrayEndDim;
		sampler->m_permutations = m_permutations;
		sampler->m_offset = m_offset;
		sampler->m_stride = m_stride;
		sampler->m_pixelPosition = m_pixelPosition;
		sampler->m_scramble = m_scramble;
		sampler->m_logHeight = m_logHeight;
		sampler->m_permutations = m_permutations;
		sampler->m_resolution = m_resolution;
		return sampler.get();
	}

	void configure() {
		Sampler::configure();
		if (m_scramble != 0) {
			/* Only create one set of permutations per address space, since this is costly
			   (taking about .5 sec and 7MB of memory on my machine) */
			LockGuard guard(m_globalPermutationsMutex);
			if (m_globalPermutations == NULL || m_globalPermutations->getScramble() != m_scramble)
				m_globalPermutations = new PermutationStorage(m_scramble);
			m_permutations = m_globalPermutations;
		}
	}

	/// Inverse integer version of the scrambled radical inverse function
	uint64_t inverseScrambledRadicalInverse(int base, uint64_t inverse, uint64_t digits, uint16_t *invPerm) {
		uint64_t index = 0;
		while (digits) {
			uint64_t digit = inverse % base;
			if (invPerm)
				digit = (uint64_t ) invPerm[digit];

			inverse /= base;
			index = index * base + digit;
			--digits;
		}
		return index;
	}

	void setFilmResolution(const Vector2i &res, bool blocked) {
		if (blocked) {
			/* Determine parameters of the space partition in the first two
			   dimensions. This is required to support blocked rendering. */
			for (int i=0; i<2; ++i)
				m_resolution[i] = std::min((uint32_t) MAX_RESOLUTION,
					math::roundToPowerOfTwo((uint32_t) res[i]));
			m_logHeight = math::log2i((uint32_t) m_resolution.y);

			m_samplesPerBatch = m_sampleCount;
			m_factor = (Float) 1.0f / (m_sampleCount *
				(size_t) m_resolution.x * (size_t) m_resolution.y);
			m_offset = 0;
			m_stride = m_resolution.y;
		} else {
			m_samplesPerBatch = m_sampleCount *
				(size_t) res.x * (size_t) res.y;
			m_factor = (Float) 1.0f / m_samplesPerBatch;
			m_resolution = Vector2i(1);
			m_offset = 0;
			m_stride = 1;
		}
		m_pixelPosition = Point2i(0);
	}

	void generate(const Point2i &pos) {
		/* Dimensions reserved to sample array requests */
		m_arrayStartDim = 5;
		m_arrayEndDim = m_arrayStartDim +
			static_cast<uint32_t>(m_req1D.size() + 2 * m_req2D.size());

		if (m_stride > 1) {
			m_pixelPosition = pos;
			m_pixelPosition.x %= MAX_RESOLUTION;
			m_pixelPosition.y %= MAX_RESOLUTION;

			m_offset = m_pixelPosition.x * m_resolution.y * m_sampleCount +
				inverseScrambledRadicalInverse(2, m_pixelPosition.y, m_logHeight,
				m_permutations.get() ? m_permutations->getInversePermutation(0) : NULL);
		}

		setSampleIndex(0);
	}

	void advance() {
		m_sampleIndex++;
		m_dimension = 0;
	}

	void setSampleIndex(size_t sampleIndex) {
		m_dimension = 0;
		m_sampleIndex = sampleIndex;
	}

	inline Float nextFloat(uint64_t idx) {
		uint32_t dim = m_dimension++;
		if (dim == 0)
			return idx * m_factor;
		else if (m_permutations != NULL)
			return scrambledRadicalInverseFast(dim-1, idx,
				m_permutations->getPermutation(dim-1));
		else
			return radicalInverseFast(dim-1, idx);
	}

	Float next1D() {
		/* Skip over dimensions that were reserved to arrays */
		if (m_dimension >= m_arrayStartDim && m_dimension < m_arrayEndDim)
			m_dimension = m_arrayEndDim;
		if (m_dimension >= primeTableSize)
			Log(EError, "Lookup dimension exceeds the prime number table size! "
				"You may have to reduce the 'maxDepth' parameter of your integrator.");
		if (m_sampleIndex >= m_samplesPerBatch)
			Log(EError, "Sample index exceeded the maximum count!");

		return nextFloat(m_offset + m_stride * m_sampleIndex);
	}

	Point2 next2D() {
		/* Skip over dimensions that were reserved to arrays */
		if (m_dimension + 1 >= m_arrayStartDim && m_dimension < m_arrayEndDim)
			m_dimension = m_arrayEndDim;
		if (m_dimension + 1 >= primeTableSize)
			Log(EError, "Lookup dimension exceeds the prime number table size! "
				"You may have to reduce the 'maxDepth' parameter of your integrator.");
		if (m_sampleIndex >= m_samplesPerBatch)
			Log(EError, "Sample index exceeded the maximum count!");

		uint64_t index = m_offset + m_stride * m_sampleIndex;

		Float value1, value2;
		if (m_dimension == 0) {
			value1 = nextFloat(index) * m_resolution.x - m_pixelPosition.x;
			value2 = nextFloat(index) * m_resolution.y - m_pixelPosition.y;
		} else {
			value1 = nextFloat(index);
			value2 = nextFloat(index);
		}

		return Point2(value1, value2);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HammersleySampler[" << endl
			<< "  sampleCount = " << m_sampleCount << "," << endl
			<< "  sampleIndex = " << m_sampleIndex << "," << endl
			<< "  scramble = " << m_scramble << endl
			<< "]";
		return oss.str();
	}

	void request1DArray(size_t size) {
		Log(EError, "request1DArray(): Not supported for the Hammersley QMC sequence!");
	}

	void request2DArray(size_t size) {
		Log(EError, "request2DArray(): Not supported for the Hammersley QMC sequence! "
			"If you used this sampler together with the 'direct' or 'ao' integrators, you "
			"will have to set the 'shadingSamples' parameter to 1 to use this sampler.");
	}

	MTS_DECLARE_CLASS()
private:
	uint32_t m_dimension;
	uint32_t m_arrayStartDim;
	uint32_t m_arrayEndDim;
	int m_scramble;
	Float m_factor;

	/* Faure permutation */
	ref<const PermutationStorage> m_permutations;
	static ref<const PermutationStorage> m_globalPermutations;
	static ref<Mutex> m_globalPermutationsMutex;

	/* Tiling-related */
	uint64_t m_offset, m_stride;
	Vector2i m_resolution;
	Point2i m_pixelPosition;
	uint32_t m_logHeight;
	size_t m_samplesPerBatch;
};

ref<Mutex> HammersleySampler::m_globalPermutationsMutex = new Mutex();
ref<const PermutationStorage> HammersleySampler::m_globalPermutations = NULL;

MTS_IMPLEMENT_CLASS_S(HammersleySampler, false, Sampler)
MTS_EXPORT_PLUGIN(HammersleySampler, "Hammersley QMC sampler");
MTS_NAMESPACE_END

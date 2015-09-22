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
   tiles a block of identical Halton sequences across the screen. This is
   important to avoid running out of single precision bits very quickly ..
   The tile size should not be too small, or visible patterns will emerge */

#define MAX_RESOLUTION 128

MTS_NAMESPACE_BEGIN

/*!\plugin{halton}{Halton QMC sampler}
 * \order{4}
 * \parameters{
 *     \parameter{sampleCount}{\Integer}{
 *        Number of samples per pixel \default{4}
 *     }
 *     \parameter{scramble}{\Integer}{
 *        This plugin can operate in one of three scrambling modes:
 *        \begin{enumerate}[(i)]
 *        \item When set to \code{0}, the implementation will provide the standard Halton sequence.
 *
 *        \item When set to \code{-1}, the implementation will compute
 *        a scrambled variant of the Halton sequence based on permutations by
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
 *     of the Faure-scrambled Halton seq. onto the first two dimensions.}{sampler_halton_scrambled_0}
 *     \unframedrendering{Projection of the first 1024 points
 *     of the Faure-scrambled Halton seq. onto the 32th and 33th dim.}{sampler_halton_scrambled_32}
 * }
 * This plugin implements a Quasi-Monte Carlo (QMC) sample generator based on the
 * Halton sequence. QMC number sequences are designed to reduce sample clumping
 * across integration dimensions, which can lead to a higher order of
 * convergence in renderings. Because of the deterministic character of the samples,
 * errors will manifest as grid or moir\'e patterns rather than random noise, but
 * these diminish as the number of samples is increased.
 *
 * The Halton sequence in particular provides a very high quality point set that unfortunately
 * becomes increasingly correlated in higher dimensions. To ameliorate this problem, the Halton
 * points are usually combined with a scrambling permutation, and this is also the default.
 * Because everything that happens inside this sampler is completely deterministic and
 * independent of operating system scheduling behavior, subsequent runs of Mitsuba will always
 * compute the same image, and this even holds when rendering with multiple threads
 * and/or machines.
  * \renderings{
 *     \unframedrendering{A projection of the first 1024 points
 *     of the \emph{original} Halton sequence onto the first two dimensions, obtained by
 *     setting \code{scramble=0}}{sampler_halton_nonscrambled_0}
 *     \unframedrendering{A projection of the first 1024 points
 *     of the \emph{original} Halton sequence onto the 32th and 33th dimensions.
 *     Note the strong correlation -- a scrambled sequence is usually preferred to avoid this problem.}{sampler_halton_nonscrambled_32}
 * }
 * \renderings{
 *     \unframedrendering{A projection of the first 1024 points
 *     of a randomly scrambled Halton sequence onto the first two dimensions (\code{scramble=1}).}{sampler_halton_rscrambled_0}
 *     \unframedrendering{A projection of the first 1024 points
 *     of a randomly scrambled Halton sequence onto the 32th and 33th dimensions.}{sampler_halton_rscrambled_32}
 * }
*
 * By default, the implementation provides a scrambled variant of the Halton sequence based
 * on permutations by Faure \cite{Faure1992Good} that has better equidistribution properties
 * in high dimensions, but this can be changed using the \code{scramble} parameter.
 * Internally, the plugin uses a table of prime numbers to provide elements
 * of the Halton sequence up to a dimension of 1024. Because of this upper bound,
 * the maximum path depth of the integrator must be limited (e.g. to 100), or
 * rendering might fail with the following error message: \emph{Lookup dimension
 * exceeds the prime number table size! You may have to reduce the 'maxDepth'
 * parameter of your integrator}.
 *
 * To support bucket-based renderings, the Halton sequence is internally enumerated
 * using a scheme proposed by Gr\"unschlo\ss\ et al. \cite{Grunschloss2010Enumerating};
 * the implementation in Mitsuba is based on a Python script by the authors of
 * this paper.
 *
 * \remarks{
 *   \item This sampler is incompatible with Metropolis Light Transport (all variants).
 *   It interoperates poorly with Bidirectional Path Tracing and Energy Redistribution
 *   Path Tracing, hence these should not be used together. The \pluginref{sobol} QMC
 *   sequence is an alternative for the latter two cases, and \pluginref{ldsampler}
 *   works as well.
 * }
 */
class HaltonSampler : public Sampler {
public:
	HaltonSampler() : Sampler(Properties()) { }

	HaltonSampler(const Properties &props) : Sampler(props) {
		/* Number of samples per pixel */
		m_sampleCount = props.getSize("sampleCount", 4);

		/* Scramble value, which can be used to break up temporally coherent
		   noise patterns when rendering the frames of an animation. */
		m_scramble = props.getInteger("scramble", -1);

		setFilmResolution(Vector2i(1), false);
		m_arrayStartDim = m_arrayEndDim = 5;
	}

	HaltonSampler(Stream *stream, InstanceManager *manager)
	 : Sampler(stream, manager) {
		m_arrayStartDim = stream->readUInt();
		m_arrayEndDim = stream->readUInt();
		m_offset = stream->readULong();
		m_stride = stream->readULong();
		m_multInverse[0] = stream->readULong();
		m_multInverse[1] = stream->readULong();
		m_scramble = stream->readInt();
		m_primePowers = Vector2i(stream);
		m_primeExponents = Vector2i(stream);
		m_pixelPosition = Point2i(0);
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Sampler::serialize(stream, manager);
		stream->writeUInt(m_arrayStartDim);
		stream->writeUInt(m_arrayEndDim);
		stream->writeULong(m_offset);
		stream->writeULong(m_stride);
		stream->writeULong(m_multInverse[0]);
		stream->writeULong(m_multInverse[1]);
		stream->writeInt(m_scramble);
		m_primePowers.serialize(stream);
		m_primeExponents.serialize(stream);
	}

	ref<Sampler> clone() {
		ref<HaltonSampler> sampler = new HaltonSampler();
		sampler->m_sampleCount = m_sampleCount;
		sampler->m_sampleIndex = m_sampleIndex;
		sampler->m_dimension = m_dimension;
		sampler->m_arrayStartDim = m_arrayStartDim;
		sampler->m_arrayEndDim = m_arrayEndDim;
		sampler->m_permutations = m_permutations;
		sampler->m_offset = m_offset;
		sampler->m_stride = m_stride;
		sampler->m_multInverse[0] = m_multInverse[0];
		sampler->m_multInverse[1] = m_multInverse[1];
		sampler->m_primePowers = m_primePowers;
		sampler->m_primeExponents = m_primeExponents;
		sampler->m_pixelPosition = m_pixelPosition;
		sampler->m_scramble = m_scramble;
		sampler->m_permutations = m_permutations;
		for (size_t i=0; i<m_req1D.size(); ++i)
			sampler->request1DArray(m_req1D[i]);
		for (size_t i=0; i<m_req2D.size(); ++i)
			sampler->request2DArray(m_req2D[i]);
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

	/**
	 * \brief Extended Euclidean algorithm
	 *
	 * 	Computes 'x' and 'y' that satisfy gcd(a, b) = ax + by,
	 * 	where x and y may be negative.
	 */
	void extendedGCD(uint64_t a, uint64_t b, int64_t &x, int64_t &y) {
		if (b == 0) {
			x = 1; y = 0;
			return;
		}

		int64_t d = a / b, x_, y_;
		extendedGCD(b, a % b, x_, y_);
		x = y_;
		y = x_ - (int64_t) (d * y_);
	}

	/**
	 * \brief Compute the multiplicative inverse of a math::modulo n,
	 * where a and n a re relative prime. The result is in the
	 * range 1, ..., n - 1.
	 */
	uint64_t multiplicativeInverse(int64_t a, int64_t n) {
		int64_t x, y;
		extendedGCD(a, n, x, y);
		return math::modulo(x, n);
	}

	void setFilmResolution(const Vector2i &res, bool blocked) {
		if (blocked) {
			/* Determine parameters of the space partition in the first two
			   dimensions. This is required to support bucketed rendering. */

			m_primeExponents = Vector2i(0, 0);
			m_stride = 1;

			for (int i=0; i<2; ++i) {
				int prime = primeTable[i],
					value = 1, exp = 0;

				while (value < std::min(res[i], MAX_RESOLUTION)) {
					value *= prime;
					++exp;
				}

				m_primePowers[i] = value;
				m_primeExponents[i] = exp;
				m_stride *= value;
			}

			m_multInverse[0] = multiplicativeInverse(m_primePowers.y, m_primePowers.x);
			m_multInverse[1] = multiplicativeInverse(m_primePowers.x, m_primePowers.y);
		} else {
			m_primeExponents[0] = m_primeExponents[1] = 0;
			m_primePowers[0] = m_primePowers[1] = 1;
			m_multInverse[0] = m_multInverse[1] = 0;
			m_stride = 1;
		}
		m_pixelPosition = Point2i(0);
		m_offset = 0;
	}

	void generate(const Point2i &pos) {
		/* Dimensions reserved to sample array requests */
		m_arrayStartDim = 5;
		m_arrayEndDim = m_arrayStartDim +
			static_cast<uint32_t>(m_req1D.size() + 2 * m_req2D.size());

		m_offset = 0;

		if (m_stride > 1) {
			for (int i=0; i<2; ++i) {
				m_pixelPosition[i] = pos[i] % MAX_RESOLUTION;

				/* Determine axis offset along each requested coordinate independently and
				   use the chinese remainder theorem to solve for a combined offset */
				uint64_t offset = inverseScrambledRadicalInverse(primeTable[i], m_pixelPosition[i],
						m_primeExponents[i], m_permutations.get()
						? m_permutations->getInversePermutation(i) : NULL);
				m_offset += offset * (m_stride / m_primePowers[i]) * m_multInverse[i];
			}
			m_offset %= m_stride;
		}

		uint32_t dim = m_arrayStartDim;
		for (size_t i=0; i<m_req1D.size(); i++) {
			if (!m_permutations.get()) {
				for (size_t j=0; j<m_sampleCount * m_req1D[i]; ++j)
					m_sampleArrays1D[i][j] = radicalInverseFast(dim, m_offset + j * m_stride);
			} else {
				uint16_t *perm = m_permutations->getPermutation(dim);
				for (size_t j=0; j<m_sampleCount * m_req1D[i]; ++j)
					m_sampleArrays1D[i][j] = scrambledRadicalInverseFast(dim, m_offset + j * m_stride, perm);
			}
			dim += 1;
		}

		for (size_t i=0; i<m_req2D.size(); i++) {
			if (!m_permutations.get()) {
				for (size_t j=0; j<m_sampleCount * m_req2D[i]; ++j) {
					m_sampleArrays2D[i][j] = Point2(
						radicalInverseFast(dim,   m_offset + j * m_stride),
						radicalInverseFast(dim+1, m_offset + j * m_stride));
				}
			} else {
				uint16_t *perm1 = m_permutations->getPermutation(dim);
				uint16_t *perm2 = m_permutations->getPermutation(dim+1);
				for (size_t j=0; j<m_sampleCount * m_req2D[i]; ++j) {
					m_sampleArrays2D[i][j] = Point2(
						scrambledRadicalInverseFast(dim,   m_offset + j * m_stride, perm1),
						scrambledRadicalInverseFast(dim+1, m_offset + j * m_stride, perm2));
				}
			}
			dim += 2;
		}

		setSampleIndex(0);
	}

	void advance() {
		m_sampleIndex++;
		m_dimension = 0;
		m_dimension1DArray = m_dimension2DArray = 0;
	}

	void setSampleIndex(size_t sampleIndex) {
		m_dimension = 0;
		m_sampleIndex = sampleIndex;
		m_dimension1DArray = m_dimension2DArray = 0;
	}

	inline Float nextFloat(uint64_t idx) {
		uint32_t dim = m_dimension++;
		if (m_permutations != NULL)
			return scrambledRadicalInverseFast(dim, idx,
				m_permutations->getPermutation(dim));
		else
			return radicalInverseFast(dim, idx);
	}

	Float next1D() {
		/* Skip over dimensions that were reserved to arrays */
		if (m_dimension >= m_arrayStartDim && m_dimension < m_arrayEndDim)
			m_dimension = m_arrayEndDim;
		if (m_dimension >= primeTableSize)
			Log(EError, "Lookup dimension exceeds the prime number table size! "
				"You may have to reduce the 'maxDepth' parameter of your integrator.");

		uint64_t index = m_offset + m_stride * m_sampleIndex;
		return nextFloat(index);
	}

	Point2 next2D() {
		/* Skip over dimensions that were reserved to arrays */
		if (m_dimension + 1 >= m_arrayStartDim && m_dimension < m_arrayEndDim)
			m_dimension = m_arrayEndDim;
		if (m_dimension + 1 >= primeTableSize)
			Log(EError, "Lookup dimension exceeds the prime number table size! "
				"You may have to reduce the 'maxDepth' parameter of your integrator.");

		uint64_t index = m_offset + m_stride * m_sampleIndex;

		Float value1, value2;
		if (m_dimension == 0) {
			value1 = nextFloat(index) * m_primePowers.x - m_pixelPosition.x;
			value2 = nextFloat(index) * m_primePowers.y - m_pixelPosition.y;
		} else {
			value1 = nextFloat(index);
			value2 = nextFloat(index);
		}

		return Point2(value1, value2);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HaltonSampler[" << endl
			<< "  sampleCount = " << m_sampleCount << "," << endl
			<< "  sampleIndex = " << m_sampleIndex << "," << endl
			<< "  stride = " << m_stride << "," << endl
			<< "  primePowers = " << m_primePowers.toString() << "," << endl
			<< "  primeExponents = " << m_primeExponents.toString() << "," << endl
			<< "  multInverse = [" << m_multInverse[0] << ", " << m_multInverse[1] << "]," << endl
			<< "  scramble = " << m_scramble << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	uint32_t m_dimension;
	uint32_t m_arrayStartDim;
	uint32_t m_arrayEndDim;
	int m_scramble;

	/* Faure permutation */
	ref<const PermutationStorage> m_permutations;
	static ref<const PermutationStorage> m_globalPermutations;
	static ref<Mutex> m_globalPermutationsMutex;

	/* Space partition */
	uint64_t m_offset, m_stride;
	uint64_t m_multInverse[2];
	Vector2i m_primePowers;
	Vector2i m_primeExponents;
	Point2i m_pixelPosition;
};

ref<Mutex> HaltonSampler::m_globalPermutationsMutex = new Mutex();
ref<const PermutationStorage> HaltonSampler::m_globalPermutations = NULL;

MTS_IMPLEMENT_CLASS_S(HaltonSampler, false, Sampler)
MTS_EXPORT_PLUGIN(HaltonSampler, "Halton QMC sampler");
MTS_NAMESPACE_END

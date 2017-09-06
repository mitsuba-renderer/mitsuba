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

#pragma once
#if !defined(__MITSUBA_CORE_QMC_H_)
#define __MITSUBA_CORE_QMC_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/** \addtogroup libcore
 *  \addtogroup libpython
 *  @{
 */

// -----------------------------------------------------------------------
//! @{ \name Elementary Quasi-Monte Carlo number sequences
//
// Based on implementations by Leonhard Gruenschloss
// -----------------------------------------------------------------------

static const size_t primeTableSize = 1024;
/// Table of the first 1024 prime numbers
extern const int MTS_EXPORT_CORE primeTable[primeTableSize];

/// Van der Corput radical inverse in base 2 with single precision
inline float radicalInverse2Single(uint32_t n, uint32_t scramble = 0U) {
    /* Efficiently reverse the bits in 'n' using binary operations */
#if (defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 2))) || defined(__clang__)
    n = __builtin_bswap32(n);
#else
    n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
#endif
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);

    // Account for the available precision and scramble
    n = (n >> (32 - 24)) ^ (scramble & ~-(1 << 24));

    return (float) n / (float) (1U << 24);
}

/// Van der Corput radical inverse in base 2 with double precision
inline double radicalInverse2Double(uint64_t n, uint64_t scramble = 0ULL) {
    /* Efficiently reverse the bits in 'n' using binary operations */
#if (defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 2))) || defined(__clang__)
    n = __builtin_bswap64(n);
#else
    n = (n << 32) | (n >> 32);
    n = ((n & 0x0000ffff0000ffffULL) << 16) | ((n & 0xffff0000ffff0000ULL) >> 16);
    n = ((n & 0x00ff00ff00ff00ffULL) << 8)  | ((n & 0xff00ff00ff00ff00ULL) >> 8);
#endif
    n = ((n & 0x0f0f0f0f0f0f0f0fULL) << 4)  | ((n & 0xf0f0f0f0f0f0f0f0ULL) >> 4);
    n = ((n & 0x3333333333333333ULL) << 2)  | ((n & 0xccccccccccccccccULL) >> 2);
    n = ((n & 0x5555555555555555ULL) << 1)  | ((n & 0xaaaaaaaaaaaaaaaaULL) >> 1);

    // Account for the available precision and scramble
    n = (n >> (64 - 53)) ^ (scramble & ~-(1LL << 53));

    return (double) n / (double) (1ULL << 53);
}

/// Sobol' radical inverse in base 2 with single precision.
inline float sobol2Single(uint32_t n, uint32_t scramble = 0U) {
    for (uint32_t v = 1U << 31; n != 0; n >>= 1, v ^= v >> 1)
        if (n & 1)
            scramble ^= v;
    return (float) scramble / (float) (1ULL << 32);
}

/// Sobol' radical inverse in base 2 with double precision.
inline double sobol2Double(uint64_t n, uint64_t scramble = 0ULL) {
    scramble &= ~-(1LL << 53);
    for (uint64_t v = 1ULL << 52; n != 0; n >>= 1, v ^= v >> 1)
        if (n & 1)
            scramble ^= v;
    return (double) scramble / (double) (1ULL << 53);
}

/// Generate an element from a (0, 2) sequence, single precision
inline Point2f sample02Single(uint32_t n, uint32_t scramble[2]) {
    return Point2f(
        radicalInverse2Single(n, scramble[0]),
        sobol2Single(n, scramble[1])
    );
}

/// Generate an element from a (0, 2) sequence, double precision version
inline Point2d sample02Double(uint64_t n, uint64_t scramble[2]) {
    return Point2d(
        radicalInverse2Double(n, scramble[0]),
        sobol2Double(n, scramble[1])
    );
}

/// Generate an element from a (0, 2) sequence (without scrambling)
inline Point2 sample02(size_t n) {
    #if defined(SINGLE_PRECISION)
        return Point2(
            radicalInverse2Single((uint32_t) n),
            sobol2Single((uint32_t) n)
        );
    #else
        return Point2(
            radicalInverse2Double((uint64_t) n),
            sobol2Double((uint64_t) n)
        );
    #endif
}

/**
 * \brief Generate fast and reasonably good pseudorandom numbers using the
 * Tiny Encryption Algorithm (TEA) by David Wheeler and Roger Needham.
 *
 * For details, refer to "GPU Random Numbers via the Tiny Encryption Algorithm"
 * by Fahad Zafar, Marc Olano, and Aaron Curtis.
 *
 * \param v0
 *     First input value to be encrypted (could be the sample index)
 * \param v1
 *     Second input value to be encrypted (e.g. the requested random number dimension)
 * \param rounds
 *     How many rounds should be executed? The default for random number
 *     generation is 4.
 * \return
 *     A uniformly distributed 64-bit integer
 */
inline uint64_t sampleTEA(uint32_t v0, uint32_t v1, int rounds = 4) {
    uint32_t sum = 0;

    for (int i=0; i<rounds; ++i) {
        sum += 0x9e3779b9;
        v0 += ((v1 << 4) + 0xA341316C) ^ (v1 + sum) ^ ((v1 >> 5) + 0xC8013EA4);
        v1 += ((v0 << 4) + 0xAD90777D) ^ (v0 + sum) ^ ((v0 >> 5) + 0x7E95761E);
    }

    return ((uint64_t) v1 << 32) + v0;
}

/**
 * \brief Generate fast and reasonably good pseudorandom numbers using the
 * Tiny Encryption Algorithm (TEA) by David Wheeler and Roger Needham.
 *
 * This function uses \ref sampleTEA to return single precision floating point
 * numbers on the interval <tt>[0, 1)</tt>
 *
 * \param v0
 *     First input value to be encrypted (could be the sample index)
 * \param v1
 *     Second input value to be encrypted (e.g. the requested random number dimension)
 * \param rounds
 *     How many rounds should be executed? The default for random number
 *     generation is 4.
 * \return
 *     A uniformly distributed floating point number on the interval <tt>[0, 1)</tt>
 */
inline float sampleTEASingle(uint32_t v0, uint32_t v1, int rounds = 4) {
    /* Trick from MTGP: generate an uniformly distributed
       single precision number in [1,2) and subtract 1. */
    union {
        uint32_t u;
        float f;
    } x;
    x.u = ((sampleTEA(v0, v1, rounds) & 0xFFFFFFFF) >> 9) | 0x3f800000UL;
    return x.f - 1.0f;
}

/**
 * \brief Generate fast and reasonably good pseudorandom numbers using the
 * Tiny Encryption Algorithm (TEA) by David Wheeler and Roger Needham.
 *
 * This function uses \ref sampleTEA to return single precision floating point
 * numbers on the interval <tt>[0, 1)</tt>
 *
 * \param v0
 *     First input value to be encrypted (could be the sample index)
 * \param v1
 *     Second input value to be encrypted (e.g. the requested random number dimension)
 * \param rounds
 *     How many rounds should be executed? The default for random number
 *     generation is 4.
 * \return
 *     A uniformly distributed floating point number on the interval <tt>[0, 1)</tt>
 */
inline double sampleTEADouble(uint32_t v0, uint32_t v1, int rounds = 4) {
    /* Trick from MTGP: generate an uniformly distributed
       single precision number in [1,2) and subtract 1. */
    union {
        uint64_t u;
        double f;
    } x;
    x.u = (sampleTEA(v0, v1, rounds) >> 12) | 0x3ff0000000000000ULL;
    return x.f - 1.0;
}

#if defined(SINGLE_PRECISION)
/// Alias to \ref sampleTEASingle or \ref sampleTEADouble based on compilation flags
inline Float sampleTEAFloat(uint32_t v0, uint32_t v1, int rounds = 4) {
    return sampleTEASingle(v0, v1, rounds);
}
#else
/// Alias to \ref sampleTEASingle or \ref sampleTEADouble based on compilation flags
inline Float sampleTEAFloat(uint32_t v0, uint32_t v1, int rounds = 4) {
    return sampleTEADouble(v0, v1, rounds);
}
#endif

/**
 * \brief Calculate the radical inverse function
 *
 * This function is used as a building block to construct Halton and
 * Hammersley sequences. Roughly, it computes a b-ary representation
 * of the input value \c index, mirrors it along the decimal
 * point, and returns the resulting fractional value.
 */
extern MTS_EXPORT_CORE Float radicalInverse(int base, uint64_t index);

/**
 * \brief Calculate a scrambled radical inverse function
 *
 * This function is used as a building block to construct permuted
 * Halton and Hammersley sequence variants. It works like the normal
 * radical inverse function \ref radicalInverse(), except that every digit
 * is run through an extra scrambling permutation specified as array
 * of size \c base.
 *
 * \remark This function is not available in the Python API
 */
extern MTS_EXPORT_CORE Float scrambledRadicalInverse(int base,
    uint64_t index, uint16_t *perm);

/**
 * \brief Incrementally calculate the next Van Der Corput sequence
 * value starting from a current entry \c x (wrt. a fixed base)
 *
 * Repeated evaluation eventually causes a loss of accuracy
 */
extern MTS_EXPORT_CORE Float radicalInverseIncremental(int base, Float x);

/**
 * \brief Calculate a radical inverse function (fast version)
 *
 * This function works similarly to \ref radicalInverse, but is potentially much
 * faster. Internally, it relies on optimized implementations of the radical inverse
 * functions for the first 1024 prime number bases. For that reason, only works for
 * such bases.
 *
 * \param baseIndex
 *    Prime number index starting at 0 (i.e. 3 would cause 7 to be
 *    used as the basis)
 * \param index
 *    Sequence index
 */
extern MTS_EXPORT_CORE Float radicalInverseFast(uint16_t baseIndex, uint64_t index);

/**
 * \brief Calculate a scrambled radical inverse function (fast version)
 *
 * This function is used as a building block to construct permuted
 * Halton and Hammersley sequence variants. It works like the fast
 * radical inverse function \ref radicalInverseFast(), except that every
 * digit is run through an extra scrambling permutation.
 *
 * \remark This function is not available in the Python API
 */
extern MTS_EXPORT_CORE Float scrambledRadicalInverseFast(uint16_t baseIndex,
        uint64_t index, uint16_t *perm);

//! @}
// -----------------------------------------------------------------------

//! @}

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_QMC_H_ */

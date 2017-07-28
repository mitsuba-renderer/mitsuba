// Copyright (c) 2012 Leonhard Gruenschloss (leonhard@gruenschloss.org)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#if !defined(__SOBOL_H)
#define __SOBOL_H

#include <mitsuba/mitsuba.h>
#include <cassert>

namespace sobol {

struct Matrices {
    static const uint32_t num_dimensions = 1024;
    static const uint32_t size = 52;
    static const uint32_t matrices32[];
    static const uint64_t matrices64[];
    static const uint64_t vdc_sobol_matrices_inv[][52];
    static const uint64_t vdc_sobol_matrices[][52];
};

// Compute one component of the Sobol'-sequence, where the component
// corresponds to the dimension parameter, and the index specifies
// the point inside the sequence. The scramble parameter can be used
// to permute elementary intervals, and might be chosen randomly to
// generate a randomized QMC sequence.
inline float sampleSingle(
    uint64_t index,
    const uint32_t dimension,
    const uint32_t scramble = 0U)
{
    assert(dimension < Matrices::num_dimensions);

    uint32_t result = scramble;
    for (uint32_t i = dimension * Matrices::size; index; index >>= 1, ++i)
    {
        if (index & 1)
            result ^= Matrices::matrices32[i];
    }

    return std::min(result * (1.0f / (1ULL << 32)), ONE_MINUS_EPS_FLT);
}

// Compute one component of the Sobol'-sequence, where the component
// corresponds to the dimension parameter, and the index specifies
// the point inside the sequence. The scramble parameter can be used
// to permute elementary intervals, and might be chosen randomly to
// generate a randomized QMC sequence. Only the Matrices::size least
// significant bits of the scramble value are used.
inline double sampleDouble(
    uint64_t index,
    const uint32_t dimension,
    const uint64_t scramble = 0ULL)
{
    assert(dimension < Matrices::num_dimensions);

    uint64_t result = scramble & ~-(1LL << Matrices::size);
    for (uint32_t i = dimension * Matrices::size; index; index >>= 1, ++i)
    {
        if (index & 1)
            result ^= Matrices::matrices64[i];
    }

    return std::min(result * (1.0 / (1ULL << Matrices::size)), ONE_MINUS_EPS_DBL);
}

// Call sampleSingle or sampleDouble depending on the compilation options
inline mitsuba::Float sample(
    const uint64_t index,
    const uint32_t dimension,
    const uint64_t scramble = 0ULL)
{
#if defined(SINGLE_PRECISION)
    return sampleSingle(index, dimension, (uint32_t) scramble);
#else
    return sampleDouble(index, dimension, (uint64_t) scramble);
#endif
}

// Return the index of the frame-th sample falling
// into the square elementary interval (px, py),
// without using look-up tables.
inline uint64_t look_up(
    const uint32_t m,
    uint32_t frame,
    const uint32_t px,
    const uint32_t py,
    uint64_t scramble)
{
    const uint32_t m2 = m << 1;
    uint64_t index = uint64_t(frame) << m2;

    // Note: the delta value only depends on frame
    // and m, thus it can be cached across multiple
    // function calls, if desired.
    uint64_t delta = 0;
    for (uint32_t c = 0; frame; frame >>= 1, ++c)
        if (frame & 1) // Add flipped column m + c + 1.
            delta ^= Matrices::vdc_sobol_matrices[m - 1][c];

#if defined(SINGLE_PRECISION)
    scramble = (scramble & 0xFFFFFFFF) >> (32 - m);
#else
    scramble = (scramble & ~-(1ULL << Matrices::size)) >> (Matrices::size - m);
#endif

    // flipped b
    uint64_t b = (((uint64_t) (px ^ scramble) << m) | (py ^ scramble)) ^ delta;

    for (uint32_t c = 0; b; b >>= 1, ++c)
        if (b & 1) // Add column 2 * m - c.
            index ^= Matrices::vdc_sobol_matrices_inv[m - 1][c];

    return index;
}


} // namespace sobol

#endif /* __SOBOL_H */

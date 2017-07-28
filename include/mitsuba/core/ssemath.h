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
#if !defined(__MITSUBA_CORE_SSEMATH_H_)
#define __MITSUBA_CORE_SSEMATH_H_

#ifdef MTS_SSE

#include <mitsuba/core/sse.h>

MTS_NAMESPACE_BEGIN

namespace math {
    /**
     * \brief SIMD (SSE2) implementation of \c log
     * \author Julien Pommier
     */
    extern MTS_EXPORT_CORE __m128 log_ps(__m128 x);

    /**
     * \brief SIMD (SSE2) implementation of \c exp
     * \author Julien Pommier
     */
    extern MTS_EXPORT_CORE __m128 exp_ps(__m128 x);

    /**
     * \brief SIMD (SSE2) implementation of \c sin
     * \author Julien Pommier
     */
    extern MTS_EXPORT_CORE __m128 sin_ps(__m128 x);

    /**
     * \brief SIMD (SSE2) implementation of \c cos
     * \author Julien Pommier
     */
    extern MTS_EXPORT_CORE __m128 cos_ps(__m128 x);

    /**
     * \brief SIMD (SSE2) implementation which simultaneously
     * computes the sine and cosine of a given value
     * \author Julien Pommier
     */
    extern MTS_EXPORT_CORE void sincos_ps(__m128 x, __m128* s, __m128* c);

    /**
     * \brief Fast SIMD (SSE2) approximation of \c log
     * which provides about 10-11 mantissa bits.
     * Inspired by the Intel Approximate Math Library.
     */
    extern MTS_EXPORT_CORE __m128 fastlog_ps(__m128 x);

    /**
     * \brief Fast SIMD (SSE2) approximation of \c pow
     * which provides about 10-11 mantissa bits.
     * Inspired by the Intel Approximate Math Library.
     */
    extern MTS_EXPORT_CORE __m128 fastpow_ps(__m128 x, __m128 y);

    /**
     * \brief The arguments <tt>row0</tt>, <tt>row1</tt>, <tt>row2</tt> and
     * <tt>row3</tt> are \c __m128 values whose elements form the corresponding
     * rows of a 4-by-4 matrix. The matrix transposition is returned in
     * arguments <tt>row0</tt>, <tt>row1</tt>, <tt>row2</tt> and <tt>row3</tt>
     * where \c row0 now holds column 0 of the original matrix, \c row1 now
     * holds column 1 of the original matrix, and so on.
     * \author Intel Intrinsics Guide for AVX2
     */
    FINLINE void transpose_ps(__m128& row0, __m128& row1,
        __m128& row2, __m128& row3) {
        __m128 tmp3, tmp2, tmp1, tmp0;
        tmp0 = _mm_unpacklo_ps(row0, row1);
        tmp2 = _mm_unpacklo_ps(row2, row3);
        tmp1 = _mm_unpackhi_ps(row0, row1);
        tmp3 = _mm_unpackhi_ps(row2, row3);

        row0 = _mm_movelh_ps(tmp0, tmp2);
        row1 = _mm_movehl_ps(tmp2, tmp0);
        row2 = _mm_movelh_ps(tmp1, tmp3);
        row3 = _mm_movehl_ps(tmp3, tmp1);
    }

    /// Component-wise clamp: <tt>max(min(x, maxVal), minVal)</tt>
    inline __m128 clamp_ps(__m128 x, __m128 minVal, __m128 maxVal) {
        return _mm_max_ps(_mm_min_ps(x, maxVal), minVal);
    }

    /// Sum of all elements in the vector
    inline float hsum_ps(__m128 vec) {
        __m128 tmp = _mm_shuffle_ps(vec, vec,  _MM_SHUFFLE(1,0,3,2));
        __m128 sum_tmp = _mm_add_ps(vec, tmp);
        tmp = _mm_shuffle_ps(sum_tmp, sum_tmp, _MM_SHUFFLE(2,3,0,1));
        sum_tmp = _mm_add_ps(sum_tmp, tmp);
        return _mm_cvtss_f32(sum_tmp);
    }

    /// Maximum across all the elements of a vector
    inline float hmax_ps(__m128 vec) {
        __m128 tmp = _mm_shuffle_ps(vec, vec,  _MM_SHUFFLE(1,0,3,2));
        __m128 tmp_max = _mm_max_ps(vec, tmp);
        tmp = _mm_shuffle_ps(tmp_max, tmp_max, _MM_SHUFFLE(2,3,0,1));
        tmp_max = _mm_max_ps(tmp_max, tmp);
        return _mm_cvtss_f32(tmp_max);
    }

    /// Minimum across all the elements of a vector
    inline float hmin_ps(__m128 vec) {
        __m128 tmp = _mm_shuffle_ps(vec, vec,  _MM_SHUFFLE(1,0,3,2));
        __m128 tmp_min = _mm_min_ps(vec, tmp);
        tmp = _mm_shuffle_ps(tmp_min, tmp_min, _MM_SHUFFLE(2,3,0,1));
        tmp_min = _mm_min_ps(tmp_min, tmp);
        return _mm_cvtss_f32(tmp_min);
    }
};

MTS_NAMESPACE_END

#endif /* MTS_SSE */

#endif /* __MITSUBA_CORE_SSEMATH_H_ */

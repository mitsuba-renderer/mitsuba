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

/*============================================================================
  HDRITools - High Dynamic Range Image Tools
  Copyright 2008-2012 Program of Computer Graphics, Cornell University

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
 -----------------------------------------------------------------------------
 Primary author:
     Edgar Velazquez-Armendariz <cs#cornell#edu - eva5>
============================================================================*/

#pragma once
#if !defined(__MITSUBA_CORE_SSEVECTOR_H_)
#define __MITSUBA_CORE_SSEVECTOR_H_

#include <mitsuba/core/platform.h>
#include <mitsuba/core/sse.h>

#if !MTS_SSE
# error "This header requires SSE support"
#endif

MTS_NAMESPACE_BEGIN

namespace math
{

// Forward declarations, required by Clang and ICL 12.1
struct SSEVector4f;
struct SSEvector4i;

template <int idx3, int idx2, int idx1, int idx0>
SSEVector4f shuffle(const SSEVector4f& low, const SSEVector4f& hi);

template <int idx3, int idx2, int idx1, int idx0>
SSEVector4f shuffle(const SSEVector4f& a);


struct SSEVector4f
{
private:
    __m128 xmm;

public:
    SSEVector4f() {}
    SSEVector4f(const SSEVector4f& other) : xmm(other.xmm) {}
    SSEVector4f(__m128 val) : xmm(val) {}
    explicit SSEVector4f(float val) : xmm(_mm_set1_ps(val)) {}
    SSEVector4f(float f3, float f2, float f1, float f0) :
    xmm(_mm_set_ps(f3, f2, f1, f0))
    {}

    inline SSEVector4f& operator= (float val) {
        xmm = _mm_set1_ps(val);
        return *this;
    }

    inline static SSEVector4f zero() {
        return _mm_setzero_ps();
    }

    operator __m128() const {
        return xmm;
    }

    friend SSEVector4f operator& (const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_and_ps(a.xmm, b.xmm);
    }
    friend SSEVector4f operator| (const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_or_ps(a.xmm, b.xmm);
    }
    friend SSEVector4f operator^ (const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_xor_ps(a.xmm, b.xmm);
    }
    /// ~a & b
    friend SSEVector4f andnot(const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_andnot_ps(a.xmm, b.xmm);
    }

    SSEVector4f& operator&= (const SSEVector4f& a) {
        xmm = _mm_and_ps(xmm, a.xmm);
        return *this;
    }
    SSEVector4f& operator|= (const SSEVector4f& a) {
        xmm = _mm_or_ps(xmm, a.xmm);
        return *this;
    }
    SSEVector4f& operator^= (const SSEVector4f& a) {
        xmm = _mm_xor_ps(xmm, a.xmm);
        return *this;
    }

    friend SSEVector4f operator+ (const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_add_ps(a.xmm, b.xmm);
    }
    friend SSEVector4f operator- (const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_sub_ps(a.xmm, b.xmm);
    }
    friend SSEVector4f operator* (const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_mul_ps(a.xmm, b.xmm);
    }
    friend SSEVector4f operator/ (const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_div_ps(a.xmm, b.xmm);
    }

    SSEVector4f& operator+= (const SSEVector4f& a) {
        xmm = _mm_add_ps(xmm, a.xmm);
        return *this;
    }
    SSEVector4f& operator-= (const SSEVector4f& a) {
        xmm = _mm_sub_ps(xmm, a.xmm);
        return *this;
    }
    SSEVector4f& operator*= (const SSEVector4f& a) {
        xmm = _mm_mul_ps(xmm, a.xmm);
        return *this;
    }
    SSEVector4f& operator/= (const SSEVector4f& a) {
        xmm = _mm_div_ps(xmm, a.xmm);
        return *this;
    }

    /**
     * \brief Newton-Rhapson Reciprocal:
     * \f[ 2 * rcp(x) - (x * rcp(x) * rcp(x)) \f]
     */
    friend inline SSEVector4f rcp_nr(const SSEVector4f& v) {
        __m128 x0 = _mm_rcp_ps(v.xmm);
        return _mm_sub_ps(_mm_add_ps(x0,x0),
            _mm_mul_ps(_mm_mul_ps(x0,v.xmm), x0));
    }

    friend inline SSEVector4f rcp(const SSEVector4f& v) {
        return _mm_rcp_ps(v.xmm);
    }

    friend SSEVector4f min(const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_min_ps(a.xmm, b.xmm);
    }
    friend SSEVector4f max(const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_max_ps(a.xmm, b.xmm);
    }

    friend SSEVector4f isnan(const SSEVector4f& a) {
        return _mm_cmpunord_ps(a.xmm, a.xmm);
    }
    friend SSEVector4f isnan(const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_cmpunord_ps(a.xmm, b.xmm);
    }

    /**
     * \brief Moves either of the values of \c low into the low 64-bits
     * of the result, and either of the values of \c high into
     * the high 64-bits of the result. Each index in the
     * template is a index in the range [0,3] to choose a value from the
     * source, 0 being the lowest and 3 the highest.
     */
    template <int idx3, int idx2, int idx1, int idx0>
    friend SSEVector4f shuffle(const SSEVector4f& low, const SSEVector4f& hi) {
        return _mm_shuffle_ps(low.xmm,hi.xmm,_MM_SHUFFLE(idx3,idx2,idx1,idx0));
    }

    /// Shuffles the elements of the given vector using the indices [0,3]
    template <int idx3, int idx2, int idx1, int idx0>
    friend SSEVector4f shuffle(const SSEVector4f& a) {
        return _mm_shuffle_ps(a.xmm, a.xmm, _MM_SHUFFLE(idx3,idx2,idx1,idx0));
    }

    /// a == b
    friend SSEVector4f cmpeq(const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_cmpeq_ps(a.xmm, b.xmm);
    }
    /// a < b
    friend SSEVector4f cmplt(const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_cmplt_ps(a.xmm, b.xmm);
    }
    /// a <= b
    friend SSEVector4f cmple(const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_cmple_ps(a.xmm, b.xmm);
    }
    /// a > b
    friend SSEVector4f cmpgt(const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_cmpgt_ps(a.xmm, b.xmm);
    }
    /// a >= b
    friend SSEVector4f cmpge(const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_cmpge_ps(a.xmm, b.xmm);
    }
    /// a != b
    friend SSEVector4f cmpneq(const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_cmpneq_ps(a.xmm, b.xmm);
    }
    /// !(a < b)
    friend SSEVector4f cmpnlt(const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_cmpnlt_ps(a.xmm, b.xmm);
    }
    /// !(a <= b)
    friend SSEVector4f cmpnle(const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_cmpnle_ps(a.xmm, b.xmm);
    }
    /// !(a > b)
    friend SSEVector4f cmpngt(const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_cmpngt_ps(a.xmm, b.xmm);
    }
    /// !(a >= b)
    friend SSEVector4f cmpnge(const SSEVector4f& a, const SSEVector4f& b) {
        return _mm_cmpnge_ps(a.xmm, b.xmm);
    }

    friend SSEVector4f operator==(const SSEVector4f& a, const SSEVector4f& b) {
        return cmpeq(a, b);
    }
    friend SSEVector4f operator!=(const SSEVector4f& a, const SSEVector4f& b) {
        return cmpneq(a, b);
    }
    friend SSEVector4f operator<(const SSEVector4f& a, const SSEVector4f& b) {
        return cmplt(a, b);
    }
    friend SSEVector4f operator<=(const SSEVector4f& a, const SSEVector4f& b) {
        return cmple(a, b);
    }
    friend SSEVector4f operator>(const SSEVector4f& a, const SSEVector4f& b) {
        return cmpgt(a, b);
    }
    friend SSEVector4f operator>=(const SSEVector4f& a, const SSEVector4f& b) {
        return cmpge(a, b);
    }

    /// Select/blend operation <tt>(mask) ? a : b</tt>
    friend inline SSEVector4f select(const SSEVector4f& mask,
        const SSEVector4f& a, const SSEVector4f& b) {
        // Alternative method by Jim Conyngham/Wikipedia MD5 page, via
        // http://markplusplus.wordpress.com/2007/03/14/fast-sse-select-operation/ [July 2012]
        return _mm_xor_ps(b.xmm, _mm_and_ps(mask.xmm, _mm_xor_ps(a.xmm, b.xmm)));
    }

    /// Round \c a towards zero
    friend inline SSEVector4f roundTruncate(const SSEVector4f& a) {
        __m128i truncated = _mm_cvttps_epi32(a.xmm);
        return _mm_cvtepi32_ps(truncated);
    }

    /// Save to \c dest without polluting the cache
    friend inline void stream(SSEVector4f* dest, const SSEVector4f& value) {
        _mm_stream_ps(reinterpret_cast<float*>(dest), value.xmm);
    }
    /// Save to \c dest without polluting the cache
    friend inline void stream(__m128* dest, const SSEVector4f& value) {
        _mm_stream_ps(reinterpret_cast<float*>(dest), value.xmm);
    }
    /// Save to \c dest without polluting the cache
    friend inline void stream(float* dest, const SSEVector4f& value) {
        _mm_stream_ps(dest, value.xmm);
    }
};



struct SSEVector4i
{
private:
    __m128i xmm;

public:
    SSEVector4i() {}
    SSEVector4i(const SSEVector4i& val) : xmm(val.xmm) {}
    SSEVector4i(__m128i val) : xmm(val) {}
    explicit SSEVector4i(int32_t val) : xmm(_mm_set1_epi32(val)) {}
    SSEVector4i(int32_t i3, int32_t i2, int32_t i1, int32_t i0) :
    xmm(_mm_set_epi32(i3, i2, i1, i0))
    {}

    SSEVector4i& operator= (int32_t val) {
        xmm = _mm_set1_epi32(val);
        return *this;
    }

    inline static SSEVector4i zero() {
        return _mm_setzero_si128();
    }

    operator __m128i() const {
        return xmm;
    }

    friend SSEVector4i operator& (const SSEVector4i& a, const SSEVector4i& b) {
        return _mm_and_si128(a.xmm, b.xmm);
    }
    friend SSEVector4i operator| (const SSEVector4i& a, const SSEVector4i& b) {
        return _mm_or_si128(a.xmm, b.xmm);
    }
    friend SSEVector4i operator^ (const SSEVector4i& a, const SSEVector4i& b) {
        return _mm_xor_si128(a.xmm, b.xmm);
    }
    /// ~a & b
    friend SSEVector4i andnot(const SSEVector4i& a, const SSEVector4i& b) {
        return _mm_andnot_si128(a.xmm, b.xmm);
    }
    SSEVector4i& operator&= (const SSEVector4i& a) {
        xmm = _mm_and_si128(xmm, a.xmm);
        return *this;
    }
    SSEVector4i& operator|= (const SSEVector4i& a) {
        xmm = _mm_or_si128(xmm, a.xmm);
        return *this;
    }
    SSEVector4i& operator^= (const SSEVector4i& a) {
        xmm = _mm_xor_si128(xmm, a.xmm);
        return *this;
    }

    friend SSEVector4i operator+ (const SSEVector4i& a, const SSEVector4i& b) {
        return _mm_add_epi32(a.xmm, b.xmm);
    }
    friend SSEVector4i operator- (const SSEVector4i& a, const SSEVector4i& b) {
        return _mm_sub_epi32(a.xmm, b.xmm);
    }
    SSEVector4i& operator+= (const SSEVector4i& a) {
        xmm = _mm_add_epi32(xmm, a.xmm);
        return *this;
    }
    SSEVector4i& operator-= (const SSEVector4i& a) {
        xmm = _mm_sub_epi32(xmm, a.xmm);
        return *this;
    }

    /// Test if all elements are zero
    inline bool isZero() const {
        const __m128i mask = _mm_cmpeq_epi32(xmm, _mm_setzero_si128());
        return _mm_movemask_epi8(mask) == 0xFFFF;
    }

    /// a == b
    friend SSEVector4i cmpeq(const SSEVector4i& a, const SSEVector4i& b) {
        return _mm_cmpeq_epi32(a.xmm, b.xmm);
    }
    /// a < b
    friend SSEVector4i cmplt(const SSEVector4i& a, const SSEVector4i& b) {
        return _mm_cmplt_epi32(a.xmm, b.xmm);
    }
    /// a > b
    friend SSEVector4i cmpgt(const SSEVector4i& a, const SSEVector4i& b) {
        return _mm_cmpgt_epi32(a.xmm, b.xmm);
    }
    friend SSEVector4i operator==(const SSEVector4i& a, const SSEVector4i& b) {
        return cmpeq(a, b);
    }
    friend SSEVector4i operator<(const SSEVector4i& a, const SSEVector4i& b) {
        return cmplt(a, b);
    }
    friend SSEVector4i operator>(const SSEVector4i& a, const SSEVector4i& b) {
        return cmpgt(a, b);
    }

    /// Select/blend: <tt>(mask) ? a : b</tt>
    friend inline SSEVector4i select(const SSEVector4i& mask,
        const SSEVector4i& a, const SSEVector4i& b) {
        // Alternative method by Jim Conyngham/Wikipedia MD5 page, via
        // http://markplusplus.wordpress.com/2007/03/14/fast-sse-select-operation/ [July 2012]
        return _mm_xor_si128(b.xmm,
            _mm_and_si128(mask.xmm, _mm_xor_si128(a.xmm, b.xmm)));
    }

    template <int32_t i3, int32_t i2, int32_t i1, int32_t i0>
    static const __m128i& constant() {
        static const union {
            int32_t i32[4];
            __m128i xmm;
        } u = {{i0, i1, i2, i3}};
        return u.xmm;
    }

    template <int32_t value>
    static const __m128i& constant() {
        static const union {
            int32_t i32[4];
            __m128i xmm;
        } u = {{value, value, value, value}};
        return u.xmm;
    }

    /// Shift right by \c count bits while shifting in zeros
    friend inline SSEVector4i srl(const SSEVector4i& a, int count) {
        return _mm_srli_epi32(a.xmm, count);
    }

    /// Shift left by \c count bits while shifting in zeros
    friend inline SSEVector4i sll(const SSEVector4i& a, int count) {
        return _mm_slli_epi32(a.xmm, count);
    }

    /// Save to \c dest without polluting the cache
    friend inline void stream(SSEVector4i* dest, const SSEVector4i& value) {
        _mm_stream_si128(&(dest->xmm), value);
    }
    /// Save to \c dest without polluting the cache
    friend inline void stream(__m128i* dest, const SSEVector4i& value) {
        _mm_stream_si128(dest, value);
    }
};

/// Reinterprets \c as a \c SSEVector4i
inline SSEVector4i castAsInt(const SSEVector4f& a) {
    return _mm_castps_si128(a);
}
/// Convert \c a to integer using truncate
inline SSEVector4i toInt(const SSEVector4f& a) {
    return _mm_cvttps_epi32(a);
}
/// Converts \c a to integer using round
inline SSEVector4i roundToInt(const SSEVector4f& a) {
    return _mm_cvtps_epi32(a);
}

/// Reinterprets \c a as a \c SSEVector4f
inline SSEVector4f castAsFloat(const SSEVector4i& a) {
    return _mm_castsi128_ps(a);
}
/// Convert \c a to floating point
inline SSEVector4f toFloat(const SSEVector4i& a) {
    return _mm_cvtepi32_ps(a);
}

/**
 * \brief The arguments <tt>row0</tt>, <tt>row1</tt>, <tt>row2</tt> and
 * <tt>row3</tt> are \c __m128 values whose elements form the corresponding
 * rows of a 4-by-4 matrix. The matrix transposition is returned in
 * arguments <tt>row0</tt>, <tt>row1</tt>, <tt>row2</tt> and <tt>row3</tt>
 * where \c row0 now holds column 0 of the original matrix, \c row1 now
 * holds column 1 of the original matrix, and so on.
 * \author Intel Intrinsics Guide for AVX2
 */
FINLINE void transpose(SSEVector4f& row0, SSEVector4f& row1,
    SSEVector4f& row2, SSEVector4f& row3) {
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

} // namespace sse

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_SSEVECTOR_H_ */

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

#if !defined(__MTS_SSE_H)
#define __MTS_SSE_H

/* Compiler hints and SSE wrappers from Radius (https://gna.org/projects/radius) */
#if defined(__LINUX__) && defined(__GNUC__)
#define __restrict             __restrict__
#define FINLINE                inline __attribute__((always_inline))
#define NOINLINE               __attribute__((noinline))
#define MM_ALIGN16             __attribute__ ((aligned (16)))
#define EXPECT_TAKEN(a)        __builtin_expect(!!(a), true)
#define EXPECT_NOT_TAKEN(a)    __builtin_expect(!!(a), false)
#define BREAKPOINT()            do { __asm__("int3"); } while(0)
#elif defined(__OSX__) && defined(__GNUC__)
#define FINLINE                inline __attribute__((always_inline))
#define NOINLINE               __attribute__((noinline))
#define MM_ALIGN16             __attribute__ ((aligned (16)))
#define EXPECT_TAKEN(a)        __builtin_expect(!!(a), true)
#define EXPECT_NOT_TAKEN(a)    __builtin_expect(!!(a), false)
#define BREAKPOINT()            do { __asm__("int3"); } while(0)
#elif defined(__MSVC__)
#define FINLINE                __forceinline
#define NOINLINE               __declspec(noinline)
#define MM_ALIGN16             __declspec(align(16))
#define EXPECT_TAKEN(a)        (a)
#define EXPECT_NOT_TAKEN(a)    (a)
#define BREAKPOINT()           do { _asm { int 3 } } while (0)
#else
#error Unsupported compiler!
#endif
#define STACK_ALIGN16(t)       reinterpret_cast<float *>((reinterpret_cast<size_t>(t)+0xF) & ~(size_t) 0xF)

/* ========= SSE intrinsics ========= */
#ifndef MTS_SSE
#define SSE_STR	"SSE2 disabled"
#define enable_fpexcept_sse()
#define query_fpexcept_sse() 0
#define disable_fpexcept_sse()
#else
/* Include SSE intrinsics header file */
#include <emmintrin.h>

#define SSE_STR	"SSE2 enabled"
#define splat_ps(ps, i)          _mm_shuffle_ps   ((ps),(ps), (i<<6) | (i<<4) | (i<<2) | i)
#define splat_epi32(ps, i)       _mm_shuffle_epi32((ps), (i<<6) | (i<<4) | (i<<2) | i)
#define mux_ps(sel, op1, op2)    _mm_or_ps   (_mm_and_ps   ((sel), (op1)), _mm_andnot_ps   ((sel), (op2)))
#define mux_epi32(sel, op1, op2) _mm_or_si128(_mm_and_si128((sel), (op1)), _mm_andnot_si128((sel), (op2)))
#define enable_fpexcept_sse()	 _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~(_MM_MASK_INVALID | _MM_MASK_DIV_ZERO))
#define query_fpexcept_sse()	 (~_MM_GET_EXCEPTION_MASK() & (_MM_MASK_INVALID | _MM_MASK_DIV_ZERO))
#define disable_fpexcept_sse()	 _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | (_MM_MASK_INVALID | _MM_MASK_DIV_ZERO))
#define load1_epi32(i)           _mm_shuffle_epi32(_mm_cvtsi32_si128(i), 0)

#ifdef __MSVC__
	union ps_pi_t {
		__m128i pi;
		__m128	ps;
		ps_pi_t(const __m128 v) : ps(v) {}
		ps_pi_t(const __m128i v) : pi(v) {}

		static __m128i from_ps(const __m128 v)  { return ps_pi_t(v).pi; }
		static __m128  from_pi(const __m128i v) { return ps_pi_t(v).ps; }

		operator __m128()	const	{ return ps; }
		operator __m128i()	const	{ return pi; }
		operator __m128()			{ return ps; }
		operator __m128i()			{ return pi; }
	};

	#define pstoepi32(ps)            ps_pi_t::from_ps(ps)
	#define epi32tops(pi)            ps_pi_t::from_pi(pi)
#else
	#define pstoepi32(ps)            _mm_castps_si128(ps)
	#define epi32tops(pi)            _mm_castsi128_ps(pi)
#endif

#ifndef SINGLE_PRECISION
#error SSE2 only supported with single precision
#endif

MTS_NAMESPACE_BEGIN

/**
 * \headerfile mitsuba/core/sse.h mitsuba/mitsuba.h
 * \brief SSE 4-vector and useful aliases
 */
union SSEVector {
	__m128 ps;
	__m128i pi;
	float f[4];
	int32_t	i[4];
	struct { float   f0,f1,f2,f3; };
	struct { int32_t i0,i1,i2,i3; };

	inline SSEVector() {
	}

	explicit SSEVector(__m128 ps)
		: ps(ps) {
	}

	explicit SSEVector(float f0, float f1, float f2, float f3) 
		: f0(f0), f1(f1), f2(f2), f3(f3) {
	}
	
	explicit SSEVector(int32_t i0, int32_t i1, int32_t i2, int32_t i3) 
		: i0(i0), i1(i1), i2(i2), i3(i3) {
	}

	inline SSEVector &operator=(const SSEVector &vec) {
		ps = vec.ps;
		return *this;
	}
};

/**
 * \brief Some useful constant values for use with SSE
 * \headerfile mitsuba/core/sse.h mitsuba/mitsuba.h
 */
class MTS_EXPORT_CORE SSEConstants {
public:
	/// (0, 0, 0, 0)
	static const MM_ALIGN16 SSEVector zero;
	/// (1, 1, 1, 1)
	static const MM_ALIGN16 SSEVector one;
	/// (flt_max, flt_max, flt_max, flt_max)
	static const MM_ALIGN16 SSEVector max;
	/// (eps, eps, eps, eps)
	static const MM_ALIGN16 SSEVector eps;
	/// (1+eps, 1+eps, 1+eps, 1+eps)
	static const MM_ALIGN16 SSEVector op_eps;
	/// (1-eps, 1-eps, 1-eps, 1-eps)
	static const MM_ALIGN16 SSEVector om_eps;
	/// (+inf, +inf, +inf, +inf)
	static const MM_ALIGN16 SSEVector p_inf;
	/// (-inf, -inf, -inf, -inf)
	static const MM_ALIGN16 SSEVector n_inf;
	/// (0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
	static const MM_ALIGN16 SSEVector ffffffff;
};

/** Four 3D vectors as SoA (structure of arrays) */
typedef SSEVector QuadVector[3];

inline void _mm_debug_ps(const char *desc, __m128 value) {
	float dest[4];
	_mm_storeu_ps(dest, value);
	printf("%s: [%f, %f, %f, %f]\n", desc, dest[0], dest[1], dest[2], dest[3]);
}

/*
inline void _mm_debug_epi32(const char *desc, __m128i value) {
	uint32_t dest[4];
	_mm_storeu_si128((__m128i * const) dest, value);
	printf("%s: [0x%x, 0x%x, 0x%x, 0x%x]\n", desc, dest[0], dest[1], dest[2], dest[3]);
}
*/

MTS_NAMESPACE_END

#endif
	
/* ====== Performance counters (not really related to SSE) ====== */

#ifdef __GNUC__
#if defined(__i386__)
static FINLINE uint64_t rdtsc(void) {
  uint64_t x;
	 __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
	 return x;
}
#elif defined(__x86_64__)
static FINLINE uint64_t rdtsc(void) {
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ((uint64_t) lo)| (((uint64_t) hi) << 32);
}
#endif
#else
#ifndef _WIN64
__declspec(naked) static FINLINE unsigned __int64 __cdecl rdtsc(void) {
	__asm {
		rdtsc
		ret
	}
}
#else
static FINLINE __int64 rdtsc(void) {
	return __rdtsc();
}
#endif
#endif

#endif /* __MTS_SSE_H */

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

#if defined(__GXX_EXPERIMENTAL_CXX0X__)
	/* Needed to prevent a segmentation fault in the Intel C++
	   compiler on Linux (as of Nov 2012) */
	#undef __GXX_EXPERIMENTAL_CXX0X__
#endif

#if MTS_SSE
#include <mitsuba/mitsuba.h>
#include <mitsuba/core/ssemath.h>
#include <mitsuba/core/ssevector.h>

MTS_NAMESPACE_BEGIN

const MM_ALIGN16 SSEVector SSEConstants::zero          = SSEVector(0.0f);
const MM_ALIGN16 SSEVector SSEConstants::one           = SSEVector(1.0f);
const MM_ALIGN16 SSEVector SSEConstants::max           = SSEVector(std::numeric_limits<float>::max());
const MM_ALIGN16 SSEVector SSEConstants::eps           = SSEVector(Epsilon);
const MM_ALIGN16 SSEVector SSEConstants::op_eps        = SSEVector(1+Epsilon);
const MM_ALIGN16 SSEVector SSEConstants::om_eps        = SSEVector(1-Epsilon);
const MM_ALIGN16 SSEVector SSEConstants::p_inf         = SSEVector( std::numeric_limits<float>::infinity());
const MM_ALIGN16 SSEVector SSEConstants::n_inf         = SSEVector(-std::numeric_limits<float>::infinity());
const MM_ALIGN16 SSEVector SSEConstants::ffffffff      = SSEVector(0xFFFFFFFF);
const MM_ALIGN16 SSEVector SSEConstants::negation_mask = SSEVector(0x80000000);

/* SIMD (SSE2) implementation of sin, cos, exp and log

   Inspired by Intel Approximate Math library, and based on the
   corresponding algorithms of the cephes math library
*/

/* Copyright (C) 2007 Julien Pommier

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  (this is the zlib license)
*/
namespace math {
namespace constants {
	static const MM_ALIGN16 SSEVector ps_1(  1.0f);
	static const MM_ALIGN16 SSEVector ps_0p5(0.5f);
	/// 0x00800000 - smallest non denormalized float number  (~1.17549435E-38)
	static const MM_ALIGN16 SSEVector min_norm_pos( 0x00800000);
	static const MM_ALIGN16 SSEVector mant_mask(    0x7f800000);
	static const MM_ALIGN16 SSEVector inv_mant_mask(~0x7f800000);

	static const MM_ALIGN16 SSEVector sign_mask(     0x80000000);
	static const MM_ALIGN16 SSEVector inv_sign_mask(~0x80000000);

	static const MM_ALIGN16 SSEVector pi32_1(1);
	static const MM_ALIGN16 SSEVector pi32_inv1(~1);
	static const MM_ALIGN16 SSEVector pi32_2(2);
	static const MM_ALIGN16 SSEVector pi32_4(4);
	static const MM_ALIGN16 SSEVector pi32_0x7f(0x7f);

	static const MM_ALIGN16 SSEVector cephes_SQRTHF(  0.707106781186547524f);
	static const MM_ALIGN16 SSEVector cephes_log_p0(  7.0376836292E-2f);
	static const MM_ALIGN16 SSEVector cephes_log_p1(- 1.1514610310E-1f);
	static const MM_ALIGN16 SSEVector cephes_log_p2(  1.1676998740E-1f);
	static const MM_ALIGN16 SSEVector cephes_log_p3(- 1.2420140846E-1f);
	static const MM_ALIGN16 SSEVector cephes_log_p4(+ 1.4249322787E-1f);
	static const MM_ALIGN16 SSEVector cephes_log_p5(- 1.6668057665E-1f);
	static const MM_ALIGN16 SSEVector cephes_log_p6(+ 2.0000714765E-1f);
	static const MM_ALIGN16 SSEVector cephes_log_p7(- 2.4999993993E-1f);
	static const MM_ALIGN16 SSEVector cephes_log_p8(+ 3.3333331174E-1f);
	static const MM_ALIGN16 SSEVector cephes_log_q1(- 2.12194440e-4f);
	static const MM_ALIGN16 SSEVector cephes_log_q2(  0.693359375f);

	static const MM_ALIGN16 SSEVector exp_hi( 88.3762626647949f);
	static const MM_ALIGN16 SSEVector exp_lo(-88.3762626647949f);

	static const MM_ALIGN16 SSEVector cephes_LOG2EF( 1.44269504088896341f);
	static const MM_ALIGN16 SSEVector cephes_exp_C1( 0.693359375f);
	static const MM_ALIGN16 SSEVector cephes_exp_C2(-2.12194440e-4f);

	static const MM_ALIGN16 SSEVector cephes_exp_p0(1.9875691500E-4f);
	static const MM_ALIGN16 SSEVector cephes_exp_p1(1.3981999507E-3f);
	static const MM_ALIGN16 SSEVector cephes_exp_p2(8.3334519073E-3f);
	static const MM_ALIGN16 SSEVector cephes_exp_p3(4.1665795894E-2f);
	static const MM_ALIGN16 SSEVector cephes_exp_p4(1.6666665459E-1f);
	static const MM_ALIGN16 SSEVector cephes_exp_p5(5.0000001201E-1f);

	static const MM_ALIGN16 SSEVector minus_cephes_DP1(-0.78515625f);
	static const MM_ALIGN16 SSEVector minus_cephes_DP2(-2.4187564849853515625e-4f);
	static const MM_ALIGN16 SSEVector minus_cephes_DP3(-3.77489497744594108e-8f);
	static const MM_ALIGN16 SSEVector sincof_p0(-1.9515295891E-4f);
	static const MM_ALIGN16 SSEVector sincof_p1( 8.3321608736E-3f);
	static const MM_ALIGN16 SSEVector sincof_p2(-1.6666654611E-1f);
	static const MM_ALIGN16 SSEVector coscof_p0( 2.443315711809948E-005f);
	static const MM_ALIGN16 SSEVector coscof_p1(-1.388731625493765E-003f);
	static const MM_ALIGN16 SSEVector coscof_p2(4.166664568298827E-002f);
	static const MM_ALIGN16 SSEVector cephes_FOPI(1.27323954473516f);

	static const MM_ALIGN16 SSEVector am_log_p0(-7.89580278884799154124e-1f);
	static const MM_ALIGN16 SSEVector am_log_p1( 1.63866645699558079767e1f);
	static const MM_ALIGN16 SSEVector am_log_p2(-6.41409952958715622951e1f);
	static const MM_ALIGN16 SSEVector am_log_q0(-3.56722798256324312549e1f);
	static const MM_ALIGN16 SSEVector am_log_q1( 3.12093766372244180303e2f);
	static const MM_ALIGN16 SSEVector am_log_q2(-7.69691943550460008604e2f);
	static const MM_ALIGN16 SSEVector am_log_c0( 0.693147180559945f);
	static const MM_ALIGN16 SSEVector am_log2_c0(1.44269504088896340735992f);

	static const MM_ALIGN16 SSEVector am_exp2_hi( 127.4999961853f);
	static const MM_ALIGN16 SSEVector am_exp2_lo(-127.4999961853f);
	static const MM_ALIGN16 SSEVector am_exp2_p0(2.30933477057345225087e-2f);
	static const MM_ALIGN16 SSEVector am_exp2_p1(2.02020656693165307700e1f);
	static const MM_ALIGN16 SSEVector am_exp2_p2(1.51390680115615096133e3f);
	static const MM_ALIGN16 SSEVector am_exp2_q0(2.33184211722314911771e2f);
	static const MM_ALIGN16 SSEVector am_exp2_q1(4.36821166879210612817e3f);
};

/* natural logarithm computed for 4 simultaneous float
   return NaN for x <= 0
*/
__m128 log_ps(__m128 x) {
    typedef __m128 v4sf;
    typedef __m128i v4si;

    v4si emm0;
    v4sf one = constants::ps_1.ps;

    v4sf invalid_mask = _mm_cmple_ps(x, _mm_setzero_ps());

    x = _mm_max_ps(x, constants::min_norm_pos.ps);  // cut off denormalized stuff

    emm0 = _mm_srli_epi32(_mm_castps_si128(x), 23);
    // keep only the fractional part
    x = _mm_and_ps(x, constants::inv_mant_mask.ps);
    x = _mm_or_ps(x,  constants::ps_0p5.ps);

    emm0 = _mm_sub_epi32(emm0, constants::pi32_0x7f.pi);
    v4sf e = _mm_cvtepi32_ps(emm0);

    e = _mm_add_ps(e, one);

    /* part2:
       if( x < SQRTHF ) {
         e -= 1;
         x = x + x - 1.0;
       } else { x = x - 1.0; }
    */
    v4sf mask = _mm_cmplt_ps(x, constants::cephes_SQRTHF.ps);
    v4sf tmp = _mm_and_ps(x, mask);
    x = _mm_sub_ps(x, one);
    e = _mm_sub_ps(e, _mm_and_ps(one, mask));
    x = _mm_add_ps(x, tmp);

    v4sf z = _mm_mul_ps(x,x);

    v4sf y = constants::cephes_log_p0.ps;
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, constants::cephes_log_p1.ps);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, constants::cephes_log_p2.ps);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, constants::cephes_log_p3.ps);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, constants::cephes_log_p4.ps);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, constants::cephes_log_p5.ps);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, constants::cephes_log_p6.ps);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, constants::cephes_log_p7.ps);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, constants::cephes_log_p8.ps);
    y = _mm_mul_ps(y, x);

    y = _mm_mul_ps(y, z);

    tmp = _mm_mul_ps(e, constants::cephes_log_q1.ps);
    y = _mm_add_ps(y, tmp);

    tmp = _mm_mul_ps(z, constants::ps_0p5.ps);
    y = _mm_sub_ps(y, tmp);

    tmp = _mm_mul_ps(e, constants::cephes_log_q2.ps);
    x = _mm_add_ps(x, y);
    x = _mm_add_ps(x, tmp);
    x = _mm_or_ps(x, invalid_mask); // negative arg will be NAN
    return x;
}

__m128 exp_ps(__m128 x) {
    typedef __m128 v4sf;
    typedef __m128i v4si;

    v4sf tmp = _mm_setzero_ps(), fx;
    v4si emm0;
    v4sf one = constants::ps_1.ps;

    x = _mm_min_ps(x, constants::exp_hi.ps);
    x = _mm_max_ps(x, constants::exp_lo.ps);

    /* express exp(x) as exp(g + n*log(2)) */
    fx = _mm_mul_ps(x,  constants::cephes_LOG2EF.ps);
    fx = _mm_add_ps(fx, constants::ps_0p5.ps);

    /* how to perform a floorf with SSE: just below */
    emm0 = _mm_cvttps_epi32(fx);
    tmp  = _mm_cvtepi32_ps(emm0);
    /* if greater, substract 1 */
    v4sf mask = _mm_cmpgt_ps(tmp, fx);
    mask = _mm_and_ps(mask, one);
    fx = _mm_sub_ps(tmp, mask);

    tmp = _mm_mul_ps(fx, constants::cephes_exp_C1.ps);
    v4sf z = _mm_mul_ps(fx, constants::cephes_exp_C2.ps);
    x = _mm_sub_ps(x, tmp);
    x = _mm_sub_ps(x, z);

    z = _mm_mul_ps(x,x);

    v4sf y = constants::cephes_exp_p0.ps;
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, constants::cephes_exp_p1.ps);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, constants::cephes_exp_p2.ps);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, constants::cephes_exp_p3.ps);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, constants::cephes_exp_p4.ps);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, constants::cephes_exp_p5.ps);
    y = _mm_mul_ps(y, z);
    y = _mm_add_ps(y, x);
    y = _mm_add_ps(y, one);

    /* build 2^n */
    emm0 = _mm_cvttps_epi32(fx);
    emm0 = _mm_add_epi32(emm0, constants::pi32_0x7f.pi);
    emm0 = _mm_slli_epi32(emm0, 23);
    v4sf pow2n = _mm_castsi128_ps(emm0);
    y = _mm_mul_ps(y, pow2n);
    return y;
}

/* evaluation of 4 sines at onces, using only SSE2.

   The code is the exact rewriting of the cephes sinf function.
   Precision is excellent as long as x < 8192 (I did not bother to
   take into account the special handling they have for greater values
   -- it does not return garbage for arguments over 8192, though, but
   the extra precision is missing).

   Note that it is such that sinf((float)M_PI) = 8.74e-8, which is the
   surprising but correct result.

   Performance is also surprisingly good, 1.33 times faster than the
   macos vsinf SSE2 function, and 1.5 times faster than the
   __vrs4_sinf of amd's ACML (which is only available in 64 bits). Not
   too bad for an SSE1 function (with no special tuning) !
   However the latter libraries probably have a much better handling of NaN,
   Inf, denormalized and other special arguments..

   On my core 1 duo, the execution of this function takes approximately 95 cycles.

   From what I have observed on the experiments with Intel AMath lib, switching to an
   SSE2 version would improve the perf by only 10%.

   Since it is based on SSE intrinsics, it has to be compiled at -O2 to
   deliver full speed.
*/
__m128 sin_ps(__m128 x) { // any x
    typedef __m128 v4sf;
    typedef __m128i v4si;

    v4sf xmm1, xmm2 = _mm_setzero_ps(), xmm3, sign_bit, y;

    v4si emm0, emm2;
    sign_bit = x;
    /* take the absolute value */
    x = _mm_and_ps(x, constants::inv_mant_mask.ps);
    /* extract the sign bit (upper one) */
    sign_bit = _mm_and_ps(sign_bit, constants::sign_mask.ps);

    /* scale by 4/Pi */
    y = _mm_mul_ps(x, constants::cephes_FOPI.ps);

    /* store the integer part of y in mm0 */
    emm2 = _mm_cvttps_epi32(y);
    /* j=(j+1) & (~1) (see the cephes sources) */
    emm2 = _mm_add_epi32(emm2, constants::pi32_1.pi);
    emm2 = _mm_and_si128(emm2, constants::pi32_inv1.pi);
    y = _mm_cvtepi32_ps(emm2);
    /* get the swap sign flag */
    emm0 = _mm_and_si128(emm2, constants::pi32_4.pi);
    emm0 = _mm_slli_epi32(emm0, 29);
    /* get the polynom selection mask
       there is one polynom for 0 <= x <= Pi/4
       and another one for Pi/4<x<=Pi/2

       Both branches will be computed.
    */
    emm2 = _mm_and_si128(emm2, constants::pi32_2.pi);
    emm2 = _mm_cmpeq_epi32(emm2, _mm_setzero_si128());

    v4sf swap_sign_bit = _mm_castsi128_ps(emm0);
    v4sf poly_mask = _mm_castsi128_ps(emm2);
    sign_bit = _mm_xor_ps(sign_bit, swap_sign_bit);

    /* The magic pass: "Extended precision modular arithmetic"
    x = ((x - y * DP1) - y * DP2) - y * DP3; */
    xmm1 = constants::minus_cephes_DP1.ps;
    xmm2 = constants::minus_cephes_DP2.ps;
    xmm3 = constants::minus_cephes_DP3.ps;
    xmm1 = _mm_mul_ps(y, xmm1);
    xmm2 = _mm_mul_ps(y, xmm2);
    xmm3 = _mm_mul_ps(y, xmm3);
    x = _mm_add_ps(x, xmm1);
    x = _mm_add_ps(x, xmm2);
    x = _mm_add_ps(x, xmm3);

    /* Evaluate the first polynom  (0 <= x <= Pi/4) */
    y = constants::coscof_p0.ps;
    v4sf z = _mm_mul_ps(x,x);

    y = _mm_mul_ps(y, z);
    y = _mm_add_ps(y, constants::coscof_p1.ps);
    y = _mm_mul_ps(y, z);
    y = _mm_add_ps(y, constants::coscof_p2.ps);
    y = _mm_mul_ps(y, z);
    y = _mm_mul_ps(y, z);
    v4sf tmp = _mm_mul_ps(z, constants::ps_0p5.ps);
    y = _mm_sub_ps(y, tmp);
    y = _mm_add_ps(y, constants::ps_1.ps);

    /* Evaluate the second polynom  (Pi/4 <= x <= 0) */

    v4sf y2 = constants::sincof_p0.ps;
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_add_ps(y2, constants::sincof_p1.ps);
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_add_ps(y2, constants::sincof_p2.ps);
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_mul_ps(y2, x);
    y2 = _mm_add_ps(y2, x);

    /* select the correct result from the two polynoms */
    xmm3 = poly_mask;
    y2 = _mm_and_ps(xmm3, y2); //, xmm3);
    y = _mm_andnot_ps(xmm3, y);
    y = _mm_add_ps(y,y2);
    /* update the sign */
    y = _mm_xor_ps(y, sign_bit);

    return y;
}

/* almost the same as sin_ps */
__m128 cos_ps(__m128 x) { // any x
    typedef __m128 v4sf;
    typedef __m128i v4si;

    v4sf xmm1, xmm2 = _mm_setzero_ps(), xmm3, y;
    v4si emm0, emm2;
    /* take the absolute value */
    x = _mm_and_ps(x, constants::inv_sign_mask.ps);

    /* scale by 4/Pi */
    y = _mm_mul_ps(x, constants::cephes_FOPI.ps);

    /* store the integer part of y in mm0 */
    emm2 = _mm_cvttps_epi32(y);
    /* j=(j+1) & (~1) (see the cephes sources) */
    emm2 = _mm_add_epi32(emm2, constants::pi32_1.pi);
    emm2 = _mm_and_si128(emm2, constants::pi32_inv1.pi);
    y = _mm_cvtepi32_ps(emm2);

    emm2 = _mm_sub_epi32(emm2, constants::pi32_2.pi);

    /* get the swap sign flag */
    emm0 = _mm_andnot_si128(emm2, constants::pi32_4.pi);
    emm0 = _mm_slli_epi32(emm0, 29);
    /* get the polynom selection mask */
    emm2 = _mm_and_si128(emm2, constants::pi32_2.pi);
    emm2 = _mm_cmpeq_epi32(emm2, _mm_setzero_si128());

    v4sf sign_bit = _mm_castsi128_ps(emm0);
    v4sf poly_mask = _mm_castsi128_ps(emm2);

    /* The magic pass: "Extended precision modular arithmetic"
    x = ((x - y * DP1) - y * DP2) - y * DP3; */
    xmm1 = constants::minus_cephes_DP1.ps;
    xmm2 = constants::minus_cephes_DP2.ps;
    xmm3 = constants::minus_cephes_DP3.ps;
    xmm1 = _mm_mul_ps(y, xmm1);
    xmm2 = _mm_mul_ps(y, xmm2);
    xmm3 = _mm_mul_ps(y, xmm3);
    x = _mm_add_ps(x, xmm1);
    x = _mm_add_ps(x, xmm2);
    x = _mm_add_ps(x, xmm3);

    /* Evaluate the first polynom  (0 <= x <= Pi/4) */
    y = constants::coscof_p0.ps;
    v4sf z = _mm_mul_ps(x,x);

    y = _mm_mul_ps(y, z);
    y = _mm_add_ps(y, constants::coscof_p1.ps);
    y = _mm_mul_ps(y, z);
    y = _mm_add_ps(y, constants::coscof_p2.ps);
    y = _mm_mul_ps(y, z);
    y = _mm_mul_ps(y, z);
    v4sf tmp = _mm_mul_ps(z, constants::ps_0p5.ps);
    y = _mm_sub_ps(y, tmp);
    y = _mm_add_ps(y, constants::ps_1.ps);

    /* Evaluate the second polynom  (Pi/4 <= x <= 0) */

    v4sf y2 = constants::sincof_p0.ps;
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_add_ps(y2, constants::sincof_p1.ps);
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_add_ps(y2, constants::sincof_p2.ps);
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_mul_ps(y2, x);
    y2 = _mm_add_ps(y2, x);

    /* select the correct result from the two polynoms */
    xmm3 = poly_mask;
    y2 = _mm_and_ps(xmm3, y2); //, xmm3);
    y = _mm_andnot_ps(xmm3, y);
    y = _mm_add_ps(y,y2);
    /* update the sign */
    y = _mm_xor_ps(y, sign_bit);

    return y;
}

/* since sin_ps and cos_ps are almost identical, sincos_ps could replace both of them..
   it is almost as fast, and gives you a free cosine with your sine */
void sincos_ps(__m128 x, __m128* s, __m128* c) {
    typedef __m128 v4sf;
    typedef __m128i v4si;

    v4sf xmm1, xmm2, xmm3 = _mm_setzero_ps(), sign_bit_sin, y;
    v4si emm0, emm2, emm4;
    sign_bit_sin = x;
    /* take the absolute value */
    x = _mm_and_ps(x, constants::inv_sign_mask.ps);
    /* extract the sign bit (upper one) */
    sign_bit_sin = _mm_and_ps(sign_bit_sin, constants::sign_mask.ps);

    /* scale by 4/Pi */
    y = _mm_mul_ps(x, constants::cephes_FOPI.ps);

    /* store the integer part of y in emm2 */
    emm2 = _mm_cvttps_epi32(y);

    /* j=(j+1) & (~1) (see the cephes sources) */
    emm2 = _mm_add_epi32(emm2, constants::pi32_1.pi);
    emm2 = _mm_and_si128(emm2, constants::pi32_inv1.pi);
    y = _mm_cvtepi32_ps(emm2);

    emm4 = emm2;

    /* get the swap sign flag for the sine */
    emm0 = _mm_and_si128(emm2, constants::pi32_4.pi);
    emm0 = _mm_slli_epi32(emm0, 29);
    v4sf swap_sign_bit_sin = _mm_castsi128_ps(emm0);

    /* get the polynom selection mask for the sine*/
    emm2 = _mm_and_si128(emm2, constants::pi32_2.pi);
    emm2 = _mm_cmpeq_epi32(emm2, _mm_setzero_si128());
    v4sf poly_mask = _mm_castsi128_ps(emm2);

    /* The magic pass: "Extended precision modular arithmetic"
    x = ((x - y * DP1) - y * DP2) - y * DP3; */
    xmm1 = constants::minus_cephes_DP1.ps;
    xmm2 = constants::minus_cephes_DP2.ps;
    xmm3 = constants::minus_cephes_DP3.ps;
    xmm1 = _mm_mul_ps(y, xmm1);
    xmm2 = _mm_mul_ps(y, xmm2);
    xmm3 = _mm_mul_ps(y, xmm3);
    x = _mm_add_ps(x, xmm1);
    x = _mm_add_ps(x, xmm2);
    x = _mm_add_ps(x, xmm3);

    emm4 = _mm_sub_epi32(emm4, constants::pi32_2.pi);
    emm4 = _mm_andnot_si128(emm4, constants::pi32_4.pi);
    emm4 = _mm_slli_epi32(emm4, 29);
    v4sf sign_bit_cos = _mm_castsi128_ps(emm4);

    sign_bit_sin = _mm_xor_ps(sign_bit_sin, swap_sign_bit_sin);


    /* Evaluate the first polynom  (0 <= x <= Pi/4) */
    v4sf z = _mm_mul_ps(x,x);
    y = constants::coscof_p0.ps;

    y = _mm_mul_ps(y, z);
    y = _mm_add_ps(y, constants::coscof_p1.ps);
    y = _mm_mul_ps(y, z);
    y = _mm_add_ps(y, constants::coscof_p2.ps);
    y = _mm_mul_ps(y, z);
    y = _mm_mul_ps(y, z);
    v4sf tmp = _mm_mul_ps(z, constants::ps_0p5.ps);
    y = _mm_sub_ps(y, tmp);
    y = _mm_add_ps(y, constants::ps_1.ps);

    /* Evaluate the second polynom  (Pi/4 <= x <= 0) */

    v4sf y2 = constants::sincof_p0.ps;
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_add_ps(y2, constants::sincof_p1.ps);
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_add_ps(y2, constants::sincof_p2.ps);
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_mul_ps(y2, x);
    y2 = _mm_add_ps(y2, x);

    /* select the correct result from the two polynoms */
    xmm3 = poly_mask;
    v4sf ysin2 = _mm_and_ps(xmm3, y2);
    v4sf ysin1 = _mm_andnot_ps(xmm3, y);
    y2 = _mm_sub_ps(y2,ysin2);
    y = _mm_sub_ps(y, ysin1);

    xmm1 = _mm_add_ps(ysin1,ysin2);
    xmm2 = _mm_add_ps(y,y2);

    /* update the sign */
    *s = _mm_xor_ps(xmm1, sign_bit_sin);
    *c = _mm_xor_ps(xmm2, sign_bit_cos);
}

__m128 fastlog_ps(__m128 x) {
	typedef SSEVector4f V4f;
	typedef SSEVector4i V4i;

	// Constants
	const V4f min_normal(constants::min_norm_pos.ps);
	const V4f inv_mantissa_mask(constants::inv_mant_mask.ps);
	const V4f const_1(constants::ps_1.ps);
	const V4i const_127(constants::pi32_0x7f.pi);

	const V4f log_p0(constants::am_log_p0.ps);
	const V4f log_p1(constants::am_log_p1.ps);
	const V4f log_p2(constants::am_log_p2.ps);

	const V4f log_q0(constants::am_log_q0.ps);
	const V4f log_q1(constants::am_log_q1.ps);
	const V4f log_q2(constants::am_log_q2.ps);

	const V4f log_c0(constants::am_log_c0.ps);


	// Kill invalid values (ignoring NaN and Inf)
	const V4f x0 = max(x, min_normal);

	// Kill the exponent and combine with the exponent of 1.0f to get the
	// actual embedded mantissa as a valid floating point value:
	// a value in the range [1.0, 2.0)
	const V4f mantissa = (x0 & inv_mantissa_mask) | const_1;

	const V4f v_min1  = mantissa - const_1;
	const V4f v_plus1 = mantissa + const_1;

	// Extract the original exponent and undo the bias
	const V4i biasedExponent = srl(castAsInt(x0), 23);
	const V4f origExponent = toFloat(biasedExponent - const_127);

	V4f vFrac = v_min1 * rcp(v_plus1); // Is it worth it to use rcp_nr?
	vFrac += vFrac;
	const V4f vFracSqr = vFrac * vFrac;

	// Evaluate the polynomial
	const V4f polyP = ((((log_p0 * vFracSqr) + log_p1) * vFracSqr)
											 + log_p2) * vFracSqr;
	const V4f polyQ =  (((log_q0 * vFracSqr) + log_q1) * vFracSqr) + log_q2;

	const V4f poly = polyP * rcp(polyQ); // Use rcp_nr?
	const V4f logApprox = poly * vFrac;

	// Scale by log(2) to get the natural logarithm of the exponent part
	const V4f logExpPart = origExponent * log_c0;

	// Combine the different parts
	const V4f result = logApprox + vFrac + logExpPart;

	return result;
}

__m128 fastpow_ps(__m128 x, __m128 y) {
	typedef SSEVector4f V4f;
	typedef SSEVector4i V4i;

	// Constants
	const V4f min_normal(constants::min_norm_pos.ps);
	const V4f inv_mantissa_mask(constants::inv_mant_mask.ps);
	const V4f const_1(constants::ps_1.ps);
	const V4i const_127(constants::pi32_0x7f.pi);

	const V4f log_p0(constants::am_log_p0.ps);
	const V4f log_p1(constants::am_log_p1.ps);
	const V4f log_p2(constants::am_log_p2.ps);

	const V4f log_q0(constants::am_log_q0.ps);
	const V4f log_q1(constants::am_log_q1.ps);
	const V4f log_q2(constants::am_log_q2.ps);

	const V4f log2_c0(constants::am_log2_c0.ps);


	// Remember negative values
	const V4f negative_mask(V4f::zero() < V4f(x));

	// Cutoff denormalized stuff (preserving NaN and Infinity)
	const V4f x0 = max(x, min_normal);

	// First step: compute log(x)

	// Kill the exponent and combine with the exponent of 1.0f to get the
	// actual embedded mantissa as a valid floating point value:
	// a value in the range [1.0, 2.0)
	const V4f mantissa = (x0 & inv_mantissa_mask) | const_1;

	const V4f v_min1  = mantissa - const_1;
	const V4f v_plus1 = mantissa + const_1;

	// Extract the original exponent and undo the bias
	const V4i biasedExponent = srl(castAsInt(x0), 23);
	const V4f origExponent = toFloat(biasedExponent - const_127);

	V4f vFrac = v_min1 * rcp(v_plus1); // Is it worth it to use rcp_nr?
	vFrac += vFrac;
	const V4f vFracSqr = vFrac * vFrac;

	// Evaluate the polynomial
	const V4f polyP = ((((log_p0 * vFracSqr) + log_p1) *
								   vFracSqr) + log_p2) * vFracSqr;
	const V4f polyQ =  (((log_q0 * vFracSqr) + log_q1) *
								   vFracSqr) + log_q2;
	const V4f logApprox = (polyP * rcp(polyQ)) * vFrac;

	// y * log2(x)
	V4f exponent = (logApprox * log2_c0) + ((vFrac * log2_c0) + origExponent);
	exponent *= y;


	// Constants for the exponential
	const V4f const_0p5(constants::ps_0p5.ps);

	const V4f exp2_hi(constants::am_exp2_hi.ps);
	const V4f exp2_lo(constants::am_exp2_lo.ps);

	const V4f exp2_p0(constants::am_exp2_p0.ps);
	const V4f exp2_p1(constants::am_exp2_p1.ps);
	const V4f exp2_p2(constants::am_exp2_p2.ps);

	const V4f exp2_q0(constants::am_exp2_q0.ps);
	const V4f exp2_q1(constants::am_exp2_q1.ps);

	// Clamp the exponent
	exponent = max(min(exponent, exp2_hi), exp2_lo);

	// More floating point tricks: normalize the mantissa to [1.0 - 1.5]
	const V4f normExponent = exponent + const_0p5;

	// Build the biased exponent
	const V4f expNegExponentMask = cmpnlt(V4f::zero(), normExponent);
	const V4f expNormalization = expNegExponentMask & const_1;
	const V4f truncExp = roundTruncate(normExponent);
	const V4f resExp = truncExp - expNormalization;
	V4i biasedExp = toInt(resExp) + const_127;
	biasedExp = sll(biasedExp, 23);
	const V4f exponentPart = castAsFloat(biasedExp) & negative_mask;

	// Get the fractional part of the exponent
	exponent -= resExp;
	const V4f exponentSqr = exponent * exponent;

	// Exp polynomial
	const V4f EPolyP = ((((exp2_p0 * exponentSqr) + exp2_p1) *
									 exponentSqr) + exp2_p2) * exponent;
	const V4f EPolyQ =   ((exp2_q0 * exponentSqr) + exp2_q1) - EPolyP;
	V4f expApprox = EPolyP * rcp(EPolyQ);
	expApprox += expApprox;
	expApprox += const_1;

	V4f result = expApprox * exponentPart;
	return result;
}

}

MTS_NAMESPACE_END

#endif /* MTS_SSE */

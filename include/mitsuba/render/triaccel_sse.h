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
#if !defined(__MITSUBA_RENDER_TRIACCEL_SSE_H_)
#define __MITSUBA_RENDER_TRIACCEL_SSE_H_

#include <mitsuba/render/trimesh.h>

MTS_NAMESPACE_BEGIN

FINLINE __m128 rayIntersectPacket(const TriAccel &tri, const RayPacket4 &packet,
	__m128 mint, __m128 maxt, __m128 inactive, Intersection4 &its) {
	static const MM_ALIGN16 int waldModulo[4] = { 1, 2, 0, 1 };
	const int ku = waldModulo[tri.k], kv = waldModulo[tri.k+1];

	/* Get the u and v components */
	const __m128
		o_u = packet.o[ku].ps, o_v = packet.o[kv].ps, o_k = packet.o[tri.k].ps,
		d_u = packet.d[ku].ps, d_v = packet.d[kv].ps, d_k = packet.d[tri.k].ps;

	/* Extract data from the first cache line */
	const __m128
		line1 = _mm_load_ps((const float *) &tri),
		n_u = splat_ps(line1, 1),
		n_v = splat_ps(line1, 2),
		n_d = splat_ps(line1, 3);

	const __m128
		ounu = _mm_mul_ps(o_u, n_u),
		ovnv = _mm_mul_ps(o_v, n_v),
		dunu = _mm_mul_ps(d_u, n_u),
		dvnv = _mm_mul_ps(d_v, n_v);

	/* Calculate the plane intersection (Typo in the thesis?) */
	const __m128
		num   = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(n_d, ounu), ovnv), o_k),
		denom = _mm_add_ps(_mm_add_ps(dunu, dvnv), d_k);

	const __m128
		t = _mm_div_ps(num, denom);

	__m128 hasIts =
		_mm_andnot_ps(inactive, _mm_and_ps(_mm_cmpgt_ps(maxt, t), _mm_cmpgt_ps(t, mint)));

	if (_mm_movemask_ps(hasIts) == 0)
		return hasIts;

	/* Extract data from the second cache line */
	const __m128
		line2 = _mm_load_ps(&tri.a_u),
		a_u   = splat_ps(line2, 0),
		a_v   = splat_ps(line2, 1),
		b_nu  = splat_ps(line2, 2),
		b_nv  = splat_ps(line2, 3);

	const __m128
		hu = _mm_add_ps(o_u, _mm_sub_ps(_mm_mul_ps(t, d_u), a_u)),
		hv = _mm_add_ps(o_v, _mm_sub_ps(_mm_mul_ps(t, d_v), a_v));

	/* Extract data from the third cache line */
	const __m128
		line3     = _mm_load_ps(&tri.c_nu),
		c_nu      = splat_ps(line3, 0),
		c_nv      = splat_ps(line3, 1);
	const __m128i
		primIndex = splat_epi32(pstoepi32(line3), 3),
		shapeIndex = splat_epi32(pstoepi32(line3), 2);

	const __m128
		u = _mm_add_ps(_mm_mul_ps(hv, b_nu), _mm_mul_ps(hu, b_nv)),
		v = _mm_add_ps(_mm_mul_ps(hu, c_nu), _mm_mul_ps(hv, c_nv));

	const __m128
		zero = _mm_setzero_ps(),
		term1 = _mm_cmpge_ps(u, zero),
		term2 = _mm_cmpge_ps(v, zero),
		term3 = _mm_add_ps(u, v);

	const __m128
		term4 = _mm_and_ps(term1, term2),
		term5 = _mm_cmpge_ps(SSEConstants::one.ps, term3);

	hasIts = _mm_and_ps(hasIts, _mm_and_ps(term4, term5));

	if (_mm_movemask_ps(hasIts) == 0)
		return hasIts;

	its.t.ps  = mux_ps(hasIts, t, its.t.ps);
	its.u.ps  = mux_ps(hasIts, u, its.u.ps);
	its.v.ps  = mux_ps(hasIts, v, its.v.ps);
	its.primIndex.pi = mux_epi32(pstoepi32(hasIts), primIndex, its.primIndex.pi);
	its.shapeIndex.pi = mux_epi32(pstoepi32(hasIts), shapeIndex, its.shapeIndex.pi);

	return hasIts;
}

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_TRIACCEL_SSE_H_ */


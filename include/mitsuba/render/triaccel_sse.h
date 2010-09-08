/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__TRIACCEL_SSE_H)
#define __TRIACCEL_SSE_H

#include <mitsuba/render/trimesh.h>

MTS_NAMESPACE_BEGIN

/**
 * Pre-computed triangle representation using Ingo Wald's TriAccel layout.
 * This is a special version, which can store and simultaneously intersect 
 * against four triangles using SSE instructions.
 */
struct TriAccel4 {
	uint8_t k;
	uint8_t nonTriFlag; // Flags any non-triangle shapes that may be referenced
	uint16_t shapeIndex[4];
	uint16_t indirectionCount;
	uint32_t indirectionIndex;
	float nu[4];
	float nv[4];
	float nd[4];
	float au[4];
	float av[4];
	float bnu[4];
	float bnv[4];
	float cnu[4];
	float cnv[4];
	uint32_t index[4];

	/// Create from vertex data. Returns the number of failures
	inline int load(const Point *A, const Point *B, const Point *C, 
		const uint32_t *shapeIndex, const uint32_t *index);

	/// Fast ray-triangle intersection test
	inline bool rayIntersect(const __m128 o, const __m128 d, float _mint, 
		float _maxt, float *_t, float *_u, float *_v, unsigned int &_shapeIndex,
		unsigned int &_index);
};

FINLINE __m128 TriAccel::rayIntersectPacket(const RayPacket4 &packet, 
	__m128 mint, __m128 maxt, __m128 inactive, Intersection4 &its) const {
	static const int waldModulo[4] = { 1, 2, 0, 1 };
	const int ku = waldModulo[k], kv = waldModulo[k+1];

	/* Get the u and v components */
	const __m128 
		o_u = packet.o[ku].ps, o_v = packet.o[kv].ps, o_k = packet.o[k].ps,
		d_u = packet.d[ku].ps, d_v = packet.d[kv].ps, d_k = packet.d[k].ps;

	/* Extract data from the first cache line */
	const __m128 
		line1 = _mm_load_ps((const float *) this),
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
		line2 = _mm_load_ps(&this->a_u),
		a_u   = splat_ps(line2, 0),
		a_v   = splat_ps(line2, 1),
		b_nu  = splat_ps(line2, 2),
		b_nv  = splat_ps(line2, 3);

	const __m128 
		hu = _mm_add_ps(o_u, _mm_sub_ps(_mm_mul_ps(t, d_u), a_u)),
		hv = _mm_add_ps(o_v, _mm_sub_ps(_mm_mul_ps(t, d_v), a_v));

	/* Extract data from the third cache line */
	const __m128 
		line3     = _mm_load_ps(&this->c_nu),
		c_nu      = splat_ps(line3, 0),
		c_nv      = splat_ps(line3, 1);
	const __m128i
		primIndex = splat_epi32(pstoepi32(line3), 2),
		shapeIndex = splat_epi32(pstoepi32(line3), 3);

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

inline int TriAccel4::load(const Point *A, const Point *B, const Point *C,
		const uint32_t *_shapeIndex, const uint32_t *_index) {
	static const int waldModulo[4] = { 1, 2, 0, 1 };
	int factor = 1, failures = 0;
	k = 0;

	for (int i=0; i<4; i++) {
		Vector b = C[i]-A[i], c = B[i]-A[i], N = cross(c, b);

		/* Determine the closest projection axis */
		int kk = 0;
		for (int j=0; j<3; j++) {
			if (std::abs(N[j]) > std::abs(N[kk]))
				kk = j;
		}

		int u = waldModulo[kk],
			v = waldModulo[kk+1];
		const Float n_k = N[kk],
			denom = b[u]*c[v] - b[v]*c[u];

		if (denom == 0) {
			if ((_shapeIndex[i] != 0 || _index[i] != 0) && _index[i] != KNoTriangleFlag)
				failures++;
		}

		/* Pre-compute intersection calculation 
		   constants */
		nu[i]   =  N[u] / n_k;
		nv[i]   =  N[v] / n_k;
		nd[i]   =  dot(A[i], N) / n_k;
		bnu[i]  =  b[u] / denom;
		bnv[i]  = -b[v] / denom;
		au[i]   =  A[i][u];
		av[i]   =  A[i][v];
		cnu[i]  =  c[v] / denom;
		cnv[i]  = -c[u] / denom;
		shapeIndex[i] = _shapeIndex[i];
		index[i] = _index[i];
		k += factor * kk;
		factor *= 3;
	}
	return failures;
}

inline bool TriAccel4::rayIntersect(const __m128 o, const __m128 d, float _mint, 
	float _maxt, float *_t, float *_u, float *_v, unsigned int &_shapeIndex, 
	unsigned int &_index) {
	__m128 o_k, o_u, o_v, d_k, d_u, d_v;
	/* Arrange the ray according to the projection axes of the
	   four packed triangles. This requires a *good* compiler
	   (read: GCC 4.2 or ICC)
	*/

#ifdef MTS_DEBUG_FP
	/* If debugging, turn off FP exceptions while testing for intersections
	   kd-tree (the code makes use of NaNs) */
	disable_fpexcept();
#endif

	switch (k) {
		case 0:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,0,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,1,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,2,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,0,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,1,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,2,2));
			break;
		case 1:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,0,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,1,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,2,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,0,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,1,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,2,0));
			break;
		case 2:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,0,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,1,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,2,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,0,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,1,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,2,1));
			break;
		case 3:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,1,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,2,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,0,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,1,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,2,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,0,2));
			break;
		case 4:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,1,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,2,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,0,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,1,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,2,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,0,0));
			break;
		case 5:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,1,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,2,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,0,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,1,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,2,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,0,1));
			break;
		case 6:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,2,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,0,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,1,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,2,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,0,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,1,2));
			break;
		case 7:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,2,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,0,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,1,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,2,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,0,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,1,0));
			break;
		case 8:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,2,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,0,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,1,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,2,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,0,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,1,1));
			break;
		case 9:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,0,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,1,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,2,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,0,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,1,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,2,2));
			break;
		case 10:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,0,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,1,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,2,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,0,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,1,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,2,0));
			break;
		case 11:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,0,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,1,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,2,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,0,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,1,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,2,1));
			break;
		case 12:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,1,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,2,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,0,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,1,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,2,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,0,2));
			break;
		case 13:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,1,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,2,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,0,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,1,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,2,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,0,0));
			break;
		case 14:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,1,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,2,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,0,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,1,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,2,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,0,1));
			break;
		case 15:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,2,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,0,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,1,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,2,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,0,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,1,2));
			break;
		case 16:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,2,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,0,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,1,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,2,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,0,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,1,0));
			break;
		case 17:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,2,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,0,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,1,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,2,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,0,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,1,1));
			break;
		case 18:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,0,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,1,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,2,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,0,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,1,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,2,2));
			break;
		case 19:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,0,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,1,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,2,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,0,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,1,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,2,0));
			break;
		case 20:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,0,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,1,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,2,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,0,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,1,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,2,1));
			break;
		case 21:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,1,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,2,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,0,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,1,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,2,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,0,2));
			break;
		case 22:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,1,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,2,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,0,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,1,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,2,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,0,0));
			break;
		case 23:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,1,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,2,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,0,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,1,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,2,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,0,1));
			break;
		case 24:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,2,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,0,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,1,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,2,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,0,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,1,2));
			break;
		case 25:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,2,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,0,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,1,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,2,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,0,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,1,0));
			break;
		case 26:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,2,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,0,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,1,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,2,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,0,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,1,1));
			break;
		case 27:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,0,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,1,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,2,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,0,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,1,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,2,2));
			break;
		case 28:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,0,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,1,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,2,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,0,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,1,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,2,0));
			break;
		case 29:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,0,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,1,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,2,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,0,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,1,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,2,1));
			break;
		case 30:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,1,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,2,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,0,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,1,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,2,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,0,2));
			break;
		case 31:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,1,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,2,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,0,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,1,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,2,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,0,0));
			break;
		case 32:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,1,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,2,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,0,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,1,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,2,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,0,1));
			break;
		case 33:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,2,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,0,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,1,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,2,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,0,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,1,2));
			break;
		case 34:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,2,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,0,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,1,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,2,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,0,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,1,0));
			break;
		case 35:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,2,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,0,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,1,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,2,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,0,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,1,1));
			break;
		case 36:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,0,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,1,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,2,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,0,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,1,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,2,2));
			break;
		case 37:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,0,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,1,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,2,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,0,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,1,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,2,0));
			break;
		case 38:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,0,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,1,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,2,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,0,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,1,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,2,1));
			break;
		case 39:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,1,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,2,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,0,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,1,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,2,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,0,2));
			break;
		case 40:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,1,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,2,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,0,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,1,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,2,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,0,0));
			break;
		case 41:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,1,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,2,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,0,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,1,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,2,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,0,1));
			break;
		case 42:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,2,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,0,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,1,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,2,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,0,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,1,2));
			break;
		case 43:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,2,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,0,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,1,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,2,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,0,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,1,0));
			break;
		case 44:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,2,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,0,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,1,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,2,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,0,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,1,1));
			break;
		case 45:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,0,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,1,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,2,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,0,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,1,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,2,2));
			break;
		case 46:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,0,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,1,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,2,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,0,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,1,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,2,0));
			break;
		case 47:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,0,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,1,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,2,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,0,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,1,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,2,1));
			break;
		case 48:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,1,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,2,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,0,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,1,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,2,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,0,2));
			break;
		case 49:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,1,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,2,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,0,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,1,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,2,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,0,0));
			break;
		case 50:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,1,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,2,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,0,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,1,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,2,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,0,1));
			break;
		case 51:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,2,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,0,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,1,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,2,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,0,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,1,2));
			break;
		case 52:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,2,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,0,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,1,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,2,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,0,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,1,0));
			break;
		case 53:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,2,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,0,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,1,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,2,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,0,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,1,1));
			break;
		case 54:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,0,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,1,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,2,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,0,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,1,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,2,2));
			break;
		case 55:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,0,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,1,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,2,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,0,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,1,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,2,0));
			break;
		case 56:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,0,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,1,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,2,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,0,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,1,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,2,1));
			break;
		case 57:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,1,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,2,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,0,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,1,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,2,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,0,2));
			break;
		case 58:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,1,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,2,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,0,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,1,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,2,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,0,0));
			break;
		case 59:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,1,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,2,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,0,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,1,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,2,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,0,1));
			break;
		case 60:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,2,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,0,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,1,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,2,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,0,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,1,2));
			break;
		case 61:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,2,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,0,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,1,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,2,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,0,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,1,0));
			break;
		case 62:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,0,2,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,1,0,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,2,1,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,0,2,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,1,0,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,2,1,1));
			break;
		case 63:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,0,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,1,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,2,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,0,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,1,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,2,2));
			break;
		case 64:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,0,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,1,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,2,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,0,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,1,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,2,0));
			break;
		case 65:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,0,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,1,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,2,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,0,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,1,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,2,1));
			break;
		case 66:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,1,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,2,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,0,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,1,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,2,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,0,2));
			break;
		case 67:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,1,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,2,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,0,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,1,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,2,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,0,0));
			break;
		case 68:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,1,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,2,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,0,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,1,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,2,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,0,1));
			break;
		case 69:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,2,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,0,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,1,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,2,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,0,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,1,2));
			break;
		case 70:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,2,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,0,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,1,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,2,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,0,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,1,0));
			break;
		case 71:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,1,2,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,2,0,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,0,1,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,1,2,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,2,0,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,0,1,1));
			break;
		case 72:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,0,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,1,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,2,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,0,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,1,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,2,2));
			break;
		case 73:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,0,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,1,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,2,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,0,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,1,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,2,0));
			break;
		case 74:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,0,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,1,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,2,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,0,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,1,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,2,1));
			break;
		case 75:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,1,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,2,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,0,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,1,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,2,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,0,2));
			break;
		case 76:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,1,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,2,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,0,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,1,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,2,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,0,0));
			break;
		case 77:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,1,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,2,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,0,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,1,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,2,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,0,1));
			break;
		case 78:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,2,0));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,0,1));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,1,2));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,2,0));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,0,1));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,1,2));
			break;
		case 79:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,2,1));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,0,2));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,1,0));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,2,1));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,0,2));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,1,0));
			break;
		case 80:
			o_k = _mm_shuffle_ps(o, o, _MM_SHUFFLE(2,2,2,2));
			o_u = _mm_shuffle_ps(o, o, _MM_SHUFFLE(0,0,0,0));
			o_v = _mm_shuffle_ps(o, o, _MM_SHUFFLE(1,1,1,1));
			d_k = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2,2,2,2));
			d_u = _mm_shuffle_ps(d, d, _MM_SHUFFLE(0,0,0,0));
			d_v = _mm_shuffle_ps(d, d, _MM_SHUFFLE(1,1,1,1));
			break;
		default:
			return false;
	}

	const __m128 
		n_d = _mm_load_ps(nd),
		n_u = _mm_load_ps(nu),
		n_v = _mm_load_ps(nv);

	const __m128 
		ounu = _mm_mul_ps(o_u, n_u),
		ovnv = _mm_mul_ps(o_v, n_v),
		dunu = _mm_mul_ps(d_u, n_u),
		dvnv = _mm_mul_ps(d_v, n_v);

	const __m128
		num   = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(n_d, ounu), ovnv), o_k),
		denom = _mm_add_ps(_mm_add_ps(dunu, dvnv), d_k);

	const __m128 
		t = _mm_div_ps(num, denom),
		mint = _mm_load1_ps(&_mint),
		maxt = _mm_load1_ps(&_maxt);

	__m128 hasIts = 
		_mm_and_ps(_mm_cmpgt_ps(maxt, t), _mm_cmpgt_ps(t, mint));

	if (_mm_movemask_ps(hasIts) == 0) {
#ifdef MTS_DEBUG_FP
		enable_fpexcept();
#endif
		return false;
	}

	const __m128 
		a_u = _mm_load_ps(au),
		a_v = _mm_load_ps(av);

	const __m128 
		hu = _mm_add_ps(o_u, _mm_sub_ps(_mm_mul_ps(t, d_u), a_u)),
		hv = _mm_add_ps(o_v, _mm_sub_ps(_mm_mul_ps(t, d_v), a_v));
	
	const __m128 
		b_nu = _mm_load_ps(bnu),
		b_nv = _mm_load_ps(bnv),
		c_nu = _mm_load_ps(cnu),
		c_nv = _mm_load_ps(cnv);

	const __m128
		u = _mm_add_ps(_mm_mul_ps(hv, b_nu), _mm_mul_ps(hu, b_nv)),
		v = _mm_add_ps(_mm_mul_ps(hu, c_nu), _mm_mul_ps(hv, c_nv));

	const __m128 
		zero = _mm_setzero_ps(),
		one  = SSEConstants::one.ps;

	hasIts = _mm_and_ps(hasIts, 
		_mm_and_ps(_mm_and_ps(_mm_cmpge_ps(u, zero), _mm_cmpge_ps(v, zero)),
			_mm_cmpge_ps(one, _mm_add_ps(u, v))));

	if (_mm_movemask_ps(hasIts) != 0) {
		_mm_store_ps(_t, _mm_and_ps(t, hasIts));
		_mm_store_ps(_u, u);
		_mm_store_ps(_v, v);
	
		int closest = 0;
		float closestValue = std::numeric_limits<float>::max();
		for (int i=0; i<4; i++) {
			if (_t[i] != 0 && _t[i] <= closestValue) {
				closest = i;
				closestValue = _t[i];
			}
		}
		_t[0] = closestValue;
		_u[0] = _u[closest];
		_v[0] = _v[closest];
		_shapeIndex = shapeIndex[closest];
		_index = index[closest];
#ifdef MTS_DEBUG_FP
		enable_fpexcept();
#endif
		return true;
	}
#ifdef MTS_DEBUG_FP
	enable_fpexcept();
#endif
	return false;
}

MTS_NAMESPACE_END

#endif /* __TRIACCEL_SSE_H */


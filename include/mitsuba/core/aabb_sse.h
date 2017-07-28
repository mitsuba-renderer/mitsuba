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
#if !defined(__MITSUBA_CORE_AABB_SSE_H_)
#define __MITSUBA_CORE_AABB_SSE_H_

#include <mitsuba/core/platform.h>
#include <mitsuba/core/sse.h>
#include <mitsuba/core/aabb.h>
#include <mitsuba/core/ray_sse.h>

MTS_NAMESPACE_BEGIN

/**
 * NaN-aware slab test using SSE by Thierry Berger-Perrin (Intersects
 * against 4 rays simultaneously). Returns false if none of the rays
 * intersect.
 */
FINLINE bool AABB::rayIntersectPacket(const RayPacket4 &ray,
                                   RayInterval4 &interval) const {
    const __m128
        xl1 = _mm_mul_ps(ray.dRcp[0].ps,
            _mm_sub_ps(_mm_set1_ps(min.x), ray.o[0].ps)),
        xl2 = _mm_mul_ps(ray.dRcp[0].ps,
            _mm_sub_ps(_mm_set1_ps(max.x), ray.o[0].ps)),
        xl1a = _mm_min_ps(xl1, SSEConstants::p_inf.ps),
        xl2a = _mm_min_ps(xl2, SSEConstants::p_inf.ps),
        xl1b = _mm_max_ps(xl1, SSEConstants::n_inf.ps),
        xl2b = _mm_max_ps(xl2, SSEConstants::n_inf.ps);

    __m128
        lmax = _mm_max_ps(xl1a, xl2a),
        lmin = _mm_min_ps(xl1b, xl2b);

    const __m128
        yl1 = _mm_mul_ps(ray.dRcp[1].ps,
            _mm_sub_ps(_mm_set1_ps(min.y), ray.o[1].ps)),
        yl2 = _mm_mul_ps(ray.dRcp[1].ps,
            _mm_sub_ps(_mm_set1_ps(max.y), ray.o[1].ps)),
        yl1a = _mm_min_ps(yl1, SSEConstants::p_inf.ps),
        yl2a = _mm_min_ps(yl2, SSEConstants::p_inf.ps),
        yl1b = _mm_max_ps(yl1, SSEConstants::n_inf.ps),
        yl2b = _mm_max_ps(yl2, SSEConstants::n_inf.ps);

    lmax = _mm_min_ps(_mm_max_ps(yl1a,yl2a), lmax);
    lmin = _mm_max_ps(_mm_min_ps(yl1b,yl2b), lmin);

    const __m128
        zl1 = _mm_mul_ps(ray.dRcp[2].ps,
            _mm_sub_ps(_mm_set1_ps(min.z), ray.o[2].ps)),
        zl2 = _mm_mul_ps(ray.dRcp[2].ps,
            _mm_sub_ps(_mm_set1_ps(max.z), ray.o[2].ps)),
        zl1a = _mm_min_ps(zl1, SSEConstants::p_inf.ps),
        zl2a = _mm_min_ps(zl2, SSEConstants::p_inf.ps),
        zl1b = _mm_max_ps(zl1, SSEConstants::n_inf.ps),
        zl2b = _mm_max_ps(zl2, SSEConstants::n_inf.ps);

    lmax = _mm_min_ps(_mm_max_ps(zl1a,zl2a), lmax);
    lmin = _mm_max_ps(_mm_min_ps(zl1b,zl2b), lmin);

    const bool hasIntersection = _mm_movemask_ps(
        _mm_and_ps(
            _mm_cmpge_ps(lmax, _mm_setzero_ps()),
            _mm_cmple_ps(lmin, lmax))) != 0;

    interval.mint.ps = lmin;
    interval.maxt.ps = lmax;
    return hasIntersection;
}

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_AABB_SSE_H_ */

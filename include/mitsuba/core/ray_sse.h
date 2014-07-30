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
#if !defined(__MITSUBA_CORE_RAY_SSE_H_)
#define __MITSUBA_CORE_RAY_SSE_H_

#if !MTS_SSE
#error "This headers requires SSE support."
#endif

#include <mitsuba/core/platform.h>
#include <mitsuba/core/sse.h>
#include <mitsuba/core/ray.h>

MTS_NAMESPACE_BEGIN

/** \brief SIMD quad-packed ray for coherent ray tracing */
struct RayPacket4 {
	QuadVector o, d;
	QuadVector dRcp;
	uint8_t signs[4][4];

	inline RayPacket4() {
	}

	inline bool load(const Ray *rays) {
		for (int i=0; i<4; i++) {
			for (int axis=0; axis<3; axis++) {
				o[axis].f[i] = rays[i].o[axis];
				d[axis].f[i] = rays[i].d[axis];
				dRcp[axis].f[i] = rays[i].dRcp[axis];
				signs[axis][i] = rays[i].d[axis] < 0 ? 1 : 0;
				if (signs[axis][i] != signs[axis][0])
					return false;
			}
		}
		return true;
	}
};

struct RayInterval4 {
	SSEVector mint;
	SSEVector maxt;

	inline RayInterval4() {
		mint = SSEConstants::eps;
		maxt = SSEConstants::p_inf;
	}

	inline RayInterval4(const Ray *rays) {
		for (int i=0; i<4; i++) {
			mint.f[i] = rays[i].mint;
			maxt.f[i] = rays[i].maxt;
		}
	}
};

struct Intersection4 {
	SSEVector t;
	SSEVector u;
	SSEVector v;
	SSEVector primIndex;
	SSEVector shapeIndex;

	inline Intersection4() {
		t          = SSEConstants::p_inf;
		u          = SSEConstants::zero;
		v          = SSEConstants::zero;
		primIndex  = SSEConstants::ffffffff;
		shapeIndex = SSEConstants::ffffffff;
	}
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_RAY_SSE_H_ */

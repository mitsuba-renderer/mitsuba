/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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
#if !defined(__MITSUBA_CORE_CONSTANTS_H)
#define __MITSUBA_CORE_CONSTANTS_H

/* Choice of precision */
#ifdef DOUBLE_PRECISION
#define Epsilon 1e-7
#define ShadowEpsilon 1e-5
#else
#define Epsilon 1e-4f
#define ShadowEpsilon 1e-3f
#endif
#define DeltaEpsilon 1e-3f

/* Assumed L1 cache line size for alignment purposes */
#if !defined(L1_CACHE_LINE_SIZE)
#define L1_CACHE_LINE_SIZE 64
#endif

#ifdef M_E
#undef M_E
#endif

#ifdef M_PI
#undef M_PI
#endif

#ifdef INFINITY
#undef INFINITY
#endif

#if defined(__WINDOWS__)
#define ONE_MINUS_EPS_FLT 0.999999940395355225f
#define ONE_MINUS_EPS_DBL 0.999999999999999888
#else
#define ONE_MINUS_EPS_FLT 0x1.fffffep-1f
#define ONE_MINUS_EPS_DBL 0x1.fffffffffffff7p-1
#endif

#ifdef SINGLE_PRECISION
#define M_E           2.71828182845904523536f
#define M_PI          3.14159265358979323846f
#define INV_PI        0.31830988618379067154f
#define INV_TWOPI     0.15915494309189533577f
#define INV_FOURPI    0.07957747154594766788f
#define SQRT_TWO      1.41421356237309504880f
#define INV_SQRT_TWO  0.70710678118654752440f
#define ONE_MINUS_EPS ONE_MINUS_EPS_FLT
#else
#define M_E           2.71828182845904523536
#define M_PI          3.14159265358979323846
#define INV_PI        0.31830988618379067154
#define INV_TWOPI     0.15915494309189533577
#define INV_FOURPI    0.07957747154594766788
#define SQRT_TWO      1.41421356237309504880
#define INV_SQRT_TWO  0.70710678118654752440
#define ONE_MINUS_EPS ONE_MINUS_EPS_DBL
#endif
#endif /* __MITSUBA_CORE_CONSTANTS_H */

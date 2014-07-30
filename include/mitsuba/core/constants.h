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
#define RCPOVERFLOW_FLT   2.93873587705571876e-39f
#define RCPOVERFLOW_DBL   5.56268464626800345e-309
#else
#define ONE_MINUS_EPS_FLT 0x1.fffffep-1f
#define ONE_MINUS_EPS_DBL 0x1.fffffffffffff7p-1
#define RCPOVERFLOW_FLT   0x1p-128f
#define RCPOVERFLOW_DBL   0x1p-1024
#endif

#define M_E_FLT           2.71828182845904523536f
#define M_PI_FLT          3.14159265358979323846f
#define INV_PI_FLT        0.31830988618379067154f
#define INV_TWOPI_FLT     0.15915494309189533577f
#define INV_FOURPI_FLT    0.07957747154594766788f
#define SQRT_TWO_FLT      1.41421356237309504880f
#define INV_SQRT_TWO_FLT  0.70710678118654752440f

#define M_E_DBL           2.71828182845904523536
#define M_PI_DBL          3.14159265358979323846
#define INV_PI_DBL        0.31830988618379067154
#define INV_TWOPI_DBL     0.15915494309189533577
#define INV_FOURPI_DBL    0.07957747154594766788
#define SQRT_TWO_DBL      1.41421356237309504880
#define INV_SQRT_TWO_DBL  0.70710678118654752440

#ifdef SINGLE_PRECISION
#define M_E               M_E_FLT
#define M_PI              M_PI_FLT
#define INV_PI            INV_PI_FLT
#define INV_TWOPI         INV_TWOPI_FLT
#define INV_FOURPI        INV_FOURPI_FLT
#define SQRT_TWO          SQRT_TWO_FLT
#define INV_SQRT_TWO      INV_SQRT_TWO_FLT
#define ONE_MINUS_EPS     ONE_MINUS_EPS_FLT
#define RCPOVERFLOW       RCPOVERFLOW_FLT
#else
#define M_E               M_E_DBL
#define M_PI              M_PI_DBL
#define INV_PI            INV_PI_DBL
#define INV_TWOPI         INV_TWOPI_DBL
#define INV_FOURPI        INV_FOURPI_DBL
#define SQRT_TWO          SQRT_TWO_DBL
#define INV_SQRT_TWO      INV_SQRT_TWO_DBL
#define ONE_MINUS_EPS     ONE_MINUS_EPS_DBL
#define RCPOVERFLOW       RCPOVERFLOW_DBL
#endif
#endif /* __MITSUBA_CORE_CONSTANTS_H */

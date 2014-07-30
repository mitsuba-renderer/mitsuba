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
#if !defined(__MITSUBA_RENDER_NOISE_H_)
#define __MITSUBA_RENDER_NOISE_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Contains a few useful noise functions
 *
 * The implementations in this class are based on PBRT
 * \ingroup librender
 * \ingroup libpython
 */
class MTS_EXPORT_RENDER Noise {
public:
	/**
	 * \brief Evaluate the Perlin noise function at \a p.
	 */
	static Float perlinNoise(const Point &p);

	/**
	 * \brief Evaluate a fractional Brownian noise function
	 * based on \ref perlinNoise() at \a p.
	 *
	 * \param dpdx Differential of p with respect to
	 *   the next horizontal pixel in screen-space.
	 * \param dpdy Differential of p with respect to
	 *   the next vertical pixel in screen-space.
	 * \param omega Controls the falloff weights applied
	 *   to higher-frequency octaves
	 * \param maxOctaves Max. number of octaves used
	 *   in the noise computation
	 */
	static Float fbm(const Point &p, const Vector &dpdx,
		const Vector &dpdy, Float omega, int maxOctaves);

	/**
	 * \brief Similar to \ref fbm, but adds first-derivative
	 * discontinuities, causing the resulting function to have
	 * infinite frequency content.
	 *
	 * \param dpdx Differential of p with respect to
	 *   the next horizontal pixel in screen-space.
	 * \param dpdy Differential of p with respect to
	 *   the next vertical pixel in screen-space.
	 * \param omega Controls the falloff weights applied
	 *   to higher-frequency octaves
	 * \param maxOctaves Max. number of octaves used
	 *   in the noise computation
	 */
	static Float turbulence(const Point &p, const Vector &dpdx,
		const Vector &dpdy, Float omega, int maxOctaves);
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_NOISE_H_ */

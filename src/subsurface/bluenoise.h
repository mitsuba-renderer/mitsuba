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

#if !defined(__BLUENOISE_H)
#define __BLUENOISE_H

#include "irrproc.h"

MTS_NAMESPACE_BEGIN

/**
 * \brief Generate a point set with blue noise properties
 *
 * Based on the paper "Parallel Poisson Disk Sampling with
 * Spectrum Analysis on Surfaces" by John Bowers, Rui Wang,
 * Li-Yi Wei and David Maletz.
 *
 * \param scene
 *    A pointer to the underlying scene
 * \param shapes
 *    A list of input shapes on which samples should be placed
 * \param radius
 *    The Poisson radius of the point set to be generated
 * \param target
 *    A position sample vector (which will be populated with the result)
 * \param sa
 *    Will be used to return the combined surface area of the shapes
 * \param aabb
 *    Will be used to store an an axis-aligned box containing all samples
 * \param data
 *    Custom pointer that will be sent along with progress messages
 *    (usually contains a pointer to the \ref RenderJob instance)
 */
extern void blueNoisePointSet(const Scene *scene,
    const std::vector<Shape *> &shapes, Float radius,
    PositionSampleVector *target, Float &sa, AABB &aabb,
    const void *data);

MTS_NAMESPACE_END

#endif /* __BLUENOISE_H */

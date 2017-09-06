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
#if !defined(__MITSUBA_BIDIR_COMMON_H_)
#define __MITSUBA_BIDIR_COMMON_H_

#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/// Activates extra debugging checks in the bidirectional layer
#define MTS_BD_DEBUG 1

#if MTS_BD_DEBUG == 1
#define BDAssert(expr) SAssert(expr)
#else
#define BDAssert(expr)
#endif

#define MAX(a, b) ((a) > (b) ? (a) : (b))

/**
 * \brief Data record associated with path endpoints (aka supernodes)
 * in the path-space framework
 *
 * This record stores the desired time value when starting a new path.
 * For sensor subpaths, it optionally specifies the desired pixel
 * position of the path.
 *
 * \ingroup libbidir
 */
struct EndpointRecord {
    /// Time value associated with the path
    Float time;

    /// Create a new endpoint record for a given time value
    inline EndpointRecord(Float time)
        : time(time) { }

    /// Create a new endpoint record for a given time value
    inline EndpointRecord(Float time,
        const Point2 &uv) : time(time) { }

    /// Return a human-readable description
    std::string toString() const;
};

/* Forward declarations */
struct PathVertex;
struct PathEdge;
struct Path;
struct PathSeed;
struct SplatList;
struct MutationRecord;
class MemoryPool;
class SeedWorkUnit;
class PathSampler;
class Mutator;
class PathSolver;
class SpecularManifold;

MTS_NAMESPACE_END

#endif /* __MITSUBA_BIDIR_COMMON_H_ */

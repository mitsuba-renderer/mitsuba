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
#if !defined(__MITSUBA_BIDIR_MEMPOOL_H_)
#define __MITSUBA_BIDIR_MEMPOOL_H_

#include <mitsuba/bidir/vertex.h>
#include <mitsuba/bidir/edge.h>
#include <mitsuba/core/mempool.h>

MTS_NAMESPACE_BEGIN

class MemoryPool {
public:
    /// Create a new memory pool with aninitial set of 128 entries
    MemoryPool(size_t nEntries = 128)
        : m_vertexPool(nEntries), m_edgePool(nEntries) { }

    /// Destruct the memory pool and release all entries
    ~MemoryPool() { }

    /// Acquire an edge
    inline PathEdge *allocEdge() {
        PathEdge *edge = m_edgePool.alloc();
        #if defined(MTS_BD_DEBUG_HEAVY)
        memset(edge, 0xFF, sizeof(PathEdge));
        #endif
        return edge;
    }

    /// Acquire an vertex
    inline PathVertex *allocVertex() {
        PathVertex *vertex = m_vertexPool.alloc();
        #if defined(MTS_BD_DEBUG_HEAVY)
        memset(vertex, 0xFF, sizeof(PathVertex));
        #endif
        return vertex;
    }

    /// Release an edge
    inline void release(PathEdge *edge) {
        m_edgePool.release(edge);
    }

    /// Release an entry
    inline void release(PathVertex *vertex) {
        m_vertexPool.release(vertex);
    }

    /// Check if every entry has been released
    bool unused() const {
        return m_vertexPool.unused() && m_edgePool.unused();
    }

    /// Return the currently allocated amount of storage for edges
    inline size_t edgeSize() {
        return m_edgePool.size();
    }

    /// Return the currently allocated amount of storage for vertices
    inline size_t vertexSize() {
        return m_vertexPool.size();
    }

    /// Return a human-readable description
    std::string toString() const {
        std::ostringstream oss;
        oss << "MemoryPool[" << endl
            << "  vertexPool = " << m_vertexPool.toString() << "," << endl
            << "  edgePool = " << m_edgePool.toString() << endl
            << "]";
        return oss.str();
    }

private:
    BasicMemoryPool<PathVertex> m_vertexPool;
    BasicMemoryPool<PathEdge> m_edgePool;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_BIDIR_MEMPOOL_H_ */

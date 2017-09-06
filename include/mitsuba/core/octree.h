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
#if !defined(__MITSUBA_CORE_OCTREE_H_)
#define __MITSUBA_CORE_OCTREE_H_

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/atomic.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Lock-free linked list data structure
 *
 * This class provides a very basic linked list data structure whose primary
 * purpose it is to efficiently service append operations from multiple parallel
 * threads. These are internally realized via atomic compare and exchange
 * operations, meaning that no lock must be acquired.
 *
 * \ingroup libcore
 */
template <typename T> class LockFreeList {
public:
    struct ListItem {
        T value;
        ListItem *next;

        inline ListItem(const T &value) :
            value(value), next(NULL) { }
    };

    inline LockFreeList() : m_head(NULL) {}

    ~LockFreeList() {
        ListItem *cur = m_head;
        while (cur) {
            ListItem *next = cur->next;
            delete cur;
            cur = next;
        }
    }

    inline const ListItem *head() const {
        return m_head;
    }

    void append(const T &value) {
        ListItem *item = new ListItem(value);
        ListItem **cur = &m_head;

        while (!atomicCompareAndExchangePtr<ListItem>(cur, item, NULL))
            cur = &((*cur)->next);
    }
private:
    ListItem *m_head;
};

/**
 * \brief Generic single-reference static octree
 *
 * This class is currently used to implement BSSRDF evaluation
 * with irradiance point clouds.
 *
 * The \c Item template parameter must implement a function
 * named <tt>getPosition()</tt> that returns a \ref Point.
 *
 * \ingroup libcore
 */
template <typename Item, typename NodeData> class StaticOctree {
public:
    struct OctreeNode {
        bool leaf : 1;
        NodeData data;

        union {
            struct {
                OctreeNode *children[8];
            };

            struct {
                uint32_t offset;
                uint32_t count;
            };
        };

        ~OctreeNode() {
            if (!leaf) {
                for (int i=0; i<8; ++i) {
                    if (children[i])
                        delete children[i];
                }
            }
        }
    };

    /**
     * \brief Create a new octree
     * \param maxDepth
     *     Maximum tree depth (24 by default)
     * \param maxItems
     *     Maximum items per interior node (8 by default)
     *
     * By default, the maximum tree depth is set to 16
     */
    inline StaticOctree(const AABB &aabb, uint32_t maxDepth = 24, uint32_t maxItems = 8) :
        m_aabb(aabb), m_maxDepth(maxDepth), m_maxItems(maxItems), m_root(NULL) { }

    /// Release all memory
    ~StaticOctree() {
        if (m_root)
            delete m_root;
    }

    void build() {
        SLog(EDebug, "Building an octree over " SIZE_T_FMT " data points (%s)..",
            m_items.size(), memString(m_items.size() * sizeof(Item)).c_str());

        ref<Timer> timer = new Timer();
        std::vector<uint32_t> perm(m_items.size()), temp(m_items.size());

        for (uint32_t i=0; i<m_items.size(); ++i)
            perm[i] = i;

        /* Build the kd-tree and compute a suitable permutation of the elements */
        m_root = build(m_aabb, 0, &perm[0], &temp[0], &perm[0], &perm[0] + m_items.size());

        /* Apply the permutation */
        permute_inplace(&m_items[0], perm);

        SLog(EDebug, "Done (took %i ms)" , timer->getMilliseconds());
    }

protected:
    struct LabelOrdering : public std::binary_function<uint32_t, uint32_t, bool> {
        LabelOrdering(const std::vector<Item> &items) : m_items(items) { }

        inline bool operator()(uint32_t a, uint32_t b) const {
            return m_items[a].label < m_items[b].label;
        }

        const std::vector<Item> &m_items;
    };

    /// Return the AABB for a child of the specified index
    inline AABB childBounds(int child, const AABB &nodeAABB, const Point &center) const {
        AABB childAABB;
        childAABB.min.x = (child & 4) ? center.x : nodeAABB.min.x;
        childAABB.max.x = (child & 4) ? nodeAABB.max.x : center.x;
        childAABB.min.y = (child & 2) ? center.y : nodeAABB.min.y;
        childAABB.max.y = (child & 2) ? nodeAABB.max.y : center.y;
        childAABB.min.z = (child & 1) ? center.z : nodeAABB.min.z;
        childAABB.max.z = (child & 1) ? nodeAABB.max.z : center.z;
        return childAABB;
    }

    OctreeNode *build(const AABB &aabb, uint32_t depth, uint32_t *base,
            uint32_t *temp, uint32_t *start, uint32_t *end) {
        if (start == end) {
            return NULL;
        } else if ((uint32_t) (end-start) < m_maxItems || depth > m_maxDepth) {
            OctreeNode *result = new OctreeNode();
            result->count = (uint32_t) (end-start);
            result->offset = (uint32_t) (start-base);
            result->leaf = true;
            return result;
        }

        Point center = aabb.getCenter();
        uint32_t nestedCounts[8];
        memset(nestedCounts, 0, sizeof(uint32_t)*8);

        /* Label all items */
        for (uint32_t *it = start; it != end; ++it) {
            Item &item = m_items[*it];
            const Point &p = item.getPosition();

            uint8_t label = 0;
            if (p.x > center.x) label |= 4;
            if (p.y > center.y) label |= 2;
            if (p.z > center.z) label |= 1;

            AABB bounds = childBounds(label, aabb, center);
            SAssert(bounds.contains(p));

            item.label = label;
            nestedCounts[label]++;
        }

        uint32_t nestedOffsets[9];
        nestedOffsets[0] = 0;
        for (int i=1; i<=8; ++i)
            nestedOffsets[i] = nestedOffsets[i-1] + nestedCounts[i-1];

        /* Sort by label */
        for (uint32_t *it = start; it != end; ++it) {
            int offset = nestedOffsets[m_items[*it].label]++;
            temp[offset] = *it;
        }
        memcpy(start, temp, (end-start) * sizeof(uint32_t));

        /* Recurse */
        OctreeNode *result = new OctreeNode();
        for (int i=0; i<8; i++) {
            AABB bounds = childBounds(i, aabb, center);

            uint32_t *it = start + nestedCounts[i];
            result->children[i] = build(bounds, depth+1, base, temp, start, it);
            start = it;
        }

        result->leaf = false;

        return result;
    }

    inline StaticOctree() : m_root(NULL) { }
protected:
    AABB m_aabb;
    std::vector<Item> m_items;
    uint32_t m_maxDepth;
    uint32_t m_maxItems;
    OctreeNode *m_root;
};

/**
 * \brief Generic multiple-reference octree with support for parallel dynamic updates
 *
 * Based on the excellent implementation in PBRT. Modifications are
 * the addition of a bounding sphere query and support for multithreading.
 *
 * This class is currently used to implement irradiance caching.
 *
 * \ingroup libcore
 */
template <typename Item> class DynamicOctree {
public:
    /**
     * \brief Create a new octree
     *
     * By default, the maximum tree depth is set to 24
     */
    inline DynamicOctree(const AABB &aabb, uint32_t maxDepth = 24)
     : m_aabb(aabb), m_maxDepth(maxDepth) {
    }

    /// Insert an item with the specified cell coverage
    inline void insert(const Item &value, const AABB &coverage) {
        insert(&m_root, m_aabb, value, coverage,
            coverage.getExtents().lengthSquared(), 0);
    }

    /// Execute <tt>functor.operator()</tt> on all records, which potentially overlap \c p
    template <typename Functor> inline void lookup(const Point &p, Functor &functor) const {
        if (!m_aabb.contains(p))
            return;
        lookup(&m_root, m_aabb, p, functor);
    }

    /// Execute <tt>functor.operator()</tt> on all records, which potentially overlap \c bsphere
    template <typename Functor> inline void searchSphere(const BSphere &sphere, Functor &functor) {
        if (!m_aabb.overlaps(sphere))
            return;
        searchSphere(&m_root, m_aabb, sphere, functor);
    }

    inline const AABB &getAABB() const { return m_aabb; }
private:
    struct OctreeNode {
    public:
        OctreeNode() {
            for (int i=0; i<8; ++i)
                children[i] = NULL;
        }

        ~OctreeNode() {
            for (int i=0; i<8; ++i) {
                if (children[i])
                    delete children[i];
            }
        }

        OctreeNode *children[8];
        LockFreeList<Item> data;
    };

    /// Return the AABB for a child of the specified index
    inline AABB childBounds(int child, const AABB &nodeAABB, const Point &center) const {
        AABB childAABB;
        childAABB.min.x = (child & 4) ? center.x : nodeAABB.min.x;
        childAABB.max.x = (child & 4) ? nodeAABB.max.x : center.x;
        childAABB.min.y = (child & 2) ? center.y : nodeAABB.min.y;
        childAABB.max.y = (child & 2) ? nodeAABB.max.y : center.y;
        childAABB.min.z = (child & 1) ? center.z : nodeAABB.min.z;
        childAABB.max.z = (child & 1) ? nodeAABB.max.z : center.z;
        return childAABB;
    }

    void insert(OctreeNode *node, const AABB &nodeAABB, const Item &value,
            const AABB &coverage, Float diag2, uint32_t depth) {
        /* Add the data item to the current octree node if the max. tree
           depth is reached or the data item's coverage area is smaller
           than the current node size */
        if (depth == m_maxDepth ||
            (nodeAABB.getExtents().lengthSquared() < diag2)) {
            node->data.append(value);
            return;
        }

        /* Otherwise: test for overlap */
        const Point center = nodeAABB.getCenter();

        /* Otherwise: test for overlap */
        bool x[2] = { coverage.min.x <= center.x, coverage.max.x > center.x };
        bool y[2] = { coverage.min.y <= center.y, coverage.max.y > center.y };
        bool z[2] = { coverage.min.z <= center.z, coverage.max.z > center.z };
        bool over[8] = { x[0] && y[0] && z[0], x[0] && y[0] && z[1],
                         x[0] && y[1] && z[0], x[0] && y[1] && z[1],
                         x[1] && y[0] && z[0], x[1] && y[0] && z[1],
                         x[1] && y[1] && z[0], x[1] && y[1] && z[1] };

        /* Recurse */
        for (int child=0; child<8; ++child) {
            if (!over[child])
                continue;
            if (!node->children[child]) {
                OctreeNode *newNode = new OctreeNode();
                if (!atomicCompareAndExchangePtr<OctreeNode>(&node->children[child], newNode, NULL))
                    delete newNode;
            }
            const AABB childAABB(childBounds(child, nodeAABB, center));
            insert(node->children[child], childAABB,
                value, coverage, diag2, depth+1);
        }
    }

    /// Internal lookup procedure - const version
    template <typename Functor> inline void lookup(const OctreeNode *node,
            const AABB &nodeAABB, const Point &p, Functor &functor) const {
        const Point center = nodeAABB.getCenter();

        const typename LockFreeList<Item>::ListItem *item = node->data.head();
        while (item) {
            functor(item->value);
            item = item->next;
        }

        int child = (p.x > center.x ? 4 : 0)
                + (p.y > center.y ? 2 : 0)
                + (p.z > center.z ? 1 : 0);

        OctreeNode *childNode = node->children[child];

        if (childNode) {
            const AABB childAABB(childBounds(child, nodeAABB, center));
            lookup(node->children[child], childAABB, p, functor);
        }
    }

    template <typename Functor> inline void searchSphere(OctreeNode *node,
            const AABB &nodeAABB, const BSphere &sphere,
            Functor &functor) {
        const Point center = nodeAABB.getCenter();

        const typename LockFreeList<Item>::ListItem *item = node->data.head();
        while (item) {
            functor(item->value);
            item = item->next;
        }

        // Potential for much optimization..
        for (int child=0; child<8; ++child) {
            if (node->children[child]) {
                const AABB childAABB(childBounds(child, nodeAABB, center));
                if (childAABB.overlaps(sphere))
                    searchSphere(node->children[child], childAABB, sphere, functor);
            }
        }
    }
private:
    OctreeNode m_root;
    AABB m_aabb;
    uint32_t m_maxDepth;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_OCTREE_H_ */

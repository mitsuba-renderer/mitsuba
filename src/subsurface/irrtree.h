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

#if !defined(__IRRTREE_H)
#define __IRRTREE_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/octree.h>
#include "irrproc.h"

MTS_NAMESPACE_BEGIN

class IrradianceOctree : public StaticOctree<IrradianceSample, IrradianceSample>, public SerializableObject {
public:
    /// Construct a new irradiance octree
    IrradianceOctree(const AABB &aabb, Float solidAngleThreshold,
        std::vector<IrradianceSample> &records);

    /// Unserialize an octree from a binary data stream
    IrradianceOctree(Stream *stream, InstanceManager *manager);

    /// Serialize an octree to a binary data stream
    void serialize(Stream *stream, InstanceManager *manager) const;

    /// Query the octree using a customizable functor, while representatives for distant nodes
    template <typename QueryType> inline void performQuery(QueryType &query) const {
        performQuery(m_aabb, m_root, query);
    }

    MTS_DECLARE_CLASS()
protected:
    /// Propagate irradiance approximations througout the tree
    void propagate(OctreeNode *node);

    /// Query the octree using a customizable functor, while representatives for distant nodes
    template <typename QueryType> void performQuery(const AABB &aabb, OctreeNode *node, QueryType &query) const {
        /* Compute the approximate solid angle subtended by samples within this node */
        Float approxSolidAngle = node->data.area / (query.p - node->data.p).lengthSquared();

        /* Use the representative if this is a distant node */
        if (!aabb.contains(query.p) && approxSolidAngle < m_solidAngleThreshold) {
            query(node->data);
        } else {
            if (node->leaf) {
                for (uint32_t i=0; i<node->count; ++i)
                    query(m_items[node->offset + i]);
            } else {
                Point center = aabb.getCenter();
                for (int i=0; i<8; i++) {
                    if (!node->children[i])
                        continue;

                    AABB childAABB = childBounds(i, aabb, center);
                    performQuery(childAABB, node->children[i], query);
                }
            }
        }
    }
private:
    Float m_solidAngleThreshold;
};

MTS_NAMESPACE_END

#endif /* __IRRTREE_H */

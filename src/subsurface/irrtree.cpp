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

#include <mitsuba/core/statistics.h>
#include "irrtree.h"

MTS_NAMESPACE_BEGIN

static StatsCounter statsNumSamples("SSS Irradiance Octree", "Created samples");
static StatsCounter statsNumNodes("SSS Irradiance Octree", "Created nodes");

IrradianceOctree::IrradianceOctree(const AABB &bounds, Float solidAngleThreshold, std::vector<IrradianceSample> &records)
    : StaticOctree<IrradianceSample, IrradianceSample>(bounds), m_solidAngleThreshold(solidAngleThreshold) {

    m_items.swap(records);

    build();
    propagate(m_root);
}

IrradianceOctree::IrradianceOctree(Stream *stream, InstanceManager *manager) {
    m_aabb = AABB(stream);
    m_maxDepth = stream->readUInt();
    m_maxItems = stream->readUInt();
    m_solidAngleThreshold = stream->readFloat();

    size_t items = stream->readSize();
    m_items.resize(items);

    for (size_t i=0; i<items; ++i)
        m_items[i] = IrradianceSample(stream);

    build();
    propagate(m_root);
}

void IrradianceOctree::serialize(Stream *stream, InstanceManager *manager) const {
    m_aabb.serialize(stream);
    stream->writeUInt(m_maxDepth);
    stream->writeUInt(m_maxItems);
    stream->writeFloat(m_solidAngleThreshold);

    stream->writeSize(m_items.size());
    for (size_t i=0; i<m_items.size(); ++i)
        m_items[i].serialize(stream);
}

void IrradianceOctree::propagate(OctreeNode *node) {
    IrradianceSample &repr = node->data;

    /* Initialize the cluster values */
    repr.E = Spectrum(0.0f);
    repr.area = 0.0f;
    repr.p = Point(0.0f, 0.0f, 0.0f);
    Float weightSum = 0.0f;

    if (node->leaf) {
        /* Inner node */
        for (uint32_t i=0; i<node->count; ++i) {
            const IrradianceSample &sample = m_items[i+node->offset];
            repr.E += sample.E * sample.area;
            repr.area += sample.area;
            Float weight = sample.E.getLuminance() * sample.area;
            repr.p += sample.p * weight;
            weightSum += weight;
        }
        statsNumSamples += node->count;
    } else {
        /* Inner node */
        for (int i=0; i<8; i++) {
            OctreeNode *child = node->children[i];
            if (!child)
                continue;
            propagate(child);
            repr.E += child->data.E * child->data.area;
            repr.area += child->data.area;
            Float weight = child->data.E.getLuminance() * child->data.area;
            repr.p += child->data.p * weight;
            weightSum += weight;
        }
    }
    if (repr.area != 0)
        repr.E /= repr.area;
    if (weightSum != 0)
        repr.p /= weightSum;

    ++statsNumNodes;
}

MTS_IMPLEMENT_CLASS_S(IrradianceOctree, false, SerializableObject)
MTS_NAMESPACE_END

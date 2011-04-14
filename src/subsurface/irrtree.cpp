/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#include "irrtree.h"
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

IrradianceOctree::IrradianceOctree(int maxDepth, Float threshold, const AABB &bounds) 
 : m_maxDepth(maxDepth), m_threshold(threshold) {
	m_root = new OctreeNode(bounds);
	m_numSamples = 0;
}

IrradianceOctree::IrradianceOctree(Stream *stream, InstanceManager *manager) : 
		SerializableObject(stream, manager) {
	m_root = new OctreeNode(AABB(stream));
	m_threshold = stream->readFloat();
	m_maxDepth = stream->readInt();
	m_numSamples = 0;
	unsigned int numSamples = stream->readUInt();

	for (unsigned int i=0; i<numSamples; ++i)
		addSample(IrradianceSample(stream));
	preprocess();
}

IrradianceOctree::~IrradianceOctree() {
	delete m_root;
}

void IrradianceOctree::serialize(Stream *stream, InstanceManager *manager) const {
	m_root->aabb.serialize(stream);
	stream->writeFloat(m_threshold);
	stream->writeInt(m_maxDepth);
	stream->writeUInt(m_numSamples);
	m_root->serialize(stream);
}
	
	
void IrradianceOctree::addSample(const IrradianceSample &sample) {
	static StatsCounter numSamples("SSS IrradianceOctree", "Number of samples");
	++numSamples;
	++m_numSamples;
	if (!sample.E.isValid())
		Log(EWarn, "Invalid sample: %s", sample.E.toString().c_str());
	else
		addSample(m_root, sample, 0);
}

void IrradianceOctree::dumpOBJ(const std::string &filename) const {
	std::ofstream os(filename.c_str());
	os << "o IrrSamples" << endl;
	m_root->dumpOBJ(os);
	/// Need to generate some fake geometry so that blender will import the points
	for (unsigned int i=3; i<=m_numSamples; i++) 
		os << "f " << i << " " << i-1 << " " << i-2 << endl;
	os.close();
}

void IrradianceOctree::addSample(OctreeNode *node, const IrradianceSample &sample, int depth) {
	static StatsCounter nodesCreated("SSS IrradianceOctree", "Number of created nodes");

	node->samples.push_back(sample);
	if (node->leaf && (m_maxDepth == depth || node->samples.size() < 8)) 
		return;

	node->leaf = false;

	for (sample_iterator it = node->samples.begin();
		it != node->samples.end(); ++it) {
		const IrradianceSample &s = *it;
		Point center = node->aabb.getCenter();
		int nodeId = (s.p.x > center.x ? 1 : 0) +
			(s.p.y > center.y ? 2 : 0) +
			(s.p.z > center.z ? 4 : 0);

		if (!node->children[nodeId]) {
			AABB childAABB(center, center);
			childAABB.expandBy(node->aabb.getCorner(nodeId));
			node->children[nodeId] = new OctreeNode(childAABB);
			++nodesCreated;
		}

		addSample(node->children[nodeId], s, depth + 1);
	}
	std::vector<IrradianceSample> empty;
	node->samples.swap(empty);
}

void IrradianceOctree::preprocess(OctreeNode *node) {
	/* Initialize the cluster values */
	node->cluster.E = Spectrum(0.0f);
	node->cluster.area = 0.0f;
	node->cluster.p = Point(0.0f, 0.0f, 0.0f);
	Float combinedWeight = 0.0f;

	if (!node->samples.empty()) {
		/* Leaf node */
		for (sample_iterator it = node->samples.begin();
			it != node->samples.end(); ++it) {
			const IrradianceSample &sample = *it;
			node->cluster.E += sample.E * sample.area;
			node->cluster.area += sample.area;

			Float pointWeight = sample.E.average() * sample.area;
			node->cluster.p += sample.p * pointWeight;
			combinedWeight += pointWeight;
		}
		node->cluster.E /= node->cluster.area;
		if (combinedWeight != 0)
			node->cluster.p /= combinedWeight;
	} else {
		int numChildren = 0;
		/* Inner node */
		for (int i=0; i<8; i++) {
			OctreeNode *child = node->children[i];
			if (!child)
				continue;
			numChildren++;
			preprocess(child);
			/* Point repulsion not used, weight everything by area */
			node->cluster.E += child->cluster.E * child->cluster.area;
			node->cluster.area += child->cluster.area;

			Float pointWeight = child->cluster.E.average() * child->cluster.area;
			node->cluster.p += child->cluster.p * pointWeight;
			combinedWeight += pointWeight;
		}
		if (combinedWeight != 0)
			node->cluster.p /= combinedWeight;
		if (node->cluster.area != 0)
			node->cluster.E /= node->cluster.area;
	}
}

MTS_IMPLEMENT_CLASS_S(IrradianceOctree, false, SerializableObject)
MTS_NAMESPACE_END

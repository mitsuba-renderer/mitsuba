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

#if !defined(__IRRTREE_H)
#define __IRRTREE_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/aabb.h>
#include <fstream>
#include "irrproc.h"

MTS_NAMESPACE_BEGIN

class IrradianceOctree : public SerializableObject {
private:
	typedef std::vector<IrradianceSample>			sample_vector;
	typedef sample_vector::const_iterator	        sample_iterator;

	/// Private octree node class
	struct OctreeNode {
		OctreeNode *children[8];
		sample_vector samples;
		IrradianceSample cluster;
		AABB aabb;
		bool leaf;

		inline OctreeNode(const AABB &_aabb) : aabb(_aabb) {
			for (int i=0; i<8; i++)
				children[i] = NULL;
			leaf = true;
		}

		void serialize(Stream *stream) {
			for (size_t i=0; i<samples.size(); ++i)
				samples[i].serialize(stream);
			for (int i=0; i<8; i++)
				if (children[i])
					children[i]->serialize(stream);
		}

		void dumpOBJ(std::ostream &os) {
			for (size_t i=0; i<samples.size(); ++i) {
				const Point &p = samples[i].p;
				os << "v " << p.x << " " << p.y << " " << p.z << endl;
			}
			for (int i=0; i<8; i++)
				if (children[i])
					children[i]->dumpOBJ(os);
		}

		~OctreeNode() {
			for (int i=0; i<8; i++) {
				if (children[i])
					delete children[i];
			}
		}
	};
public:
	/// Construct an empty octree
	IrradianceOctree(int maxDepth, Float threshold, const AABB &bounds); 

	/// Unserialize an octree from a binary data stream
	IrradianceOctree(Stream *stream, InstanceManager *manager);

	/// Serialize an octree to a binary data stream
	void serialize(Stream *stream, InstanceManager *manager) const;

	/// Add an irradiance sample to the octree
	void addSample(const IrradianceSample &sample);

	/// Query the octree using a customizable implementation
	template <typename QueryType> inline void execute(QueryType &query) const {
		execute(m_root, query);
	}

	/// Write the samples to an OBJ file for debugging
	void dumpOBJ(const std::string &filename) const;

	/// Pre-process step (clusters the samples)
	inline void preprocess() { 
		Log(EDebug, "Sub-surface integrator - clustering %i samples..", m_numSamples);
		preprocess(m_root);
	}

	MTS_DECLARE_CLASS()
protected:
	/// Add an irradiance sample to the octree
	void addSample(OctreeNode *node, const IrradianceSample &sample, int depth);

	/// Pre-process step (clusters the samples)
	void preprocess(OctreeNode *node);

	/// Query the octree using a customizable implementation
	template <typename QueryType> void execute(OctreeNode *node, QueryType &query) const {
		if (node->leaf) {
			for (sample_iterator it = node->samples.begin();
				it != node->samples.end(); ++it)
				query(*it);
		} else {
			Float approxSolidAngle = node->cluster.area / (query.p - node->cluster.p).lengthSquared();
			if (!node->aabb.contains(query.p) && approxSolidAngle < m_threshold) {
				query(node->cluster);
			} else {
				for (int i=0; i<8; i++)
					if (node->children[i])
						execute(node->children[i], query);
			}
		}
	}

	/// Destruct the tree
	virtual ~IrradianceOctree();
private:
	OctreeNode *m_root;
	int m_maxDepth;
	unsigned int m_numSamples;
	Float m_threshold;
};

MTS_NAMESPACE_END

#endif /* __IRRTREE_H */

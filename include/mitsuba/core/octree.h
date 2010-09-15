/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__OCTREE_H)
#define __OCTREE_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Generic multiple-reference octree.
 *
 * Based on the excellent implementation in PBRT. Modifications are 
 * the addition of a bounding sphere query and support for multithreading.
 */
template <typename T> class Octree {
public:
	inline Octree(const AABB &aabb, int maxDepth = 16) 
	 : m_aabb(aabb), m_maxDepth(maxDepth) {
	}

	inline void insert(const T &value, const AABB &coverage) {
		insert(&m_root, m_aabb, value, coverage,
			coverage.getExtents().lengthSquared(), 0);
	}

	/// Execute operator() of <tt>functor</tt> on all records, which potentially overlap <tt>p</tt>
	template <typename Functor> inline void lookup(const Point &p, Functor &functor) const {
		if (!m_aabb.contains(p))
			return;
		lookup(&m_root, m_aabb, p, functor);
	}

	/// Execute operator() of <tt>functor</tt> on all records, which potentially overlap <tt>bsphere</tt>
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
			pthread_rwlock_init(&lock, NULL);
			for (int i=0; i<8; ++i)
				children[i] = NULL;
		}

		~OctreeNode() {
			pthread_rwlock_destroy(&lock);
			for (int i=0; i<8; ++i) {
				if (children[i])
					delete children[i];
			}
		}

		inline void readLock() const { pthread_rwlock_rdlock(&lock); }
		inline void readUnlock() const { pthread_rwlock_unlock(&lock); }
		inline void writeLock() { pthread_rwlock_wrlock(&lock); }
		inline void writeUnlock() { pthread_rwlock_unlock(&lock); }

		OctreeNode *children[8];
		mutable pthread_rwlock_t lock;
		std::vector<T> data;
	};

	void insert(OctreeNode *node, const AABB &nodeAABB, const T &value, 
			const AABB &coverage, Float diag2, int depth) {
		/* Add the data item to the current octree node if the max. tree
		   depth is reached or the data item's coverage area is smaller
		   than the current node size */
		if (depth == m_maxDepth || 
			(nodeAABB.getExtents().lengthSquared() < diag2)) {
			node->writeLock();
			node->data.push_back(value);
			node->writeUnlock();
			return;
		}

		/* Otherwise: test for overlap */
		const Point center = nodeAABB.getCenter();
		bool over[8];
		AABB childAABB;

		over[0] = over[1] = over[2] = over[3] = (coverage.min.x <= center.x);
		over[4] = over[5] = over[6] = over[7] = (coverage.max.x  > center.x);
		over[0] &= (coverage.min.y <= center.y);
		over[1] &= (coverage.min.y <= center.y);
		over[4] &= (coverage.min.y <= center.y);
		over[5] &= (coverage.min.y <= center.y);
		over[2] &= (coverage.max.y  > center.y);
		over[3] &= (coverage.max.y  > center.y);
		over[6] &= (coverage.max.y  > center.y);
		over[7] &= (coverage.max.y  > center.y);
		over[0] &= (coverage.min.z <= center.z);
		over[2] &= (coverage.min.z <= center.z);
		over[4] &= (coverage.min.z <= center.z);
		over[6] &= (coverage.min.z <= center.z);
		over[1] &= (coverage.max.z  > center.z);
		over[3] &= (coverage.max.z  > center.z);
		over[5] &= (coverage.max.z  > center.z);
		over[7] &= (coverage.max.z  > center.z);

		/* Recurse */
		for (int child=0; child<8; ++child) {
			if (!over[child])
				continue;
			if (!node->children[child]) {
				node->writeLock();
				node->children[child] = new OctreeNode();
				node->writeUnlock();
			}
			childAABB.min.x = (child & 4) ? center.x : nodeAABB.min.x;
			childAABB.max.x = (child & 4) ? nodeAABB.max.x : center.x;
			childAABB.min.y = (child & 2) ? center.y : nodeAABB.min.y;
			childAABB.max.y = (child & 2) ? nodeAABB.max.y : center.y;
			childAABB.min.z = (child & 1) ? center.z : nodeAABB.min.z;
			childAABB.max.z = (child & 1) ? nodeAABB.max.z : center.z;
			insert(node->children[child], childAABB,
				value, coverage, diag2, depth+1);
		}
	}

	/// Internal lookup procedure - const version
	template <typename Functor> inline void lookup(const OctreeNode *node, 
			const AABB &nodeAABB, const Point &p, Functor &functor) const {
		const Point center = nodeAABB.getCenter();

		node->readLock();
		for (size_t i=0; i<node->data.size(); ++i)
			functor(node->data[i]);

		int child = (p.x > center.x ? 4 : 0)
				+ (p.y > center.y ? 2 : 0) 
				+ (p.z > center.z ? 1 : 0);

		OctreeNode *childNode = node->children[child];
		node->readUnlock();

		if (childNode) {
			AABB childAABB;
			childAABB.min.x = (child & 4) ? center.x : nodeAABB.min.x;
			childAABB.max.x = (child & 4) ? nodeAABB.max.x : center.x;
			childAABB.min.y = (child & 2) ? center.y : nodeAABB.min.y;
			childAABB.max.y = (child & 2) ? nodeAABB.max.y : center.y;
			childAABB.min.z = (child & 1) ? center.z : nodeAABB.min.z;
			childAABB.max.z = (child & 1) ? nodeAABB.max.z : center.z;
			lookup(node->children[child], childAABB, p, functor);
		}
	}

	template <typename Functor> inline void searchSphere(OctreeNode *node, 
			const AABB &nodeAABB, const BSphere &sphere, 
			Functor &functor) {
		const Point center = nodeAABB.getCenter();

		node->readLock();
		for (size_t i=0; i<node->data.size(); ++i)
			functor(node->data[i]);
		node->readUnlock();

		// Potential for much optimization..
		for (int child=0; child<8; ++child) { 
			if (node->children[child]) {
				AABB childAABB;
				childAABB.min.x = (child & 4) ? center.x : nodeAABB.min.x;
				childAABB.max.x = (child & 4) ? nodeAABB.max.x : center.x;
				childAABB.min.y = (child & 2) ? center.y : nodeAABB.min.y;
				childAABB.max.y = (child & 2) ? nodeAABB.max.y : center.y;
				childAABB.min.z = (child & 1) ? center.z : nodeAABB.min.z;
				childAABB.max.z = (child & 1) ? nodeAABB.max.z : center.z;

				if (childAABB.overlaps(sphere))
					searchSphere(node->children[child], childAABB, sphere, functor);
			}
		}
	}
private:
	OctreeNode m_root;
	AABB m_aabb;
	int m_maxDepth;
};

MTS_NAMESPACE_END

#endif /* __OCTREE_H */

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

#if !defined(__SAH_KDTREE4_H)
#define __SAH_KDTREE4_H

#include <mitsuba/render/sahkdtree3.h>

MTS_NAMESPACE_BEGIN

typedef TAABB<Point4> AABB4;

/**
 * \brief Implements the four-dimensional surface area heuristic for use
 * by the \ref GenericKDTree construction algorithm.
 */
class SurfaceAreaHeuristic4 {
public:
	/**
	 * \brief Initialize the surface area heuristic with the bounds of
	 * a parent node
	 *
	 * Precomputes some information so that traversal probabilities
	 * of potential split planes can be evaluated efficiently
	 */
	inline SurfaceAreaHeuristic4(const AABB4 &aabb) : m_aabb(aabb) {
		const Vector4 extents(aabb.getExtents());
		const Float temp = 1.0f / (extents.x * extents.y
				+ extents.y*extents.z + extents.x*extents.z);

		m_temp0 = Vector4(
			extents.y * extents.z * temp,
			extents.x * extents.z * temp,
			extents.x * extents.y * temp,
			0.0f);

		m_temp1 = Vector4(
			(extents.y + extents.z) * temp,
			(extents.x + extents.z) * temp,
			(extents.x + extents.y) * temp,
			1.0f / extents.w);
	}

	/**
	 * Given a split on axis \a axis that produces children having extents
	 * \a leftWidth and \a rightWidth along \a axis, compute the probability
	 * of traversing the left and right child during a typical query
	 * operation.
	 */
	inline std::pair<Float, Float> operator()(int axis, Float leftWidth, Float rightWidth) const {
		if (axis == 3 && m_temp1.w == std::numeric_limits<Float>::infinity()) {
			return std::pair<Float, Float>(
				std::numeric_limits<Float>::infinity(),
				std::numeric_limits<Float>::infinity()
			);
		}

		return std::pair<Float, Float>(
			m_temp0[axis] + m_temp1[axis] * leftWidth,
			m_temp0[axis] + m_temp1[axis] * rightWidth);
	}

	/**
	 * Compute the underlying quantity used by the tree construction
	 * heuristic. This is used to compute the final cost of a kd-tree.
	 */
	inline static Float getQuantity(const AABB4 &aabb) {
		const Vector4 extents(aabb.getExtents());
		Float result = 2 * (extents[0] * extents[1] + extents[1] * extents[2]
				  + extents[2] * extents[0]);
		if (extents[3] != 0)
			result *= extents[3];
		return result;
	}
private:
	Vector4 m_temp0, m_temp1;
	AABB4 m_aabb;
};

/**
 * This class specializes \ref GenericKDTree to a four-dimensional
 * tree to be used for spacetime ray tracing. One additional function call
 * must be implemented by subclasses:
 *
 * /// Check whether a primitive is intersected by the given ray.
 * /// Some temporary space is supplied, which can be used to cache
 * /// information about the intersection
 * bool intersect(const Ray &ray, IndexType idx,
 *     Float mint, Float maxt, Float &t, void *tmp);
 *
 * This class implements an epsilon-free version of the optimized ray
 * traversal algorithm (TA^B_{rec}), which is explained in Vlastimil
 * Havran's PhD thesis "Heuristic Ray Shooting Algorithms".
 *
 * \author Wenzel Jakob
 */
template <typename Derived>
	class SAHKDTree4D : public GenericKDTree<AABB4, SurfaceAreaHeuristic4, Derived> {
public:
	typedef typename KDTreeBase<AABB4>::SizeType  SizeType;
	typedef typename KDTreeBase<AABB4>::IndexType IndexType;
	typedef typename KDTreeBase<AABB4>::KDNode     KDNode;

protected:
	void buildInternal() {
		SizeType primCount = this->cast()->getPrimitiveCount();
		KDLog(EInfo, "Constructing a 4D SAH kd-tree (%i primitives) ..", primCount);
		GenericKDTree<AABB4, SurfaceAreaHeuristic4, Derived>::buildInternal();
	}

	/**
	 * \brief Hashed mailbox implementation
	 */
	struct HashedMailbox {
		inline HashedMailbox() {
			memset(entries, 0xFF, sizeof(IndexType)*MTS_KD_MAILBOX_SIZE);
		}

		inline void put(IndexType primIndex) {
			entries[primIndex & MTS_KD_MAILBOX_MASK] = primIndex;
		}

		inline bool contains(IndexType primIndex) const {
			return entries[primIndex & MTS_KD_MAILBOX_MASK] == primIndex;
		}

		IndexType entries[MTS_KD_MAILBOX_SIZE];
	};

	/// Ray traversal stack entry for Havran-style incoherent ray tracing
	struct KDStackEntryHavran {
		/* Pointer to the far child */
		const KDNode * __restrict node;
		/* Distance traveled along the ray (entry or exit) */
		Float t;
		/* Previous stack item */
		uint32_t prev;
		/* Associated point */
		Point p;
	};

	/**
	 * \brief Ray tracing kd-tree traversal loop (Havran variant)
	 *
	 * This is generally the most robust and fastest traversal routine
	 * of the methods implemented in this class.
	 */
	template<bool shadowRay> FINLINE
			bool rayIntersectHavran(const Ray &ray, Float mint, Float maxt,
			Float &t, void *temp) const {
		KDStackEntryHavran stack[MTS_KD_MAXDEPTH];
		#if 0
		static const int prevAxisTable[] = { 2, 0, 1 };
		static const int nextAxisTable[] = { 1, 2, 0 };
		#endif

		#if defined(MTS_KD_MAILBOX_ENABLED)
		HashedMailbox mailbox;
		#endif

		/* Set up the entry point */
		uint32_t enPt = 0;
		stack[enPt].t = mint;
		stack[enPt].p = ray(mint);

		/* Set up the exit point */
		uint32_t exPt = 1;
		stack[exPt].t = maxt;
		stack[exPt].p = ray(maxt);
		stack[exPt].node = NULL;

		bool foundIntersection = false;
		const KDNode * __restrict currNode = this->m_nodes;
		while (currNode != NULL) {
			while (EXPECT_TAKEN(!currNode->isLeaf())) {
				const Float splitVal = (Float) currNode->getSplit();
				const int axis = currNode->getAxis();
				const KDNode * __restrict farChild;

				if (axis == 3) {
					if (ray.time <= splitVal)
						currNode = currNode->getLeft();
					else
						currNode = currNode->getRight();
					continue;
				} else if (stack[enPt].p[axis] <= splitVal) {
					if (stack[exPt].p[axis] <= splitVal) {
						/* Cases N1, N2, N3, P5, Z2 and Z3 (see thesis) */
						currNode = currNode->getLeft();
						continue;
					}

					/* Typo in Havran's thesis:
					   (it specifies "stack[exPt].p == splitVal", which
					    is clearly incorrect) */
					if (stack[enPt].p[axis] == splitVal) {
						/* Case Z1 */
						currNode = currNode->getRight();
						continue;
					}

					/* Case N4 */
					currNode = currNode->getLeft();
					farChild = currNode + 1; // getRight()
				} else { /* stack[enPt].p[axis] > splitVal */
					if (splitVal < stack[exPt].p[axis]) {
						/* Cases P1, P2, P3 and N5 */
						currNode = currNode->getRight();
						continue;
					}
					/* Case P4 */
					farChild = currNode->getLeft();
					currNode = farChild + 1; // getRight()
				}

				/* Cases P4 and N4 -- calculate the distance to the split plane */
				Float distToSplit = (splitVal - ray.o[axis]) * ray.dRcp[axis];

				/* Set up a new exit point */
				const uint32_t tmp = exPt++;
				if (exPt == enPt) /* Do not overwrite the entry point */
					++exPt;

				KDAssert(exPt < MTS_KD_MAXDEPTH);
				stack[exPt].prev = tmp;
				stack[exPt].t = distToSplit;
				stack[exPt].node = farChild;

				#if 1
				/* Intrestingly, this appears to be faster than the
				   original code with the prevAxis & nextAxis table */
				stack[exPt].p = ray(distToSplit);
				stack[exPt].p[axis] = splitVal;
				#else
				const int nextAxis = nextAxisTable[axis];
				const int prevAxis = prevAxisTable[axis];
				stack[exPt].p[axis] = splitVal;
				stack[exPt].p[nextAxis] = ray.o[nextAxis]
					+ distToSplit*ray.d[nextAxis];
				stack[exPt].p[prevAxis] = ray.o[prevAxis]
					+ distToSplit*ray.d[prevAxis];
				#endif

			}

			/* Reached a leaf node */
			for (IndexType entry=currNode->getPrimStart(),
					last = currNode->getPrimEnd(); entry != last; entry++) {
				const IndexType primIdx = this->m_indices[entry];

				#if defined(MTS_KD_MAILBOX_ENABLED)
				if (mailbox.contains(primIdx))
					continue;
				#endif

				bool result;
				if (!shadowRay)
					result = this->cast()->intersect(ray, primIdx, mint, maxt, t, temp);
				else
					result = this->cast()->intersect(ray, primIdx, mint, maxt);

				if (result) {
					if (shadowRay)
						return true;
					maxt = t;
					foundIntersection = true;
				}

				#if defined(MTS_KD_MAILBOX_ENABLED)
				mailbox.put(primIdx);
				#endif
			}

			if (stack[exPt].t > maxt)
				break;

			/* Pop from the stack and advance to the next node on the interval */
			enPt = exPt;
			currNode = stack[exPt].node;
			exPt = stack[enPt].prev;
		}

		return foundIntersection;
	}
};

MTS_NAMESPACE_END

#endif /* __SAH_KDTREE4_H */

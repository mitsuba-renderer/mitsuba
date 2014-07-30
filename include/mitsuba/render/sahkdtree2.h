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
#if !defined(__MITSUBA_RENDER_SAHKDTREE2_H_)
#define __MITSUBA_RENDER_SAHKDTREE2_H_

#include <mitsuba/core/aabb.h>
#include <mitsuba/render/gkdtree.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Implements the 2D surface area heuristic for use
 * by the \ref GenericKDTree construction algorithm.
 * \ingroup librender
 */
class SurfaceAreaHeuristic2 {
public:
	/**
	 * \brief Initialize the surface area heuristic with the bounds of
	 * a parent node
	 *
	 * Precomputes some information so that traversal probabilities
	 * of potential split planes can be evaluated efficiently
	 */
	inline SurfaceAreaHeuristic2(const AABB2 &aabb) {
		m_extents = aabb.getExtents();
		m_normalization = 1.0f / (m_extents.x + m_extents.y);
	}

	/**
	 * Given a split on axis \a axis that produces children having extents
	 * \a leftWidth and \a rightWidth along \a axis, compute the probability
	 * of traversing the left and right child during a typical query
	 * operation.
	 */
	inline std::pair<Float, Float> operator()(int axis, Float leftWidth, Float rightWidth) const {
		return std::pair<Float, Float>(
			(m_extents[1-axis] + leftWidth) * m_normalization,
			(m_extents[1-axis] + rightWidth) * m_normalization);
	}

	/**
	 * Compute the underlying quantity used by the tree construction
	 * heuristic. This is used to compute the final cost of a kd-tree.
	 */
	inline static Float getQuantity(const AABB2 &aabb) {
		return aabb.getSurfaceArea();
	}
private:
	Float m_normalization;
	Vector2 m_extents;
};

/**
 * \brief Specializes \ref GenericKDTree to a two-dimensional
 * tree to be used for flatland ray tracing.
 *
 * One additional function call must be implemented by subclasses:
 * \code
 * /// Check whether a primitive is intersected by the given ray.
 * /// Some temporary space is supplied, which can be used to cache
 * /// information about the intersection
 * bool intersect(const Ray2 &ray, IndexType idx,
 *     Float mint, Float maxt, Float &t, void *tmp);
 * \endcode
 *
 * This class implements an epsilon-free version of the optimized ray
 * traversal algorithm (TA^B_{rec}), which is explained in Vlastimil
 * Havran's PhD thesis "Heuristic Ray Shooting Algorithms".
 *
 * \author Wenzel Jakob
 * \ingroup librender
 */
template <typename Derived>
	class SAHKDTree2D : public GenericKDTree<AABB2, SurfaceAreaHeuristic2, Derived> {
public:
	typedef GenericKDTree<AABB2, SurfaceAreaHeuristic2, Derived> Parent;
	typedef typename KDTreeBase<AABB2>::SizeType                 SizeType;
	typedef typename KDTreeBase<AABB2>::IndexType                IndexType;
	typedef typename KDTreeBase<AABB2>::KDNode                   KDNode;

	using Parent::m_nodes;
	using Parent::m_aabb;
	using Parent::m_indices;

protected:
	void buildInternal() {
		SizeType primCount = cast()->getPrimitiveCount();
		KDLog(EInfo, "Constructing a SAH kd-tree (%i primitives) ..", primCount);
		GenericKDTree<AABB2, SurfaceAreaHeuristic2, Derived>::buildInternal();
	}

	/// Cast to the derived class
	inline Derived *cast() {
		return static_cast<Derived *>(this);
	}

	/// Cast to the derived class (const version)
	inline const Derived *cast() const {
		return static_cast<const Derived *>(this);
	}

	/// Ray traversal stack entry for Wald-style incoherent ray tracing
	struct KDStackEntry {
		const KDNode * __restrict node;
		Float mint, maxt;
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
		Point2 p;
	};

	/**
	 * \brief Ray tracing kd-tree traversal loop (Havran variant)
	 *
	 * This is generally the most robust and fastest traversal routine
	 * of the methods implemented in this class.
	 */
	FINLINE bool rayIntersectHavran(const Ray2 &ray, Float mint, Float maxt,
			Float &t, void *temp) const {
		KDStackEntryHavran stack[MTS_KD_MAXDEPTH];

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
		const KDNode * __restrict currNode = m_nodes;
		while (currNode != NULL) {
			while (EXPECT_TAKEN(!currNode->isLeaf())) {
				const Float splitVal = (Float) currNode->getSplit();
				const int axis = currNode->getAxis();
				const KDNode * __restrict farChild;

				if (stack[enPt].p[axis] <= splitVal) {
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

				/* Intrestingly, this appears to be faster than the
				   original code with the prevAxis & nextAxis table */
				stack[exPt].p = ray(distToSplit);
				stack[exPt].p[axis] = splitVal;
			}

			/* Reached a leaf node */
			for (IndexType entry=currNode->getPrimStart(),
					last = currNode->getPrimEnd(); entry != last; entry++) {
				const IndexType primIdx = m_indices[entry];

				bool result = cast()->intersect(ray, primIdx, mint, maxt, t, temp);

				if (result) {
					maxt = t;
					foundIntersection = true;
				}
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

#endif /* __MITSUBA_RENDER_SAHKDTREE2_H_ */

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
#if !defined(__MITSUBA_RENDER_SAHKDTREE3_H_)
#define __MITSUBA_RENDER_SAHKDTREE3_H_

#include <mitsuba/render/gkdtree.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/random.h>

MTS_NAMESPACE_BEGIN

/// Use a simple hashed 8-entry mailbox per thread
#define MTS_KD_MAILBOX_ENABLED 1
#define MTS_KD_MAILBOX_SIZE 8
#define MTS_KD_MAILBOX_MASK (MTS_KD_MAILBOX_SIZE-1)

/**
 * \brief Implements the 3D surface area heuristic for use
 * by the \ref GenericKDTree construction algorithm.
 * \ingroup librender
 */
class SurfaceAreaHeuristic3 {
public:
	/**
	 * \brief Initialize the surface area heuristic with the bounds of
	 * a parent node
	 *
	 * Precomputes some information so that traversal probabilities
	 * of potential split planes can be evaluated efficiently
	 */
	inline SurfaceAreaHeuristic3(const AABB &aabb) {
		const Vector extents(aabb.getExtents());
		const Float temp = 1.0f / (extents.x * extents.y
				+ extents.y*extents.z + extents.x*extents.z);
		m_temp0 = Vector(
			extents[1] * extents[2],
			extents[0] * extents[2],
			extents[0] * extents[1]) * temp;
		m_temp1 = Vector(
			extents[1] + extents[2],
			extents[0] + extents[2],
			extents[0] + extents[1]) * temp;
	}

	/**
	 * Given a split on axis \a axis that produces children having extents
	 * \a leftWidth and \a rightWidth along \a axis, compute the probability
	 * of traversing the left and right child during a typical query
	 * operation. In the case of the surface area heuristic, this is simply
	 * the ratio of surface areas.
	 */
	inline std::pair<Float, Float> operator()(int axis, Float leftWidth, Float rightWidth) const {
		return std::pair<Float, Float>(
			m_temp0[axis] + m_temp1[axis] * leftWidth,
			m_temp0[axis] + m_temp1[axis] * rightWidth);
	}

	/**
	 * Compute the underlying quantity used by the tree construction
	 * heuristic. This is used to compute the final cost of a kd-tree.
	 */
	inline static Float getQuantity(const AABB &aabb) {
		return aabb.getSurfaceArea();
	}
private:
	Vector m_temp0, m_temp1;
};

/**
 * \brief Specializes \ref GenericKDTree to a three-dimensional
 * tree to be used for ray tracing.
 *
 * One additional function call must be implemented by subclasses:
 * \code
 * /// Check whether a primitive is intersected by the given ray.
 * /// Some temporary space is supplied, which can be used to cache
 * /// information about the intersection
 * bool intersect(const Ray &ray, IndexType idx,
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
	class SAHKDTree3D : public GenericKDTree<AABB, SurfaceAreaHeuristic3, Derived> {
public:
	typedef GenericKDTree<AABB, SurfaceAreaHeuristic3, Derived> Parent;
	typedef typename KDTreeBase<AABB>::SizeType                 SizeType;
	typedef typename KDTreeBase<AABB>::IndexType                IndexType;
	typedef typename KDTreeBase<AABB>::KDNode                   KDNode;

	using Parent::m_nodes;
	using Parent::m_aabb;
	using Parent::m_indices;

protected:
	void buildInternal() {
		SizeType primCount = cast()->getPrimitiveCount();
		KDLog(EInfo, "Constructing a SAH kd-tree (%i primitives) ..", primCount);
		GenericKDTree<AABB, SurfaceAreaHeuristic3, Derived>::buildInternal();
	}

	/// Cast to the derived class
	inline Derived *cast() {
		return static_cast<Derived *>(this);
	}

	/// Cast to the derived class (const version)
	inline const Derived *cast() const {
		return static_cast<const Derived *>(this);
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
				const IndexType primIdx = m_indices[entry];

				#if defined(MTS_KD_MAILBOX_ENABLED)
				if (mailbox.contains(primIdx))
					continue;
				#endif

				bool result;
				if (!shadowRay)
					result = cast()->intersect(ray, primIdx, mint, maxt, t, temp);
				else
					result = cast()->intersect(ray, primIdx, mint, maxt);

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

	struct RayStatistics {
		bool foundIntersection;
		uint32_t numTraversals;
		uint32_t numIntersections;
		uint64_t time;

		RayStatistics(bool foundIntersection, uint32_t numTraversals,
			uint32_t numIntersections, uint64_t time) :
			foundIntersection(foundIntersection), numTraversals(numTraversals),
			numIntersections(numIntersections), time(time) { }
	};

	/**
	 * \brief Internal kd-tree traversal implementation (Havran variant)
	 *
	 * This method is almost identical to \ref rayIntersectHavran, except
	 * that it additionally returns statistics on the number of traversed
	 * nodes, intersected shapes, as well as the time taken to do this
	 * (measured using rtdsc).
	 */
	FINLINE RayStatistics rayIntersectHavranCollectStatistics(
			const Ray &ray, Float mint, Float maxt, Float &t, void *temp) const {
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

		uint32_t numTraversals = 0;
		uint32_t numIntersections = 0;
		uint64_t timer = rdtsc();
		bool foundIntersection = false;

		const KDNode * __restrict currNode = m_nodes;
		while (currNode != NULL) {
			while (EXPECT_TAKEN(!currNode->isLeaf())) {
				const Float splitVal = (Float) currNode->getSplit();
				const int axis = currNode->getAxis();
				const KDNode * __restrict farChild;

				++numTraversals;
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
				t = (splitVal - ray.o[axis]) * ray.dRcp[axis];

				/* Set up a new exit point */
				const uint32_t tmp = exPt++;
				if (exPt == enPt) /* Do not overwrite the entry point */
					++exPt;

				KDAssert(exPt < MTS_KD_MAXDEPTH);
				stack[exPt].prev = tmp;
				stack[exPt].t = t;
				stack[exPt].node = farChild;
				stack[exPt].p = ray(t);
				stack[exPt].p[axis] = splitVal;
			}

			/* Reached a leaf node */
			for (unsigned int entry=currNode->getPrimStart(),
					last = currNode->getPrimEnd(); entry != last; entry++) {
				const IndexType primIdx = m_indices[entry];

				++numIntersections;
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

		return RayStatistics(foundIntersection, numTraversals,
				numIntersections, rdtsc() - timer);
	}

	/**
	 * \brief Ray tracing kd-tree traversal loop (PBRT variant)
	 */
	template<bool shadowRay> FINLINE bool rayIntersectPBRT(const Ray &ray,
			Float mint_, Float maxt_, Float &t, void *temp) const {
		KDStackEntry stack[MTS_KD_MAXDEPTH];
		int stackPos = 0;
		Float mint = mint_, maxt=maxt_;
		const KDNode *node = m_nodes;
		bool foundIntersection = false;

		while (node != NULL) {
			if (maxt_ < mint)
				break;

			if (EXPECT_TAKEN(!node->isLeaf())) {
				const Float split = (Float) node->getSplit();
				const int axis = node->getAxis();
				const float tPlane = (split - ray.o[axis]) * ray.dRcp[axis];
				bool leftOfSplit = (ray.o[axis] < split)
					|| (ray.o[axis] == split && ray.d[axis] <= 0);

				const KDNode * __restrict left = node->getLeft();
				const KDNode * __restrict right = left + 1;
				const KDNode * __restrict first  = leftOfSplit ? left : right;
				const KDNode * __restrict second = leftOfSplit ? right : left;

				if (tPlane > maxt || tPlane <= 0) {
					node = first;
				} else if (tPlane < mint) {
					node = second;
				} else {
					stack[stackPos].node = second;
					stack[stackPos].mint = tPlane;
					stack[stackPos].maxt = maxt;
					++stackPos;
					node = first;
					maxt = tPlane;
				}
			} else {
				for (unsigned int entry=node->getPrimStart(),
						last = node->getPrimEnd(); entry != last; entry++) {
					const IndexType primIdx = m_indices[entry];

					bool result;
					if (!shadowRay)
						result = cast()->intersect(ray, primIdx, mint, maxt, t, temp);
					else
						result = cast()->intersect(ray, primIdx, mint, maxt);

					if (result) {
						if (shadowRay)
							return true;
						maxt_ = t;
						foundIntersection = true;
					}
				}

				if (stackPos > 0) {
					--stackPos;
					node = stack[stackPos].node;
					mint = stack[stackPos].mint;
					maxt = stack[stackPos].maxt;
				} else {
					break;
				}
			}
		}
		return foundIntersection;
	}
public:
	/**
	 * \brief Empirically find the best traversal and intersection
	 * cost values
	 *
	 * This is done by running the traversal code on random rays
	 * and fitting the SAH cost model to the collected statistics.
	 */
	void findCosts(Float &traversalCost, Float &intersectionCost) {
		ref<Random> random = new Random();
		uint8_t temp[128];
		BSphere bsphere = m_aabb.getBSphere();
		int nRays = 10000000, warmup = nRays/4;
		Vector *A = new Vector[nRays-warmup];
		Float *b = new Float[nRays-warmup];
		int nIntersections = 0, idx = 0;

		for (int i=0; i<nRays; ++i) {
			Point2 sample1(random->nextFloat(), random->nextFloat()),
				sample2(random->nextFloat(), random->nextFloat());
			Point p1 = bsphere.center + warp::squareToUniformSphere(sample1) * bsphere.radius;
			Point p2 = bsphere.center + warp::squareToUniformSphere(sample2) * bsphere.radius;
			Ray ray(p1, normalize(p2-p1), 0.0f);
			Float mint, maxt, t;
			if (m_aabb.rayIntersect(ray, mint, maxt)) {
				if (ray.mint > mint) mint = ray.mint;
				if (ray.maxt < maxt) maxt = ray.maxt;
				if (EXPECT_TAKEN(maxt > mint)) {
					RayStatistics statistics =
						rayIntersectHavranCollectStatistics(ray, mint, maxt, t, temp);
					if (statistics.foundIntersection)
						nIntersections++;
					if (i > warmup) {
						A[idx].x = 1;
						A[idx].y = (Float) statistics.numTraversals;
						A[idx].z = (Float) statistics.numIntersections;
						b[idx]   = (Float) statistics.time;
						idx++;
					}
				}
			}
		}

		KDLog(EDebug, "Fitting to " SIZE_T_FMT " samples (" SIZE_T_FMT
				" intersections)", idx, nIntersections);

		/* Solve using normal equations */
		Matrix4x4 M(0.0f), Minv;
		Vector4 rhs(0.0f), x;

		for (int i=0; i<3; ++i) {
			for (int j=0; j<3; ++j)
				for (int k=0; k<idx; ++k)
					M.m[i][j] += A[k][i]*A[k][j];
			for (int k=0; k<idx; ++k)
				rhs[i] += A[k][i]*b[k];
		}
		M.m[3][3] = 1.0f;
		bool success = M.invert(Minv);
		SAssert(success);

		Transform(Minv, M)(rhs, x);

		Float avgRdtsc = 0, avgResidual = 0;
		for (int i=0; i<idx; ++i) {
			avgRdtsc += b[i];
			Float model = x[0] * A[i][0]
				+ x[1] * A[i][1]
				+ x[2] * A[i][2];
			avgResidual += std::abs(b[i] - model);
		}
		avgRdtsc /= idx;
		avgResidual /= idx;

		for (int k=0; k<idx; ++k)
			avgRdtsc += b[k];
		avgRdtsc /= idx;

		KDLog(EDebug, "Least squares fit:");
		KDLog(EDebug, "   Constant overhead    = %.2f", x[0]);
		KDLog(EDebug, "   Traversal cost       = %.2f", x[1]);
		KDLog(EDebug, "   Intersection cost    = %.2f", x[2]);
		KDLog(EDebug, "   Average rdtsc value  = %.2f", avgRdtsc);
		KDLog(EDebug, "   Avg. residual        = %.2f", avgResidual);
		x *= 10/x[1];
		KDLog(EDebug, "Re-scaled:");
		KDLog(EDebug, "   Constant overhead    = %.2f", x[0]);
		KDLog(EDebug, "   Traversal cost       = %.2f", x[1]);
		KDLog(EDebug, "   Intersection cost    = %.2f", x[2]);

		delete[] A;
		delete[] b;
		traversalCost = x[1];
		intersectionCost = x[2];
	}
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_SAHKDTREE3_H_ */

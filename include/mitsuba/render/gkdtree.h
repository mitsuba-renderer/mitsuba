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

#if !defined(__KDTREE_GENERIC_H)
#define __KDTREE_GENERIC_H

#include <mitsuba/core/timer.h>
#include <boost/static_assert.hpp>

#define MTS_KD_MAX_DEPTH 48   ///< Compile-time KD-tree depth limit
#define MTS_KD_STATISTICS 1   ///< Collect statistics during building/traversal 
#define MTS_KD_MINMAX_BINS 32

MTS_NAMESPACE_BEGIN

/**
 * \brief Generic kd-tree for data structure used to accelerate 
 * ray intersections against large amounts of three-dimensional
 * shapes.
 * 
 */
template <typename ShapeType> class GenericKDTree : public Object {
public:
	/// Index number format (max 2^32 prims)
	typedef uint32_t index_type;
	/// Size number format
	typedef uint32_t size_type;

	/**
	 * \brief Data type for split candidates computed by 
	 * the SAH optimization routines.
	 * */
	struct SplitCandidate {
		Float cost;
		float pos;
		int axis;
		int numLeft, numRight;
		bool planarLeft;
	};

	GenericKDTree() {
		m_traversalCost = 15;
		m_intersectionCost = 20;
		m_emptySpaceBonus = 0.8f;
	}

	/**
	 * \brief Build a KD-tree over the supplied axis-aligned bounding
	 * boxes.
	 */
	template <typename Functor>
	void build(const Functor &functor) {
		/* Establish an ad-hoc depth cutoff value (Formula from PBRT) */
		m_maxDepth = std::min((int) (8 + 1.3f * log2i(functor.getPrimitiveCount())),
			MTS_KD_MAX_DEPTH);

		ref<Timer> timer = new Timer();
		m_aabb.reset();
		for (size_type i=0; i<functor.getPrimitiveCount(); ++i)
			m_aabb.expandBy(functor.getAABB(i));
		Log(EInfo, "Computed scene bounds in %i ms", timer->getMilliseconds());

		Log(EDebug, "kd-tree configuration:");
		Log(EDebug, "   Traversal cost           : %.2f", m_traversalCost);
		Log(EDebug, "   Intersection cost        : %.2f", m_intersectionCost);
		Log(EDebug, "   Empty space bonus        : %.2f", m_emptySpaceBonus);
		Log(EDebug, "   Max. tree depth          : %i", m_maxDepth);
		Log(EDebug, "   Scene bounding box (min) : %s", m_aabb.min.toString().c_str());
		Log(EDebug, "   Scene bounding box (max) : %s", m_aabb.max.toString().c_str());
	
		Log(EInfo, "Constructing a SAH kd-tree (%i primitives) ..", 
				functor.getPrimitiveCount());

		timer->reset();

		Log(EInfo, "Performing initial clustering ..");
		MinMaxBins<MTS_KD_MINMAX_BINS> bins(m_aabb);
		bins.bin(functor);
		SplitCandidate candidate = bins.maximizeSAH(m_traversalCost,
			m_intersectionCost, m_emptySpaceBonus);

		AABB leftAABB, rightAABB;
		bins.partition(functor, candidate, leftAABB, rightAABB,
				m_traversalCost, m_intersectionCost);

		cout << "originalCost = " << functor.getPrimitiveCount() * m_intersectionCost << endl;

		cout << "bestCost = " 
			<< candidate.cost
			<< ", bestPos = "
			<< candidate.pos
			<< ", bestAxis = "
			<< candidate.axis
			<< ", numLeft = "
			<< candidate.numLeft
			<< ", numRight = "
			<< candidate.numRight
			<< endl;
		Log(EInfo, "Initial clustering took %i ms", timer->getMilliseconds());
	}
protected:
	/**
	 * \brief Describes the beginning or end of a primitive
	 * when projected onto a certain dimension.
	 */
	struct EdgeEvent {
		/// Possible event types
		enum EEventType {
			EEdgeEnd = 0,
			EEdgePlanar = 1,
			EEdgeStart = 2
		};
	
		/// Dummy constructor
		inline EdgeEvent() { }

		/// Create a new edge event
		inline EdgeEvent(int type, float t, index_type index)
		 : t(t), index(index), type(type) { }

		/// Plane position
		float t;
		/// Primitive index
		index_type index;
		/// Event type: end/planar/start
		int type;
	};

	BOOST_STATIC_ASSERT(sizeof(EdgeEvent) == 12);

	typedef std::vector<EdgeEvent> edge_event_vector;
	typedef edge_event_vector edge_event_vector3[3];

	/// Edge event comparison functor
	struct EdgeEventOrdering : public std::binary_function<EdgeEvent, EdgeEvent, bool> {
		inline bool operator()(const EdgeEvent &a, const EdgeEvent &b) const {
			if (a.t != b.t)
				return a.t < b.t;
			return a.type < b.type;
		}
	};

	/// KD-tree node in 64bit
	struct KDNode {
		union {
			/* Inner node */
			struct {
				/* Bit layout:
				   31   : False (inner node)
				   30-3 : Offset to the right child
				   3-0  : Split axis
				*/
				uint32_t combined;

				/// Split plane coordinate
				float split;
			} inner;

			/* Leaf node */
			struct {
				/* Bit layout:
				   31   : True (leaf node)
				   30-0 : Offset to the node's primitive list
				*/
				uint32_t combined;

				/// End offset of the primitive list
				uint32_t end;
			} leaf;
		};

		enum EMask {
			ETypeMask = 1 << 31,
			ELeafOffsetMask = ~ETypeMask,
			EInnerAxisMask = 0x3,
			EInnerOffsetMask = ~EInnerAxisMask
		};

		/// Initialize a leaf kd-Tree node
		inline void setLeaf(unsigned int offset, unsigned int numPrims) {
			leaf.combined = ETypeMask | offset;
			leaf.end = offset + numPrims;
		}

		/// Initialize an interior kd-Tree node
		inline void setInner(int axis, unsigned int offset, Float split) {
			inner.combined = axis | (offset << 2);
			inner.split = (float) split;
		}

		/// Is this a leaf node?
		FINLINE bool isLeaf() const {
			return leaf.combined & ETypeMask;
		}

		/// Assuming this is a leaf node, return the first primitive index
		FINLINE index_type getPrimStart() const {
			return leaf.combined & ELeafOffsetMask;
		}

		/// Assuming this is a leaf node, return the last primitive index
		FINLINE index_type getPrimEnd() const {
			return leaf.end;
		}
		
		/// Return the sibling of this node
		FINLINE const KDNode * __restrict getSibling() const {
			return (const KDNode *) ((ptrdiff_t) this ^ (ptrdiff_t) 8);
		}

		/// Return the left child (assuming that this is an interior node)
		FINLINE const KDNode * __restrict getLeft() const {
			return this + 
				((inner.combined & EInnerOffsetMask) >> 2);
		}
		
		/// Return the right child (assuming that this is an interior node)
		FINLINE const KDNode * __restrict getRight() const {
			return getLeft() + 1;
		}

		/// Return the split plane location (assuming that this is an interior node)
		FINLINE float getSplit() const {
			return inner.split;
		}

		/// Return the split axis (assuming that this is an interior node)
		FINLINE int getAxis() const {
			return inner.combined & EInnerAxisMask;
		}
	};

	BOOST_STATIC_ASSERT(sizeof(KDNode) == 8);

	class ChunkAllocator {
	public:
		ChunkAllocator() {
		}
	};

	/**
	 * \brief Min-max binning as described in
	 * "Highly Parallel Fast KD-tree Construction for Interactive
	 *  Ray Tracing of Dynamic Scenes"
	 * by M. Shevtsov, A. Soupikov and A. Kapustin
	 *
	 * @tparam BinCount Number of bins to be allocated
	 */
	template <int BinCount = 32> struct MinMaxBins {
		MinMaxBins(AABB &aabb) : m_aabb(aabb) {
			m_binSize = m_aabb.getExtents() / BinCount;
		}

		/**
		 * \brief Run min-max binning
		 *
		 * \param func Functor to be used to determine the AABB for
		 *     a given list of primitives
		 * \param count Number of primitives
		 */
		template <typename Functor> void bin(
				const Functor &functor) {
			m_primCount = functor.getPrimitiveCount();
			memset(m_minBins, 0, sizeof(size_type) * 3 * BinCount);
			memset(m_maxBins, 0, sizeof(size_type) * 3 * BinCount);
			Vector invBinSize;

			for (int axis=0; axis<3; ++axis) 
				invBinSize[axis] = 1/m_binSize[axis];

			for (size_type i=0; i<m_primCount; ++i) {
				const AABB aabb = functor.getAABB(i);
				for (int axis=0; axis<3; ++axis) {
					int minIdx = (int) ((aabb.min[axis] - m_aabb.min[axis]) 
							* invBinSize[axis]);
					int maxIdx = (int) ((aabb.max[axis] - m_aabb.min[axis]) 
							* invBinSize[axis]);
					Assert(minIdx >= 0 && maxIdx >= 0);
					m_maxBins[axis * BinCount + std::min(maxIdx, BinCount-1)]++;
					m_minBins[axis * BinCount + std::min(minIdx, BinCount-1)]++;
				}
			}
		}

		/**
		 * \brief Add the bin counts of another \ref MinMaxBins instance.
		 *    This is used when distributing min-max-binning over multiple
		 *    processors in a SMP machine.
		 */
		MinMaxBins& operator+=(const MinMaxBins &otherBins) {
			m_primCount += otherBins.primCount;
			for (int i=0; i<3*BinCount; ++i) {
				m_minBins[i] += otherBins.m_minBins[i];
				m_maxBins[i] += otherBins.m_maxBins[i];
			}
			return *this;
		}

		/**
		 * \brief Evaluate the surface area heuristic at each bin boundary
		 * and return the maximizer for the given cost constants.
		 */
		SplitCandidate maximizeSAH(Float traversalCost,
				Float intersectionCost, Float emptySpaceBonus) {
			SplitCandidate candidate;
			Float normalization = 2.0f / m_aabb.getSurfaceArea();
			int binIdx = 0, leftBin = 0;
			candidate.planarLeft = false;
			candidate.cost = std::numeric_limits<Float>::infinity();

			for (int axis=0; axis<3; ++axis) {
				Vector extents = m_aabb.getExtents();
				size_type numLeft = 0, numRight = m_primCount;
				Float leftWidth = 0, rightWidth = extents[axis];
				const Float binSize = m_binSize[axis];

				for (int i=0; i<BinCount-1; ++i) {
					numLeft  += m_minBins[binIdx];
					numRight -= m_maxBins[binIdx];
					leftWidth += binSize;
					rightWidth -= binSize;

					extents[axis] = leftWidth;
					Float pLeft = normalization * (extents.x*extents.y 
							+ extents.x*extents.z + extents.y*extents.z);

					extents[axis] = rightWidth;
					Float pRight = normalization * (extents.x*extents.y 
							+ extents.x*extents.z + extents.y*extents.z);

					Float cost = traversalCost + intersectionCost 
						* (pLeft * numLeft + pRight * numRight);

					if (numLeft == 0 || numRight == 0)
						cost *= emptySpaceBonus;

					if (cost < candidate.cost) {
						candidate.cost = cost;
						candidate.axis = axis;
						candidate.numLeft = numLeft;
						candidate.numRight = numRight;
						leftBin = i;
					}

					binIdx++;
				}
				binIdx++;
			}
			
			Assert(candidate.cost != std::numeric_limits<Float>::infinity());

			const int axis = candidate.axis;
			const float min = m_aabb.min[axis];

			/* This part is ensures that the returned split plane is consistent
			 * with the floating point calculations done by the binning code 
			 * in \ref bin(). Since reciprocals and various floating point 
			 * roundoff errors are involved, simply setting
			 *
			 * candidate.pos = m_aabb.min[axis] + (leftBin+1) * m_binSize[axis];
			 *
			 * will potentially lead to a different number primitives being
			 * classified to the left and right compared to the numbers stored
			 * in candidate.numLeft and candidate.numRight. We can't have that,
			 * however, since the partitioning code assumes that these 
			 * numbers are correct. This removes the need for an extra sweep
			 * through the whole primitive list.
			 */
			float invBinSize = 1/m_binSize[axis],
				  split = min + (leftBin + 1) * m_binSize[axis];
			float splitNext = nextafterf(split, 
				  std::numeric_limits<float>::max());
			int idx     = (int) ((split - min) * invBinSize);
			int idxNext = (int) ((splitNext - min) * invBinSize);

			/**
			 * The split plane should be along the last discrete floating
			 * floating position, which would still be classified into
			 * the left bin.
			 */
			if (!(idx == leftBin && idxNext > leftBin)) {
				float direction;
				
				/* First, determine the search direction */
				if (idx > leftBin)
					direction = -std::numeric_limits<float>::max();
				else
					direction = std::numeric_limits<float>::max();

				while (true) {
					split     = nextafterf(split, direction);
					splitNext = nextafterf(split, 
							std::numeric_limits<float>::max());
					idx     = (int) ((split - min) * invBinSize);
					idxNext = (int) ((splitNext - min) * invBinSize);
					if (idx == leftBin && idxNext > leftBin)
						break;
				}
			}

			candidate.pos = split;

			return candidate;
		}

		/**
		 * \brief Given a suitable split candiate, compute tight bounding
		 * boxes for the left and right subtrees and return associated
		 * primitive lists.
		 */
		template <typename Functor> void partition(
				const Functor &functor, SplitCandidate &split, 
				AABB &left, AABB &right, Float traversalCost,
				Float intersectionCost) {
			const float splitPos = split.pos;
			const int axis = split.axis;
			int numLeft = 0, numRight = 0;
			
			Assert(m_primCount == functor.getPrimitiveCount());

			for (size_type i=0; i<m_primCount; ++i) {
				const AABB aabb = functor.getAABB(i);

				if (aabb.max[axis] <= splitPos) {
					left.expandBy(aabb);
					numLeft++;
				} else if (aabb.min[axis] >= splitPos) {
					right.expandBy(aabb);
					numRight++;
				} else {
					left.expandBy(aabb);
					right.expandBy(aabb);
					numLeft++;
					numRight++;
				}
			}

			Assert(numLeft == split.numLeft);
			Assert(numRight == split.numRight);

			left.max[axis] = std::min(left.max[axis], (Float) splitPos);
			right.min[axis] = std::max(right.min[axis], (Float) splitPos);

			if (left.max[axis] != right.min[axis]) {
				/* There is some space between the child nodes -- move the
				   split plane onto one of the AABBs to minimize the SAH */
				Float normalization = 2.0f / m_aabb.getSurfaceArea();
				Vector extents = m_aabb.getExtents();

				cout << "There is some space (split=" << splitPos << ") :" << endl;
				cout << m_aabb.toString() << endl;
				cout << left.toString() << endl;
				cout << right.toString() << endl;

				extents[axis] = left.max[axis] - m_aabb.min[axis];
				Float pLeft1 = normalization * (extents.x*extents.y 
							+ extents.x*extents.z + extents.y*extents.z);
				extents[axis] = m_aabb.max[axis] - left.max[axis];
				Float pRight1 = normalization * (extents.x*extents.y 
							+ extents.x*extents.z + extents.y*extents.z);
				Float sahCost1 = traversalCost + intersectionCost 
					* (pLeft1 * numLeft + pRight1 * numRight);

				extents[axis] = right.min[axis] - m_aabb.min[axis];
				Float pLeft2 = normalization * (extents.x*extents.y 
							+ extents.x*extents.z + extents.y*extents.z);
				extents[axis] = m_aabb.max[axis] - right.min[axis];
				Float pRight2 = normalization * (extents.x*extents.y 
							+ extents.x*extents.z + extents.y*extents.z);
				Float sahCost2 = traversalCost + intersectionCost 
					* (pLeft2 * numLeft + pRight2 * numRight);

				if (sahCost1 <= sahCost2) {
					cout << "Adjusting cost from " << split.cost << " to " << sahCost1 << endl;
					split.cost = sahCost1;
					split.pos = left.max[axis];
				} else {
					cout << "Adjusting cost from " << split.cost << " to " << sahCost2 << endl;
					split.cost = sahCost2;
					split.pos = right.min[axis];
				}
			}
		}
	private:
		size_type m_minBins[3*BinCount], m_maxBins[3*BinCount];
		size_type m_primCount;
		AABB m_aabb;
		Vector m_binSize;
	};
private:
	Float m_traversalCost, m_intersectionCost, m_emptySpaceBonus;
	AABB m_aabb;
	int m_maxDepth;
};

MTS_NAMESPACE_END

#endif /* __KDTREE_GENERIC_H */

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
#include <boost/tuple/tuple.hpp>

#define MTS_KD_MAX_DEPTH 48   ///< Compile-time KD-tree depth limit
#define MTS_KD_STATISTICS 1   ///< Collect statistics during building/traversal 
#define MTS_KD_MINMAX_BINS 32 ///< Min-max bins resolution
#define MTS_KD_MINMAX_DEPTH 4 ///< Use min-max binning for the first 4 levels
#define MTS_KD_MIN_ALLOC 128  ///< Allocate memory in 128 KB chunks

MTS_NAMESPACE_BEGIN

/**
 * \brief Special "ordered" memory allocator
 *
 * During kd-tree construction, large amounts of memory are required 
 * to temporarily hold index and edge event lists. When not implemented
 * properly, these allocations can become a critical bottleneck.
 * The class \ref OrderedChunkAllocator provides a specialized
 * memory allocator, which reserves memory in chunks of at least
 * 128KiB. An important assumption made by the allocator is that
 * memory will be released in the same order, in which it is 
 * allocated. This makes it possible to create an implementation
 * with a very low memory overhead. Note that no locking is done, 
 * hence each thread will need its own allocator.
 */
class OrderedChunkAllocator {
public:
	inline OrderedChunkAllocator(const std::string &name) : m_name(name) {
		m_chunks.reserve(16);
	}

	/**
	 * \brief Release all memory used by the allocator
	 */
	inline void cleanup() {
		for (std::vector<Chunk>::iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it)
			freeAligned((*it).start);
	}

	/**
	 * \brief Request a block of memory from the allocator
	 *
	 * Walks through the list of chunks to find one with enough
	 * free memory. If no chunk could be found, a new one is created.
	 */
	template <typename T> T *allocate(size_t size) {
		cout << m_name << ":" << " alloc(" << size << ")" << endl;

		for (std::vector<Chunk>::iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it) {
			Chunk &chunk = *it;
			if (chunk.getRemainder() >= size) {
				T* result = reinterpret_cast<T *>(chunk.cur);
				chunk.cur += size;
				return result;
			}
		}

		/* No chunk had enough free memory */
		size_t allocSize = std::max(size, 
			(size_t) (MTS_KD_MIN_ALLOC * 1024));

		Chunk chunk;
		chunk.start = (uint8_t *) allocAligned(allocSize);
		chunk.cur = chunk.start + size;
		chunk.size = allocSize;
		m_chunks.push_back(chunk);

		return reinterpret_cast<T *>(chunk.start);
	}

	template <typename T> void release(T *ptr) {
		for (std::vector<Chunk>::iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it) {
			Chunk &chunk = *it;
			if ((uint8_t *) ptr >= chunk.start && 
				(uint8_t *) ptr < chunk.start + chunk.size) {
				cout << m_name << ":" << " release(" << (chunk.cur - (uint8_t *) ptr) << ")" << endl;
				chunk.cur = (uint8_t *) ptr;
				return;
			}
		}
		SLog(EError, "OrderedChunkAllocator: Internal error while"
			" releasing memory");
	}

	/**
	 * \brief Shrink the size of the last allocated chunk
	 */
	template <typename T> void shrinkAllocation(T *ptr, size_t newSize) {
		for (std::vector<Chunk>::iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it) {
			Chunk &chunk = *it;
			if ((uint8_t *) ptr >= chunk.start &&
				(uint8_t *) ptr < chunk.start + chunk.size) {
				cout << m_name << ":" << " shrink(" << (chunk.cur - (uint8_t *) ptr) << " -> " << newSize << ")" << endl;
				chunk.cur = (uint8_t *) ptr + newSize;
				return;
			}
		}
		SLog(EError, "OrderedChunkAllocator: Internal error in shrinkLast");
	}

	inline size_t getChunkCount() const { return m_chunks.size(); }

	/**
	 * \brief Return the total amount of chunk memory in bytes
	 */
	inline size_t getSize() const {
		size_t result = 0;
		for (std::vector<Chunk>::const_iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it)
			result += (*it).size;
		return result;
	}
	
	/**
	 * \brief Return the total amount of used memory in bytes
	 */
	inline size_t getUsed() const {
		size_t result = 0;
		for (std::vector<Chunk>::const_iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it)
			result += (*it).getUsed();
		return result;
	}

private:
	struct Chunk {
		size_t size;
		uint8_t *start, *cur;

		inline size_t getUsed() const {
			return cur - start;
		}

		inline size_t getRemainder() const {
			return size - getUsed();
		}
	};

	std::string m_name;
	std::vector<Chunk> m_chunks;
};

/**
 * \brief Generic kd-tree for data structure used to accelerate 
 * ray intersections against large amounts of three-dimensional
 * shapes.
 */
template <typename Derived> class GenericKDTree : public Object {
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

	struct BuildContext {
		OrderedChunkAllocator leftAlloc, rightAlloc;
		OrderedChunkAllocator tempAlloc, nodeAlloc;

		inline BuildContext() : leftAlloc("left"), rightAlloc("right"),
			tempAlloc("temp"), nodeAlloc("node") { }
	};

	GenericKDTree() {
		m_traversalCost = 15;
		m_intersectionCost = 20;
		m_emptySpaceBonus = 0.8f;
		m_clip = true;
		m_stopPrims = 4;
	}

	/// Cast to the derived class
	inline Derived *downCast() {
		return static_cast<Derived *>(this);
	}

	/// Cast to the derived class (const version)
	inline const Derived *downCast() const {
		return static_cast<Derived *>(this);
	}

	/**
	 * \brief Build a KD-tree over the supplied axis-aligned bounding
	 * boxes.
	 */
	void build() {
		BuildContext ctx;
		size_type primCount = downCast()->getPrimitiveCount();

		/* Establish an ad-hoc depth cutoff value (Formula from PBRT) */
		m_maxDepth = std::min((int) (8 + 1.3f * log2i(primCount)), MTS_KD_MAX_DEPTH);

		Log(EDebug, "Creating a preliminary index list (%.2f KiB)", 
				primCount * sizeof(index_type) / 1024.0f);

		OrderedChunkAllocator &leftAlloc = ctx.leftAlloc;
		index_type *indices = leftAlloc.allocate<index_type>(
				primCount * sizeof(index_type));

		ref<Timer> timer = new Timer();
		m_aabb.reset();
		for (index_type i=0; i<primCount; ++i) {
			m_aabb.expandBy(downCast()->getAABB(i));
			indices[i] = i;
		}

		Log(EDebug, "Computed scene bounds in %i ms", timer->getMilliseconds());

		Log(EDebug, "kd-tree configuration:");
		Log(EDebug, "   Traversal cost           : %.2f", m_traversalCost);
		Log(EDebug, "   Intersection cost        : %.2f", m_intersectionCost);
		Log(EDebug, "   Empty space bonus        : %.2f", m_emptySpaceBonus);
		Log(EDebug, "   Max. tree depth          : %i", m_maxDepth);
		Log(EDebug, "   Stopping primitive count : %i", m_stopPrims);
		Log(EDebug, "   Scene bounding box (min) : %s", m_aabb.min.toString().c_str());
		Log(EDebug, "   Scene bounding box (max) : %s", m_aabb.max.toString().c_str());
		Log(EDebug, "   Min-max bins             : %i", MTS_KD_MINMAX_BINS);
		Log(EDebug, "   Exact SAH eval. depth    : %i", MTS_KD_MINMAX_DEPTH);
		Log(EDebug, "   Perfect splits           : %s", m_clip ? "yes" : "no");

		Log(EInfo, "Constructing a SAH kd-tree (%i primitives) ..", primCount); 

		buildTree(ctx, 1, m_aabb, indices, primCount, true);
		ctx.leftAlloc.release(indices);

		Assert(ctx.leftAlloc.getUsed() == 0);
		Assert(ctx.rightAlloc.getUsed() == 0);
		Assert(ctx.tempAlloc.getUsed() == 0);

		Log(EInfo, "Finished. Construction took %i ms.", timer->getMilliseconds());

		Log(EDebug, "Chunk allocator statistics");
		Log(EDebug, "   Left:  " SIZE_T_FMT " chunks (%.2f KiB)",
				ctx.leftAlloc.getChunkCount(), ctx.leftAlloc.getSize() / 1024.0f);
		Log(EDebug, "   Right: " SIZE_T_FMT " chunks (%.2f KiB)",
				ctx.rightAlloc.getChunkCount(), ctx.rightAlloc.getSize() / 1024.0f);
		Log(EDebug, "   Temp:  " SIZE_T_FMT " chunks (%.2f KiB)",
				ctx.tempAlloc.getChunkCount(), ctx.tempAlloc.getSize() / 1024.0f);
		Log(EDebug, "   Nodes: " SIZE_T_FMT " chunks (%.2f KiB)",
				ctx.nodeAlloc.getChunkCount(), ctx.nodeAlloc.getSize() / 1024.0f);

		ctx.leftAlloc.cleanup();
		ctx.rightAlloc.cleanup();
		ctx.tempAlloc.cleanup();
	}

	void buildTree(BuildContext &ctx, unsigned int depth, const AABB &aabb, 
			index_type *indices, size_type primCount, bool isLeftChild) {
		if (primCount <= m_stopPrims || depth >= m_maxDepth) 
			return;

		MinMaxBins<MTS_KD_MINMAX_BINS> bins(aabb);
		bins.bin(downCast(), indices, primCount);
		SplitCandidate candidate = bins.maximizeSAH(m_traversalCost,
			m_intersectionCost, m_emptySpaceBonus);

		const Float originalCost = primCount * m_intersectionCost;

		cout << "originalCost = " << originalCost << endl;

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

		boost::tuple<AABB, index_type *, AABB, index_type *> partition = 
			bins.partition(ctx, downCast(), indices, candidate, isLeftChild, 
			m_traversalCost, m_intersectionCost);

		buildTree(ctx, depth+1, boost::get<0>(partition), boost::get<1>(partition),
				candidate.numLeft, true);
		buildTree(ctx, depth+1, boost::get<2>(partition), boost::get<3>(partition),
				candidate.numLeft, false);

		ctx.rightAlloc.release(boost::get<3>(partition));
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

	/// Edge event comparison functor
	struct EdgeEventOrdering : public std::binary_function<EdgeEvent, EdgeEvent, bool> {
		inline bool operator()(const EdgeEvent &a, const EdgeEvent &b) const {
			if (a.t != b.t)
				return a.t < b.t;
			return a.type < b.type;
		}
	};

	/**
	 * \brief KD-tree node in 8 bytes. 
	 */
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

	/**
	 * \brief Min-max binning as described in
	 * "Highly Parallel Fast KD-tree Construction for Interactive
	 *  Ray Tracing of Dynamic Scenes"
	 * by M. Shevtsov, A. Soupikov and A. Kapustin
	 *
	 * @tparam BinCount Number of bins to be allocated
	 */
	template <int BinCount> struct MinMaxBins {
		MinMaxBins(const AABB &aabb) : m_aabb(aabb) {
			m_binSize = m_aabb.getExtents() / BinCount;
		}

		/**
		 * \brief Run min-max binning
		 *
		 * \param derived Derived class to be used to determine the AABB for
		 *     a given list of primitives
		 * \param indices Primitive indirection list
		 * \param primCount Specifies the length of \a indices
		 */
		void bin(const Derived *derived, index_type *indices, size_type primCount) {
			m_primCount = primCount;
			memset(m_minBins, 0, sizeof(size_type) * 3 * BinCount);
			memset(m_maxBins, 0, sizeof(size_type) * 3 * BinCount);
			Vector invBinSize;

			for (int axis=0; axis<3; ++axis) 
				invBinSize[axis] = 1/m_binSize[axis];

			for (size_type i=0; i<m_primCount; ++i) {
				const AABB aabb = derived->getAABB(indices[i]);
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
		boost::tuple<AABB, index_type *, AABB, index_type *> partition(
				BuildContext &ctx, const Derived *derived, index_type *primIndices,
				SplitCandidate &split, bool isLeftChild, Float traversalCost, Float intersectionCost) {
			const float splitPos = split.pos;
			const int axis = split.axis;
			int numLeft = 0, numRight = 0;
			AABB leftBounds, rightBounds;

			index_type *leftIndices, *rightIndices;
			if (isLeftChild) {
				OrderedChunkAllocator &rightAlloc = ctx.rightAlloc;
				leftIndices = primIndices;
				rightIndices = rightAlloc.allocate<index_type>(
					sizeof(index_type) * split.numRight);
			} else {
				OrderedChunkAllocator &leftAlloc = ctx.leftAlloc;
				leftIndices = leftAlloc.allocate<index_type>(
					sizeof(index_type) * split.numLeft);
				rightIndices = primIndices;
			}

			for (size_type i=0; i<m_primCount; ++i) {
				const index_type primIndex = primIndices[i];
				const AABB aabb = derived->getAABB(primIndex);

				if (aabb.max[axis] <= splitPos) {
					Assert(numLeft < split.numLeft);
					leftBounds.expandBy(aabb);
					leftIndices[numLeft++] = primIndex;
				} else if (aabb.min[axis] >= splitPos) {
					Assert(numRight < split.numRight);
					rightBounds.expandBy(aabb);
					rightIndices[numRight++] = primIndex;
				} else {
					leftBounds.expandBy(aabb);
					rightBounds.expandBy(aabb);
					Assert(numLeft < split.numLeft);
					Assert(numRight < split.numRight);
					leftIndices[numLeft++] = primIndex;
					rightIndices[numRight++] = primIndex;
				}
			}

			Assert(numLeft == split.numLeft);
			Assert(numRight == split.numRight);

			/// Release the unused memory regions
			if (isLeftChild)
				ctx.leftAlloc.shrinkAllocation(leftIndices, 
						split.numLeft * sizeof(index_type));
			else
				ctx.rightAlloc.shrinkAllocation(rightIndices, 
						split.numRight * sizeof(index_type));

			leftBounds.max[axis] = std::min(leftBounds.max[axis], (Float) splitPos);
			rightBounds.min[axis] = std::max(rightBounds.min[axis], (Float) splitPos);

			if (leftBounds.max[axis] != rightBounds.min[axis] && false) {
				/* There is some space between the child nodes -- move
				   the split plane onto one of the AABBs so that the
				   surface area heuristic is minimized */
				Float normalization = 2.0f / m_aabb.getSurfaceArea();
				Vector extents = m_aabb.getExtents();

				cout << "There is some space (split=" << splitPos << ") :" << endl;
				cout << m_aabb.toString() << endl;
				cout << leftBounds.toString() << endl;
				cout << rightBounds.toString() << endl;

				extents[axis] = leftBounds.max[axis] - m_aabb.min[axis];
				Float pLeft1 = normalization * (extents.x*extents.y 
							+ extents.x*extents.z + extents.y*extents.z);
				extents[axis] = m_aabb.max[axis] - leftBounds.max[axis];
				Float pRight1 = normalization * (extents.x*extents.y 
							+ extents.x*extents.z + extents.y*extents.z);
				Float sahCost1 = traversalCost + intersectionCost 
					* (pLeft1 * numLeft + pRight1 * numRight);

				extents[axis] = rightBounds.min[axis] - m_aabb.min[axis];
				Float pLeft2 = normalization * (extents.x*extents.y 
							+ extents.x*extents.z + extents.y*extents.z);
				extents[axis] = m_aabb.max[axis] - rightBounds.min[axis];
				Float pRight2 = normalization * (extents.x*extents.y 
							+ extents.x*extents.z + extents.y*extents.z);
				Float sahCost2 = traversalCost + intersectionCost 
					* (pLeft2 * numLeft + pRight2 * numRight);

				if (sahCost1 <= sahCost2) {
					cout << "Adjusting cost from " << split.cost << " to " << sahCost1 << endl;
					split.cost = sahCost1;
					split.pos = leftBounds.max[axis];
				} else {
					cout << "Adjusting cost from " << split.cost << " to " << sahCost2 << endl;
					split.cost = sahCost2;
					split.pos = rightBounds.min[axis];
				}
			}

			return boost::make_tuple(leftBounds, leftIndices,
					rightBounds, rightIndices);
		}
	private:
		size_type m_minBins[3*BinCount], m_maxBins[3*BinCount];
		size_type m_primCount;
		AABB m_aabb;
		Vector m_binSize;
	};

private:
	Float m_traversalCost;
	Float m_intersectionCost;
	Float m_emptySpaceBonus;
	bool m_clip;
	AABB m_aabb;
	size_type m_maxDepth;
	size_type m_stopPrims;
};

MTS_NAMESPACE_END

#endif /* __KDTREE_GENERIC_H */

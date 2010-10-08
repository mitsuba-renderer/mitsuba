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
#define MTS_KD_MINMAX_BINS 32 ///< Min-max bin count
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
	inline OrderedChunkAllocator() {
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
		size *= sizeof(T);
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
		newSize *= sizeof(T);
		for (std::vector<Chunk>::iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it) {
			Chunk &chunk = *it;
			if ((uint8_t *) ptr >= chunk.start &&
				(uint8_t *) ptr < chunk.start + chunk.size) {
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

	std::vector<Chunk> m_chunks;
};

/**
 * \brief Generic kd-tree for data structure used to accelerate 
 * ray intersections against large amounts of three-dimensional
 * shapes.
 */
template <typename Derived> class GenericKDTree : public Object {
protected:
	struct KDNode;
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
	};

	GenericKDTree() : m_root(NULL) {
		m_traversalCost = 15;
		m_intersectionCost = 20;
		m_emptySpaceBonus = 0.8f;
		m_clip = true;
		m_stopPrims = 2;
		m_badSplits = 0;
		m_leafNodeCount = 0;
		m_innerNodeCount = 0;
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
	 * \brief Build a KD-tree over supplied geometry
	 */
	void build() {
		if (m_root != NULL)
			Log(EError, "KD-tree has already been built!");

		BuildContext ctx;
		size_type primCount = downCast()->getPrimitiveCount();

		/* Establish an ad-hoc depth cutoff value (Formula from PBRT) */
		m_maxDepth = std::min((int) (8 + 1.3f * log2i(primCount)), MTS_KD_MAX_DEPTH);

		Log(EDebug, "Creating a preliminary index list (%.2f KiB)", 
			primCount * sizeof(index_type) / 1024.0f);

		OrderedChunkAllocator &leftAlloc = ctx.leftAlloc, &nodeAlloc = ctx.nodeAlloc;
		index_type *indices = leftAlloc.allocate<index_type>(primCount);

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

		m_root = nodeAlloc.allocate<KDNode>(1);
		Float finalSAHCost
			= buildTree(ctx, 1, m_root, m_aabb, m_aabb, indices, primCount, true);

		ctx.leftAlloc.release(indices);

		Assert(ctx.leftAlloc.getUsed() == 0);
		Assert(ctx.rightAlloc.getUsed() == 0);
		Assert(ctx.tempAlloc.getUsed() == 0);

		Log(EInfo, "Finished -- took %i ms.", timer->getMilliseconds());

		Log(EDebug, "Chunk allocator statistics");
		Log(EDebug, "   Left:  " SIZE_T_FMT " chunks (%.2f KiB)",
				ctx.leftAlloc.getChunkCount(), ctx.leftAlloc.getSize() / 1024.0f);
		Log(EDebug, "   Right: " SIZE_T_FMT " chunks (%.2f KiB)",
				ctx.rightAlloc.getChunkCount(), ctx.rightAlloc.getSize() / 1024.0f);
		Log(EDebug, "   Temp:  " SIZE_T_FMT " chunks (%.2f KiB)",
				ctx.tempAlloc.getChunkCount(), ctx.tempAlloc.getSize() / 1024.0f);
		Log(EDebug, "   Nodes: " SIZE_T_FMT " chunks (%.2f KiB)",
				ctx.nodeAlloc.getChunkCount(), ctx.nodeAlloc.getSize() / 1024.0f);
		Log(EDebug, "Detailed kd-tree statistics:");
		Log(EDebug, "   Final SAH cost    : %.2f", finalSAHCost);
		Log(EDebug, "   # of Leaf nodes   : %i", m_leafNodeCount);
		Log(EDebug, "   # of Inner nodes  : %i", m_innerNodeCount);
		Log(EDebug, "   # of bad splits   : %i", m_badSplits);
		Log(EDebug, "   Indirection table : " SIZE_T_FMT " entries",
				m_indirectionTable.size());

		ctx.leftAlloc.cleanup();
		ctx.rightAlloc.cleanup();
		ctx.tempAlloc.cleanup();
	}

	/**
	 * \brief Leaf node creation helper function
	 *
	 * \param ctx 
	 *     Thread-specific build context containing allocators etc.
	 * \param node
	 *     KD-tree node entry to be filled
	 * \param primCount
	 *     Total primitive count for the current node
	 *
	 * \returns 
	 */
	Float createLeaf(BuildContext &ctx, KDNode *node, size_type primCount) {
		node->initLeafNode(0, primCount);
		m_leafNodeCount++;
		return m_intersectionCost * primCount;
	}

	/**
	 * \brief Build helper function
	 *
	 * \param ctx 
	 *     Thread-specific build context containing allocators etc.
	 * \param depth 
	 *     Current tree depth (1 == root node)
	 * \param node
	 *     KD-tree node entry to be filled
	 * \param nodeAABB
	 *     Axis-aligned bounding box of the current node
	 * \param tightAABB
	 *     Tight bounding box of the contained geometry (for min-max binning)
	 * \param indices
	 *     Index list of all triangles in the current node (for min-max binning)
	 * \param primCount
	 *     Total primitive count for the current node
	 * \param isLeftChild
	 *     Is this node the left child of its parent? This is important for
	 *     memory management.
	 * \returns 
	 *     A tuple specifying the final SAH cost and a pointer to the node.
	 */
	Float buildTree(BuildContext &ctx, unsigned int depth, KDNode *node, 
			const AABB &nodeAABB, const AABB &tightAABB, index_type *indices,
			size_type primCount, bool isLeftChild) {
		Assert(nodeAABB.contains(tightAABB));

		if (primCount <= m_stopPrims || depth >= m_maxDepth) {
			return createLeaf(ctx, node, primCount);
		}

		MinMaxBins<MTS_KD_MINMAX_BINS> bins(tightAABB);
		bins.bin(downCast(), indices, primCount);
		SplitCandidate split = bins.maximizeSAH(m_traversalCost,
			m_intersectionCost, m_emptySpaceBonus);

		boost::tuple<AABB, index_type *, AABB, index_type *> partition = 
			bins.partition(ctx, downCast(), indices, split, isLeftChild, 
			m_traversalCost, m_intersectionCost);

		OrderedChunkAllocator &nodeAlloc = ctx.nodeAlloc;
		KDNode *children = nodeAlloc.allocate<KDNode>(2);

		size_t initialIndirectionTableSize = m_indirectionTable.size();
		if (!node->initInnerNode(split.axis, split.pos, children-node)) {
			/* Unable to store relative offset -- create an indirection
			   table entry */
			node->initIndirectionNode(split.axis, split.pos, 
					initialIndirectionTableSize);
			m_indirectionTable.push_back(children);
		}
		m_innerNodeCount++;

		AABB childAABB(nodeAABB);
		childAABB.max[split.axis] = split.pos;
		Float saLeft = childAABB.getSurfaceArea();

		Float leftSAHCost = buildTree(ctx, depth+1, children,
				childAABB, boost::get<0>(partition), boost::get<1>(partition), 
				split.numLeft, true);

		childAABB.min[split.axis] = split.pos;
		childAABB.max[split.axis] = nodeAABB.max[split.axis];
		Float saRight = childAABB.getSurfaceArea();

		Float rightSAHCost = buildTree(ctx, depth+1, children + 1,
				childAABB, boost::get<2>(partition), boost::get<3>(partition), 
				split.numRight, false);

		/* Compute the final SAH cost given the updated cost 
		   values received from the children */
		Float finalSAHCost = m_traversalCost + 
			(saLeft * leftSAHCost + saRight * rightSAHCost)
			/ nodeAABB.getSurfaceArea();

		/* Release the index lists not needed by the children anymore */
		if (isLeftChild)
			ctx.rightAlloc.release(boost::get<3>(partition));
		else
			ctx.leftAlloc.release(boost::get<1>(partition));	

		if (finalSAHCost < primCount * m_intersectionCost) {
			return finalSAHCost;
		} else {
			/* In the end, splitting didn't help to reduce the SAH cost.
			   Tear up everything below this node and create a leaf */
			if (m_indirectionTable.size() > initialIndirectionTableSize)
				m_indirectionTable.resize(initialIndirectionTableSize);
			tearUp(ctx, node);
			m_badSplits++;
			return createLeaf(ctx, node, primCount);
		}
	}

	/**
	 * \brief Tear up a subtree after a split did not reduce the SAH cost
	 */
	void tearUp(BuildContext &ctx, KDNode *node) {
		if (node->isLeaf()) {
			/// XXX Create primitive list for leaf
		} else {
			KDNode *left;
			if (EXPECT_TAKEN(!node->isIndirection()))
				left = node->getLeft();
			else
				left = m_indirectionTable[node->getIndirectionIndex()];

			tearUp(ctx, left);
			tearUp(ctx, left+1);

			ctx.nodeAlloc.release(left);
		}
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
				   30   : Indirection node flag
				   29-3 : Offset to the left child
				   2-0  : Split axis
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
			EIndirectionMask = 1 << 30,
			ELeafOffsetMask = ~ETypeMask,
			EInnerAxisMask = 0x3,
			EInnerOffsetMask = ~(EInnerAxisMask + EIndirectionMask),
			ERelOffsetLimit = (1<<28) - 1
		};

		/// Initialize a leaf kd-Tree node
		inline void initLeafNode(unsigned int offset, unsigned int numPrims) {
			leaf.combined = ETypeMask | offset;
			leaf.end = offset + numPrims;
		}

		/**
		 * Initialize an interior kd-Tree node. Reports a failure if the
		 * relative offset to the left child node is too large.
		 */
		inline bool initInnerNode(int axis, float split, ptrdiff_t relOffset) {
			if (relOffset < 0 || relOffset > ERelOffsetLimit)
				return false;
			inner.combined = axis | ((uint32_t) relOffset << 2);
			inner.split = split;
			return true;
		}

		/**
		 * \brief Initialize an interior indirection node.
		 *
		 * Indirections are necessary whenever the children cannot be 
		 * referenced using a relative pointer, which can happen when 
		 * they lie in different memory chunks. In this case, the node
		 * stores an index into a globally shared pointer list.
		 */
		inline void initIndirectionNode(int axis, float split, uint32_t indirectionEntry) {
			inner.combined = EIndirectionMask | axis | ((uint32_t) indirectionEntry << 2);
			inner.split = split;
		}

		/// Is this a leaf node?
		FINLINE bool isLeaf() const {
			return leaf.combined & ETypeMask;
		}

		/// Is this an indirection node?
		FINLINE bool isIndirection() const {
			return leaf.combined & EIndirectionMask;
		}

		/// Assuming this is a leaf node, return the first primitive index
		FINLINE index_type getPrimStart() const {
			return leaf.combined & ELeafOffsetMask;
		}

		/// Assuming this is a leaf node, return the last primitive index
		FINLINE index_type getPrimEnd() const {
			return leaf.end;
		}
				
		/// Return the index of an indirection node
		FINLINE index_type getIndirectionIndex() const {
			return(inner.combined & EInnerOffsetMask) >> 2;
		}

		/// Return the left child (assuming that this is an interior node)
		FINLINE const KDNode * __restrict getLeft() const {
			return this + 
				((inner.combined & EInnerOffsetMask) >> 2);
		}

		/// Return the left child (assuming that this is an interior node)
		FINLINE KDNode * __restrict getLeft() {
			return this + 
				((inner.combined & EInnerOffsetMask) >> 2);
		}

		/// Return the left child (assuming that this is an interior node)
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
					m_maxBins[axis * BinCount + std::max(0, std::min(maxIdx, BinCount-1))]++;
					m_minBins[axis * BinCount + std::max(0, std::min(minIdx, BinCount-1))]++;
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

			Assert(split > min && split < m_aabb.max[axis]);

			if (split <= min || split >= m_aabb.max[axis]) {
				cout << "Ran into some problems!" << endl;
				cout << "AABB = " << m_aabb.toString() << endl;
				cout << "candidate axis = " << candidate.axis << endl;
				cout << "candidate cost = " << candidate.cost << endl;
				cout << "prims = " << m_primCount << endl;
				cout << "invBinSize = " << invBinSize << endl;
				cout << "split pos = " << split << endl;
				cout << "left bin = " << leftBin << endl;
				cout << "idx = " << idx << ", " << idxNext << endl;
				exit(-1);
			}

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
				int it = 0;

				while (true) {
					split     = nextafterf(split, direction);
					splitNext = nextafterf(split, 
							std::numeric_limits<float>::max());
					idx     = (int) ((split - min) * invBinSize);
					idxNext = (int) ((splitNext - min) * invBinSize);
					if (idx == leftBin && idxNext > leftBin)
						break;
					++it;
					if (it % 100 == 0)
						cout << "In iteration " << it << endl;
				}

			}

			Assert(split > m_aabb.min[axis] && split < m_aabb.max[axis]);

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
				rightIndices = rightAlloc.allocate<index_type>(split.numRight);
			} else {
				OrderedChunkAllocator &leftAlloc = ctx.leftAlloc;
				leftIndices = leftAlloc.allocate<index_type>(split.numLeft);
				rightIndices = primIndices;
			}

			for (size_type i=0; i<m_primCount; ++i) {
				const index_type primIndex = primIndices[i];
				const AABB aabb = derived->getAABB(primIndex);

				if (aabb.max[axis] <= splitPos) {
					Assert(numLeft < split.numLeft);
					leftBounds.expandBy(aabb);
					leftIndices[numLeft++] = primIndex;
				} else if (aabb.min[axis] > splitPos) {
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

			leftBounds.clip(m_aabb);
			rightBounds.clip(m_aabb);

			Assert(numLeft == split.numLeft);
			Assert(numRight == split.numRight);

			/// Release the unused memory regions
			if (isLeftChild)
				ctx.leftAlloc.shrinkAllocation(leftIndices, split.numLeft);
			else
				ctx.rightAlloc.shrinkAllocation(rightIndices, split.numRight);

			leftBounds.max[axis] = std::min(leftBounds.max[axis], (Float) splitPos);
			rightBounds.min[axis] = std::max(rightBounds.min[axis], (Float) splitPos);

			if (leftBounds.max[axis] != rightBounds.min[axis]) {
				/* There is some space between the child nodes -- move
				   the split plane onto one of the AABBs so that the
				   surface area heuristic is minimized */
				Float normalization = 2.0f / m_aabb.getSurfaceArea();
				Vector extents = m_aabb.getExtents();

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
					if (sahCost1 > split.cost) {
						cout << "Original splitting cost = " << split.cost << endl;
						cout << "New splitting cost      = " << sahCost1 << endl;
					}
///					Assert(sahCost1 <= split.cost);
					split.cost = sahCost1;
					split.pos = leftBounds.max[axis];
				} else {
					if (sahCost2 > split.cost) {
						cout << "Original splitting cost = " << split.cost << endl;
						cout << "New splitting cost      = " << sahCost2 << endl;
					}
//					Assert(sahCost2 <= split.cost);
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
	KDNode *m_root;
	std::vector<KDNode *> m_indirectionTable;
	Float m_traversalCost;
	Float m_intersectionCost;
	Float m_emptySpaceBonus;
	bool m_clip;
	AABB m_aabb;
	size_type m_maxDepth;
	size_type m_stopPrims;
	size_type m_leafNodeCount;
	size_type m_innerNodeCount;
	size_type m_badSplits;
};

MTS_NAMESPACE_END

#endif /* __KDTREE_GENERIC_H */

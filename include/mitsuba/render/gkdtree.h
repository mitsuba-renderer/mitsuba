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
#if !defined(__MITSUBA_RENDER_GKDTREE_H_)
#define __MITSUBA_RENDER_GKDTREE_H_

#include <mitsuba/core/timer.h>
#include <mitsuba/core/lock.h>
#include <boost/static_assert.hpp>
#include <stack>

#if defined(__LINUX__)
#include <malloc.h>
#endif

/// Activate lots of extra checks
//#define MTS_KD_DEBUG 1

/** Compile-time KD-tree depth limit. Allows to put certain
    data structures on the stack */
#define MTS_KD_MAXDEPTH 48

/// OrderedChunkAllocator: don't create chunks smaller than 512 KiB
#define MTS_KD_MIN_ALLOC 512*1024

/// Allocate nodes & index lists in blocks of 512 KiB
#define MTS_KD_BLOCKSIZE_KD  (512*1024/sizeof(KDNode))
#define MTS_KD_BLOCKSIZE_IDX (512*1024/sizeof(uint32_t))

/**
 * \brief To avoid numerical issues, the size of the scene
 * bounding box is increased by this amount
 */
#define MTS_KD_AABB_EPSILON 1e-3f

#if defined(MTS_KD_DEBUG)
#define KDAssert(expr) SAssert(expr)
#define KDAssertEx(expr, text) SAssertEx(expr, text)
#else
#define KDAssert(expr)
#define KDAssertEx(expr, text)
#endif

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
 * memory will be released in the exact same order, in which it was
 * previously allocated. This makes it possible to create an
 * implementation with a very low memory overhead. Note that no locking
 * is done, hence each thread will need its own allocator.
 *
 * \author Wenzel Jakob
 */
class OrderedChunkAllocator {
public:
	inline OrderedChunkAllocator(size_t minAllocation = MTS_KD_MIN_ALLOC)
			: m_minAllocation(minAllocation) {
		m_chunks.reserve(16);
	}

	~OrderedChunkAllocator() {
		cleanup();
	}

	/**
	 * \brief Release all memory used by the allocator
	 */
	void cleanup() {
		for (std::vector<Chunk>::iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it)
			freeAligned((*it).start);
		m_chunks.clear();
	}

	/**
	 * \brief Merge the chunks of another allocator into this one
	 */
	void merge(const OrderedChunkAllocator &other) {
		m_chunks.reserve(m_chunks.size() + other.m_chunks.size());
		m_chunks.insert(m_chunks.end(), other.m_chunks.begin(),
				other.m_chunks.end());
	}

	/**
	 * \brief Forget about all chunks without actually freeing them.
	 * This is useful when the chunks have been merged into another
	 * allocator.
	 */
	void forget() {
		m_chunks.clear();
	}

	/**
	 * \brief Request a block of memory from the allocator
	 *
	 * Walks through the list of chunks to find one with enough
	 * free memory. If no chunk could be found, a new one is created.
	 */
	template <typename T> T * __restrict allocate(size_t size) {
		size *= sizeof(T);
		for (std::vector<Chunk>::iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it) {
			Chunk &chunk = *it;
			if (chunk.remainder() >= size) {
				T* result = reinterpret_cast<T *>(chunk.cur);
				chunk.cur += size;
				return result;
			}
		}

		/* No chunk had enough free memory */
		size_t allocSize = std::max(size,
			m_minAllocation);

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
#if defined(MTS_KD_DEBUG)
		/* Uh oh, allocation could not be found. Check if it has size==0 */
		for (std::vector<Chunk>::iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it) {
			const Chunk &chunk = *it;
			if ((uint8_t *) ptr == chunk.start + chunk.size)
				return;
		}
		SLog(EError, "OrderedChunkAllocator: Internal error while"
			" releasing memory");
#endif
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
#if defined(MTS_KD_DEBUG)
		/* Uh oh, allocation could not be found. Check if it has size == 0 */
		if (newSize == 0) {
			for (std::vector<Chunk>::iterator it = m_chunks.begin();
					it != m_chunks.end(); ++it) {
				const Chunk &chunk = *it;
				if ((uint8_t *) ptr == chunk.start + chunk.size)
					return;
			}
		}
		SLog(EError, "OrderedChunkAllocator: Internal error while"
			" releasing memory");
#endif
	}

	inline size_t getChunkCount() const { return m_chunks.size(); }

	/**
	 * \brief Return the total amount of chunk memory in bytes
	 */
	size_t size() const {
		size_t result = 0;
		for (std::vector<Chunk>::const_iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it)
			result += (*it).size;
		return result;
	}

	/**
	 * \brief Return the total amount of used memory in bytes
	 */
	size_t used() const {
		size_t result = 0;
		for (std::vector<Chunk>::const_iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it)
			result += (*it).used();
		return result;
	}

	/**
	 * \brief Return a string representation of the chunks
	 */
	std::string toString() const {
		std::ostringstream oss;
		oss << "OrderedChunkAllocator[" << endl;
		for (size_t i=0; i<m_chunks.size(); ++i)
			oss << "    Chunk " << i << ": " << m_chunks[i].toString() << endl;
		oss << "]";
		return oss.str();
	}

private:
	struct Chunk {
		size_t size;
		uint8_t *start, *cur;

		inline size_t used() const {
			return cur - start;
		}

		inline size_t remainder() const {
			return size - used();
		}

		std::string toString() const {
			return formatString("0x%llx-0x%llx (size=" SIZE_T_FMT
				", used=" SIZE_T_FMT ")", start, start+size,
				size, used());
		}
	};

	size_t m_minAllocation;
	std::vector<Chunk> m_chunks;
};

/**
 * \brief Basic vector implementation, which stores all data
 * in a list of fixed-sized blocks.
 *
 * This leads to a more conservative memory usage when the
 * final size of a (possibly very large) growing vector is
 * unknown. Also, frequent reallocations & copies are avoided.
 *
 * \author Wenzel Jakob
 */
template <typename T, size_t BlockSize> class BlockedVector {
public:
	BlockedVector() : m_pos(0) {}

	~BlockedVector() {
		clear();
	}

	/**
	 * \brief Append an element to the end
	 */
	inline void push_back(const T &value) {
		size_t blockIdx = m_pos / BlockSize;
		size_t offset = m_pos % BlockSize;
		if (blockIdx == m_blocks.size())
			m_blocks.push_back(new T[BlockSize]);
		m_blocks[blockIdx][offset] = value;
		m_pos++;
	}

	/**
	 * \brief Allocate a certain number of elements and
	 * return a pointer to the first one.
	 *
	 * The implementation will ensure that they lie
	 * contiguous in memory -- note that this can potentially
	 * create unused elements in the previous block if a new
	 * one has to be allocated.
	 */
	inline T * __restrict allocate(size_t size) {
#if defined(MTS_KD_DEBUG)
		SAssert(size <= BlockSize);
#endif
		size_t blockIdx = m_pos / BlockSize;
		size_t offset = m_pos % BlockSize;
		T *result;
		if (EXPECT_TAKEN(offset + size <= BlockSize)) {
			if (blockIdx == m_blocks.size())
				m_blocks.push_back(new T[BlockSize]);
			result = m_blocks[blockIdx] + offset;
			m_pos += size;
		} else {
			++blockIdx;
			if (blockIdx == m_blocks.size())
				m_blocks.push_back(new T[BlockSize]);
			result = m_blocks[blockIdx];
			m_pos += BlockSize - offset + size;
		}
		return result;
	}

	inline T &operator[](size_t index) {
		return *(m_blocks[index / BlockSize] +
			(index % BlockSize));
	}

	inline const T &operator[](size_t index) const {
		return *(m_blocks[index / BlockSize] +
			(index % BlockSize));
	}


	/**
	 * \brief Return the currently used number of items
	 */
	inline size_t size() const {
		return m_pos;
	}

	/**
	 * \brief Return the number of allocated blocks
	 */
	inline size_t blockCount() const {
		return m_blocks.size();
	}

	/**
	 * \brief Return the total capacity
	 */
	inline size_t capacity() const {
		return m_blocks.size() * BlockSize;
	}

	/**
	 * \brief Resize the vector to the given size.
	 *
	 * Note: this implementation doesn't support
	 * enlarging the vector and simply changes the
	 * last item pointer.
	 */
	inline void resize(size_t pos) {
#if defined(MTS_KD_DEBUG)
		SAssert(pos <= capacity());
#endif
		m_pos = pos;
	}

	/**
	 * \brief Release all memory
	 */
	void clear() {
		for (typename std::vector<T *>::iterator it = m_blocks.begin();
				it != m_blocks.end(); ++it)
			delete[] *it;
		m_blocks.clear();
		m_pos = 0;
	}
private:
	std::vector<T *> m_blocks;
	size_t m_pos;
};

/**
 * \brief Compact storage for primitive classifcation
 *
 * When classifying primitives with respect to a split plane,
 * a data structure is needed to hold the tertiary result of
 * this operation. This class implements a compact storage
 * (2 bits per entry) in the spirit of the std::vector<bool>
 * specialization.
 *
 * \author Wenzel Jakob
 */
class ClassificationStorage {
public:
	ClassificationStorage()
		: m_buffer(NULL), m_bufferSize(0) { }

	~ClassificationStorage() {
		if (m_buffer)
			delete[] m_buffer;
	}

	void setPrimitiveCount(size_t size) {
		if (m_buffer)
			delete[] m_buffer;
		if (size > 0) {
			m_bufferSize = size/4 + ((size % 4) > 0 ? 1 : 0);
			m_buffer = new uint8_t[m_bufferSize];
		} else {
			m_buffer = NULL;
		}
	}

	inline void set(uint32_t index, int value) {
		uint8_t *ptr = m_buffer + (index >> 2);
		uint8_t shift = (index & 3) << 1;
		*ptr = (*ptr & ~(3 << shift)) | (value << shift);
	}

	inline int get(uint32_t index) const {
		uint8_t *ptr = m_buffer + (index >> 2);
		uint8_t shift = (index & 3) << 1;
		return (*ptr >> shift) & 3;
	}

	inline size_t size() const {
		return m_bufferSize;
	}
private:
	uint8_t *m_buffer;
	size_t m_bufferSize;
};

/**
 * \brief Base class of all kd-trees.
 *
 * This class defines the byte layout for KD-tree nodes and
 * provides methods for querying the tree structure.
 */
template <typename AABBType> class KDTreeBase : public Object {
public:
	/// Index number format (max 2^32 prims)
	typedef uint32_t IndexType;

	/// Size number format
	typedef uint32_t SizeType;

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
				          or indirection table entry
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
			leaf.combined = (uint32_t) ETypeMask | offset;
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
		inline void initIndirectionNode(int axis, float split,
				uint32_t indirectionEntry) {
			inner.combined = EIndirectionMask
				| ((uint32_t) indirectionEntry << 2)
				| axis;
			inner.split = split;
		}

		/// Is this a leaf node?
		FINLINE bool isLeaf() const {
			return leaf.combined & (uint32_t) ETypeMask;
		}

		/// Is this an indirection node?
		FINLINE bool isIndirection() const {
			return leaf.combined & (uint32_t) EIndirectionMask;
		}

		/// Assuming this is a leaf node, return the first primitive index
		FINLINE IndexType getPrimStart() const {
			return leaf.combined & (uint32_t) ELeafOffsetMask;
		}

		/// Assuming this is a leaf node, return the last primitive index
		FINLINE IndexType getPrimEnd() const {
			return leaf.end;
		}

		/// Return the index of an indirection node
		FINLINE IndexType getIndirectionIndex() const {
			return(inner.combined & (uint32_t) EInnerOffsetMask) >> 2;
		}

		/// Return the left child (assuming that this is an interior node)
		FINLINE const KDNode * __restrict getLeft() const {
			return this +
				((inner.combined & (uint32_t) EInnerOffsetMask) >> 2);
		}

		/// Return the sibling of the current node
		FINLINE const KDNode * __restrict getSibling() const {
			return (const KDNode *) ((ptrdiff_t) this ^ (ptrdiff_t) 8);
		}

		/// Return the left child (assuming that this is an interior node)
		FINLINE KDNode * __restrict getLeft() {
			return this +
				((inner.combined & (uint32_t) EInnerOffsetMask) >> 2);
		}

		/// Return the left child (assuming that this is an interior node)
		FINLINE const KDNode * __restrict getRight() const {
			return getLeft() + 1;
		}

		/**
		 * \brief Return the split plane location (assuming that this
		 * is an interior node)
		 */
		FINLINE float getSplit() const {
			return inner.split;
		}

		/// Return the split axis (assuming that this is an interior node)
		FINLINE int getAxis() const {
			return inner.combined & (uint32_t) EInnerAxisMask;
		}

		/// Return a string representation
		std::string toString() const {
			std::ostringstream oss;
			if (isLeaf()) {
				oss << "KDNode[leaf, primStart=" << getPrimStart()
					<< ", primCount=" << getPrimEnd()-getPrimStart() << "]";
			} else {
				oss << "KDNode[interior, axis=" << getAxis()
					<< ", split=" << getAxis()
					<< ", leftOffset="
					<< ((inner.combined & EInnerOffsetMask) >> 2)
					<< "]";
			}
			return oss.str();
		}
	};
	BOOST_STATIC_ASSERT(sizeof(KDNode) == 8);

	/// Return the log level of kd-tree status messages
	inline ELogLevel getLogLevel() const { return m_logLevel; }

	/// Return the log level of kd-tree status messages
	inline void setLogLevel(ELogLevel level) { m_logLevel = level; }

	/// Return the root node of the kd-tree
	inline const KDNode *getRoot() const {
		return m_nodes;
	}

	/// Return whether or not the kd-tree has been built
	inline bool isBuilt() const {
		return m_nodes != NULL;
	}

	/// Return a (slightly enlarged) axis-aligned bounding box containing all primitives
	inline const AABBType &getAABB() const { return m_aabb; }

	/// Return a tight axis-aligned bounding box containing all primitives
	inline const AABBType &getTightAABB() const { return m_tightAABB;}

	MTS_DECLARE_CLASS()
protected:
	virtual ~KDTreeBase() { }
protected:
	KDNode *m_nodes;
	ELogLevel m_logLevel;
	AABBType m_aabb, m_tightAABB;
};

#if defined(_MSC_VER) && !defined(__INTELLISENSE__)
/* Use strict IEEE 754 floating point computations
   for the following kd-tree building code */
MTS_NAMESPACE_END
#pragma float_control(precise, on)
MTS_NAMESPACE_BEGIN
#endif

#define KDLog(level, fmt, ...) Thread::getThread()->getLogger()->log(\
	level, KDTreeBase<AABB>::m_theClass, __FILE__, __LINE__, \
		fmt, ## __VA_ARGS__)

/**
 * \brief Optimized KD-tree acceleration data structure for
 * n-dimensional (n<=4) shapes and various queries on them.
 *
 * Note that this class mainly concerns itself with data that cover a
 * region of space. For point data, other implementations will be more
 * suitable. The most important application in Mitsuba is the fast
 * construction of high-quality trees for ray tracing. See the class
 * \ref SAHKDTree for this specialization.
 *
 * The code in this class is a fully generic kd-tree implementation, which
 * can theoretically support any kind of shape. However, subclasses still
 * need to provide the following signatures for a functional implementation:
 *
 * \code
 * /// Return the total number of primitives
 * inline SizeType getPrimitiveCount() const;
 *
 * /// Return the axis-aligned bounding box of a certain primitive
 * inline AABB getAABB(IndexType primIdx) const;
 *
 * /// Return the AABB of a primitive when clipped to another AABB
 * inline AABB getClippedAABB(IndexType primIdx, const AABBType &aabb) const;
 * \endcode
 *
 * This class follows the "Curiously recurring template" design pattern
 * so that the above functions can be inlined (in particular, no virtual
 * calls will be necessary!).
 *
 * When the kd-tree is initially built, this class optimizes a cost
 * heuristic every time a split plane has to be chosen. For ray tracing,
 * the heuristic is usually the surface area heuristic (SAH), but other
 * choices are possible as well. The tree construction heuristic must be
 * passed as a template argument, which can use a supplied AABB and
 * split candidate to compute approximate probabilities of recursing into
 * the left and right subrees during a typical kd-tree query operation.
 * See \ref SurfaceAreaHeuristic for an example of the interface that
 * must be implemented.
 *
 * The kd-tree construction algorithm creates 'perfect split' trees as
 * outlined in the paper "On Building fast kd-Trees for Ray Tracing, and on
 * doing that in O(N log N)" by Ingo Wald and Vlastimil Havran. This works
 * even when the tree is not meant to be used for ray tracing.
 * For polygonal meshes, the involved Sutherland-Hodgman iterations can be
 * quite expensive in terms of the overall construction time. The
 * \ref setClip method can be used to deactivate perfect splits at the
 * cost of a lower-quality tree.
 *
 * Because the O(N log N) construction algorithm tends to cause many
 * incoherent memory accesses, a fast approximate technique (Min-max
 * binning) is used near the top of the tree, which significantly reduces
 * cache misses. Once the input data has been narrowed down to a
 * reasonable amount, the implementation switches over to the O(N log N)
 * builder. When multiple processors are available, the build process runs
 * in parallel.
 *
 * \author Wenzel Jakob
 * \ingroup librender
 */
template <typename AABBType, typename TreeConstructionHeuristic, typename Derived>
	class GenericKDTree : public KDTreeBase<AABBType> {
protected:
	// Some forward declarations
	struct MinMaxBins;
	struct EdgeEvent;
	struct EdgeEventOrdering;

public:
	typedef KDTreeBase<AABBType>            Parent;
	typedef typename Parent::SizeType       SizeType;
	typedef typename Parent::IndexType      IndexType;
	typedef typename Parent::KDNode         KDNode;
	typedef typename AABBType::Scalar       Scalar;
	typedef typename AABBType::PointType    PointType;
	typedef typename AABBType::VectorType   VectorType;

	using Parent::m_nodes;
	using Parent::m_aabb;
	using Parent::m_tightAABB;
	using Parent::m_logLevel;
	using Parent::isBuilt;

	/**
	 * \brief Create a new kd-tree instance initialized with
	 * the default parameters.
	 */
	GenericKDTree() : m_indices(NULL) {
		m_nodes = NULL;
		m_traversalCost = 15;
		m_queryCost = 20;
		m_emptySpaceBonus = 0.9f;
		m_clip = true;
		m_stopPrims = 6;
		m_maxBadRefines = 3;
		m_exactPrimThreshold = 65536;
		m_maxDepth = 0;
		m_retract = true;
		m_parallelBuild = true;
		m_minMaxBins = 128;
		m_logLevel = EDebug;
	}

	/**
	 * \brief Release all memory
	 */
	virtual ~GenericKDTree() {
		if (m_indices)
			delete[] m_indices;
		if (m_nodes)
			freeAligned(m_nodes-1); // undo alignment shift
	}

	/**
	 * \brief Set the traversal cost used by the tree construction heuristic
	 */
	inline void setTraversalCost(Float traversalCost) {
		m_traversalCost = traversalCost;
	}

	/**
	 * \brief Returns the underlying kd-tree index buffer
	 */
	inline IndexType *getIndices() const { return m_indices; }

	/**
	 * \brief Return the traversal cost used by the tree construction heuristic
	 */
	inline Float getTraversalCost() const {
		return m_traversalCost;
	}

	/**
	 * \brief Set the query cost used by the tree construction heuristic
	 * (This is the average cost for testing a contained shape against
	 *  a kd-tree search query)
	 */
	inline void setQueryCost(Float queryCost) {
		m_queryCost = queryCost;
	}

	/**
	 * \brief Return the query cost used by the tree construction heuristic
	 * (This is the average cost for testing a contained shape against
	 *  a kd-tree search query)
	 */
	inline Float getQueryCost() const {
		return m_queryCost;
	}

	/**
	 * \brief Set the bonus factor for empty space used by the
	 * tree construction heuristic
	 */
	inline void setEmptySpaceBonus(Float emptySpaceBonus) {
		m_emptySpaceBonus = emptySpaceBonus;
	}

	/**
	 * \brief Return the bonus factor for empty space used by the
 	 * tree construction heuristic
	 */
	inline Float getEmptySpaceBonus() const {
		return m_emptySpaceBonus;
	}

	/**
	 * \brief Set the maximum tree depth (0 = use heuristic)
	 */
	inline void setMaxDepth(SizeType maxDepth) {
		m_maxDepth = maxDepth;
	}

	/**
	 * \brief Set the number of bins used for Min-Max binning
	 */
	inline void setMinMaxBins(SizeType minMaxBins) {
		m_minMaxBins = minMaxBins;
	}

	/**
	 * \brief Return the number of bins used for Min-Max binning
	 */
	inline SizeType getMinMaxBins() const {
		return m_minMaxBins;
	}

	/**
	 * \brief Return maximum tree depth (0 = use heuristic)
	 */
	inline SizeType getMaxDepth() const {
		return m_maxDepth;
	}

	/**
	 * \brief Specify whether or not to use primitive clipping will
	 * be used in the tree construction.
	 */
	inline void setClip(bool clip) {
		m_clip = clip;
	}

	/**
	 * \brief Return whether or not to use primitive clipping will
	 * be used in the tree construction.
	 */
	inline bool getClip() const {
		return m_clip;
	}

	/**
	 * \brief Specify whether or not bad splits can be "retracted".
	 */
	inline void setRetract(bool retract) {
		m_retract = retract;
	}

	/**
	 * \brief Return whether or not bad splits can be "retracted".
	 */
	inline bool getRetract() const {
		return m_retract;
	}

	/**
	 * \brief Set the number of bad refines allowed to happen
	 * in succession before a leaf node will be created.
	 */
	inline void setMaxBadRefines(SizeType maxBadRefines) {
		m_maxBadRefines = maxBadRefines;
	}

	/**
	 * \brief Return the number of bad refines allowed to happen
	 * in succession before a leaf node will be created.
	 */
	inline SizeType getMaxBadRefines() const {
		return m_maxBadRefines;
	}

	/**
	 * \brief Set the number of primitives, at which recursion will
	 * stop when building the tree.
	 */
	inline void setStopPrims(SizeType stopPrims) {
		m_stopPrims = stopPrims;
	}

	/**
	 * \brief Return the number of primitives, at which recursion will
	 * stop when building the tree.
	 */
	inline SizeType getStopPrims() const {
		return m_stopPrims;
	}

	/**
	 * \brief Specify whether or not tree construction
	 * should run in parallel.
	 */
	inline void setParallelBuild(bool parallel) {
		m_parallelBuild = parallel;
	}

	/**
	 * \brief Return whether or not tree construction
	 * will run in parallel.
	 */
	inline bool getParallelBuild() const {
		return m_parallelBuild;
	}

	/**
	 * \brief Specify the number of primitives, at which the builder will
	 * switch from (approximate) Min-Max binning to the accurate
	 * O(n log n) optimization method.
	 */
	inline void setExactPrimitiveThreshold(SizeType exactPrimThreshold) {
		m_exactPrimThreshold = exactPrimThreshold;
	}

	/**
	 * \brief Return the number of primitives, at which the builder will
	 * switch from (approximate) Min-Max binning to the accurate
	 * O(n log n) optimization method.
	 */
	inline SizeType getExactPrimitiveThreshold() const {
		return m_exactPrimThreshold;
	}
protected:
	/**
	 * \brief Once the tree has been constructed, it is rewritten into
	 * a more convenient binary storage format.
	 *
	 * This data structure is used to store temporary information about
	 * this process.
	 */
	struct RewriteItem {
		const KDNode *node;
		KDNode *target;
		const typename GenericKDTree::BuildContext *context;
		AABBType aabb;

		RewriteItem(const KDNode *node, KDNode *target,
			const typename GenericKDTree::BuildContext *context, const AABBType &aabb) :
			node(node), target(target), context(context), aabb(aabb) { }
	};

	/**
	 * \brief Build a KD-tree over the supplied geometry
	 *
	 * To be called by the subclass.
	 */
	void buildInternal() {
		/* Some samity checks */
		if (isBuilt())
			KDLog(EError, "The kd-tree has already been built!");
		if (m_traversalCost <= 0)
			KDLog(EError, "The traveral cost must be > 0");
		if (m_queryCost <= 0)
			KDLog(EError, "The query cost must be > 0");
		if (m_emptySpaceBonus <= 0 || m_emptySpaceBonus > 1)
			KDLog(EError, "The empty space bonus must be in [0, 1]");
		if (m_minMaxBins <= 1)
			KDLog(EError, "The number of min-max bins must be > 2");

		SizeType primCount = cast()->getPrimitiveCount();
		if (primCount == 0) {
			KDLog(EWarn, "kd-tree contains no geometry!");
			// +1 shift is for alignment purposes (see KDNode::getSibling)
			m_nodes = static_cast<KDNode *>(allocAligned(sizeof(KDNode) * 2))+1;
			m_nodes[0].initLeafNode(0, 0);
			return;
		}

		if (primCount <= m_exactPrimThreshold)
			m_parallelBuild = false;

		BuildContext ctx(primCount, m_minMaxBins);

		/* Establish an ad-hoc depth cutoff value (Formula from PBRT) */
		if (m_maxDepth == 0)
			m_maxDepth = (int) (8 + 1.3f * math::log2i(primCount));
		m_maxDepth = std::min(m_maxDepth, (SizeType) MTS_KD_MAXDEPTH);

		KDLog(m_logLevel, "Creating a preliminary index list (%s)",
			memString(primCount * sizeof(IndexType)).c_str());

		OrderedChunkAllocator &leftAlloc = ctx.leftAlloc;
		IndexType *indices = leftAlloc.allocate<IndexType>(primCount);

		ref<Timer> timer = new Timer();
		AABBType &aabb = m_aabb;
		aabb.reset();
		for (IndexType i=0; i<primCount; ++i) {
			aabb.expandBy(cast()->getAABB(i));
			indices[i] = i;
		}

		#if defined(DOUBLE_PRECISION)
			for (int i=0; i<3; ++i) {
				aabb.min[i] = math::castflt_down(aabb.min[i]);
				aabb.max[i] = math::castflt_up(aabb.max[i]);
			}
		#endif

		KDLog(m_logLevel, "Computed scene bounds in %i ms",
				timer->getMilliseconds());
		KDLog(m_logLevel, "");

		KDLog(m_logLevel, "kd-tree configuration:");
		KDLog(m_logLevel, "   Traversal cost           : %.2f", m_traversalCost);
		KDLog(m_logLevel, "   Query cost               : %.2f", m_queryCost);
		KDLog(m_logLevel, "   Empty space bonus        : %.2f", m_emptySpaceBonus);
		KDLog(m_logLevel, "   Max. tree depth          : %i", m_maxDepth);
		KDLog(m_logLevel, "   Scene bounding box (min) : %s",
				aabb.min.toString().c_str());
		KDLog(m_logLevel, "   Scene bounding box (max) : %s",
				aabb.max.toString().c_str());
		KDLog(m_logLevel, "   Min-max bins             : %i", m_minMaxBins);
		KDLog(m_logLevel, "   O(n log n) method        : use for <= %i primitives",
				m_exactPrimThreshold);
		KDLog(m_logLevel, "   Perfect splits           : %s", m_clip ? "yes" : "no");
		KDLog(m_logLevel, "   Retract bad splits       : %s",
				m_retract ? "yes" : "no");
		KDLog(m_logLevel, "   Stopping primitive count : %i", m_stopPrims);
		KDLog(m_logLevel, "   Build tree in parallel   : %s",
				m_parallelBuild ? "yes" : "no");
		KDLog(m_logLevel, "");

		SizeType procCount = getCoreCount();
		if (procCount == 1)
			m_parallelBuild = false;

		if (m_parallelBuild) {
			m_builders.resize(procCount);
			for (SizeType i=0; i<procCount; ++i) {
				m_builders[i] = new TreeBuilder(i, this);
				m_builders[i]->incRef();
				m_builders[i]->start();
			}
		}

		m_indirectionLock = new Mutex();
		KDNode *prelimRoot = ctx.nodes.allocate(1);
		buildTreeMinMax(ctx, 1, prelimRoot, aabb, aabb,
				indices, primCount, true, 0);
		ctx.leftAlloc.release(indices);

		KDAssert(ctx.leftAlloc.used() == 0);
		KDAssert(ctx.rightAlloc.used() == 0);

		if (m_parallelBuild) {
			UniqueLock lock(m_interface.mutex);
			m_interface.done = true;
			m_interface.cond->broadcast();
			lock.unlock();
			for (SizeType i=0; i<m_builders.size(); ++i)
				m_builders[i]->join();
		}

		KDLog(EInfo, "Finished -- took %i ms.", timer->getMilliseconds());
		KDLog(m_logLevel, "");

		KDLog(m_logLevel, "Temporary memory statistics:");
		KDLog(m_logLevel, "   Classification storage : %s",
				memString((ctx.classStorage.size() * (1+procCount))).c_str());
		KDLog(m_logLevel, "   Indirection entries    : " SIZE_T_FMT " (%s)",
				m_indirections.size(), memString(m_indirections.capacity()
				* sizeof(KDNode *)).c_str());

		KDLog(m_logLevel, "   Main thread:");
		ctx.printStats(m_logLevel);
		size_t totalUsage = m_indirections.capacity()
			* sizeof(KDNode *) + ctx.size();

		/// Clean up event lists and print statistics
		ctx.leftAlloc.cleanup();
		ctx.rightAlloc.cleanup();
		for (SizeType i=0; i<m_builders.size(); ++i) {
			KDLog(m_logLevel, "   Worker thread %i:", i+1);
			BuildContext &subCtx = m_builders[i]->getContext();
			subCtx.printStats(m_logLevel);
			totalUsage += subCtx.size();
			subCtx.leftAlloc.cleanup();
			subCtx.rightAlloc.cleanup();
			ctx.accumulateStatisticsFrom(subCtx);
		}
		KDLog(m_logLevel, "   Total: %s", memString(totalUsage).c_str());

		KDLog(m_logLevel, "");
		timer->reset();
		KDLog(m_logLevel, "Optimizing memory layout ..");

		Float expTraversalSteps = 0;
		Float expLeavesVisited = 0;
		Float expPrimitivesIntersected = 0;
		Float heuristicCost = 0;

		SizeType nodePtr = 0, indexPtr = 0;
		SizeType maxPrimsInLeaf = 0;
		const SizeType primBucketCount = 16;
		SizeType primBuckets[primBucketCount];
		memset(primBuckets, 0, sizeof(SizeType)*primBucketCount);
		m_nodeCount = ctx.innerNodeCount + ctx.leafNodeCount;
		m_indexCount = ctx.primIndexCount;

		// +1 shift is for alignment purposes (see KDNode::getSibling)
		m_nodes = static_cast<KDNode *> (allocAligned(
				sizeof(KDNode) * (m_nodeCount+1)))+1;
		m_indices = new IndexType[m_indexCount];

		/* The following code rewrites all tree nodes with proper relative
		   indices. It also computes the final tree cost and some other
		   useful heuristics */
		std::stack<RewriteItem> stack;
		stack.push(RewriteItem(prelimRoot, &m_nodes[nodePtr++], &ctx, aabb));
		while (!stack.empty()) {
			RewriteItem item = stack.top();
			stack.pop();

			typename std::map<const KDNode *, IndexType>::const_iterator it
				= m_interface.threadMap.find(item.node);
			// Check if we're switching to a subtree built by a worker thread
			if (it != m_interface.threadMap.end())
				item.context = &m_builders[(*it).second]->getContext();

			if (item.node->isLeaf()) {
				SizeType primStart = item.node->getPrimStart(),
						  primEnd = item.node->getPrimEnd(),
						  primsInLeaf = primEnd-primStart;
				item.target->initLeafNode(indexPtr, primsInLeaf);

				Float quantity = TreeConstructionHeuristic::getQuantity(item.aabb),
					  weightedQuantity = quantity * primsInLeaf;
				expLeavesVisited += quantity;
				expPrimitivesIntersected += weightedQuantity;
				heuristicCost += weightedQuantity * m_queryCost;
				if (primsInLeaf < primBucketCount)
					primBuckets[primsInLeaf]++;
				if (primsInLeaf > maxPrimsInLeaf)
					maxPrimsInLeaf = primsInLeaf;

				const BlockedVector<IndexType, MTS_KD_BLOCKSIZE_IDX> &indices
					= item.context->indices;
				for (SizeType idx = primStart; idx<primEnd; ++idx) {
					KDAssert(indices[idx] >= 0 && indices[idx] < primCount);
					m_indices[indexPtr++] = indices[idx];
				}
			} else {
				Float quantity = TreeConstructionHeuristic::getQuantity(item.aabb);
				expTraversalSteps += quantity;
				heuristicCost += quantity * m_traversalCost;

				const KDNode *left;
				if (EXPECT_TAKEN(!item.node->isIndirection()))
					left = item.node->getLeft();
				else
					left = m_indirections[item.node->getIndirectionIndex()];

				KDNode *children = &m_nodes[nodePtr];
				nodePtr += 2;
				int axis = item.node->getAxis();
				float split = item.node->getSplit();
				bool result = item.target->initInnerNode(axis, split, children - item.target);
				if (!result)
					KDLog(EError, "Cannot represent relative pointer -- "
						"too many primitives?");

				Float tmp = item.aabb.min[axis];
				item.aabb.min[axis] = split;
				stack.push(RewriteItem(left+1, children+1, item.context, item.aabb));
				item.aabb.min[axis] = tmp;
				item.aabb.max[axis] = split;
				stack.push(RewriteItem(left, children, item.context, item.aabb));
			}
		}

		KDAssert(nodePtr == ctx.innerNodeCount + ctx.leafNodeCount);
		KDAssert(indexPtr == m_indexCount);

		KDLog(m_logLevel, "Finished -- took %i ms.", timer->getMilliseconds());

		/* Free some more memory */
		ctx.nodes.clear();
		ctx.indices.clear();
		for (SizeType i=0; i<m_builders.size(); ++i) {
			BuildContext &subCtx = m_builders[i]->getContext();
			subCtx.nodes.clear();
			subCtx.indices.clear();
		}
		m_indirectionLock = NULL;
		std::vector<KDNode *>().swap(m_indirections);

		if (m_builders.size() > 0) {
			for (SizeType i=0; i<m_builders.size(); ++i)
				m_builders[i]->decRef();
			m_builders.clear();
		}

		KDLog(m_logLevel, "");

		Float rootQuantity = TreeConstructionHeuristic::getQuantity(aabb);
		expTraversalSteps /= rootQuantity;
		expLeavesVisited /= rootQuantity;
		expPrimitivesIntersected /= rootQuantity;
		heuristicCost /= rootQuantity;

		/* Slightly enlarge the bounding box
		   (necessary e.g. when the scene is planar) */
		m_tightAABB = aabb;

		const Float eps = MTS_KD_AABB_EPSILON;

		aabb.min -= (aabb.max-aabb.min) * eps + VectorType(eps);
		aabb.max += (aabb.max-aabb.min) * eps + VectorType(eps);

		KDLog(m_logLevel, "Structural kd-tree statistics:");
		KDLog(m_logLevel, "   Parallel work units         : " SIZE_T_FMT,
				m_interface.threadMap.size());
		KDLog(m_logLevel, "   Node storage cost           : %s",
				memString(nodePtr * sizeof(KDNode)).c_str());
		KDLog(m_logLevel, "   Index storage cost          : %s",
				memString(indexPtr * sizeof(IndexType)).c_str());
		KDLog(m_logLevel, "   Inner nodes                 : %i", ctx.innerNodeCount);
		KDLog(m_logLevel, "   Leaf nodes                  : %i", ctx.leafNodeCount);
		KDLog(m_logLevel, "   Nonempty leaf nodes         : %i",
				ctx.nonemptyLeafNodeCount);
		std::ostringstream oss;
		oss << "   Leaf node histogram         : ";
		for (SizeType i=0; i<primBucketCount; i++) {
			oss << i << "(" << primBuckets[i] << ") ";
			if ((i+1)%4==0 && i+1<primBucketCount) {
				KDLog(m_logLevel, "%s", oss.str().c_str());
				oss.str("");
				oss << "                                 ";
			}
		}
		KDLog(m_logLevel, "%s", oss.str().c_str());
		KDLog(m_logLevel, "");
		KDLog(m_logLevel, "Qualitative kd-tree statistics:");
		KDLog(m_logLevel, "   Retracted splits            : %i", ctx.retractedSplits);
		KDLog(m_logLevel, "   Pruned primitives           : %i", ctx.pruned);
		KDLog(m_logLevel, "   Largest leaf node           : %i primitives",
				maxPrimsInLeaf);
		KDLog(m_logLevel, "   Avg. prims/nonempty leaf    : %.2f",
				ctx.primIndexCount / (Float) ctx.nonemptyLeafNodeCount);
		KDLog(m_logLevel, "   Expected traversals/query   : %.2f", expTraversalSteps);
		KDLog(m_logLevel, "   Expected leaf visits/query  : %.2f", expLeavesVisited);
		KDLog(m_logLevel, "   Expected prim. visits/query : %.2f",
				expPrimitivesIntersected);
		KDLog(m_logLevel, "   Final cost                  : %.2f", heuristicCost);
		KDLog(m_logLevel, "");

		#if defined(__LINUX__)
			/* Forcefully release Heap memory back to the OS */
			malloc_trim(0);
		#endif
	}

protected:
	/// Primitive classification during tree-construction
	enum EClassificationResult {
		/// Straddling primitive
		EBothSides = 0,
		/// Primitive is entirely on the left side of the split
		ELeftSide = 1,
		/// Primitive is entirely on the right side of the split
		ERightSide = 2,
		/// Edge events have been generated for the straddling primitive
		EBothSidesProcessed = 3
	};

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
		inline EdgeEvent(uint16_t type, int axis, float pos, IndexType index)
		 : pos(pos), index(index), type(type), axis(axis) { }

		/// Return a string representation
		std::string toString() const {
			std::ostringstream oss;
			oss << "EdgeEvent[" << endl
				<< "  pos = " << pos << "," << endl
				<< "  index = " << index << "," << endl
				<< "  type = ";
			if (type == EEdgeEnd)
				oss << "end";
			else if (type == EEdgePlanar)
				oss << "planar";
			else if (type == EEdgeStart)
				oss << "start";
			else
				oss << "unknown!";
			oss << "," << endl
				<< "  axis = " << axis << endl
				<<"]";
			return oss.str();
		}

		/// Plane position
		float pos;
		/// Primitive index
		IndexType index;
		/// Event type: end/planar/start
		unsigned int type:2;
		/// Event axis
		unsigned int axis:2;
	};

	BOOST_STATIC_ASSERT(sizeof(EdgeEvent) == 12);

	/// Edge event comparison functor
	struct EdgeEventOrdering : public std::binary_function<EdgeEvent, EdgeEvent, bool> {
		inline bool operator()(const EdgeEvent &a, const EdgeEvent &b) const {
			if (a.axis != b.axis)
				return a.axis < b.axis;
			if (a.pos != b.pos)
				return a.pos < b.pos;
			if (a.type != b.type)
				return a.type < b.type;
			return a.index < b.index;
		}
	};

	/**
	 * \brief Data type for split candidates computed by
	 * the O(n log n) greedy optimization method.
	 * */
	struct SplitCandidate {
		Float cost;
		float pos;
		int axis;
		SizeType numLeft, numRight;
		bool planarLeft;
		int leftBin;

		inline SplitCandidate() :
			cost(std::numeric_limits<Float>::infinity()),
			pos(0), axis(0), numLeft(0), numRight(0), planarLeft(false), leftBin(-1) {
		}

		std::string toString() const {
			std::ostringstream oss;
			oss << "SplitCandidate[" << endl
				<< "  cost=" << cost << "," << endl
				<< "  pos=" << pos << "," << endl
				<< "  axis=" << axis << "," << endl
				<< "  numLeft=" << numLeft << "," << endl
				<< "  numRight=" << numRight << "," << endl
				<< "  leftBin=" << leftBin << "," << endl
				<< "  planarLeft=" << (planarLeft ? "yes" : "no") << endl
				<< "]";
			return oss.str();
		}
	};

	/**
	 * \brief Per-thread context used to manage memory allocations,
	 * also records some useful statistics.
	 */
	struct BuildContext {
		OrderedChunkAllocator leftAlloc, rightAlloc;
		BlockedVector<KDNode, MTS_KD_BLOCKSIZE_KD> nodes;
		BlockedVector<IndexType, MTS_KD_BLOCKSIZE_IDX> indices;
		ClassificationStorage classStorage;
		MinMaxBins minMaxBins;

		SizeType leafNodeCount;
		SizeType nonemptyLeafNodeCount;
		SizeType innerNodeCount;
		SizeType primIndexCount;
		SizeType retractedSplits;
		SizeType pruned;

		BuildContext(SizeType primCount, SizeType binCount)
				: minMaxBins(binCount) {
			classStorage.setPrimitiveCount(primCount);
			leafNodeCount = 0;
			nonemptyLeafNodeCount = 0;
			innerNodeCount = 0;
			primIndexCount = 0;
			retractedSplits = 0;
			pruned = 0;
		}

		size_t size() {
			return leftAlloc.size() + rightAlloc.size()
				+ nodes.capacity() * sizeof(KDNode)
				+ indices.capacity() * sizeof(IndexType)
				+ classStorage.size();
		}

		void printStats(ELogLevel level) {
			KDLog(level, "      Left events   : " SIZE_T_FMT " chunks (%s)",
					leftAlloc.getChunkCount(),
					memString(leftAlloc.size()).c_str());
			KDLog(level, "      Right events  : " SIZE_T_FMT " chunks (%s)",
					rightAlloc.getChunkCount(),
					memString(rightAlloc.size()).c_str());
			KDLog(level, "      kd-tree nodes : " SIZE_T_FMT " entries, "
					SIZE_T_FMT " blocks (%s)", nodes.size(), nodes.blockCount(),
					memString(nodes.capacity() * sizeof(KDNode)).c_str());
			KDLog(level, "      Indices       : " SIZE_T_FMT " entries, "
					SIZE_T_FMT " blocks (%s)", indices.size(),
					indices.blockCount(), memString(indices.capacity()
					* sizeof(IndexType)).c_str());
		}

		void accumulateStatisticsFrom(const BuildContext &ctx) {
			leafNodeCount += ctx.leafNodeCount;
			nonemptyLeafNodeCount += ctx.nonemptyLeafNodeCount;
			innerNodeCount += ctx.innerNodeCount;
			primIndexCount += ctx.primIndexCount;
			retractedSplits += ctx.retractedSplits;
			pruned += ctx.pruned;
		}
	};

	/**
	 * \brief Communication data structure used to pass jobs to
	 * kd-tree builder threads
	 */
	struct BuildInterface {
		/* Communcation */
		ref<Mutex> mutex;
		ref<ConditionVariable> cond, condJobTaken;
		std::map<const KDNode *, IndexType> threadMap;
		bool done;

		/* Job description for building a subtree */
		int depth;
		KDNode *node;
		AABBType nodeAABB;
		EdgeEvent *eventStart, *eventEnd;
		SizeType primCount;
		int badRefines;

		inline BuildInterface() {
			mutex = new Mutex();
			cond = new ConditionVariable(mutex);
			condJobTaken = new ConditionVariable(mutex);
			node = NULL;
			done = false;
		}
	};

	/**
	 * \brief kd-tree builder thread
	 */
	class TreeBuilder : public Thread {
	public:
		TreeBuilder(IndexType id, GenericKDTree *parent)
			: Thread(formatString("bld%i", id)),
			m_id(id),
			m_parent(parent),
			m_context(parent->cast()->getPrimitiveCount(),
					  parent->getMinMaxBins()),
			m_interface(parent->m_interface) {
			setCritical(true);
		}

		~TreeBuilder() {
			KDAssert(m_context.leftAlloc.used() == 0);
			KDAssert(m_context.rightAlloc.used() == 0);
		}

		void run() {
			OrderedChunkAllocator &leftAlloc = m_context.leftAlloc;
			while (true) {
				UniqueLock lock(m_interface.mutex);
				while (!m_interface.done && !m_interface.node)
					m_interface.cond->wait();
				if (m_interface.done) {
					break;
				}
				int depth = m_interface.depth;
				KDNode *node = m_interface.node;
				AABBType nodeAABB = m_interface.nodeAABB;
				size_t eventCount = m_interface.eventEnd - m_interface.eventStart;
				SizeType primCount = m_interface.primCount;
				int badRefines = m_interface.badRefines;
				EdgeEvent *eventStart = leftAlloc.allocate<EdgeEvent>(eventCount),
						  *eventEnd = eventStart + eventCount;
				memcpy(eventStart, m_interface.eventStart,
						eventCount * sizeof(EdgeEvent));
				m_interface.threadMap[node] = m_id;
				m_interface.node = NULL;
				m_interface.condJobTaken->signal();
				lock.unlock();

				std::sort(eventStart, eventEnd, EdgeEventOrdering());
				m_parent->buildTree(m_context, depth, node,
					nodeAABB, eventStart, eventEnd, primCount, true, badRefines);
				leftAlloc.release(eventStart);
			}
		}

		inline BuildContext &getContext() {
			return m_context;
		}

	private:
		IndexType m_id;
		GenericKDTree *m_parent;
		BuildContext m_context;
		BuildInterface &m_interface;
	};

	/// Cast to the derived class
	inline Derived *cast() {
		return static_cast<Derived *>(this);
	}

	/// Cast to the derived class (const version)
	inline const Derived *cast() const {
		return static_cast<const Derived *>(this);
	}

	struct EventList {
		EdgeEvent *start, *end;
		SizeType primCount;

		EventList(EdgeEvent *start, EdgeEvent *end, SizeType primCount)
			: start(start), end(end), primCount(primCount) { }
	};

	/**
	 * \brief Create an edge event list for a given list of primitives.
	 *
	 * This is necessary when passing from Min-Max binning to the more
	 * accurate O(n log n) optimizier.
	 */
	EventList createEventList(
			OrderedChunkAllocator &alloc, const AABBType &nodeAABB,
			IndexType *prims, SizeType primCount) {
		SizeType initialSize = primCount * 2 * PointType::dim, actualPrimCount = 0;
		EdgeEvent *eventStart = alloc.allocate<EdgeEvent>(initialSize);
		EdgeEvent *eventEnd = eventStart;

		for (SizeType i=0; i<primCount; ++i) {
			IndexType index = prims[i];
			AABBType aabb;
			if (m_clip) {
				aabb = cast()->getClippedAABB(index, nodeAABB);
				if (!aabb.isValid() || aabb.getSurfaceArea() == 0)
					continue;
			} else {
				aabb = cast()->getAABB(index);
			}

			for (int axis=0; axis<PointType::dim; ++axis) {
				float min = math::castflt_down(aabb.min[axis]),
				      max = math::castflt_up(aabb.max[axis]);

				if (min == max) {
					*eventEnd++ = EdgeEvent(EdgeEvent::EEdgePlanar, axis,
							min, index);
				} else {
					*eventEnd++ = EdgeEvent(EdgeEvent::EEdgeStart, axis,
							min, index);
					*eventEnd++ = EdgeEvent(EdgeEvent::EEdgeEnd, axis,
							max, index);
				}
			}
			++actualPrimCount;
		}

		SizeType newSize = (SizeType) (eventEnd - eventStart);
		if (newSize != initialSize)
			alloc.shrinkAllocation<EdgeEvent>(eventStart, newSize);

		return EventList(eventStart, eventEnd, actualPrimCount);
	}

	/**
	 * \brief Leaf node creation helper function
	 *
	 * \param ctx
	 *     Thread-specific build context containing allocators etc.
	 * \param node
	 *     KD-tree node entry to be filled
	 * \param eventStart
	 *     Start pointer of an edge event list
	 * \param eventEnd
	 *     End pointer of an edge event list
	 * \param primCount
	 *     Total primitive count for the current node
	 */
	void createLeaf(BuildContext &ctx, KDNode *node, EdgeEvent *eventStart,
			EdgeEvent *eventEnd, SizeType primCount) {
		node->initLeafNode((SizeType) ctx.indices.size(), primCount);
		if (primCount > 0) {
			SizeType seenPrims = 0;
			ctx.nonemptyLeafNodeCount++;

			for (EdgeEvent *event = eventStart; event != eventEnd
					&& event->axis == 0; ++event) {
				if (event->type == EdgeEvent::EEdgeStart ||
					event->type == EdgeEvent::EEdgePlanar) {
					ctx.indices.push_back(event->index);
					seenPrims++;
				}
			}
			KDAssert(seenPrims == primCount);
			ctx.primIndexCount += primCount;
		}
		ctx.leafNodeCount++;
	}

	/**
	 * \brief Leaf node creation helper function
	 *
	 * \param ctx
	 *     Thread-specific build context containing allocators etc.
	 * \param node
	 *     KD-tree node entry to be filled
	 * \param indices
	 *     Start pointer of an index list
	 * \param primCount
	 *     Total primitive count for the current node
	 */
	void createLeaf(BuildContext &ctx, KDNode *node, SizeType *indices,
			SizeType primCount) {
		node->initLeafNode((SizeType) ctx.indices.size(), primCount);
		if (primCount > 0) {
			ctx.nonemptyLeafNodeCount++;
			for (SizeType i=0; i<primCount; ++i)
				ctx.indices.push_back(indices[i]);
			ctx.primIndexCount += primCount;
		}
		ctx.leafNodeCount++;
	}

	/**
	 * \brief Leaf node creation helper function.
	 *
	 * Creates a unique index list by collapsing
	 * a subtree with a bad cost.
	 *
	 * \param ctx
	 *     Thread-specific build context containing allocators etc.
	 * \param node
	 *     KD-tree node entry to be filled
	 * \param start
	 *     Start pointer of the subtree indices
	 */
	void createLeafAfterRetraction(BuildContext &ctx, KDNode *node, SizeType start) {
		SizeType indexCount = (SizeType) (ctx.indices.size() - start);
		SAssert(indexCount > 0);

		OrderedChunkAllocator &alloc = ctx.leftAlloc;

		/* A temporary list is allocated to do the sorting (the indices
		   are not guaranteed to be contiguous in memory) */
		IndexType *tempStart = alloc.allocate<IndexType>(indexCount),
				   *tempEnd = tempStart + indexCount,
				   *ptr = tempStart;

		for (SizeType i=start, end = start + indexCount; i<end; ++i)
			*ptr++ = ctx.indices[i];

		/* Generate an index list without duplicate entries */
		std::sort(tempStart, tempEnd, std::less<IndexType>());
		ptr = tempStart;

		int idx = start;
		while (ptr < tempEnd) {
			ctx.indices[idx] = *ptr++;
			while (ptr < tempEnd && *ptr == ctx.indices[idx])
				++ptr;
			idx++;
		}

		int nSeen = idx-start;
		ctx.primIndexCount = ctx.primIndexCount - indexCount + nSeen;
		ctx.indices.resize(idx);
		alloc.release(tempStart);
		node->initLeafNode(start, nSeen);
		ctx.nonemptyLeafNodeCount++;
		ctx.leafNodeCount++;
	}

	/**
	 * \brief Implements the transition from min-max-binning to the
	 * O(n log n) optimization.
	 *
	 * \param ctx
	 *     Thread-specific build context containing allocators etc.
	 * \param depth
	 *     Current tree depth (1 == root node)
	 * \param node
	 *     KD-tree node entry to be filled
	 * \param nodeAABB
	 *     Axis-aligned bounding box of the current node
	 * \param indices
	 *     Index list of all triangles in the current node (for min-max binning)
	 * \param primCount
	 *     Total primitive count for the current node
	 * \param isLeftChild
	 *     Is this node the left child of its parent? This is important for
	 *     memory management using the \ref OrderedChunkAllocator.
	 * \param badRefines
	 *     Number of "probable bad refines" further up the tree. This makes
	 *     it possible to split along an initially bad-looking candidate in
	 *     the hope that the cost was significantly overestimated. The
	 *     counter makes sure that only a limited number of such splits can
	 *     happen in succession.
	 * \returns
	 *     Final cost of the node
	 */
	inline Float transitionToNLogN(BuildContext &ctx, unsigned int depth, KDNode *node,
			const AABBType &nodeAABB, IndexType *indices,
			SizeType primCount, bool isLeftChild, SizeType badRefines) {
		OrderedChunkAllocator &alloc = isLeftChild
				? ctx.leftAlloc : ctx.rightAlloc;
		EventList events = createEventList(alloc, nodeAABB, indices, primCount);
		Float cost;
		if (m_parallelBuild) {
			LockGuard lock(m_interface.mutex);
			m_interface.depth = depth;
			m_interface.node = node;
			m_interface.nodeAABB = nodeAABB;
			m_interface.eventStart = events.start;
			m_interface.eventEnd = events.end;
			m_interface.primCount = events.primCount;
			m_interface.badRefines = badRefines;
			m_interface.cond->signal();

			/* Wait for a worker thread to take this job */
			while (m_interface.node)
				m_interface.condJobTaken->wait();

			// Never tear down this subtree (return a cost of -infinity)
			cost = -std::numeric_limits<Float>::infinity();
		} else {
			std::sort(events.start, events.end, EdgeEventOrdering());

			cost = buildTree(ctx, depth, node, nodeAABB, events.start,
				events.end, events.primCount, isLeftChild, badRefines);
		}
		alloc.release(events.start);
		return cost;
	}

	/**
	 * \brief Build helper function (min-max binning)
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
	 *     memory management using the \ref OrderedChunkAllocator.
	 * \param badRefines
	 *     Number of "probable bad refines" further up the tree. This makes
	 *     it possible to split along an initially bad-looking candidate in
	 *     the hope that the cost was significantly overestimated. The
	 *     counter makes sure that only a limited number of such splits can
	 *     happen in succession.
	 * \returns
	 *     Final cost of the node
	 */
	Float buildTreeMinMax(BuildContext &ctx, unsigned int depth, KDNode *node,
			const AABBType &nodeAABB, const AABBType &tightAABB, IndexType *indices,
			SizeType primCount, bool isLeftChild, SizeType badRefines) {
		KDAssert(nodeAABB.contains(tightAABB));

		Float leafCost = primCount * m_queryCost;
		if (primCount <= m_stopPrims || depth >= m_maxDepth) {
			createLeaf(ctx, node, indices, primCount);
			return leafCost;
		}

		if (primCount <= m_exactPrimThreshold)
			return transitionToNLogN(ctx, depth, node, nodeAABB, indices,
				primCount, isLeftChild, badRefines);

		/* ==================================================================== */
	    /*                              Binning                                 */
	    /* ==================================================================== */

		ctx.minMaxBins.setAABB(tightAABB);
		ctx.minMaxBins.bin(cast(), indices, primCount);

		/* ==================================================================== */
	    /*                        Split candidate search                        */
    	/* ==================================================================== */
		SplitCandidate bestSplit = ctx.minMaxBins.minimizeCost(m_traversalCost,
				m_queryCost);

		if (bestSplit.cost == std::numeric_limits<Float>::infinity()) {
			/* This is bad: we have either run out of floating point precision to
			   accurately represent split planes (e.g. 'tightAABB' is almost collapsed
			   along an axis), or the compiler made overly liberal use of floating point
			   optimizations, causing the two stages of the min-max binning code to
			   become inconsistent. The two ways to proceed at this point are to
			   either create a leaf (bad) or switch over to the O(n log n) greedy
			   optimization, which is done below */
			KDLog(EWarn, "Min-max binning was unable to split %i primitives with %s "
				"-- retrying with the O(n log n) greedy optimization",
				primCount, tightAABB.toString().c_str());
			return transitionToNLogN(ctx, depth, node, nodeAABB, indices,
				primCount, isLeftChild, badRefines);
		}

		/* "Bad refines" heuristic from PBRT */
		if (bestSplit.cost >= leafCost) {
			if ((bestSplit.cost > 4 * leafCost && primCount < 16)
				|| badRefines >= m_maxBadRefines) {
				createLeaf(ctx, node, indices, primCount);
				return leafCost;
			}
			++badRefines;
		}

		/* ==================================================================== */
	    /*                            Partitioning                              */
	    /* ==================================================================== */

		typename MinMaxBins::Partition partition =
			ctx.minMaxBins.partition(ctx, cast(), indices, bestSplit,
			isLeftChild, m_traversalCost, m_queryCost);

		/* ==================================================================== */
	    /*                              Recursion                               */
	    /* ==================================================================== */

		KDNode *children = ctx.nodes.allocate(2);

		SizeType nodePosBeforeSplit = (SizeType) ctx.nodes.size();
		SizeType indexPosBeforeSplit = (SizeType) ctx.indices.size();
		SizeType leafNodeCountBeforeSplit = ctx.leafNodeCount;
		SizeType nonemptyLeafNodeCountBeforeSplit = ctx.nonemptyLeafNodeCount;
		SizeType innerNodeCountBeforeSplit = ctx.innerNodeCount;

		if (!node->initInnerNode(bestSplit.axis, bestSplit.pos, children-node)) {
			LockGuard lock(m_indirectionLock);
			SizeType indirectionIdx = (SizeType) m_indirections.size();
			m_indirections.push_back(children);
			/* Unable to store relative offset -- create an indirection
			   table entry */
			node->initIndirectionNode(bestSplit.axis, bestSplit.pos,
					indirectionIdx);
		}
		ctx.innerNodeCount++;

		AABBType childAABB(nodeAABB);
		childAABB.max[bestSplit.axis] = bestSplit.pos;

		Float leftCost = buildTreeMinMax(ctx, depth+1, children,
				childAABB, partition.left, partition.leftIndices,
				bestSplit.numLeft, true, badRefines);

		childAABB.min[bestSplit.axis] = bestSplit.pos;
		childAABB.max[bestSplit.axis] = nodeAABB.max[bestSplit.axis];

		Float rightCost = buildTreeMinMax(ctx, depth+1, children + 1,
				childAABB, partition.right, partition.rightIndices,
				bestSplit.numRight, false, badRefines);

		TreeConstructionHeuristic tch(nodeAABB);
		std::pair<Float, Float> prob = tch(bestSplit.axis,
			bestSplit.pos - nodeAABB.min[bestSplit.axis],
			nodeAABB.max[bestSplit.axis] - bestSplit.pos);

		/* Compute the final cost given the updated cost
		   values received from the children */
		Float finalCost = m_traversalCost +
			(prob.first * leftCost + prob.second * rightCost);

		/* Release the index lists not needed by the children anymore */
		if (isLeftChild)
			ctx.rightAlloc.release(partition.rightIndices);
		else
			ctx.leftAlloc.release(partition.leftIndices);

		/* ==================================================================== */
	    /*                           Final decision                             */
	    /* ==================================================================== */

		if (!m_retract || finalCost < primCount * m_queryCost) {
			return finalCost;
		} else {
			/* In the end, splitting didn't help to reduce the cost.
			   Tear up everything below this node and create a leaf */
			ctx.nodes.resize(nodePosBeforeSplit);
			ctx.retractedSplits++;
			ctx.leafNodeCount = leafNodeCountBeforeSplit;
			ctx.nonemptyLeafNodeCount = nonemptyLeafNodeCountBeforeSplit;
			ctx.innerNodeCount = innerNodeCountBeforeSplit;
			createLeafAfterRetraction(ctx, node, indexPosBeforeSplit);
			return leafCost;
		}
	}

	/*
	 * \brief Build helper function (greedy O(n log n) optimization)
	 *
	 * \param ctx
	 *     Thread-specific build context containing allocators etc.
	 * \param depth
	 *     Current tree depth (1 == root node)
	 * \param node
	 *     KD-tree node entry to be filled
	 * \param nodeAABB
	 *     Axis-aligned bounding box of the current node
	 * \param eventStart
	 *     Pointer to the beginning of a sorted edge event list
	 * \param eventEnd
	 *     Pointer to the end of a sorted edge event list
	 * \param primCount
	 *     Total primitive count for the current node
	 * \param isLeftChild
	 *     Is this node the left child of its parent? This is important for
	 *     memory management using the \ref OrderedChunkAllocator.
	 * \param badRefines
	 *     Number of "probable bad refines" further up the tree. This makes
	 *     it possible to split along an initially bad-looking candidate in
	 *     the hope that the cost was significantly overestimated. The
	 *     counter makes sure that only a limited number of such splits can
	 *     happen in succession.
	 * \returns
	 *     Final cost of the node
	 */
	Float buildTree(BuildContext &ctx, unsigned int depth, KDNode *node,
		const AABBType &nodeAABB, EdgeEvent *eventStart, EdgeEvent *eventEnd,
		SizeType primCount, bool isLeftChild, SizeType badRefines) {

		Float leafCost = primCount * m_queryCost;
		if (primCount <= m_stopPrims || depth >= m_maxDepth) {
			createLeaf(ctx, node, eventStart, eventEnd, primCount);
			return leafCost;
		}

		SplitCandidate bestSplit;

		/* ==================================================================== */
	    /*                        Split candidate search                        */
    	/* ==================================================================== */

		/* First, find the optimal splitting plane according to the
		   tree construction heuristic. To do this in O(n), the search is
		   implemented as a sweep over the edge events */

		/* Initially, the split plane is placed left of the scene
		   and thus all geometry is on its right side */
		SizeType numLeft[PointType::dim],
				  numRight[PointType::dim];

		for (int i=0; i<PointType::dim; ++i) {
			numLeft[i] = 0;
			numRight[i] = primCount;
		}

		EdgeEvent *eventsByAxis[PointType::dim];
		int eventsByAxisCtr = 1;
		eventsByAxis[0] = eventStart;
		TreeConstructionHeuristic tch(nodeAABB);

		/* Iterate over all events on the current axis */
		for (EdgeEvent *event = eventStart; event < eventEnd; ) {
			/* Record the current position and count the number
			   and type of remaining events, which are also here.
			   Due to the sort ordering, there is no need to worry
			   about an axis change in the loops below */
			int axis = event->axis;
			float pos = event->pos;
			SizeType numStart = 0, numEnd = 0, numPlanar = 0;

			/* Count "end" events */
			while (event < eventEnd && event->pos == pos
				&& event->axis == axis
				&& event->type == EdgeEvent::EEdgeEnd) {
				++numEnd; ++event;
			}

			/* Count "planar" events */
			while (event < eventEnd && event->pos == pos
				&& event->axis == axis
				&& event->type == EdgeEvent::EEdgePlanar) {
				++numPlanar; ++event;
			}

			/* Count "start" events */
			while (event < eventEnd && event->pos == pos
				&& event->axis == axis
				&& event->type == EdgeEvent::EEdgeStart) {
				++numStart; ++event;
			}

			/* Keep track of the beginning of dimensions */
			if (event < eventEnd && event->axis != axis) {
				KDAssert(eventsByAxisCtr < PointType::dim);
				eventsByAxis[eventsByAxisCtr++] = event;
			}

			/* The split plane can now be moved onto 't'. Accordingly, all planar
			   and ending primitives are removed from the right side */
			numRight[axis] -= numPlanar + numEnd;

			/* Calculate a score using the tree construction heuristic */
			if (EXPECT_TAKEN(pos > nodeAABB.min[axis] && pos < nodeAABB.max[axis])) {
				const SizeType nL = numLeft[axis], nR = numRight[axis];
				const Float nLF = (Float) nL, nRF = (Float) nR;

				std::pair<Float, Float> prob = tch(axis,
						pos - nodeAABB.min[axis],
						nodeAABB.max[axis] - pos);

				if (numPlanar == 0) {
					Float cost = m_traversalCost + m_queryCost
						* (prob.first * nLF + prob.second * nRF);
					if (nL == 0 || nR == 0)
						cost *= m_emptySpaceBonus;
					if (cost < bestSplit.cost) {
						bestSplit.pos = pos;
						bestSplit.axis = axis;
						bestSplit.cost = cost;
						bestSplit.numLeft = nL;
						bestSplit.numRight = nR;
					}
				} else {
					Float costPlanarLeft  = m_traversalCost
						+ m_queryCost * (prob.first * (nL+numPlanar) + prob.second * nRF);
					Float costPlanarRight = m_traversalCost
						+ m_queryCost * (prob.first * nLF + prob.second * (nR+numPlanar));

					if (nL + numPlanar == 0 || nR == 0)
						costPlanarLeft *= m_emptySpaceBonus;
					if (nL == 0 || nR + numPlanar == 0)
						costPlanarRight *= m_emptySpaceBonus;

					if (costPlanarLeft < bestSplit.cost ||
						costPlanarRight < bestSplit.cost) {
						bestSplit.pos = pos;
						bestSplit.axis = axis;

						if (costPlanarLeft < costPlanarRight) {
							bestSplit.cost = costPlanarLeft;
							bestSplit.numLeft = nL + numPlanar;
							bestSplit.numRight = nR;
							bestSplit.planarLeft = true;
						} else {
							bestSplit.cost = costPlanarRight;
							bestSplit.numLeft = nL;
							bestSplit.numRight = nR + numPlanar;
							bestSplit.planarLeft = false;
						}
					}
				}
			} else {
				#if defined(MTS_KD_DEBUG)
				if (m_clip && (pos < nodeAABB.min[axis]
							|| pos > nodeAABB.max[axis])) {
					/* When primitive clipping is active, this should  never happen! */
					KDLog(EError, "Internal error: edge event is out of bounds");
				}
				#endif
			}

			/* The split plane is moved past 't'. All prims,
				which were planar on 't', are moved to the left
				side. Also, starting prims are now also left of
				the split plane. */
			numLeft[axis] += numStart + numPlanar;
		}

#if defined(MTS_KD_DEBUG)
		/* Sanity checks. Everything should now be left of the split plane */
		for (int i=0; i<PointType::dim; ++i)
			KDAssert(numRight[i] == 0 && numLeft[i] == primCount);
		for (int i=1; i<PointType::dim; ++i)
			KDAssert(eventsByAxis[i]->axis == i && (eventsByAxis[i]-1)->axis == i-1);
#endif

		/* "Bad refines" heuristic from PBRT */
		if (bestSplit.cost >= leafCost) {
			if ((bestSplit.cost > 4 * leafCost && primCount < 16)
				|| badRefines >= m_maxBadRefines
				|| bestSplit.cost == std::numeric_limits<Float>::infinity()) {
				createLeaf(ctx, node, eventStart, eventEnd, primCount);
				return leafCost;
			}
			++badRefines;
		}

		/* ==================================================================== */
		/*                      Primitive Classification                        */
		/* ==================================================================== */

		ClassificationStorage &storage = ctx.classStorage;

		/* Initially mark all prims as being located on both sides */
		for (EdgeEvent *event = eventsByAxis[bestSplit.axis];
			 event < eventEnd && event->axis == bestSplit.axis; ++event)
			storage.set(event->index, EBothSides);

		SizeType primsLeft = 0, primsRight = 0, primsBoth = primCount;
		/* Sweep over all edge events and classify the primitives wrt. the split */
		for (EdgeEvent *event = eventsByAxis[bestSplit.axis];
			 event < eventEnd && event->axis == bestSplit.axis; ++event) {
			if (event->type == EdgeEvent::EEdgeEnd && event->pos <= bestSplit.pos) {
				/* The primitive's interval ends before or on the split plane
				   -> classify to the left side */
				KDAssert(storage.get(event->index) == EBothSides);
				storage.set(event->index, ELeftSide);
				primsBoth--;
				primsLeft++;
			} else if (event->type == EdgeEvent::EEdgeStart
					&& event->pos >= bestSplit.pos) {
				/* The primitive's interval starts after or on the split plane
				   -> classify to the right side */
				KDAssert(storage.get(event->index) == EBothSides);
				storage.set(event->index, ERightSide);
				primsBoth--;
				primsRight++;
			} else if (event->type == EdgeEvent::EEdgePlanar) {
				/* If the planar primitive is not on the split plane, the
				   classification is easy. Otherwise, place it on the side with
				   the lower cost */
				KDAssert(storage.get(event->index) == EBothSides);
				if (event->pos < bestSplit.pos || (event->pos == bestSplit.pos
						&& bestSplit.planarLeft)) {
					storage.set(event->index, ELeftSide);
					primsBoth--;
					primsLeft++;
				} else if (event->pos > bestSplit.pos
					   || (event->pos == bestSplit.pos && !bestSplit.planarLeft)) {
					storage.set(event->index, ERightSide);
					primsBoth--;
					primsRight++;
				} else {
					KDAssertEx(false, "Internal error!");
				}
			}
		}

		/* Some sanity checks */
		KDAssert(primsLeft + primsRight + primsBoth == primCount);
		KDAssert(primsLeft + primsBoth == bestSplit.numLeft);
		KDAssert(primsRight + primsBoth == bestSplit.numRight);

		OrderedChunkAllocator &leftAlloc = ctx.leftAlloc,
			&rightAlloc = ctx.rightAlloc;

		EdgeEvent *leftEventsStart, *rightEventsStart;
		if (isLeftChild) {
			leftEventsStart = eventStart;
			rightEventsStart = rightAlloc.allocate<EdgeEvent>(bestSplit.numRight * 2 * PointType::dim);
		} else {
			leftEventsStart = leftAlloc.allocate<EdgeEvent>(bestSplit.numLeft * 2 * PointType::dim);
			rightEventsStart = eventStart;
		}

		EdgeEvent *leftEventsEnd = leftEventsStart, *rightEventsEnd = rightEventsStart;

		AABBType leftNodeAABB = nodeAABB, rightNodeAABB = nodeAABB;
		leftNodeAABB.max[bestSplit.axis] = bestSplit.pos;
		rightNodeAABB.min[bestSplit.axis] = bestSplit.pos;

		SizeType prunedLeft = 0, prunedRight = 0;

		/* ==================================================================== */
		/*                            Partitioning                              */
		/* ==================================================================== */

		if (m_clip) {
			EdgeEvent
			  *leftEventsTempStart = leftAlloc.allocate<EdgeEvent>(primsLeft * 2 * PointType::dim),
			  *rightEventsTempStart = rightAlloc.allocate<EdgeEvent>(primsRight * 2 * PointType::dim),
			  *newEventsLeftStart = leftAlloc.allocate<EdgeEvent>(primsBoth * 2 * PointType::dim),
			  *newEventsRightStart = rightAlloc.allocate<EdgeEvent>(primsBoth * 2 * PointType::dim);

			EdgeEvent *leftEventsTempEnd = leftEventsTempStart,
					*rightEventsTempEnd = rightEventsTempStart,
					*newEventsLeftEnd = newEventsLeftStart,
					*newEventsRightEnd = newEventsRightStart;

			for (EdgeEvent *event = eventStart; event<eventEnd; ++event) {
				int classification = storage.get(event->index);

				if (classification == ELeftSide) {
					/* Left-only primitive. Move to the left list and advance */
					*leftEventsTempEnd++ = *event;
				} else if (classification == ERightSide) {
					/* Right-only primitive. Move to the right list and advance */
					*rightEventsTempEnd++ = *event;
				} else if (classification == EBothSides) {
					/* The primitive overlaps the split plane. Re-clip and
					   generate new events for each side */
					const IndexType index = event->index;

					AABBType clippedLeft = cast()->getClippedAABB(index, leftNodeAABB);
					AABBType clippedRight = cast()->getClippedAABB(index, rightNodeAABB);

					KDAssert(leftNodeAABB.contains(clippedLeft));
					KDAssert(rightNodeAABB.contains(clippedRight));

					if (clippedLeft.isValid() && clippedLeft.getSurfaceArea() > 0) {
						for (int axis=0; axis<PointType::dim; ++axis) {
							float min = clippedLeft.min[axis],
								  max = clippedLeft.max[axis];

							if (min == max) {
								*newEventsLeftEnd++ = EdgeEvent(
										EdgeEvent::EEdgePlanar,
										axis, min, index);
							} else {
								*newEventsLeftEnd++ = EdgeEvent(
										EdgeEvent::EEdgeStart,
										axis, min, index);
								*newEventsLeftEnd++ = EdgeEvent(
										EdgeEvent::EEdgeEnd,
										axis, max, index);
							}
						}
					} else {
						prunedLeft++;
					}

					if (clippedRight.isValid() && clippedRight.getSurfaceArea() > 0) {
						for (int axis=0; axis<PointType::dim; ++axis) {
							float min = clippedRight.min[axis],
								  max = clippedRight.max[axis];

							if (min == max) {
								*newEventsRightEnd++ = EdgeEvent(
										EdgeEvent::EEdgePlanar,
										axis, min, index);
							} else {
								*newEventsRightEnd++ = EdgeEvent(
										EdgeEvent::EEdgeStart,
										axis, min, index);
								*newEventsRightEnd++ = EdgeEvent(
										EdgeEvent::EEdgeEnd,
										axis, max, index);
							}
						}
					} else {
						prunedRight++;
					}

					/* Mark this primitive as processed so that clipping
						is only done once */
					storage.set(index, EBothSidesProcessed);
				}
			}

			KDAssert((SizeType) (leftEventsTempEnd - leftEventsTempStart) <= primsLeft * 2 * PointType::dim);
			KDAssert((SizeType) (rightEventsTempEnd - rightEventsTempStart) <= primsRight * 2 * PointType::dim);
			KDAssert((SizeType) (newEventsLeftEnd - newEventsLeftStart) <= primsBoth * 2 * PointType::dim);
			KDAssert((SizeType) (newEventsRightEnd - newEventsRightStart) <= primsBoth * 2 * PointType::dim);
			ctx.pruned += prunedLeft + prunedRight;

			/* Sort the events from overlapping prims */
			std::sort(newEventsLeftStart, newEventsLeftEnd, EdgeEventOrdering());
			std::sort(newEventsRightStart, newEventsRightEnd, EdgeEventOrdering());

			/* Merge the left list */
			leftEventsEnd = std::merge(leftEventsTempStart,
				leftEventsTempEnd, newEventsLeftStart, newEventsLeftEnd,
				leftEventsStart, EdgeEventOrdering());

			/* Merge the right list */
			rightEventsEnd = std::merge(rightEventsTempStart,
				rightEventsTempEnd, newEventsRightStart, newEventsRightEnd,
				rightEventsStart, EdgeEventOrdering());

			/* Release temporary memory */
			leftAlloc.release(newEventsLeftStart);
			leftAlloc.release(leftEventsTempStart);
			rightAlloc.release(newEventsRightStart);
			rightAlloc.release(rightEventsTempStart);
		} else {
			for (EdgeEvent *event = eventStart; event < eventEnd; ++event) {
				int classification = storage.get(event->index);

				if (classification == ELeftSide) {
					/* Left-only primitive. Move to the left list and advance */
					*leftEventsEnd++ = *event;
				} else if (classification == ERightSide) {
					/* Right-only primitive. Move to the right list and advance */
					*rightEventsEnd++ = *event;
				} else if (classification == EBothSides) {
					/* The primitive overlaps the split plane. Its edge events
					   must be added to both lists. */
					*leftEventsEnd++ = *event;
					*rightEventsEnd++ = *event;
				}
			}
			KDAssert((SizeType) (leftEventsEnd - leftEventsStart) <= bestSplit.numLeft * 2 * PointType::dim);
			KDAssert((SizeType) (rightEventsEnd - rightEventsStart) <= bestSplit.numRight * 2 * PointType::dim);
		}

		/* Shrink the edge event storage now that we know exactly how
		   many are on each side */
		ctx.leftAlloc.shrinkAllocation(leftEventsStart,
				leftEventsEnd - leftEventsStart);

		ctx.rightAlloc.shrinkAllocation(rightEventsStart,
				rightEventsEnd - rightEventsStart);

		/* ==================================================================== */
		/*                              Recursion                               */
		/* ==================================================================== */

		KDNode *children = ctx.nodes.allocate(2);

		SizeType nodePosBeforeSplit = (SizeType) ctx.nodes.size();
		SizeType indexPosBeforeSplit = (SizeType) ctx.indices.size();
		SizeType leafNodeCountBeforeSplit = ctx.leafNodeCount;
		SizeType nonemptyLeafNodeCountBeforeSplit = ctx.nonemptyLeafNodeCount;
		SizeType innerNodeCountBeforeSplit = ctx.innerNodeCount;

		if (!node->initInnerNode(bestSplit.axis, bestSplit.pos, children-node)) {
			LockGuard lock(m_indirectionLock);
			SizeType indirectionIdx = (SizeType) m_indirections.size();
			m_indirections.push_back(children);
			/* Unable to store relative offset -- create an indirection
			   table entry */
			node->initIndirectionNode(bestSplit.axis, bestSplit.pos,
					indirectionIdx);
		}
		ctx.innerNodeCount++;

		Float leftCost = buildTree(ctx, depth+1, children,
				leftNodeAABB, leftEventsStart, leftEventsEnd,
				bestSplit.numLeft - prunedLeft, true, badRefines);

		Float rightCost = buildTree(ctx, depth+1, children+1,
				rightNodeAABB, rightEventsStart, rightEventsEnd,
				bestSplit.numRight - prunedRight, false, badRefines);

		std::pair<Float, Float> prob = tch(bestSplit.axis,
			bestSplit.pos - nodeAABB.min[bestSplit.axis],
			nodeAABB.max[bestSplit.axis] - bestSplit.pos);

		/* Compute the final cost given the updated cost
		   values received from the children */
		Float finalCost = m_traversalCost +
			(prob.first * leftCost + prob.second * rightCost);

		/* Release the index lists not needed by the children anymore */
		if (isLeftChild)
			ctx.rightAlloc.release(rightEventsStart);
		else
			ctx.leftAlloc.release(leftEventsStart);

		/* ==================================================================== */
	    /*                           Final decision                             */
	    /* ==================================================================== */

		if (!m_retract || finalCost < primCount * m_queryCost) {
			return finalCost;
		} else {
			/* In the end, splitting didn't help to reduce the SAH cost.
			   Tear up everything below this node and create a leaf */
			ctx.nodes.resize(nodePosBeforeSplit);
			ctx.retractedSplits++;
			ctx.leafNodeCount = leafNodeCountBeforeSplit;
			ctx.nonemptyLeafNodeCount = nonemptyLeafNodeCountBeforeSplit;
			ctx.innerNodeCount = innerNodeCountBeforeSplit;
			createLeafAfterRetraction(ctx, node, indexPosBeforeSplit);
			return leafCost;
		}

		return bestSplit.cost;
	}

	/**
	 * \brief Min-max binning as described in
	 * "Highly Parallel Fast KD-tree Construction for Interactive
	 *  Ray Tracing of Dynamic Scenes"
	 * by M. Shevtsov, A. Soupikov and A. Kapustin
	 */
	struct MinMaxBins {
		MinMaxBins(SizeType nBins) : m_binCount(nBins) {
			m_minBins = new SizeType[m_binCount*PointType::dim];
			m_maxBins = new SizeType[m_binCount*PointType::dim];
		}

		~MinMaxBins() {
			delete[] m_minBins;
			delete[] m_maxBins;
		}

		/**
		 * \brief Prepare to bin for the specified bounds
		 */
		void setAABB(const AABBType &aabb) {
			for (int axis=0; axis<PointType::dim; ++axis) {
				float min = math::castflt_down(aabb.min[axis]);
				float max = math::castflt_up(aabb.max[axis]);
				m_min[axis] = min;
				m_aabb.min[axis] = min;
				m_aabb.max[axis] = max;
				m_binSize[axis] = (max-min) / m_binCount;
				m_invBinSize[axis] = 1 / m_binSize[axis];
			}
		}

		/// Compute the bin location for a given position and axis
		inline IndexType computeIndex(float pos, int axis) {
			return (IndexType) std::min((float) (m_binCount-1), std::max(0.0f, (pos - m_min[axis]) * m_invBinSize[axis]));
		}

		/**
		 * \brief Run min-max binning
		 *
		 * \param derived Derived class to be used to determine the AABB for
		 *     a given list of primitives
		 * \param indices Primitive indirection list
		 * \param primCount Specifies the length of \a indices
		 */
		void bin(const Derived *derived, IndexType *indices,
				SizeType primCount) {
			m_primCount = primCount;
			memset(m_minBins, 0, sizeof(SizeType) * PointType::dim * m_binCount);
			memset(m_maxBins, 0, sizeof(SizeType) * PointType::dim * m_binCount);

			for (SizeType i=0; i<m_primCount; ++i) {
				const AABBType aabb = derived->getAABB(indices[i]);
				for (int axis=0; axis<PointType::dim; ++axis) {
					m_minBins[axis * m_binCount + computeIndex(math::castflt_down(aabb.min[axis]), axis)]++;
					m_maxBins[axis * m_binCount + computeIndex(math::castflt_up  (aabb.max[axis]), axis)]++;
				}
			}
		}

		/**
		 * \brief Evaluate the tree construction heuristic at each bin boundary
		 * and return the minimizer for the given cost constants. Min-max
		 * binning uses no "empty space bonus" since it cannot create such
		 * splits.
		 */
		SplitCandidate minimizeCost(Float traversalCost, Float queryCost) {
			TreeConstructionHeuristic tch(m_aabb);
			SplitCandidate candidate;
			int binIdx = 0;

			for (int axis=0; axis<PointType::dim; ++axis) {
				SizeType numLeft = 0, numRight = m_primCount;
				Float leftWidth = 0, rightWidth = m_aabb.max[axis] - m_aabb.min[axis];
				const float binSize = m_binSize[axis];

				for (int i=0; i<m_binCount-1; ++i) {
					numLeft  += m_minBins[binIdx];
					numRight -= m_maxBins[binIdx];
					leftWidth += binSize;
					rightWidth -= binSize;
					std::pair<Float, Float> prob =
						tch(axis, leftWidth, rightWidth);

					Float cost = traversalCost + queryCost
						* (prob.first * numLeft + prob.second * numRight);

					if (cost < candidate.cost) {
						candidate.cost = cost;
						candidate.axis = axis;
						candidate.numLeft = numLeft;
						candidate.numRight = numRight;
						candidate.leftBin = i;
					}

					binIdx++;
				}
				binIdx++;
			}

			KDAssert(candidate.cost != std::numeric_limits<Float>::infinity() &&
			         candidate.leftBin >= 0);

			return candidate;
		}

		struct Partition {
			AABBType left;
			IndexType *leftIndices;
			AABBType right;
			IndexType *rightIndices;

			Partition(const AABBType &left, IndexType *leftIndices,
				const AABBType &right, IndexType *rightIndices) :
				left(left), leftIndices(leftIndices), right(right),
				rightIndices(rightIndices) { }
		};

		/**
		 * \brief Given a suitable split candiate, compute tight bounding
		 * boxes for the left and right subtrees and return associated
		 * primitive lists.
		 */
		Partition partition(
				BuildContext &ctx, const Derived *derived, IndexType *primIndices,
				SplitCandidate &split, bool isLeftChild, Float traversalCost,
				Float queryCost) {
			SizeType numLeft = 0, numRight = 0;
			AABBType leftBounds, rightBounds;
			const int axis = split.axis;

			IndexType *leftIndices, *rightIndices;
			if (isLeftChild) {
				OrderedChunkAllocator &rightAlloc = ctx.rightAlloc;
				leftIndices = primIndices;
				rightIndices = rightAlloc.allocate<IndexType>(split.numRight);
			} else {
				OrderedChunkAllocator &leftAlloc = ctx.leftAlloc;
				leftIndices = leftAlloc.allocate<IndexType>(split.numLeft);
				rightIndices = primIndices;
			}

			for (SizeType i=0; i<m_primCount; ++i) {
				const IndexType primIndex = primIndices[i];
				const AABBType aabb = derived->getAABB(primIndex);
				int startIdx = computeIndex(math::castflt_down(aabb.min[axis]), axis);
				int endIdx   = computeIndex(math::castflt_up  (aabb.max[axis]), axis);

				if (endIdx <= split.leftBin) {
					KDAssert(numLeft < split.numLeft);
					leftBounds.expandBy(aabb);
					leftIndices[numLeft++] = primIndex;
				} else if (startIdx > split.leftBin) {
					KDAssert(numRight < split.numRight);
					rightBounds.expandBy(aabb);
					rightIndices[numRight++] = primIndex;
				} else {
					leftBounds.expandBy(aabb);
					rightBounds.expandBy(aabb);
					KDAssert(numLeft < split.numLeft);
					KDAssert(numRight < split.numRight);
					leftIndices[numLeft++] = primIndex;
					rightIndices[numRight++] = primIndex;
				}
			}
			leftBounds.clip(m_aabb);
			rightBounds.clip(m_aabb);
			split.pos = m_min[axis] + m_binSize[axis] * (split.leftBin + 1);
			leftBounds.max[axis] = std::min(leftBounds.max[axis], (Float) split.pos);
			rightBounds.min[axis] = std::max(rightBounds.min[axis], (Float) split.pos);

			KDAssert(numLeft == split.numLeft);
			KDAssert(numRight == split.numRight);

			/// Release the unused memory regions
			if (isLeftChild)
				ctx.leftAlloc.shrinkAllocation(leftIndices, split.numLeft);
			else
				ctx.rightAlloc.shrinkAllocation(rightIndices, split.numRight);

			return Partition(leftBounds, leftIndices,
				rightBounds, rightIndices);
		}
	private:
		SizeType *m_minBins;
		SizeType *m_maxBins;
		SizeType m_primCount;
		int m_binCount;
		float m_min[PointType::dim];
		float m_binSize[PointType::dim];
		float m_invBinSize[PointType::dim];
		AABBType m_aabb;
	};

protected:
	IndexType *m_indices;
	Float m_traversalCost;
	Float m_queryCost;
	Float m_emptySpaceBonus;
	bool m_clip, m_retract, m_parallelBuild;
	SizeType m_maxDepth;
	SizeType m_stopPrims;
	SizeType m_maxBadRefines;
	SizeType m_exactPrimThreshold;
	SizeType m_minMaxBins;
	SizeType m_nodeCount;
	SizeType m_indexCount;
	std::vector<TreeBuilder *> m_builders;
	std::vector<KDNode *> m_indirections;
	ref<Mutex> m_indirectionLock;
	BuildInterface m_interface;
};

#if defined(_MSC_VER)
/* Revert back to fast / non-strict IEEE 754
   floating point computations */
MTS_NAMESPACE_END
#pragma float_control(precise, off)
MTS_NAMESPACE_BEGIN
#endif

template <typename AABBType>
	Class *KDTreeBase<AABBType>::m_theClass
		= new Class("KDTreeBase", true, "Object");

template <typename AABBType>
	const Class *KDTreeBase<AABBType>::getClass() const {
	return m_theClass;
}


// Explicit instantiation to avoid linking error outside the library
#ifdef _MSC_VER
template class MTS_EXPORT_RENDER KDTreeBase<AABB>;
#endif

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_GKDTREE_H_ */

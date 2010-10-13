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
#include <mitsuba/core/random.h>
#include <mitsuba/core/lock.h>
#include <boost/static_assert.hpp>
#include <boost/tuple/tuple.hpp>
#include <stack>

/// Activate lots of extra checks
#define MTS_KD_DEBUG 1

/** Compile-time KD-tree depth limit. Allows to put certain
    data structures on the stack */
#define MTS_KD_MAXDEPTH 48

/// Min-max bin count
#define MTS_KD_MINMAX_BINS 128

/// OrderedChunkAllocator: don't create chunks smaller than 512 KiB
#define MTS_KD_MIN_ALLOC 512*1024

/// Allocate nodes & index lists in blocks of 512 KiB
#define MTS_KD_BLOCKSIZE_KD  (512*1024/sizeof(KDNode))
#define MTS_KD_BLOCKSIZE_IDX (512*1024/sizeof(uint32_t))

/// 64 byte temporary storage for intersection computations 
#define MTS_KD_INTERSECTION_TEMP 64

/// Use a simple hashed 8-entry mailbox per thread
//#define MTS_KD_MAILBOX_ENABLED 1
//#define MTS_KD_MAILBOX_SIZE 8
//#define MTS_KD_MAILBOX_MASK (MTS_KD_MAILBOX_SIZE-1)

#if defined(MTS_KD_DEBUG)
#define KDAssert(expr) Assert(expr)
#define KDAssertEx(expr, text) AssertEx(expr, text)
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
 * memory will be released in the same order, in which it is 
 * allocated. This makes it possible to create an implementation
 * with a very low memory overhead. Note that no locking is done, 
 * hence each thread will need its own allocator.
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
	template <typename T> T * __restrict__ allocate(size_t size) {
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
	inline T * __restrict__ allocate(size_t size) {
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
	ClassificationStorage(size_t size = 0) : m_buffer(NULL), m_bufferSize(0) { }

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
 * \brief SAH KD-tree acceleration data structure for fast ray-object 
 * intersection computations.
 *
 * The code in this class is a fully generic kd-tree implementation, which
 * can theoretically support any kind of shape. However, subclasses still
 * need to provide the following signatures for a functional implementation:
 *
 * /// Return the total number of primitives
 * inline size_type getPrimitiveCount() const;
 *
 * /// Return the axis-aligned bounding box of a certain primitive
 * inline AABB getAABB(index_type primIdx) const;
 *
 * /// Return the AABB of a primitive when clipped to another AABB
 * inline AABB getClippedAABB(index_type primIdx, const AABB &aabb) const;
 *
 * /// Check whether a primitive is intersected by the given ray. Some
 * /// temporary space is supplied, which can be used to store extra 
 * /// information about the intersection
 * EIntersectionResult intersect(const Ray &ray, index_type idx, Float mint, 
 *		Float maxt, Float &t, void *tmp);
 *
 * This class follows the "Curiously recurring template" design pattern so that
 * the above functions can be inlined (no virtual calls will be necessary!).
 *
 * The kd-tree construction algorithm creates 'perfect split' trees as 
 * outlined in the paper "On Building fast kd-Trees for Ray Tracing, and on
 * doing that in O(N log N)" by Ingo Wald and Vlastimil Havran. 
 * For polygonal meshes, the involved Sutherland-Hodgman iterations can be 
 * quite expensive in terms of the overall construction time. The \ref setClip
 * method can be used to deactivate perfect splits at the cost of a 
 * lower-quality tree.
 *
 * Because the O(N log N) construction algorithm tends to cause many
 * incoherent memory accesses, a fast approximate technique (Min-max binning)
 * is used near the top of the tree, which significantly reduces cache misses.
 * Once the input data has been narrowed down to a reasonable amount, the
 * implementation switches over to the O(N log N) builder. When multiple
 * processors are available, the build process runs in parallel.
 *
 * This implementation also provides an optimized ray traversal algorithm 
 * (TA^B_{rec}), which is explained in Vlastimil Havran's PhD thesis 
 * "Heuristic Ray Shooting Algorithms". 
 *
 * \author Wenzel Jakob
 */
template <typename Derived> class GenericKDTree : public Object {
protected:
	struct KDNode;
	struct EdgeEvent;
	struct EdgeEventOrdering;

public:
	/// Index number format (max 2^32 prims)
	typedef uint32_t index_type;

	/// Size number format
	typedef uint32_t size_type;

	/// Documents the possible outcomes of a single ray-primitive intersection
	enum EIntersectionResult {
		/// An intersection was found on the specified interval
		EYes = 0,

		/** While an intersection could not be found on the specified 
		    interval, the primitive might still intersect the ray if a
			larger interval was considered */
		ENo,

		/// The primitive will never intersect the ray in question
		ENever
	};

	/**
	 * \brief Create a new kd-tree instance initialized with 
	 * the default parameters.
	 */
	GenericKDTree() : m_nodes(NULL), m_indices(NULL) {
		m_traversalCost = 10;
		m_intersectionCost = 17;
		m_emptySpaceBonus = 0.9f;
		m_clip = true;
		m_stopPrims = 4;
		m_maxBadRefines = 3;
		m_exactPrimThreshold = 65536;
		m_maxDepth = 0;
		m_retract = true;
		m_parallelBuild = true;
	}

	/**
	 * \brief Release all memory
	 */
	virtual ~GenericKDTree() {
		if (m_indices)
			delete[] m_indices;
		if (m_nodes)
			freeAligned(m_nodes);
	}

	/**
	 * \brief Set the traversal cost used by the surface area heuristic
	 */
	inline void setTraversalCost(Float traversalCost) {
		m_traversalCost = traversalCost;
	}

	/**
	 * \brief Return the traversal cost used by the surface area heuristic
	 */
	inline Float getTraversalCost() const {
		return m_traversalCost;
	}

	/**
	 * \brief Set the intersection cost used by the surface area heuristic
	 */
	inline void setIntersectionCost(Float intersectionCost) {
		m_intersectionCost = intersectionCost;
	}

	/**
	 * \brief Return the intersection cost used by the surface area heuristic
	 */
	inline Float getIntersectionCost() const {
		return m_intersectionCost;
	}

	/**
	 * \brief Set the bonus factor for empty space used by the
	 * surface area heuristic
	 */
	inline void setEmptySpaceBonus(Float emptySpaceBonus) {
		m_emptySpaceBonus = emptySpaceBonus;
	}

	/**
	 * \brief Return the bonus factor for empty space used by the
 	 * surface area heuristic
	 */
	inline Float getEmptySpaceBonus() const {
		return m_emptySpaceBonus;
	}

	/**
	 * \brief Set the maximum tree depth (0 = use heuristic)
	 */
	inline void setMaxDepth(size_type maxDepth) {
		m_maxDepth = maxDepth;
	}

	/**
	 * \brief Return maximum tree depth (0 = use heuristic)
	 */
	inline size_type getMaxDepth() const {
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
	inline void setMaxBadRefines(size_type maxBadRefines) {
		m_maxBadRefines = maxBadRefines;
	}

	/**
	 * \brief Return the number of bad refines allowed to happen
	 * in succession before a leaf node will be created.
	 */
	inline size_type getMaxBadRefines() const {
		return m_maxBadRefines;
	}

	/**
	 * \brief Set the number of primitives, at which recursion will
	 * stop when building the tree.
	 */
	inline void setStopPrims(size_type stopPrims) {
		m_stopPrims = stopPrims;
	}

	/**
	 * \brief Return the number of primitives, at which recursion will
	 * stop when building the tree.
	 */
	inline size_type getStopPrims() const {
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
	 * \brief Specify the number of primitives, at which the builder will switch
	 * from (approximate) Min-Max binning to the accurate O(n log n) SAH-based 
	 * optimization method.
	 */
	inline void setExactPrimitiveThreshold(size_type exactPrimThreshold) {
		m_exactPrimThreshold = exactPrimThreshold;
	}

	/**
	 * \brief Return the number of primitives, at which the builder will switch
	 * from (approximate) Min-Max binning to the accurate O(n log n) SAH-based 
	 * optimization method.
	 */
	inline size_type getExactPrimitiveThreshold() const {
		return m_exactPrimThreshold;
	}

	/**
	 * \brief Return whether or not the kd-tree has been built
	 */
	inline bool isBuilt() const {
		return m_nodes != NULL;
	}

	/**
	 * \brief Return an axis-aligned bounding box containing all primitives
	 */
	inline const AABB &getAABB() const { return m_aabb; }

	/**
	 * \brief Return an bounding sphere containing all primitives
	 */
	inline const BSphere &getBSphere() const { return m_bsphere; }

	/**
	 * \brief Empirically find the best traversal and intersection cost values
	 *
	 * This is done by running the traversal code on random rays and fitting 
	 * the SAH cost model to the collected statistics.
	 */
	void findCosts(Float &traversalCost, Float &intersectionCost);

	MTS_DECLARE_CLASS()
protected:
	/**
	 * \brief Build a KD-tree over the supplied geometry
	 *
	 * To be called by the subclass.
	 */
	void buildInternal() {
		if (isBuilt()) 
			Log(EError, "The kd-tree has already been built!");

		size_type primCount = cast()->getPrimitiveCount();
		BuildContext ctx(primCount);

		/* Establish an ad-hoc depth cutoff value (Formula from PBRT) */
		if (m_maxDepth == 0)
			m_maxDepth = (int) (8 + 1.3f * log2i(primCount));
		m_maxDepth = std::min(m_maxDepth, (size_type) MTS_KD_MAXDEPTH);

		Log(EDebug, "Creating a preliminary index list (%s)", 
			memString(primCount * sizeof(index_type)).c_str());

		OrderedChunkAllocator &leftAlloc = ctx.leftAlloc;
		index_type *indices = leftAlloc.allocate<index_type>(primCount);

		ref<Timer> timer = new Timer();
		m_aabb.reset();
		for (index_type i=0; i<primCount; ++i) {
			m_aabb.expandBy(cast()->getAABB(i));
			indices[i] = i;
		}

		Log(EDebug, "Computed scene bounds in %i ms", timer->getMilliseconds());
		Log(EDebug, "");

		Log(EDebug, "kd-tree configuration:");
		Log(EDebug, "   Traversal cost           : %.2f", m_traversalCost);
		Log(EDebug, "   Intersection cost        : %.2f", m_intersectionCost);
		Log(EDebug, "   Empty space bonus        : %.2f", m_emptySpaceBonus);
		Log(EDebug, "   Max. tree depth          : %i", m_maxDepth);
		Log(EDebug, "   Scene bounding box (min) : %s", m_aabb.min.toString().c_str());
		Log(EDebug, "   Scene bounding box (max) : %s", m_aabb.max.toString().c_str());
		Log(EDebug, "   Min-max bins             : %i", MTS_KD_MINMAX_BINS);
		Log(EDebug, "   Greedy SAH optimization  : use for <= %i primitives", 
				m_exactPrimThreshold);
		Log(EDebug, "   Perfect splits           : %s", m_clip ? "yes" : "no");
		Log(EDebug, "   Retract bad splits       : %s", m_retract ? "yes" : "no");
		Log(EDebug, "   Stopping primitive count : %i", m_stopPrims);
		Log(EDebug, "");

		size_type procCount = getProcessorCount();
		if (procCount == 1)
			m_parallelBuild = false;

		if (m_parallelBuild) {
			m_builders.resize(procCount);
			for (size_type i=0; i<procCount; ++i) {
				m_builders[i] = new SAHTreeBuilder(i, this);
				m_builders[i]->incRef();
				m_builders[i]->start();
			}
		}

		Log(EInfo, "Constructing a SAH kd-tree (%i primitives) ..", primCount);

		m_indirectionLock = new Mutex();
		KDNode *prelimRoot = ctx.nodes.allocate(1);
		buildTreeMinMax(ctx, 1, prelimRoot, m_aabb, m_aabb, 
				indices, primCount, true, 0);
		ctx.leftAlloc.release(indices);

		KDAssert(ctx.leftAlloc.used() == 0);
		KDAssert(ctx.rightAlloc.used() == 0);

		if (m_parallelBuild) {
			m_interface.done = true;
			m_interface.cond->broadcast();
			for (size_type i=0; i<m_builders.size(); ++i) 
				m_builders[i]->join();
		}

		Log(EInfo, "Finished -- took %i ms.", timer->getMilliseconds());
		Log(EDebug, "");

		Log(EDebug, "Temporary memory statistics:");
		Log(EDebug, "   Classification storage : %s", 
				memString((ctx.classStorage.size() * (1+procCount))).c_str());
		Log(EDebug, "   Indirection entries    : " SIZE_T_FMT " (%s)", 
				m_indirections.size(), memString(m_indirections.capacity()
				* sizeof(KDNode *)).c_str());

		Log(EDebug, "   Main thread:");
		ctx.printStats();
		size_t totalUsage = m_indirections.capacity() 
			* sizeof(KDNode *) + ctx.size();

		/// Clean up event lists and print statistics
		ctx.leftAlloc.cleanup();
		ctx.rightAlloc.cleanup();
		for (size_type i=0; i<m_builders.size(); ++i) {
			Log(EDebug, "   Worker thread %i:", i+1);
			BuildContext &subCtx = m_builders[i]->getContext();
			subCtx.printStats();
			totalUsage += subCtx.size();
			subCtx.leftAlloc.cleanup();
			subCtx.rightAlloc.cleanup();
			ctx.accumulateStatisticsFrom(subCtx);
		}
		Log(EDebug, "   Total: %s", memString(totalUsage).c_str());

		Log(EDebug, "");
		timer->reset();
		Log(EDebug, "Optimizing memory layout ..");

		std::stack<boost::tuple<const KDNode *, KDNode *, 
				const BuildContext *, AABB> > stack;

		Float expTraversalSteps = 0;
		Float expLeavesVisited = 0;
		Float expPrimitivesIntersected = 0;
		Float sahCost = 0;
		size_type nodePtr = 0, indexPtr = 0;
		size_type maxPrimsInLeaf = 0;
		const size_type primBucketCount = 16;
		size_type primBuckets[primBucketCount];
		memset(primBuckets, 0, sizeof(size_type)*primBucketCount);
		m_nodeCount = ctx.innerNodeCount + ctx.leafNodeCount;

		m_nodes = static_cast<KDNode *> (allocAligned(
				sizeof(KDNode) * m_nodeCount));
		m_indices = new index_type[ctx.primIndexCount];

		stack.push(boost::make_tuple(prelimRoot, &m_nodes[nodePtr++], &ctx, m_aabb));
		while (!stack.empty()) {
			const KDNode *node = boost::get<0>(stack.top());
			KDNode *target = boost::get<1>(stack.top());
			const BuildContext *context = boost::get<2>(stack.top());
			AABB aabb = boost::get<3>(stack.top());
			stack.pop();

			if (node->isLeaf()) {
				size_type primStart = node->getPrimStart(),
						  primEnd = node->getPrimEnd(),
						  primCount = primEnd-primStart;
				target->initLeafNode(indexPtr, primCount);

				Float sa = aabb.getSurfaceArea(), weightedSA = sa * primCount;
				expLeavesVisited += aabb.getSurfaceArea();
				expPrimitivesIntersected += weightedSA;
				sahCost += weightedSA * m_intersectionCost;
				if (primCount < primBucketCount)
					primBuckets[primCount]++;
				if (primCount > maxPrimsInLeaf)
					maxPrimsInLeaf = primCount;

				const BlockedVector<index_type, MTS_KD_BLOCKSIZE_IDX> &indices 
					= context->indices;
				for (size_type idx = primStart; idx<primEnd; ++idx)
					m_indices[indexPtr++] = indices[idx];
			} else {
				typename std::map<const KDNode *, index_type>::const_iterator it 
					= m_interface.threadMap.find(node);
				// Check if we're switching to a subtree built by a worker thread
				if (it != m_interface.threadMap.end()) 
					context = &m_builders[(*it).second]->getContext();

				Float sa = aabb.getSurfaceArea();
				expTraversalSteps += sa;
				sahCost += sa * m_traversalCost;

				const KDNode *left;
				if (EXPECT_TAKEN(!node->isIndirection()))
					left = node->getLeft();
				else 
					left = m_indirections[node->getIndirectionIndex()];

				KDNode *children = &m_nodes[nodePtr];
				nodePtr += 2;
				int axis = node->getAxis();
				float split = node->getSplit();
				bool result = target->initInnerNode(axis, split, children - target);
				if (!result)
					Log(EError, "Cannot represent relative pointer -- "
						"too many primitives?");

				Float tmp = aabb.min[axis];
				aabb.min[axis] = split;
				stack.push(boost::make_tuple(left+1, children+1, context, aabb));
				aabb.min[axis] = tmp;
				aabb.max[axis] = split;
				stack.push(boost::make_tuple(left, children, context, aabb));
			}
		}

		KDAssert(nodePtr == ctx.innerNodeCount + ctx.leafNodeCount);
		KDAssert(indexPtr == ctx.primIndexCount);

		Log(EDebug, "Finished -- took %i ms.", timer->getMilliseconds());

		/* Free some more memory */
		ctx.nodes.clear();
		ctx.indices.clear();
		for (size_type i=0; i<m_builders.size(); ++i) {
			BuildContext &subCtx = m_builders[i]->getContext();
			subCtx.nodes.clear();
			subCtx.indices.clear();
		}
		m_indirectionLock = NULL;
		std::vector<KDNode *>().swap(m_indirections);

		if (m_builders.size() > 0) {
			for (size_type i=0; i<m_builders.size(); ++i)
				m_builders[i]->decRef();
			m_builders.clear();
		}

		Log(EDebug, "");

		Float rootSA = m_aabb.getSurfaceArea();
		expTraversalSteps /= rootSA;
		expLeavesVisited /= rootSA;
		expPrimitivesIntersected /= rootSA;
		sahCost /= rootSA;

		/* Slightly enlarge the bounding box 
		   (necessary e.g. when the scene is planar) */
		m_aabb.min -= (m_aabb.max-m_aabb.min) * Epsilon
			+ Vector(Epsilon, Epsilon, Epsilon);
		m_aabb.max += (m_aabb.max-m_aabb.min) * Epsilon
			+ Vector(Epsilon, Epsilon, Epsilon);
		m_bsphere = m_aabb.getBSphere();

		Log(EDebug, "Structural kd-tree statistics:");
		Log(EDebug, "   Parallel work units       : " SIZE_T_FMT, 
				m_interface.threadMap.size());
		Log(EDebug, "   Node storage cost         : %s", 
				memString(nodePtr * sizeof(KDNode)).c_str());
		Log(EDebug, "   Index storage cost        : %s", 
				memString(indexPtr * sizeof(index_type)).c_str());
		Log(EDebug, "   Inner nodes               : %i", ctx.innerNodeCount);
		Log(EDebug, "   Leaf nodes                : %i", ctx.leafNodeCount);
		Log(EDebug, "   Nonempty leaf nodes       : %i", ctx.nonemptyLeafNodeCount);
		std::ostringstream oss;
		oss << "   Leaf node histogram       : ";
		for (size_type i=0; i<primBucketCount; i++) {
			oss << i << "(" << primBuckets[i] << ") ";
			if ((i+1)%4==0 && i+1<primBucketCount) {
				Log(EDebug, "%s", oss.str().c_str());
				oss.str("");
				oss << "                               ";
			}
		}
		Log(EDebug, "%s", oss.str().c_str());
		Log(EDebug, "");
		Log(EDebug, "Qualitative kd-tree statistics:");
		Log(EDebug, "   Retracted splits          : %i", ctx.retractedSplits);
		Log(EDebug, "   Pruned primitives         : %i", ctx.pruned);
		Log(EDebug, "   Largest leaf node         : %i primitives",
				maxPrimsInLeaf);
		Log(EDebug, "   Avg. prims/nonempty leaf  : %.2f", 
				ctx.primIndexCount / (Float) ctx.nonemptyLeafNodeCount);
		Log(EDebug, "   Expected traversals/ray   : %.2f", expTraversalSteps);
		Log(EDebug, "   Expected leaf visits/ray  : %.2f", expLeavesVisited);
		Log(EDebug, "   Expected prim. visits/ray : %.2f", expPrimitivesIntersected);
		Log(EDebug, "   Final SAH cost            : %.2f", sahCost);
		Log(EDebug, "");
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
		inline EdgeEvent(uint16_t type, int axis, float pos, index_type index)
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
		index_type index;
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
			return a.type < b.type;
		}
	};

	/**
	 * \brief Data type for split candidates computed by 
	 * the SAH optimization routines.
	 * */
	struct SplitCandidate {
		Float sahCost;
		float pos;
		int axis;
		size_type numLeft, numRight;
		bool planarLeft;

		inline SplitCandidate() :
			sahCost(std::numeric_limits<Float>::infinity()),
			pos(0), axis(0), numLeft(0), numRight(0), planarLeft(false) {
		}

		std::string toString() const {
			std::ostringstream oss;
			oss << "SplitCandidate[" << endl
				<< "  sahCost=" << sahCost << "," << endl
				<< "  pos=" << pos << "," << endl
				<< "  axis=" << axis << "," << endl
				<< "  numLeft=" << numLeft << "," << endl
				<< "  numRight=" << numRight << "," << endl
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
		BlockedVector<index_type, MTS_KD_BLOCKSIZE_IDX> indices;
		ClassificationStorage classStorage;

		size_type leafNodeCount;
		size_type nonemptyLeafNodeCount;
		size_type innerNodeCount;
		size_type primIndexCount;
		size_type retractedSplits;
		size_type pruned;

		BuildContext(size_type primCount) : classStorage(primCount) {
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
				+ indices.capacity() * sizeof(index_type)
				+ classStorage.size();
		}

		void printStats() {
			Log(EDebug, "      Left events   : " SIZE_T_FMT " chunks (%s)",
					leftAlloc.getChunkCount(), memString(leftAlloc.size()).c_str());
			Log(EDebug, "      Right events  : " SIZE_T_FMT " chunks (%s)",
					rightAlloc.getChunkCount(), memString(rightAlloc.size()).c_str());
			Log(EDebug, "      kd-tree nodes : " SIZE_T_FMT " entries, " SIZE_T_FMT 
					" blocks (%s)", nodes.size(), nodes.blockCount(), 
					memString(nodes.capacity() * sizeof(KDNode)).c_str());
			Log(EDebug, "      Indices       : " SIZE_T_FMT " entries, " SIZE_T_FMT 
					" blocks (%s)", indices.size(), indices.blockCount(), 
					memString(indices.capacity() * sizeof(index_type)).c_str());
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
	 * SAH kd-tree builder threads
	 */
	struct BuildInterface {
		/* Communcation */
		ref<Mutex> mutex;
		ref<ConditionVariable> cond, condJobTaken;
		std::map<const KDNode *, index_type> threadMap;
		bool done;

		/* Job description for building a subtree */
		int depth;
		KDNode *node;
		AABB nodeAABB;
		EdgeEvent *eventStart, *eventEnd;
		size_type primCount;
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
		inline void initIndirectionNode(int axis, float split, 
				uint32_t indirectionEntry) {
			inner.combined = EIndirectionMask 
				| ((uint32_t) indirectionEntry << 2)
				| axis;
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

		/// Return a string representation
		std::string toString() const {
			std::ostringstream oss;
			if (isLeaf()) {
				oss << "KDNode[leaf, primStart=" << getPrimStart() 
					<< ", primCount=" << getPrimEnd()-getPrimStart() << "]";
			} else {
				oss << "KDNode[interior, axis=" << getAxis() 
					<< ", split=" << getAxis() 
					<< ", leftOffset=" << ((inner.combined & EInnerOffsetMask) >> 2)
					<< "]";
			}
			return oss.str();
		}
	};

	BOOST_STATIC_ASSERT(sizeof(KDNode) == 8);

#if defined(MTS_KD_MAILBOX_ENABLED)
	/**
	 * \brief Hashed mailbox implementation
	 */
	struct HashedMailbox {
		inline HashedMailbox() {
			memset(entries, 0xFF, sizeof(index_type)*MTS_KD_MAILBOX_SIZE);
		}

		inline void put(index_type primIndex) {
			entries[primIndex & MTS_KD_MAILBOX_MASK] = primIndex;
		}

		inline bool contains(index_type primIndex) const {
			return entries[primIndex & MTS_KD_MAILBOX_MASK] == primIndex;
		}

		index_type entries[MTS_KD_MAILBOX_SIZE];
	};
#endif

	/**
	 * \brief SAH kd-tree builder thread
	 */
	class SAHTreeBuilder : public Thread {
	public:
		SAHTreeBuilder(index_type id, GenericKDTree *parent) 
			: Thread(formatString("bld%i", id)),
			m_id(id),
			m_parent(parent),
			m_context(parent->cast()->getPrimitiveCount()),
			m_interface(parent->m_interface) {
			setCritical(true);
		}

		~SAHTreeBuilder() {
			KDAssert(m_context.leftAlloc.used() == 0);
			KDAssert(m_context.rightAlloc.used() == 0);
		}

		void run() {
			OrderedChunkAllocator &leftAlloc = m_context.leftAlloc;
			while (true) {
				m_interface.mutex->lock();
				while (!m_interface.done && !m_interface.node)
					m_interface.cond->wait();
				if (m_interface.done) {
					m_interface.mutex->unlock();
					break;
				}
				int depth = m_interface.depth;
				KDNode *node = m_interface.node;
				AABB nodeAABB = m_interface.nodeAABB;
				size_t eventCount = m_interface.eventEnd - m_interface.eventStart;
				size_type primCount = m_interface.primCount;
				int badRefines = m_interface.badRefines;
				EdgeEvent *eventStart = leftAlloc.allocate<EdgeEvent>(eventCount),
						  *eventEnd = eventStart + eventCount;
				memcpy(eventStart, m_interface.eventStart, 
						eventCount * sizeof(EdgeEvent));
				m_interface.threadMap[node] = m_id;
				m_interface.node = NULL;
				m_interface.condJobTaken->signal();
				m_interface.mutex->unlock();

				std::sort(eventStart, eventEnd, EdgeEventOrdering());
				m_parent->buildTreeSAH(m_context, depth, node,
					nodeAABB, eventStart, eventEnd, primCount, true, badRefines);
				leftAlloc.release(eventStart);
			}
		}

		inline BuildContext &getContext() {
			return m_context;
		}

	private:
		index_type m_id;
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

	/**
	 * \brief Create an edge event list for a given list of primitives. 
	 *
	 * This is necessary when passing from Min-Max binning to the more 
	 * accurate SAH-based optimizier.
	 */
	boost::tuple<EdgeEvent *, EdgeEvent *, size_type> createEventList(
			OrderedChunkAllocator &alloc, const AABB &nodeAABB, 
			index_type *prims, size_type primCount) {
		size_type initialSize = primCount * 6, actualPrimCount = 0;
		EdgeEvent *eventStart = alloc.allocate<EdgeEvent>(initialSize);
		EdgeEvent *eventEnd = eventStart;

		for (size_type i=0; i<primCount; ++i) {
			index_type index = prims[i];
			AABB aabb;
			if (m_clip) {
				aabb = cast()->getClippedAABB(index, nodeAABB);
				if (!aabb.isValid() || aabb.getSurfaceArea() == 0)
					continue;
			} else {
				aabb = cast()->getAABB(index);
			}

			for (int axis=0; axis<3; ++axis) {
				float min = (float) aabb.min[axis], max = (float) aabb.max[axis];

				if (min == max) {
					*eventEnd++ = EdgeEvent(EdgeEvent::EEdgePlanar, axis, min, index);
				} else {
					*eventEnd++ = EdgeEvent(EdgeEvent::EEdgeStart, axis, min, index);
					*eventEnd++ = EdgeEvent(EdgeEvent::EEdgeEnd, axis, max, index);
				}
			}
			++actualPrimCount;
		}

		size_type newSize = eventEnd - eventStart;
		if (newSize != initialSize)
			alloc.shrinkAllocation<EdgeEvent>(eventStart, newSize);

		return boost::make_tuple(eventStart, eventEnd, actualPrimCount);
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
			EdgeEvent *eventEnd, size_type primCount) {
		node->initLeafNode(ctx.indices.size(), primCount);
		if (primCount > 0) {
			size_type seenPrims = 0;
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
	void createLeaf(BuildContext &ctx, KDNode *node, size_type *indices,
			size_type primCount) {
		node->initLeafNode(ctx.indices.size(), primCount);
		if (primCount > 0) {
			ctx.nonemptyLeafNodeCount++;
			for (size_type i=0; i<primCount; ++i) 
				ctx.indices.push_back(indices[i]);
			ctx.primIndexCount += primCount;
		}
		ctx.leafNodeCount++;
	}

	/**
	 * \brief Leaf node creation helper function. 
	 *
	 * Creates a unique index list by collapsing
	 * a subtree with a bad SAH cost.
	 *
	 * \param ctx 
	 *     Thread-specific build context containing allocators etc.
	 * \param node
	 *     KD-tree node entry to be filled
	 * \param start
	 *     Start pointer of the subtree indices
	 * \param primCount
	 *     Total primitive count for the current node
	 */
	void createLeafAfterRetraction(BuildContext &ctx, KDNode *node, 
			size_type start, size_type primCount) {
		node->initLeafNode(start, primCount);
		size_t actualCount = ctx.indices.size() - start;

		if (primCount > 0)
			ctx.nonemptyLeafNodeCount++;

		if (actualCount != primCount) {
			KDAssert(primCount > 0);
			OrderedChunkAllocator &alloc = ctx.leftAlloc;

			/* A temporary list is allocated to do the sorting (the indices
			   are not guaranteed to be contiguous in memory) */
			index_type *tempStart = alloc.allocate<index_type>(actualCount);
			index_type *tempEnd = tempStart, *ptr = tempStart;

			for (size_type i=start, end = start + actualCount; i<end; ++i)
				*tempEnd++ = ctx.indices[i];

			std::sort(tempStart, tempEnd, std::less<index_type>());

			for (size_type i=start, end = start + primCount; i<end; ++i) {
				ctx.indices[i] = *ptr++;
				while (ptr < tempEnd && *ptr == ctx.indices[i])
					++ptr;
			}

			ctx.indices.resize(start + primCount);
			ctx.primIndexCount = ctx.primIndexCount - actualCount + primCount;
			alloc.release(tempStart);
		}

		ctx.leafNodeCount++;
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
	 *     the hope that the SAH cost was significantly overestimated. The
	 *     counter makes sure that only a limited number of such splits can
	 *     happen in succession.
	 * \returns 
	 *     Final SAH cost of the node
	 */
	Float buildTreeMinMax(BuildContext &ctx, unsigned int depth, KDNode *node, 
			const AABB &nodeAABB, const AABB &tightAABB, index_type *indices,
			size_type primCount, bool isLeftChild, size_type badRefines) {
		KDAssert(nodeAABB.contains(tightAABB));

		Float leafCost = primCount * m_intersectionCost;
		if (primCount <= m_stopPrims || depth >= m_maxDepth) {
			createLeaf(ctx, node, indices, primCount);
			return leafCost;
		}

		if (primCount <= m_exactPrimThreshold) {
			OrderedChunkAllocator &alloc = isLeftChild 
					? ctx.leftAlloc : ctx.rightAlloc;
			boost::tuple<EdgeEvent *, EdgeEvent *, size_type> events  
					= createEventList(alloc, nodeAABB, indices, primCount);
			Float sahCost;
			if (m_parallelBuild) {
				m_interface.mutex->lock();
				m_interface.depth = depth;
				m_interface.node = node;
				m_interface.nodeAABB = nodeAABB;
				m_interface.eventStart = boost::get<0>(events);
				m_interface.eventEnd = boost::get<1>(events);
				m_interface.primCount = boost::get<2>(events);
				badRefines = badRefines;
				m_interface.cond->signal();

				/* Wait for a worker thread to take this job */
				while (m_interface.node)
					m_interface.condJobTaken->wait();
				m_interface.mutex->unlock();


				// Never tear down this subtree (return a SAH cost of -infinity)
				sahCost = -std::numeric_limits<Float>::infinity();
			} else {
				std::sort(boost::get<0>(events), boost::get<1>(events), 
						EdgeEventOrdering());

				sahCost = buildTreeSAH(ctx, depth, node, nodeAABB,
					boost::get<0>(events), boost::get<1>(events), boost::get<2>(events),
					isLeftChild, badRefines);
			}
			alloc.release(boost::get<0>(events));
			return sahCost;
		}

		/* ==================================================================== */
	    /*                              Binning                                 */
	    /* ==================================================================== */

		MinMaxBins<MTS_KD_MINMAX_BINS> bins(tightAABB);
		bins.bin(cast(), indices, primCount);

		/* ==================================================================== */
	    /*                        Split candidate search                        */
    	/* ==================================================================== */
		SplitCandidate bestSplit = bins.maximizeSAH(m_traversalCost,
			m_intersectionCost);

		/* "Bad refines" heuristic from PBRT */
		if (bestSplit.sahCost >= leafCost) {
			if ((bestSplit.sahCost > 4 * leafCost && primCount < 16)
				|| bestSplit.sahCost == std::numeric_limits<Float>::infinity()
				|| badRefines >= m_maxBadRefines) {
				createLeaf(ctx, node, indices, primCount);
				return leafCost;
			}
			++badRefines;
		}

		/* ==================================================================== */
	    /*                            Partitioning                              */
	    /* ==================================================================== */

		boost::tuple<AABB, index_type *, AABB, index_type *> partition = 
			bins.partition(ctx, cast(), indices, bestSplit, isLeftChild, 
			m_traversalCost, m_intersectionCost);

		/* ==================================================================== */
	    /*                              Recursion                               */
	    /* ==================================================================== */

		KDNode *children = ctx.nodes.allocate(2);

		size_type nodePosBeforeSplit = ctx.nodes.size();
		size_type indexPosBeforeSplit = ctx.indices.size();
		size_type leafNodeCountBeforeSplit = ctx.leafNodeCount;
		size_type nonemptyLeafNodeCountBeforeSplit = ctx.nonemptyLeafNodeCount;
		size_type innerNodeCountBeforeSplit = ctx.innerNodeCount;

		if (!node->initInnerNode(bestSplit.axis, bestSplit.pos, children-node)) {
			m_indirectionLock->lock();
			size_t indirectionIdx = m_indirections.size();
			m_indirections.push_back(children);
			/* Unable to store relative offset -- create an indirection
			   table entry */
			node->initIndirectionNode(bestSplit.axis, bestSplit.pos, indirectionIdx);
			m_indirectionLock->unlock();
		}
		ctx.innerNodeCount++;

		AABB childAABB(nodeAABB);
		childAABB.max[bestSplit.axis] = bestSplit.pos;
		Float saLeft = childAABB.getSurfaceArea();

		Float leftSAHCost = buildTreeMinMax(ctx, depth+1, children,
				childAABB, boost::get<0>(partition), boost::get<1>(partition), 
				bestSplit.numLeft, true, badRefines);

		childAABB.min[bestSplit.axis] = bestSplit.pos;
		childAABB.max[bestSplit.axis] = nodeAABB.max[bestSplit.axis];
		Float saRight = childAABB.getSurfaceArea();

		Float rightSAHCost = buildTreeMinMax(ctx, depth+1, children + 1,
				childAABB, boost::get<2>(partition), boost::get<3>(partition), 
				bestSplit.numRight, false, badRefines);

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

		/* ==================================================================== */
	    /*                           Final decision                             */
	    /* ==================================================================== */

		if (!m_retract || finalSAHCost < primCount * m_intersectionCost) {
			return finalSAHCost;
		} else {
			/* In the end, splitting didn't help to reduce the SAH cost.
			   Tear up everything below this node and create a leaf */
			ctx.nodes.resize(nodePosBeforeSplit);
			ctx.retractedSplits++;
			ctx.leafNodeCount = leafNodeCountBeforeSplit;
			ctx.nonemptyLeafNodeCount = nonemptyLeafNodeCountBeforeSplit;
			ctx.innerNodeCount = innerNodeCountBeforeSplit;
			createLeafAfterRetraction(ctx, node, indexPosBeforeSplit, primCount);
			return leafCost;
		}
	}

	/*
	 * \brief Build helper function (greedy SAH-based optimization)
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
	 *     the hope that the SAH cost was significantly overestimated. The
	 *     counter makes sure that only a limited number of such splits can
	 *     happen in succession.
	 * \returns 
	 *     Final SAH cost of the node
	 */
	Float buildTreeSAH(BuildContext &ctx, unsigned int depth, KDNode *node,
		const AABB &nodeAABB, EdgeEvent *eventStart, EdgeEvent *eventEnd, 
		size_type primCount, bool isLeftChild, size_type badRefines) {

		Float leafCost = primCount * m_intersectionCost;
		if (primCount <= m_stopPrims || depth >= m_maxDepth) {
			createLeaf(ctx, node, eventStart, eventEnd, primCount);
			return leafCost;
		}

		SplitCandidate bestSplit;

		/* ==================================================================== */
	    /*                        Split candidate search                        */
    	/* ==================================================================== */

		/* First, find the optimal splitting plane according to the
		   surface area heuristic. To do this in O(n), the search is
		   implemented as a sweep over the edge events */

		/* Initially, the split plane is placed left of the scene
		   and thus all geometry is on its right side */
		size_type numLeft[3]  = { 0, 0, 0 },
				  numRight[3] = { primCount, primCount, primCount };

		EdgeEvent *eventsByAxis[3];
		int eventsByAxisCtr = 1;
		eventsByAxis[0] = eventStart;

		const Vector extents(nodeAABB.getExtents());
		const Float invSA = 0.5f / (extents.x * extents.y 
				+ extents.y*extents.z + extents.x*extents.z);
		const Vector temp0 = Vector(
			(extents[1] * extents[2]),
			(extents[0] * extents[2]),
			(extents[0] * extents[1])) * 2 * invSA;

		const Vector temp1 = Vector(
			(extents[1] + extents[2]),
			(extents[0] + extents[2]),
			(extents[0] + extents[1])) * 2 * invSA;

		/* Iterate over all events on the current axis */
		for (EdgeEvent *event = eventStart; event < eventEnd; ) {
			/* Record the current position and count the number
			   and type of remaining events, which are also here.
			   Due to the sort ordering, there is no need to worry 
			   about an axis change in the loops below */
			int axis = event->axis;
			float pos = event->pos;
			size_type numStart = 0, numEnd = 0, numPlanar = 0;

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
				KDAssert(eventsByAxisCtr < 3);
				eventsByAxis[eventsByAxisCtr++] = event;
			}

			/* The split plane can now be moved onto 't'. Accordingly, all planar 
			   and ending primitives are removed from the right side */
			numRight[axis] -= numPlanar + numEnd;

			/* Calculate a score using the surface area heuristic */
			if (EXPECT_TAKEN(pos > nodeAABB.min[axis] && pos < nodeAABB.max[axis])) {
				const size_type nL = numLeft[axis], nR = numRight[axis];
				const Float nLF = (Float) nL, nRF = (Float) nR;

				Float pLeft  = temp0[axis] + temp1[axis] * (pos - nodeAABB.min[axis]);
				Float pRight = temp0[axis] + temp1[axis] * (nodeAABB.max[axis] - pos);

				if (numPlanar == 0) {
					Float sahCost = m_intersectionCost + m_traversalCost 
						* (pLeft * nLF + pRight * nRF);
					if (nL == 0 || nR == 0)
						sahCost *= m_emptySpaceBonus;
					if (sahCost < bestSplit.sahCost) {
						bestSplit.pos = pos;
						bestSplit.axis = axis;
						bestSplit.sahCost = sahCost;
						bestSplit.numLeft = nL;
						bestSplit.numRight = nR;
					}
				} else {
					Float sahCostPlanarLeft  = m_intersectionCost + m_traversalCost 
						* (pLeft * (nL+numPlanar) + pRight * nRF);
					Float sahCostPlanarRight = m_intersectionCost + m_traversalCost 
						* (pLeft * nLF + pRight * (nR+numPlanar));

					if (nL + numPlanar == 0 || nR == 0)
						sahCostPlanarLeft *= m_emptySpaceBonus;
					if (nL == 0 || nR + numPlanar == 0)
						sahCostPlanarRight *= m_emptySpaceBonus;

					if (sahCostPlanarLeft  < bestSplit.sahCost ||
						sahCostPlanarRight < bestSplit.sahCost) {
						bestSplit.pos = pos;
						bestSplit.axis = axis;

						if (sahCostPlanarLeft < sahCostPlanarRight) {
							bestSplit.sahCost = sahCostPlanarLeft;
							bestSplit.numLeft = nL + numPlanar;
							bestSplit.numRight = nR;
							bestSplit.planarLeft = true;
						} else {
							bestSplit.sahCost = sahCostPlanarRight;
							bestSplit.numLeft = nL;
							bestSplit.numRight = nR + numPlanar;
							bestSplit.planarLeft = false;
						}
					}
				}
			} else {
				#if defined(MTS_KD_DEBUG)
				if (m_clip && (pos < nodeAABB.min[axis] || pos > nodeAABB.max[axis])) {
					/* When primitive clipping is active, this should  never happen! */
					Log(EError, "Internal error: edge event is out of bounds");
				}
				#endif
			}

			/* The split plane is moved past 't'. All prims,
				which were planar on 't', are moved to the left
				side. Also, starting prims are now also left of
				the split plane. */
			numLeft[axis] += numStart + numPlanar;
		}

		/* Sanity checks. Everything should now be left of the split plane */
		KDAssert(numRight[0] == 0 && numLeft[0] == primCount &&
			   numRight[1] == 0 && numLeft[1] == primCount &&
			   numRight[2] == 0 && numLeft[2] == primCount);

		KDAssert(eventsByAxis[1]->axis == 1 && (eventsByAxis[1]-1)->axis == 0);
		KDAssert(eventsByAxis[2]->axis == 2 && (eventsByAxis[2]-1)->axis == 1);

		/* "Bad refines" heuristic from PBRT */
		if (bestSplit.sahCost >= leafCost) {
			if ((bestSplit.sahCost > 4 * leafCost && primCount < 16)
				|| badRefines >= m_maxBadRefines
				|| bestSplit.sahCost == std::numeric_limits<Float>::infinity()) {
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

		size_type primsLeft = 0, primsRight = 0, primsBoth = primCount;
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
				   the better SAH score */
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
			rightEventsStart = rightAlloc.allocate<EdgeEvent>(bestSplit.numRight * 6);
		} else {
			leftEventsStart = leftAlloc.allocate<EdgeEvent>(bestSplit.numLeft * 6);
			rightEventsStart = eventStart;
		}

		EdgeEvent *leftEventsEnd = leftEventsStart, *rightEventsEnd = rightEventsStart;

		AABB leftNodeAABB = nodeAABB, rightNodeAABB = nodeAABB;
		leftNodeAABB.max[bestSplit.axis] = bestSplit.pos;
		rightNodeAABB.min[bestSplit.axis] = bestSplit.pos;

		size_type prunedLeft = 0, prunedRight = 0;

		/* ==================================================================== */
		/*                            Partitioning                              */
		/* ==================================================================== */

		if (m_clip) {
			EdgeEvent
			  *leftEventsTempStart = leftAlloc.allocate<EdgeEvent>(primsLeft * 6),
			  *rightEventsTempStart = rightAlloc.allocate<EdgeEvent>(primsRight * 6),
			  *newEventsLeftStart = leftAlloc.allocate<EdgeEvent>(primsBoth * 6),
			  *newEventsRightStart = rightAlloc.allocate<EdgeEvent>(primsBoth * 6);

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
					const index_type index = event->index;

					AABB clippedLeft = cast()->getClippedAABB(index, leftNodeAABB);
					AABB clippedRight = cast()->getClippedAABB(index, rightNodeAABB);

					KDAssert(leftNodeAABB.contains(clippedLeft));
					KDAssert(rightNodeAABB.contains(clippedRight));

					if (clippedLeft.isValid() && clippedLeft.getSurfaceArea() > 0) {
						for (int axis=0; axis<3; ++axis) {
							float min = (float) clippedLeft.min[axis],
								  max = (float) clippedLeft.max[axis];

							if (min == max) {
								*newEventsLeftEnd++ = EdgeEvent(EdgeEvent::EEdgePlanar, 
										axis, min, index);
							} else {
								*newEventsLeftEnd++ = EdgeEvent(EdgeEvent::EEdgeStart, 
										axis, min, index);
								*newEventsLeftEnd++ = EdgeEvent(EdgeEvent::EEdgeEnd, 
										axis, max, index);
							}
						}
					} else {
						prunedLeft++;
					}

					if (clippedRight.isValid() && clippedRight.getSurfaceArea() > 0) {
						for (int axis=0; axis<3; ++axis) {
							float min = (float) clippedRight.min[axis],
								  max = (float) clippedRight.max[axis];

							if (min == max) {
								*newEventsRightEnd++ = EdgeEvent(EdgeEvent::EEdgePlanar,
										axis, min, index);
							} else {
								*newEventsRightEnd++ = EdgeEvent(EdgeEvent::EEdgeStart, 
										axis, min, index);
								*newEventsRightEnd++ = EdgeEvent(EdgeEvent::EEdgeEnd,
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

			KDAssert(leftEventsTempEnd - leftEventsTempStart <= primsLeft * 6);
			KDAssert(rightEventsTempEnd - rightEventsTempStart <= primsRight * 6);
			KDAssert(newEventsLeftEnd - newEventsLeftStart <= primsBoth * 6);
			KDAssert(newEventsRightEnd - newEventsRightStart <= primsBoth * 6);
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
			KDAssert(leftEventsEnd - leftEventsStart <= bestSplit.numLeft * 6);
			KDAssert(rightEventsEnd - rightEventsStart <= bestSplit.numRight * 6);
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

		size_type nodePosBeforeSplit = ctx.nodes.size();
		size_type indexPosBeforeSplit = ctx.indices.size();
		size_type leafNodeCountBeforeSplit = ctx.leafNodeCount;
		size_type nonemptyLeafNodeCountBeforeSplit = ctx.nonemptyLeafNodeCount;
		size_type innerNodeCountBeforeSplit = ctx.innerNodeCount;

		if (!node->initInnerNode(bestSplit.axis, bestSplit.pos, children-node)) {
			m_indirectionLock->lock();
			size_t indirectionIdx = m_indirections.size();
			m_indirections.push_back(children);
			/* Unable to store relative offset -- create an indirection
			   table entry */
			node->initIndirectionNode(bestSplit.axis, bestSplit.pos, indirectionIdx);
			m_indirectionLock->unlock();
		}
		ctx.innerNodeCount++;

		Float leftSAHCost = buildTreeSAH(ctx, depth+1, children,
				leftNodeAABB, leftEventsStart, leftEventsEnd,
				bestSplit.numLeft - prunedLeft, true, badRefines);

		Float rightSAHCost = buildTreeSAH(ctx, depth+1, children+1,
				rightNodeAABB, rightEventsStart, rightEventsEnd,
				bestSplit.numRight - prunedRight, false, badRefines);

		Float saLeft = leftNodeAABB.getSurfaceArea();
		Float saRight = rightNodeAABB.getSurfaceArea();

		/* Compute the final SAH cost given the updated cost 
		   values received from the children */
		Float finalSAHCost = m_traversalCost + 
			(saLeft * leftSAHCost + saRight * rightSAHCost) * invSA;

		/* Release the index lists not needed by the children anymore */
		if (isLeftChild)
			ctx.rightAlloc.release(rightEventsStart);
		else
			ctx.leftAlloc.release(leftEventsStart);

		/* ==================================================================== */
	    /*                           Final decision                             */
	    /* ==================================================================== */

		if (!m_retract || finalSAHCost < primCount * m_intersectionCost) {
			return finalSAHCost;
		} else {
			/* In the end, splitting didn't help to reduce the SAH cost.
			   Tear up everything below this node and create a leaf */
			ctx.nodes.resize(nodePosBeforeSplit);
			ctx.retractedSplits++;
			ctx.leafNodeCount = leafNodeCountBeforeSplit;
			ctx.nonemptyLeafNodeCount = nonemptyLeafNodeCountBeforeSplit;
			ctx.innerNodeCount = innerNodeCountBeforeSplit;
			createLeafAfterRetraction(ctx, node, indexPosBeforeSplit, primCount);
			return leafCost;
		}

		return bestSplit.sahCost;
	}

	/**
	 * \brief Min-max binning as described in
	 * "Highly Parallel Fast KD-tree Construction for Interactive
	 *  Ray Tracing of Dynamic Scenes"
	 * by M. Shevtsov, A. Soupikov and A. Kapustin
	 *
	 * \tparam BinCount Number of bins to be allocated
	 */
	template <int BinCount> struct MinMaxBins {
		MinMaxBins(const AABB &aabb) : m_aabb(aabb) {
			m_binSize = m_aabb.getExtents() / BinCount;
			for (int axis=0; axis<3; ++axis) 
				m_invBinSize[axis] = 1/m_binSize[axis];
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

			for (size_type i=0; i<m_primCount; ++i) {
				const AABB aabb = derived->getAABB(indices[i]);
				for (int axis=0; axis<3; ++axis) {
					int minIdx = (int) ((aabb.min[axis] - m_aabb.min[axis]) 
							* m_invBinSize[axis]);
					int maxIdx = (int) ((aabb.max[axis] - m_aabb.min[axis]) 
							* m_invBinSize[axis]);
					m_maxBins[axis * BinCount 
						+ std::max(0, std::min(maxIdx, BinCount-1))]++;
					m_minBins[axis * BinCount 
						+ std::max(0, std::min(minIdx, BinCount-1))]++;
				}
			}
		}

		/**
		 * \brief Evaluate the surface area heuristic at each bin boundary
		 * and return the maximizer for the given cost constants. Min-max
		 * binning uses no "empty space bonus" since it cannot create such
		 * splits.
		 */
		SplitCandidate maximizeSAH(Float traversalCost, Float intersectionCost) {
			SplitCandidate candidate;
			Float normalization = 2.0f / m_aabb.getSurfaceArea();
			int binIdx = 0, leftBin = 0;

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

					Float sahCost = traversalCost + intersectionCost 
						* (pLeft * numLeft + pRight * numRight);

					if (sahCost < candidate.sahCost) {
						candidate.sahCost = sahCost;
						candidate.axis = axis;
						candidate.numLeft = numLeft;
						candidate.numRight = numRight;
						leftBin = i;
					}

					binIdx++;
				}
				binIdx++;
			}

			KDAssert(candidate.sahCost != std::numeric_limits<Float>::infinity());

			const int axis = candidate.axis;
			const Float min = m_aabb.min[axis];

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
			Float invBinSize = m_invBinSize[axis];
			float split = min + (leftBin + 1) * m_binSize[axis];
			float splitNext = nextafterf(split, 
				  std::numeric_limits<float>::max());
			int idx     = (int) ((split - min) * invBinSize);
			int idxNext = (int) ((splitNext - min) * invBinSize);

			/**
			 * The split plane should be along the last discrete floating
			 * floating position, which would still be classified into
			 * the left bin.
			 */
			if (!(idx <= leftBin && idxNext > leftBin)) {
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
					if (idx < leftBin && idxNext > leftBin) {
						/* Insufficient floating point resolution
						   -> a leaf will be created. */
						candidate.sahCost = std::numeric_limits<Float>::infinity();
						break;
					}

					++it;
				}
			}

			if (split <= m_aabb.min[axis] || split > m_aabb.max[axis]) {
				/* Insufficient floating point resolution 
				   -> a leaf will be created. */
				candidate.sahCost = std::numeric_limits<Float>::infinity();
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
				SplitCandidate &split, bool isLeftChild, Float traversalCost, 
				Float intersectionCost) {
			const float splitPos = split.pos;
			const int axis = split.axis;
			size_type numLeft = 0, numRight = 0;
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
					KDAssert(numLeft < split.numLeft);
					leftBounds.expandBy(aabb);
					leftIndices[numLeft++] = primIndex;
				} else if (aabb.min[axis] > splitPos) {
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

			KDAssert(numLeft == split.numLeft);
			KDAssert(numRight == split.numRight);

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
					split.sahCost = sahCost1;
					split.pos = leftBounds.max[axis];
				} else {
					split.sahCost = sahCost2;
					split.pos = rightBounds.min[axis];
				}

				leftBounds.max[axis] = std::min(leftBounds.max[axis], (Float) split.pos);
				rightBounds.min[axis] = std::max(rightBounds.min[axis], (Float) split.pos);
			}

			return boost::make_tuple(leftBounds, leftIndices,
					rightBounds, rightIndices);
		}
	private:
		size_type m_minBins[3*BinCount], m_maxBins[3*BinCount];
		size_type m_primCount;
		AABB m_aabb;
		Vector m_binSize, m_invBinSize;
	};

	/// Ray traversal stack entry for incoherent ray tracing
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
	 * \brief Internal kd-tree traversal implementation (Havran variant)
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
				/* Faster than the original code with the 
				   prevAxis & nextAxis table */
				stack[exPt].p = ray(distToSplit);
				stack[exPt].p[axis] = splitVal;
				#else
				const int nextAxis = nextAxisTable[axis];
				const int prevAxis = prevAxisTable[axis];
				stack[exPt].p[axis] = splitVal;
				stack[exPt].p[nextAxis] = ray.o[nextAxis] + distToSplit*ray.d[nextAxis];
				stack[exPt].p[prevAxis] = ray.o[prevAxis] + distToSplit*ray.d[prevAxis];
				#endif

			}
			/* Reached a leaf node */

			/* Floating-point arithmetic.. - use both absolute and relative 
			   epsilons when looking for intersections in the subinterval */
			#if defined(SINGLE_PRECISION)
			const Float eps = 1e-3;
			#else
			const Float eps = 1e-5;
			#endif
			const Float m_eps = 1-eps, p_eps = 1+eps;

			const Float searchStart = std::max(mint, stack[enPt].t * m_eps - eps);
			Float searchEnd   = std::min(maxt, stack[exPt].t * p_eps + eps);

			bool foundIntersection = false;

			for (index_type entry=currNode->getPrimStart(),
					last = currNode->getPrimEnd(); entry != last; entry++) {
				const index_type primIdx = m_indices[entry];

				#if defined(MTS_KD_MAILBOX_ENABLED)
				if (mailbox.contains(primIdx)) 
					continue;
				#endif

				EIntersectionResult result = cast()->intersect(ray, 
						primIdx, searchStart, searchEnd, t, temp);

				if (result == EYes) {
					if (shadowRay)
						return true;
					searchEnd = t;
					foundIntersection = true;
				}

				#if defined(MTS_KD_MAILBOX_ENABLED)
				else if (result == ENever) 
					mailbox.put(primIdx);
				#endif
			}

			if (foundIntersection) 
				return true;

			/* Pop from the stack and advance to the next node on the interval */
			enPt = exPt;
			currNode = stack[exPt].node;
			exPt = stack[enPt].prev;
		}

		return false;
	}

	/**
	 * \brief Internal kd-tree traversal implementation (Havran variant)
	 *
	 * This method is almost identical to \ref rayIntersectHavran, except
	 * that it additionally returns statistics on the number of traversed
	 * nodes, intersected shapes, as well as the time taken to do this
	 * (measured using rtdsc).
	 */
	FINLINE boost::tuple<bool, uint32_t, uint32_t, uint64_t> 
			rayIntersectHavranCollectStatistics(const Ray &ray, 
			Float mint, Float maxt, Float &t, void *temp) const {
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

			#if defined(SINGLE_PRECISION)
				const Float eps = 1e-3;
			#else
				const Float eps = 1e-5;
			#endif
			const Float m_eps = 1-eps, p_eps = 1+eps;

			/* Floating-point arithmetic.. - use both absolute and relative 
			   epsilons when looking for intersections in the subinterval */
			const Float searchStart = std::max(mint, stack[enPt].t * m_eps - eps);
			Float searchEnd = std::min(maxt, stack[exPt].t * p_eps + eps);

			/* Reached a leaf node */
			bool foundIntersection = false;
			for (unsigned int entry=currNode->getPrimStart(),
					last = currNode->getPrimEnd(); entry != last; entry++) {
				const index_type primIdx = m_indices[entry];

				++numIntersections;
				EIntersectionResult result = cast()->intersect(ray, 
						primIdx, searchStart, searchEnd, t, temp);

				if (result == EYes) {
					searchEnd = t;
					foundIntersection = true;
				}
			}

			if (foundIntersection)
				return boost::make_tuple(true, numTraversals, 
						numIntersections, rdtsc() - timer);


			/* Pop from the stack and advance to the next node on the interval */
			enPt = exPt;
			currNode = stack[exPt].node;
			exPt = stack[enPt].prev;
		}

		return boost::make_tuple(false, numTraversals, 
				numIntersections, rdtsc() - timer);
	}

	struct KDStackEntry {
		const KDNode * __restrict node;
		Float mint, maxt;
	};

	/**
	 * \brief Internal kd-tree traversal implementation (Plain variant)
	 */
	template<bool shadowRay> FINLINE bool rayIntersectPlain(const Ray &ray, 
			Float mint_, Float maxt_, Float &t, void *temp) const {
		KDStackEntry stack[MTS_KD_MAXDEPTH];
		int stackPos = 0;
		Float mint = mint_, maxt = maxt_;
		const KDNode *node = m_nodes;

		#if defined(MTS_KD_MAILBOX_ENABLED)
		HashedMailbox mailbox;
		#endif

		const int signs[3] = {
			ray.d[0] >= 0 ? 1 : 0,
			ray.d[1] >= 0 ? 1 : 0,
			ray.d[2] >= 0 ? 1 : 0
		};

		while (node != NULL) {
			while (EXPECT_TAKEN(!node->isLeaf())) {
				const Float split = (Float) node->getSplit();
				const int axis = node->getAxis();
				const float t = (split - ray.o[axis]) * ray.dRcp[axis];

				int sign = signs[axis];

				const KDNode * __restrict base = node->getLeft();
				const KDNode * __restrict first = base + (sign^1);
				const KDNode * __restrict second = base + sign;

				node = first;

				bool front = t > maxt,
					 back = t < mint;

				if (front) {
					;
				} else if (back) {
					node = second;
				} else {
					stack[stackPos].node = second;
					stack[stackPos].mint = t;
					stack[stackPos].maxt = maxt;
					maxt = t;
					++stackPos;
				}
			}
			bool foundIntersection = false;
			for (unsigned int entry=node->getPrimStart(),
					last = node->getPrimEnd(); entry != last; entry++) {
				const index_type primIdx = m_indices[entry];

				#if defined(MTS_KD_MAILBOX_ENABLED)
				if (mailbox.contains(primIdx)) 
					continue;
				#endif

				EIntersectionResult result = cast()->intersect(ray, 
						primIdx, mint, maxt, t, temp);

				if (result == EYes) {
					if (shadowRay)
						return true;
					maxt = t;
					foundIntersection = true;
				}

				#if defined(MTS_KD_MAILBOX_ENABLED)
				else if (result == ENever) 
					mailbox.put(primIdx);
				#endif
			}

			if (foundIntersection) 
				return true;

			if (stackPos > 0) {
				--stackPos;
				node = stack[stackPos].node;
				mint = stack[stackPos].mint;
				maxt = stack[stackPos].maxt;
			} else {
				break;
			}
		}
		return false;
	}

	/**
	 * \brief Internal kd-tree traversal implementation (PBRT variant)
	 */
	template<bool shadowRay> FINLINE bool rayIntersectPBRT(const Ray &ray, 
			Float mint_, Float maxt_, Float &t, void *temp) const {
		KDStackEntry stack[MTS_KD_MAXDEPTH];
		int stackPos = 0;
		Float mint = mint_, maxt=maxt_;
		const KDNode *node = m_nodes;
		bool foundIntersection = false;

		#if defined(MTS_KD_MAILBOX_ENABLED)
		HashedMailbox mailbox;
		#endif

		while (node != NULL) {
			if (t < mint)
				break;

			if (EXPECT_TAKEN(!node->isLeaf())) {
				const Float split = (Float) node->getSplit();
				const int axis = node->getAxis();
				const float tPlane = (split - ray.o[axis]) * ray.dRcp[axis];
				bool leftOfSplit = (ray.o[axis] < split) 
					|| (ray.o[axis] == split && ray.d[axis] >= 0);

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
					const index_type primIdx = m_indices[entry];

					#if defined(MTS_KD_MAILBOX_ENABLED)
					if (mailbox.contains(primIdx)) 
						continue;
					#endif

					EIntersectionResult result = cast()->intersect(ray, 
							primIdx, mint_, maxt_, t, temp);

					if (result == EYes) {
						if (shadowRay)
							return true;
						maxt_ = t;
						foundIntersection = true;
					}

					#if defined(MTS_KD_MAILBOX_ENABLED)
					else if (result == ENever) 
						mailbox.put(primIdx);
					#endif
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

protected:
	KDNode *m_nodes;
	index_type *m_indices;
	Float m_traversalCost;
	Float m_intersectionCost;
	Float m_emptySpaceBonus;
	bool m_clip, m_retract, m_parallelBuild;
	AABB m_aabb;
	BSphere m_bsphere;
	size_type m_maxDepth;
	size_type m_stopPrims;
	size_type m_maxBadRefines;
	size_type m_exactPrimThreshold;
	size_type m_nodeCount;
	std::vector<SAHTreeBuilder *> m_builders;
	std::vector<KDNode *> m_indirections;
	ref<Mutex> m_indirectionLock;
	BuildInterface m_interface;
};

template <typename Derived> void GenericKDTree<Derived>::findCosts(
		Float &traversalCost, Float &intersectionCost) {
	ref<Random> random = new Random();
	uint8_t temp[MTS_KD_INTERSECTION_TEMP];
	BSphere bsphere = m_aabb.getBSphere();
	const size_t nRays = 10000000, warmup = nRays/4;
	Vector *A = new Vector[nRays-warmup];
	Float *b = new Float[nRays-warmup];
	size_t nIntersections = 0, idx = 0;

	for (size_t i=0; i<nRays; ++i) {
		Point2 sample1(random->nextFloat(), random->nextFloat()),
			sample2(random->nextFloat(), random->nextFloat());
		Point p1 = bsphere.center + squareToSphere(sample1) * bsphere.radius;
		Point p2 = bsphere.center + squareToSphere(sample2) * bsphere.radius;
		Ray ray(p1, normalize(p2-p1));
		Float mint, maxt, t;
		if (ray.mint > mint) mint = ray.mint;
		if (ray.maxt < maxt) maxt = ray.maxt;
		if (m_aabb.rayIntersect(ray, mint, maxt)) {
			if (EXPECT_TAKEN(maxt > mint)) {
				boost::tuple<bool, uint32_t, uint32_t, uint64_t> statistics =
					rayIntersectHavranCollectStatistics(ray, mint, maxt, t, temp);
				if (boost::get<0>(statistics))
					nIntersections++;
				if (i > warmup) {
					A[idx].x = 1;
					A[idx].y = boost::get<1>(statistics);
					A[idx].z = boost::get<2>(statistics);
					b[idx]   = boost::get<3>(statistics);
					idx++;
				}
			}
		}
	}

	Log(EDebug, "Fitting to " SIZE_T_FMT " samples (" SIZE_T_FMT 
			" intersections)", idx, nIntersections);

	/* Solve using normal equations */
	ref<Matrix4x4> M = new Matrix4x4(0.0f);
	Vector4 rhs(0.0f), x;

	for (size_t i=0; i<3; ++i) {
		for (size_t j=0; j<3; ++j) 
			for (size_t k=0; k<idx; ++k) 
				M->m[i][j] += A[k][i]*A[k][j];
		for (size_t k=0; k<idx; ++k) 
			rhs[i] += A[k][i]*b[k];
	}
	M->m[3][3] = 1.0f;

	Transform(M->inverse())(rhs, x);

	Float avgRdtsc = 0, avgResidual = 0;
	for (size_t i=0; i<idx; ++i) {
		avgRdtsc += b[i];
		Float model = x[0] * A[i][0]
			+ x[1] * A[i][1]
			+ x[2] * A[i][2];
		//cout << b[i] << " vs " << (int) model << endl;
		avgResidual += std::abs(b[i] - model);
	}
	avgRdtsc /= idx;
	avgResidual /= idx;

	for (size_t k=0; k<idx; ++k) 
		avgRdtsc += b[k];
	avgRdtsc /= idx;

	Log(EDebug, "Least squares fit:");
	Log(EDebug, "   Constant overhead    = %.2f", x[0]);
	Log(EDebug, "   Traversal cost       = %.2f", x[1]);
	Log(EDebug, "   Intersection cost    = %.2f", x[2]);
	Log(EDebug, "   Average rdtsc value  = %.2f", avgRdtsc);
	Log(EDebug, "   Avg. residual        = %.2f", avgResidual);
	x *= 10/x[1];
	Log(EDebug, "Re-scaled:");
	Log(EDebug, "   Constant overhead    = %.2f", x[0]);
	Log(EDebug, "   Traversal cost       = %.2f", x[1]);
	Log(EDebug, "   Intersection cost    = %.2f", x[2]);

	delete[] A;
	delete[] b;
	traversalCost = x[1];
	intersectionCost = x[2];
}

template <typename Derived> Class *GenericKDTree<Derived>::m_theClass 
	= new Class("GenericKDTree", true, "Object");

template <typename Derived> const Class *GenericKDTree<Derived>::getClass() const {
	return m_theClass;
}

MTS_NAMESPACE_END

#endif /* __KDTREE_GENERIC_H */

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

#define MTS_KD_MAX_DEPTH 48     ///< Compile-time KD-tree depth limit
#define MTS_KD_STATISTICS 1     ///< Collect statistics during building/traversal 
#define MTS_KD_MINMAX_BINS 32   ///< Min-max bin count
#define MTS_KD_MIN_ALLOC 128    ///< Allocate memory in 128 KB chunks

#if 1
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
 */
class OrderedChunkAllocator {
public:
	inline OrderedChunkAllocator() {
		m_chunks.reserve(16);
	}

	/**
	 * \brief Release all memory used by the allocator
	 */
	void cleanup() {
		for (std::vector<Chunk>::iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it)
			freeAligned((*it).start);
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
		/* Uh oh, allocation could not be found. Check if it has size==0 */
		for (std::vector<Chunk>::iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it) {
			const Chunk &chunk = *it;
			if ((uint8_t *) ptr == chunk.start + chunk.size) 
				return;
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
		/* Uh oh, allocation could not be found. Check if it has size==0 */
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
	}

	inline size_t getChunkCount() const { return m_chunks.size(); }

	/**
	 * \brief Return the total amount of chunk memory in bytes
	 */
	size_t getSize() const {
		size_t result = 0;
		for (std::vector<Chunk>::const_iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it)
			result += (*it).size;
		return result;
	}
	
	/**
	 * \brief Return the total amount of used memory in bytes
	 */
	size_t getUsed() const {
		size_t result = 0;
		for (std::vector<Chunk>::const_iterator it = m_chunks.begin();
				it != m_chunks.end(); ++it)
			result += (*it).getUsed();
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

		inline size_t getUsed() const {
			return cur - start;
		}

		inline size_t getRemainder() const {
			return size - getUsed();
		}

		std::string toString() const {
			return formatString("0x%llx-0x%llx (size=" SIZE_T_FMT 
					", used=" SIZE_T_FMT ")", start, start+size, 
					size, getUsed());
		}
	};

	std::vector<Chunk> m_chunks;
};

/**
 * \brief Compact storage for primitive classifcation
 *
 * When classifying primitives with respect to a split plane,
 * a data structure is needed to hold the tertiary result of
 * this operation. This class implements a compact storage
 * (2 bits per entry) in the spirit of the std::vector<bool> 
 * specialization.
 */
class ClassificationStorage {
public:
	inline ClassificationStorage() : m_buffer(NULL),
		m_bufferSize(0) {
	}

	inline ClassificationStorage(size_t size) {
		m_bufferSize = size/4 + ((size % 4) > 0 ? 1 : 0);
		m_buffer = new uint8_t[m_bufferSize];
	}

	inline ~ClassificationStorage() {
		if (m_buffer)
			delete[] m_buffer;
	}

	inline void set(uint32_t index, uint8_t value) {
		uint8_t *ptr = m_buffer + (index >> 2);
		uint8_t shift = (index & 3) << 1;
		*ptr = (*ptr & ~(3 << shift)) | (value << shift);
	}

	inline uint8_t get(uint32_t index) const {
		uint8_t *ptr = m_buffer + (index >> 2);
		uint8_t shift = (index & 3) << 1;
		return (*ptr >> shift) & 3;
	}

	inline size_t getSize() const {
		return m_bufferSize;
	}
private:
	uint8_t *m_buffer;
	size_t m_bufferSize;
};

/**
 * \brief Generic kd-tree for data structure used to accelerate 
 * ray intersections against large amounts of three-dimensional
 * shapes.
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

	/**
	 * \brief Create a new kd-tree instance initialized with 
	 * the default parameters.
	 */
	GenericKDTree() : m_root(NULL) {
		m_traversalCost = 15;
		m_intersectionCost = 20;
		m_emptySpaceBonus = 0.9f;
		m_clip = true;
		m_stopPrims = 1;
		m_maxBadRefines = 2;
		m_exactPrimThreshold = 1024;
		m_maxDepth = 0;
		m_retract = true;
	}

	/**
	 * \brief Build a KD-tree over supplied geometry
	 */
	void build() {
		if (m_root != NULL)
			Log(EError, "The kd-tree has already been built!");

		size_type primCount = downCast()->getPrimitiveCount();
		BuildContext ctx(primCount);

		/* Establish an ad-hoc depth cutoff value (Formula from PBRT) */
		if (m_maxDepth == 0)
			m_maxDepth = (int) (8 + 1.3f * log2i(primCount));
		m_maxDepth = std::min(m_maxDepth, (size_type) MTS_KD_MAX_DEPTH);

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
		Log(EDebug, "");

		Log(EDebug, "kd-tree configuration:");
		Log(EDebug, "   Traversal cost           : %.2f", m_traversalCost);
		Log(EDebug, "   Intersection cost        : %.2f", m_intersectionCost);
		Log(EDebug, "   Empty space bonus        : %.2f", m_emptySpaceBonus);
		Log(EDebug, "   Max. tree depth          : %i", m_maxDepth);
		Log(EDebug, "   Stopping primitive count : %i", m_stopPrims);
		Log(EDebug, "   Scene bounding box (min) : %s", m_aabb.min.toString().c_str());
		Log(EDebug, "   Scene bounding box (max) : %s", m_aabb.max.toString().c_str());
		Log(EDebug, "   Min-max bins             : %i", MTS_KD_MINMAX_BINS);
		Log(EDebug, "   Greedy SAH optimization  : <= %i primitives", m_exactPrimThreshold);
		Log(EDebug, "   Perfect splits           : %s", m_clip ? "yes" : "no");
		Log(EDebug, "   Retract bad splits       : %s", m_retract ? "yes" : "no");
		Log(EDebug, "");

		size_type procCount = getProcessorCount();
		m_builders.resize(procCount);
		for (size_type i=0; i<procCount; ++i) {
			m_builders[i] = new SAHTreeBuilder(i+1, primCount, m_interface);
			m_builders[i]->incRef();
			m_builders[i]->start();
		}

		Log(EInfo, "Constructing a SAH kd-tree (%i primitives) ..", primCount);

		m_root = nodeAlloc.allocate<KDNode>(1);
		Float finalSAHCost = buildTreeMinMax(ctx, 1, m_root, 
			m_aabb, m_aabb, indices, primCount, true, 0);
		ctx.leftAlloc.release(indices);

		KDAssert(ctx.leftAlloc.getUsed() == 0);
		KDAssert(ctx.rightAlloc.getUsed() == 0);

		m_interface.done = true;
		m_interface.cond->broadcast();
		for (size_type i=0; i<procCount; ++i) 
			m_builders[i]->join();
		Log(EDebug, "");
		Log(EInfo, "Finished -- took %i ms.", timer->getMilliseconds());

		Log(EDebug, "Memory allocation statistics:");
		Log(EDebug, "   Classification storage : %.2f KiB", 
				(ctx.storage.getSize() * (1+procCount)) / 1024.0f);

		Log(EDebug, "   Main:");
		ctx.printStats();

		for (size_type i=0; i<procCount; ++i) {
			Log(EDebug, "   Thread %i:", i+1);
			BuildContext &subCtx = m_builders[i]->getContext();
			subCtx.printStats();
			ctx.accumulateStatistics(subCtx);
			ctx.nodeAlloc.merge(subCtx.nodeAlloc);
			m_builders[i]->decRef();
		}
		m_builders.clear();
		Log(EDebug, "");

		Float rootSA = m_aabb.getSurfaceArea();
		ctx.expTraversalSteps /= rootSA;
		ctx.expLeavesVisited /= rootSA;
		ctx.expPrimitivesIntersected /= rootSA;

		Log(EDebug, "Detailed kd-tree statistics:");
		Log(EDebug, "   Final SAH cost      : %.2f", finalSAHCost);
		Log(EDebug, "   Inner nodes         : %i", ctx.innerNodeCount);
		Log(EDebug, "   Leaf nodes          : %i", ctx.leafNodeCount);
		Log(EDebug, "   Nonempty leaf nodes : %i", ctx.nonemptyLeafNodeCount);
		Log(EDebug, "   Retracted splits    : %i", ctx.retractedSplits);
		Log(EDebug, "   Pruned primitives   : %i", ctx.pruned);
		Log(EDebug, "   Exp. traversals     : %.2f", ctx.expTraversalSteps);
		Log(EDebug, "   Exp. leaf visits    : %.2f", ctx.expLeavesVisited);
		Log(EDebug, "   Exp. intersections  : %.2f", ctx.expPrimitivesIntersected);
		Log(EDebug, "   Indirection table   : " SIZE_T_FMT " entries",
				m_indirectionTable.size());
		Log(EDebug, "");

		ctx.leftAlloc.cleanup();
		ctx.rightAlloc.cleanup();
		m_aabb.getSurfaceArea();
	}
protected:
	/// Primitive classification during tree-construction
	enum EClassificationResult {
		EBothSides = 0,
		ELeftSide = 1,
		ERightSide = 2,
		EBothSidesProcessed = 3 //< Used to indicate that edge events have already been generated for a straddling primitive
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
		inline EdgeEvent(uint16_t type, uint16_t axis, float pos, index_type index)
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
		uint16_t type;
		/// Event axis
		uint16_t axis;
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
		OrderedChunkAllocator leftAlloc, rightAlloc, nodeAlloc;
		ClassificationStorage storage;
	
		size_type leafNodeCount;
		size_type nonemptyLeafNodeCount;
		size_type innerNodeCount;
		size_type retractedSplits;
		size_type pruned;
		Float expTraversalSteps;
		Float expLeavesVisited;
		Float expPrimitivesIntersected;

		inline BuildContext(size_type primCount) : storage(primCount) {
			retractedSplits = 0;
			leafNodeCount = 0;
			nonemptyLeafNodeCount = 0;
			innerNodeCount = 0;
			pruned = 0;
			expTraversalSteps = 0;
			expLeavesVisited = 0;
			expPrimitivesIntersected = 0;
		}

		void printStats() {
			Log(EDebug, "      Left:  " SIZE_T_FMT " chunks (%.2f KiB)",
					leftAlloc.getChunkCount(), leftAlloc.getSize() / 1024.0f);
			Log(EDebug, "      Right: " SIZE_T_FMT " chunks (%.2f KiB)",
					rightAlloc.getChunkCount(), rightAlloc.getSize() / 1024.0f);
			Log(EDebug, "      Nodes: " SIZE_T_FMT " chunks (%.2f KiB)",
					nodeAlloc.getChunkCount(), nodeAlloc.getSize() / 1024.0f);
		}

		void accumulateStatistics(const BuildContext &ctx) {
			leafNodeCount += ctx.leafNodeCount;
			nonemptyLeafNodeCount += ctx.nonemptyLeafNodeCount;
			innerNodeCount += ctx.innerNodeCount;
			retractedSplits += ctx.retractedSplits;
			pruned += ctx.pruned;
			expTraversalSteps += ctx.expTraversalSteps;
			expLeavesVisited += ctx.expLeavesVisited;
			expPrimitivesIntersected += ctx.expPrimitivesIntersected;
		}
	};

	/**
	 * \brief Communication data structure used to pass jobs to
	 * SAH kd-tree builder threads
	 */
	struct BuildInterface {
		/* Communcation */
		ref<Mutex> mutex;
		ref<ConditionVariable> cond;
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
	 * \brief SAH kd-tree builder thread
	 */
	class SAHTreeBuilder : public Thread {
	public:
		SAHTreeBuilder(size_type idx, size_type primCount, BuildInterface &interface) 
			: Thread(formatString("bld%i", idx)), m_context(primCount), m_interface(interface) { }

		void run() {
			m_interface.mutex->lock();
			while (!m_interface.done) {
				m_interface.cond->wait();
			}
			m_interface.mutex->unlock();
		}

		inline BuildContext &getContext() {
			return m_context;
		}

	private:
		BuildContext m_context;
		BuildInterface &m_interface;
	};

	/// Cast to the derived class
	inline Derived *downCast() {
		return static_cast<Derived *>(this);
	}

	/// Cast to the derived class (const version)
	inline const Derived *downCast() const {
		return static_cast<Derived *>(this);
	}

	/**
	 * \brief Create an edge event list for a given list of primitives. 
	 *
	 * This is necessary when passing from Min-Max binning to the more 
	 * accurate SAH-based optimizier.
	 */
	boost::tuple<EdgeEvent *, EdgeEvent *> createEventList(
			OrderedChunkAllocator &alloc, const AABB &nodeAABB, index_type *prims, size_type primCount) {
		size_type initialSize = primCount * 6;
		EdgeEvent *eventStart = alloc.allocate<EdgeEvent>(initialSize);
		EdgeEvent *eventEnd = eventStart;

		for (size_type i=0; i<primCount; ++i) {
			index_type index = prims[i];
			AABB aabb = downCast()->getClippedAABB(index, nodeAABB);

			for (int axis=0; axis<3; ++axis) {
				float min = (float) aabb.min[axis], max = (float) aabb.max[axis];

				if (min == max) {
					*eventEnd++ = EdgeEvent(EdgeEvent::EEdgePlanar, axis, min, index);
				} else {
					*eventEnd++ = EdgeEvent(EdgeEvent::EEdgeStart, axis, min, index);
					*eventEnd++ = EdgeEvent(EdgeEvent::EEdgeEnd, axis, max, index);
				}
			}
		}

		size_type newSize = eventEnd - eventStart;
		if (newSize != initialSize)
			alloc.shrinkAllocation<EdgeEvent>(eventStart, newSize);

		return boost::make_tuple(eventStart, eventEnd);
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
	 */
	void createLeaf(BuildContext &ctx, KDNode *node, const AABB &nodeAABB, size_type primCount) {
		node->initLeafNode(0, primCount);
		ctx.leafNodeCount++;
		ctx.expLeavesVisited += nodeAABB.getSurfaceArea();
		ctx.expPrimitivesIntersected += primCount * nodeAABB.getSurfaceArea();
		if (primCount > 0)
			ctx.nonemptyLeafNodeCount++;
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
			createLeaf(ctx, node, nodeAABB, primCount);
			return leafCost;
		}

		if (primCount <= m_exactPrimThreshold) {
			OrderedChunkAllocator &alloc = isLeftChild ? ctx.leftAlloc : ctx.rightAlloc;
			boost::tuple<EdgeEvent *, EdgeEvent *> events  
					= createEventList(alloc, nodeAABB, indices, primCount);
		
			std::sort(boost::get<0>(events), boost::get<1>(events), EdgeEventOrdering());

			Float sahCost = buildTreeSAH(ctx, depth, node, nodeAABB,
					boost::get<0>(events), boost::get<1>(events), primCount, 
					isLeftChild, badRefines);

			alloc.release(boost::get<0>(events));
			return sahCost;
		}

		/* ==================================================================== */
	    /*                              Binning                                 */
	    /* ==================================================================== */

		MinMaxBins<MTS_KD_MINMAX_BINS> bins(tightAABB);
		bins.bin(downCast(), indices, primCount);

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
				createLeaf(ctx, node, nodeAABB, primCount);
				return leafCost;
			}
			++badRefines;
		}

		/* ==================================================================== */
	    /*                            Partitioning                              */
	    /* ==================================================================== */

		boost::tuple<AABB, index_type *, AABB, index_type *> partition = 
			bins.partition(ctx, downCast(), indices, bestSplit, isLeftChild, 
			m_traversalCost, m_intersectionCost);

		/* ==================================================================== */
	    /*                              Recursion                               */
	    /* ==================================================================== */

		OrderedChunkAllocator &nodeAlloc = ctx.nodeAlloc;
		KDNode *children = nodeAlloc.allocate<KDNode>(2);

		size_t initialIndirectionTableSize = m_indirectionTable.size();
		if (!node->initInnerNode(bestSplit.axis, bestSplit.pos, children-node)) {
			/* Unable to store relative offset -- create an indirection
			   table entry */
			node->initIndirectionNode(bestSplit.axis, bestSplit.pos, 
					initialIndirectionTableSize);
			m_indirectionTable.push_back(children);
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
			ctx.expTraversalSteps += nodeAABB.getSurfaceArea();
			return finalSAHCost;
		} else {
			/* In the end, splitting didn't help to reduce the SAH cost.
			   Tear up everything below this node and create a leaf */
			if (m_indirectionTable.size() > initialIndirectionTableSize)
				m_indirectionTable.resize(initialIndirectionTableSize);
			tearUp(ctx, node);
			ctx.retractedSplits++;
			createLeaf(ctx, node, nodeAABB, primCount);
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
			createLeaf(ctx, node, nodeAABB, primCount);
			return leafCost;
		}

		SplitCandidate bestSplit;
		Float invSA = 1.0f / nodeAABB.getSurfaceArea();

		/* ==================================================================== */
	    /*                        Split candidate search                        */
    	/* ==================================================================== */

		/* First, find the optimal splitting plane according to the
		   surface area heuristic. To do this in O(n), the search is
		   implemented as a sweep over the edge events */
			
		/* Initially, the split plane is placed left of the scene
		   and thus all geometry is on its right side */
		size_type numLeft[3], numRight[3];
		AABB aabb(nodeAABB);
		for (int i=0; i<3; ++i) {
			numLeft[i] = 0;
			numRight[i] = primCount;
		}
		EdgeEvent *eventsByAxis[3];
		int eventsByAxisCtr = 1;
		eventsByAxis[0] = eventStart;

		/* Iterate over all events on the current axis */
		for (EdgeEvent *event = eventStart; event < eventEnd; ) {
			/* Record the current position and count all 
				other events, which are also here */
			uint16_t axis = event->axis;
			float pos = event->pos;
			size_type numStart = 0, numEnd = 0, numPlanar = 0;

			/* Count "end" events */
			while (event < eventEnd && event->pos == pos && event->axis == axis 
				&& event->type == EdgeEvent::EEdgeEnd) {
				++numEnd; ++event;
			}

			/* Count "planar" events */
			while (event < eventEnd && event->pos == pos && event->axis == axis 
				&& event->type == EdgeEvent::EEdgePlanar) {
				++numPlanar; ++event;
			}

			/* Count "start" events */
			while (event < eventEnd && event->pos == pos && event->axis == axis 
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
			if (EXPECT_TAKEN(pos >= nodeAABB.min[axis] && pos <= nodeAABB.max[axis])) {
				size_type nL = numLeft[axis], nR = numRight[axis];
				Float tmp = nodeAABB.max[axis];
				aabb.max[axis] = pos;
				Float pLeft = invSA * aabb.getSurfaceArea();
				aabb.max[axis] = tmp;
				tmp = aabb.min[axis];
				aabb.min[axis] = pos;
				Float pRight = invSA * aabb.getSurfaceArea();
				aabb.min[axis] = tmp;
				Float sahCostPlanarLeft = m_traversalCost + m_intersectionCost 
					* (pLeft * (nL + numPlanar) + pRight * nR);
				Float sahCostPlanarRight = m_traversalCost + m_intersectionCost 
					* (pLeft * nL + pRight * (nR + numPlanar));
				if (nL + numPlanar == 0 || nR == 0)
					sahCostPlanarLeft *= m_emptySpaceBonus;
				if (nL == 0 || nR + numPlanar == 0)
					sahCostPlanarRight *= m_emptySpaceBonus;

				if (sahCostPlanarLeft < bestSplit.sahCost || sahCostPlanarRight < bestSplit.sahCost) {
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
			} else {
				/* When primitive clipping is active, this should 
					never happen! */
				KDAssertEx(!m_clip, "Internal error: edge event is out of bounds");
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

		KDAssert(bestSplit.sahCost != std::numeric_limits<Float>::infinity());

		/* "Bad refines" heuristic from PBRT */
		if (bestSplit.sahCost >= leafCost) {
			if ((bestSplit.sahCost > 4 * leafCost && primCount < 16)
				|| badRefines >= m_maxBadRefines) {
				createLeaf(ctx, node, nodeAABB, primCount);
				return leafCost;
			}
			++badRefines;
		}

		/* ==================================================================== */
		/*                      Primitive Classification                        */
		/* ==================================================================== */

		ClassificationStorage &storage = ctx.storage;

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
				} else if (event->pos > bestSplit.pos || (event->pos == bestSplit.pos && 
					!bestSplit.planarLeft)) {
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
			EdgeEvent *leftEventsTempStart = leftAlloc.allocate<EdgeEvent>(primsLeft * 6),
					  *rightEventsTempStart = rightAlloc.allocate<EdgeEvent>(primsRight * 6),
					  *newEventsLeftStart = leftAlloc.allocate<EdgeEvent>(primsBoth * 6),
					  *newEventsRightStart = rightAlloc.allocate<EdgeEvent>(primsBoth * 6);

			EdgeEvent *leftEventsTempEnd = leftEventsTempStart, 
					*rightEventsTempEnd = rightEventsTempStart,
					*newEventsLeftEnd = newEventsLeftStart,
					*newEventsRightEnd = newEventsRightStart;

			for (EdgeEvent *event = eventStart; event<eventEnd; ++event) {
				uint8_t classification = storage.get(event->index);
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

					AABB clippedLeft = downCast()->getClippedAABB(index, leftNodeAABB);
					AABB clippedRight = downCast()->getClippedAABB(index, rightNodeAABB);

					KDAssert(leftNodeAABB.contains(clippedLeft));
					KDAssert(rightNodeAABB.contains(clippedRight));

					if (clippedLeft.isValid() && clippedLeft.getSurfaceArea() > 0) {
						for (int axis=0; axis<3; ++axis) {
							float min = (float) clippedLeft.min[axis], max = (float) clippedLeft.max[axis];

							if (min == max) {
								*newEventsLeftEnd++ = EdgeEvent(EdgeEvent::EEdgePlanar, axis, min, index);
							} else {
								*newEventsLeftEnd++ = EdgeEvent(EdgeEvent::EEdgeStart, axis, min, index);
								*newEventsLeftEnd++ = EdgeEvent(EdgeEvent::EEdgeEnd, axis, max, index);
							}
						}
					} else {
						prunedLeft++;
					}

					if (clippedRight.isValid() && clippedRight.getSurfaceArea() > 0) {
						for (int axis=0; axis<3; ++axis) {
							float min = (float) clippedRight.min[axis], max = (float) clippedRight.max[axis];

							if (min == max) {
								*newEventsRightEnd++ = EdgeEvent(EdgeEvent::EEdgePlanar, axis, min, index);
							} else {
								*newEventsRightEnd++ = EdgeEvent(EdgeEvent::EEdgeStart, axis, min, index);
								*newEventsRightEnd++ = EdgeEvent(EdgeEvent::EEdgeEnd, axis, max, index);
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
				uint8_t classification = storage.get(event->index);
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

		OrderedChunkAllocator &nodeAlloc = ctx.nodeAlloc;
		KDNode *children = nodeAlloc.allocate<KDNode>(2);

		size_t initialIndirectionTableSize = m_indirectionTable.size();
		if (!node->initInnerNode(bestSplit.axis, bestSplit.pos, children-node)) {
			/* Unable to store relative offset -- create an indirection
			   table entry */
			node->initIndirectionNode(bestSplit.axis, bestSplit.pos, 
					initialIndirectionTableSize);
			m_indirectionTable.push_back(children);
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
			ctx.expTraversalSteps += nodeAABB.getSurfaceArea();
			return finalSAHCost;
		} else {
			/* In the end, splitting didn't help to reduce the SAH cost.
			   Tear up everything below this node and create a leaf */
			if (m_indirectionTable.size() > initialIndirectionTableSize)
				m_indirectionTable.resize(initialIndirectionTableSize);
			tearUp(ctx, node);
			ctx.retractedSplits++;
			createLeaf(ctx, node, nodeAABB, primCount);
			return leafCost;
		}

		return bestSplit.sahCost;
	}

	/**
	 * \brief Tear up a subtree after a split did not reduce the SAH cost
	 */
	void tearUp(BuildContext &ctx, KDNode *node) {
		if (node->isLeaf()) {
			if (node->getPrimStart() != node->getPrimEnd())
				ctx.nonemptyLeafNodeCount--;
			ctx.leafNodeCount--;
			/// XXX Create primitive list for leaf
		} else {
			KDNode *left;
			ctx.innerNodeCount--;
			if (EXPECT_TAKEN(!node->isIndirection()))
				left = node->getLeft();
			else
				left = m_indirectionTable[node->getIndirectionIndex()];

			tearUp(ctx, left);
			tearUp(ctx, left+1);

			ctx.nodeAlloc.release(left);
		}
	}

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
						/* Insufficient floating point resolution -- a leaf will be created. */
						candidate.sahCost = std::numeric_limits<Float>::infinity();
						break;
					}

					++it;
				}
			}

			if (split <= m_aabb.min[axis] || split > m_aabb.max[axis]) {
				/* Insufficient floating point resolution -- a leaf will be created. */
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
				SplitCandidate &split, bool isLeftChild, Float traversalCost, Float intersectionCost) {
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
	bool m_clip, m_retract;
	AABB m_aabb;
	size_type m_maxDepth;
	size_type m_stopPrims;
	size_type m_maxBadRefines;
	size_type m_exactPrimThreshold;
	std::vector<SAHTreeBuilder *> m_builders;
	BuildInterface m_interface;
};

MTS_NAMESPACE_END

#endif /* __KDTREE_GENERIC_H */

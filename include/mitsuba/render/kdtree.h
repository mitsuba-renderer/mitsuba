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

#if !defined(__KD_TREE_H)
#define __KD_TREE_H

#include <mitsuba/render/shape.h>

/**
 * First, some default configuration definitions:
 *
 * If SSE is available, use the quad-packed TriAccel4 intersection method. 
 * Otherwise, fall back to the scalar TriAccel implementation. To use the 
 * scalar Moeller-Trumbore intersection algorithm (slower), remove this
 * whole following block:
 */
#ifdef MTS_SSE
#define MTS_USE_TRIACCEL4 1
#define MTS_USE_TRIACCEL 1
#else
#define MTS_USE_TRIACCEL 1
#endif

#if defined(MTS_HAS_COHERENT_RT)
/* Coherent ray tracing needs the plain TriAccel data stucture.
   Note that both TriAccel and TriAccel4 can be pre-computed
   at the same time in order to have a performance gain with
   both coherent and non-coherent methods (at the cost of higher 
   memory usage) */
#define MTS_USE_TRIACCEL 1
#if !defined(MTS_SSE)
#error Coherent ray tracing requires SSE
#endif
#endif

#if !defined(MTS_USE_TRIACCEL4) && !defined(MTS_USE_TRIACCEL)
/* Switch to Moeller-Trumbore if neither TriAccel4 nor TriAccel is selected */
#define MTS_USE_MT 1
#endif

/**
 * Pre-defined max. stack size for the ray traversal algorithm
 */
#define MTS_KD_MAXDEPTH 35

MTS_NAMESPACE_BEGIN

/**
 * SAH KD-tree acceleration data structure for fast ray-triangle intersections.
 * Implements the construction algorithm for 'perfect split' trees as outlined 
 * in the paper "On Bulding fast kd-Trees for Ray Tracing, and on doing that in 
 * O(N log N)" by Ingo Wald and Vlastimil Havran. Non-triangle shapes are also
 * supported, but most optimizations here target large triangle meshes.
 *
 * This class offers a choice of several different triangle intersection algorithms:
 * By default, intersections are computed using the "TriAccel" projection with 
 * pre-computation method from Ingo Wald's PhD thesis "Realtime Ray Tracing 
 * and Interactive Global Illumination". When SSE is available on the target system, 
 * leaf triangles are stored in an explicit, packed format so that up to four 
 * intersections can be performed simultaneously. While about 15% faster, this also 
 * tends to duplicate a lot of geometry, which results in high memory usage. 
 * Therefore, this behavior can optionally be turned off with a #define. The third 
 * choice is the Moeller-Trumbore intersection test, which requires the least amount
 * of memory, but is also the slowest.
 *
 * When SSE is enabled, packets of 4 rays can efficiently be traced to make 
 * use of any ray coherence. This requires "TriAccel" or "TriAccel4" to be
 * enabled.
 *
 * During the kd-tree construction, this class uses a technique named 
 * "primitive clipping" to significantly improve the quality of the resulting
 * trees. However, the involved Sutherland-Hodgman iterations are expensive
 * and can lead to long construction times. The setClip method can be used
 * to deactivate primitive clipping at the cost of slower intersections.
 *
 * Finally, this class also uses an optimized ray traversal algorithm 
 * (TA^B_{rec}), which is explained in Vlastimil Havran's PhD thesis 
 * "Heuristic Ray Shooting Algorithms". 
 *
 * @author Wenzel Jakob
 */
class MTS_EXPORT_RENDER KDTree : public Object {
public:
	/// Construct a new kd-tree in an unbuilt state
	KDTree();

	/// Add geometry to the kd-tree
	void addShape(const Shape *shape);
	
	/// Return the list of stored shapes
	inline const std::vector<const Shape *> &getShapes() const { return m_shapes; }

	/// Set the relative cost of an intersection operation
	inline void setIntersectionCost(Float pCost) { m_intersectionCost = pCost; }

	/// Return the relative cost of an intersection operation
	inline Float getIntersectionCost() const { return m_intersectionCost; }
	
	/// Set the relative cost of a traversal operation
	inline void setTraversalCost(Float pCost) { m_traversalCost = pCost; }

	/// Return the relative cost of a traversal operation
	inline Float getTraversalCost() const { return m_traversalCost; }

	/// Set a bonus factor for cutting away empty space
	inline void setEmptyBonus(Float pCost) { m_emptyBonus = pCost; }

	/// Return the bonus factor for cutting away empty space
	inline Float getEmptyBonus() const { return m_emptyBonus; }

	/// Set the min. number of primitives, which will never be split.
	inline void setStopPrims(int pCost) { m_stopPrims = pCost; }

	/// Return the min. number of primitives, which will never be split.
	inline int getStopPrims() const { return m_stopPrims; }

	/// Enable or disable primitive clipping
	inline void setClip(bool pClip) { m_clip = pClip; }

	/// Return whether primitive clipping is enabled
	inline bool getClip() const { return m_clip; }

	/// Has the kd-tree been built?
	inline bool isBuilt() const { return m_built; }

	/// Return the axis-aligned bounding box containing all primitives
	inline const AABB &getAABB() const { return m_rootBounds; }

	/// Return the bounding sphere containing all primitives
	inline const BSphere &getBSphere() const { return m_bsphere; }

	/// Build the kd-tree
	void build();

	/**
	 * Intersect a ray with the stored triangle meshes and only
	 * check for intersections. This is the fastest intersection test.
	 */
	bool rayIntersect(const Ray &ray) const;

	/**
	 * Intersect a ray with the stored triangle meshes and return
	 * a detailed intersection information record
	 */
	bool rayIntersect(const Ray &ray, Intersection &its) const;

	/**
	 * Intersect four rays with the stored triangle meshes while making
	 * use of ray coherence to do this very efficiently. If the coherent
	 * ray tracing #define is missing, this function simply does four
	 * separate mono-ray traversals.
	 */
	void rayIntersectPacket(const Ray *rays, Intersection *its) const;

#if defined(MTS_HAS_COHERENT_RT)
	/**
	 * Intersect four rays with the stored triangle meshes while making
	 * use of ray coherence to do this very efficiently. Requires SSE.
	 */
	void rayIntersectPacket(const RayPacket4 &packet, 
		const RayInterval4 &interval, Intersection4 &its) const;

	/**
	 * Fallback for incoherent rays
	 */
	void rayIntersectPacketIncoherent(const RayPacket4 &packet, 
		const RayInterval4 &interval, Intersection4 &its) const;
#endif

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~KDTree();

	/// \cond

	/// KD-tree node in 8 bytes
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
				   30-0 : Offset to the node's triangle list
				*/
				uint32_t combined;

				/// End offset of the triangle list
				uint32_t end;
			} leaf;
		};

		enum EMask {
			ETypeMask = 1 << 31,
			ELeafOffsetMask = ~ETypeMask,
			EInnerAxisMask = 3,
			EInnerOffsetMask = ~EInnerAxisMask
		};

		inline void setLeaf(unsigned int offset, unsigned int numPrims) {
			leaf.combined = ETypeMask | offset;
			leaf.end = offset + numPrims;
		}

		inline void setInner(int axis, unsigned int offset, Float split) {
			inner.combined = axis | (offset << 2);
			inner.split = (float) split;
		}

		FINLINE bool isLeaf() const {
			return leaf.combined & ETypeMask;
		}

		FINLINE uint32_t getPrimStart() const {
			return leaf.combined & ELeafOffsetMask;
		}

		FINLINE uint32_t getPrimEnd() const {
			return leaf.end;
		}

		FINLINE const KDNode * __restrict getLeft() const {
			return this + 
				((inner.combined & EInnerOffsetMask) >> 2);
		}

		FINLINE const KDNode * __restrict getOtherChild() const {
			return (const KDNode *) ((ptrdiff_t) this ^ (ptrdiff_t) 8);
		}

		FINLINE const KDNode * __restrict getRight() const {
			return getLeft() + 1;
		}

		inline float getSplit() const {
			return inner.split;
		}

		inline int getAxis() const {
			return inner.combined & EInnerAxisMask;
		}
	};

	/// Primitive classification during tree-construction
	enum EClassificationResult {
		EBothSides = 1,
		ELeftSide = 2,
		ERightSide = 3,
		EBothSidesProcessed = 4
	};

#if defined(MTS_USE_TRIACCEL)
	typedef TriAccel KDTriangle;
#else
	struct KDTriangle {
		uint32_t k;
		uint32_t index;
		uint32_t shapeIndex;
	};
#endif

	/// AABB edge event data structure
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
		inline EdgeEvent(uint8_t type, Float t, int index)
		 : t(t), index(index), type(type) {
		}

		/* Plane position */
		Float t;
		/* Triangle index */
		int index;
		/* Event type: end/planar/start */
		uint8_t type;
	};
	
	typedef std::vector<EdgeEvent> EdgeEventVec;
	typedef EdgeEventVec EdgeEventVec3[3];

	/// Edge event comparison functor
	struct EdgeEventSorter : public std::binary_function<EdgeEvent, EdgeEvent, bool> {
		inline bool operator()(const EdgeEvent &a, const EdgeEvent &b) const {
			if (a.t != b.t)
				return a.t < b.t;
			return a.type < b.type;
		}
	};

	/// Score for a split candidate
	struct Score {
		enum EPlanarSide {
			EPlanarLeft = 0,
			EPlanarRight = 1
		};

		/// Create an upper bound score
		inline Score() : score(std::numeric_limits<Float>::infinity()) {
		}

		/// Create a new store record
		inline Score(Float score, Float t, int nLeft, int nRight, 
			uint8_t axis, uint8_t planarSide) : score(score), t(t), 
			nLeft(nLeft), nRight(nRight), axis(axis), planarSide(planarSide) {
		}

		/// Return a string representation
		std::string toString() const {
			std::ostringstream oss;
			oss << "Score[value=" << score << ", t=" << t
				<< ", axis=" << (int) axis << ", planarSide="
				<< (planarSide == EPlanarLeft ? "left" : "right")
				<< ", nLeft=" << nLeft << ", nRight=" << nRight
				<< "]";
			return oss.str();
		}

		/// Numerical score value
		Float score;
		/// Split position
		Float t;
		/// Primitive counts on the left and right side
		int nLeft, nRight;
		/// Split axis
		uint8_t axis;
		/// Should planar prims be placed left or right?
		uint8_t planarSide;
		/// Score comparison operator
		inline bool operator<(const Score &b) const {
			return score < b.score;
		}
	};

	/// Surface area heuristic (SAH) cost function
	inline Float SAH(Float prLeft, Float prRight, int numLeft, int numRight) {
		Float cost = m_traversalCost + m_intersectionCost 
			* (prLeft * numLeft + prRight * numRight);

		/* Favor splits, which cut off empty regions of space */
		if (numLeft == 0 || numRight == 0)
			cost *= m_emptyBonus;
		return cost;
	}

	/// Surface area heuristic
	inline Score SAH(int axis, Float invArea, AABB aabb, Float t, int numLeft, 
			int numRight, int numPlanar) {
		/* Generate the left+right node bounding boxes */
		AABB leftBounds = aabb, rightBounds = aabb;
		leftBounds.max[axis] = t; rightBounds.min[axis] = t;

		if (std::abs(aabb.min[axis]-t) < Epsilon ||
			std::abs(aabb.max[axis]-t) < Epsilon) {
			/* Do not allow tiny splits */
			return Score(std::numeric_limits<Float>::infinity(), 
				t, numLeft, numRight+numPlanar,  axis, Score::EPlanarRight);
		}

		/* Determinate approximate intersection probabilities
		   for uniformly distributed rays */
		Float prLeft = leftBounds.getSurfaceArea() * invArea,
		      prRight = rightBounds.getSurfaceArea() * invArea;

		/* Two costs are calculated depending on whether planar
		   primitives are put on the left or right side */
		Float costPlanarLeft = SAH(prLeft, prRight,
			numLeft + numPlanar, numRight);
		Float costPlanarRight = SAH(prLeft, prRight,
			numLeft, numRight + numPlanar);

		if (costPlanarLeft < costPlanarRight)
			return Score(costPlanarLeft, t, numLeft+numPlanar, numRight, 
				axis, Score::EPlanarLeft);
		else
			return Score(costPlanarRight, t, numLeft, numRight+numPlanar, 
				axis, Score::EPlanarRight);
	}

	/// Ray traversal stack entry for incoherent ray tracing
	struct KDStackEntry {
		/* Pointer to the far child */
		const KDNode * __restrict node;
		/* Distance traveled along the ray (entry or exit) */
		Float t;
		/* Previous stack item */
		uint32_t prev;
		/* Associated point */
		Point pb;
	};
	
#if defined(MTS_HAS_COHERENT_RT)
	/// Ray traversal stack entry for uncoherent ray tracing
	struct CoherentKDStackEntry {
		/* Current ray interval */
		RayInterval4 MM_ALIGN16 interval;
		/* Pointer to the far child */
		const KDNode * __restrict node;
	};
#endif
	/// \endcond

	/**
	 * Intersection method from Vlastimil Havran's 
	 * PhD thesis (algorithm TA^B_{rec})
	 */
	bool rayIntersect(const Ray &ray, Intersection &its, Float mint, Float maxt, 
		bool shadowRay, unsigned int &shapeIndex, unsigned int &primIndex) const;
	
	/// Recursive tree-building algorithm
	void buildTree(int nodeIndex, int depth, int badRefines,
		int numPrims, const AABB &aabb, EdgeEventVec3 &allEvents);

	/// Create a leaf kd-tree node
	void createLeaf(int nodeIndex, int depth, int numPrims,
		EdgeEventVec3 &allEvents);

	/**
	 * Fallback for incoherent rays
	 */
	inline void rayIntersectPacketIncoherent(const Ray *rays, Intersection *its) const {
		for (int i=0; i<4; i++) {
			if (!rayIntersect(rays[i], its[i]))
				its[i].t = std::numeric_limits<float>::infinity();
		}
	}
protected:
	/* Has the kd-tree been built yet? */
	bool m_built;
	/* Axis-aligned bounding box of the root node */
	AABB m_rootBounds;
	/* Bounding sphere of the root node */
	BSphere m_bsphere;
	/* Pointers to all contained shapes */
	std::vector<const Shape *> m_shapes;
	/// Storage for kd-tree nodes
	std::vector<KDNode, aligned_allocator<KDNode> > m_nodes;
#if defined(MTS_USE_TRIACCEL4)
	/// Explicitly store per-leaf quad-packed triangles
	TriAccel4 *m_packedTriangles;
	std::vector<TriAccel4> m_tempPackedTriangles;
	unsigned int m_packedTriangleCount;
#endif
#if defined(MTS_USE_TRIACCEL) || defined(MTS_USE_MT)
	/// Storage for triangle redirection indices
	std::vector<unsigned int> m_indices;
#endif
	/// Total number of triangles/non-triangles in the tree
	unsigned int m_triangleCount, m_nonTriangleCount, m_primitiveCount;
	/// Geometry storage
	KDTriangle* m_triangles;
	/// Cost values for the surface are heuristic
	Float m_traversalCost, m_intersectionCost, m_emptyBonus;
	/// Ad-hoc depth cutoff value when building the tree
	int m_maxDepth;
	/// Allowed number of bad refines before creating a leaf
	int m_maxBadRefines;
	/// Minimal number of primitives per node
	int m_stopPrims;
	/// Use primitive clipping?
	bool m_clip;
	/// Storage for edge events
	EdgeEventVec3 m_events, *m_rightEvents;
	/// Short-term storage for edge events completely located left or right
	EdgeEventVec m_eventsL, m_eventsR;
	/// Short-term storage for edge events from by prims overlapping the split
	EdgeEventVec m_sEventsL, m_sEventsR;
	/// Statistics
	int m_leafNodes, m_innerNodes;
	int m_minLeafPrims, m_maxLeafPrims;
	int m_totalLeafDepth, m_actualMaxDepth;
	int m_leafPrims, m_clippedAway;
	int m_badRefines, m_bucketCount;
	int m_failures, *m_triBuckets;
};

MTS_NAMESPACE_END

#endif /* __KD_TREE_H */

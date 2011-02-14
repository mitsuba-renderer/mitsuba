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

#if !defined(__KDTREE_H)
#define __KDTREE_H

#include <mitsuba/core/aabb.h>
#include <boost/foreach.hpp>

MTS_NAMESPACE_BEGIN

/**
 * \brief Basic point node for use with \ref TKDTree.
 *
 * \tparam PointType Underlying point data type (e.g. \ref TPoint3<float>)
 * \tparam DataRecord Custom payload to be attached to each node
 */
template <typename PointType, typename DataRecord> struct BasicKDNode {
	typedef PointType point_type;

	PointType position;
	uint32_t right;
	uint16_t flags;
	uint16_t axis;
	DataRecord value;

	/// Initialize a KD-tree node
	inline BasicKDNode() : position((Float) 0), 
		right(0), flags(0), axis(0), value() { }
	/// Initialize a KD-tree node with the given data record
	inline BasicKDNode(const DataRecord &value) : position((Float) 0), 
		right(0), flags(0), axis(0), value(value) { }

	/// Given the current node's index, return the index of the right child 
	inline uint32_t getRightIndex(uint32_t curIndex) { return right; }
	/// Given the current node's index, return the index of the right child (const version)
	inline const uint32_t getRightIndex(uint32_t curIndex) const { return right; }
	/// Given the current node's index, set the right child index
	inline void setRightIndex(uint32_t curIndex, uint32_t value) { right = value; }

	/// Given the current node's index, return the index of the left child 
	inline uint32_t getLeftIndex(uint32_t curIndex) { return curIndex + 1; }
	/// Given the current node's index, return the index of the left child (const version)
	inline const uint32_t getLeftIndex(uint32_t curIndex) const { return curIndex + 1; }
	/// Given the current node's index, set the left child index
	inline void setLeftIndex(uint32_t curIndex, uint32_t value) {
		if (value != curIndex+1) SLog(EError, "Not supported!"); }

	/// Check whether this is a leaf node
	inline bool isLeaf() const { return flags & 1; }
	/// Specify whether this is a leaf node
	inline void setLeaf(bool value) { if (value) flags |= 1; else flags &= ~1; }

	/// Return the split axis associated with this node
	inline uint16_t getAxis() const { return axis; }
	/// Set the split axis associated with this node
	inline void setAxis(uint16_t value) { axis = value; }

	/// Return the position associated with this node
	inline const PointType &getPosition() const { return position; }
	/// Set the position associated with this node
	inline void setPosition(const PointType &value) { position = value; }

	/// Return the data record associated with this node
	inline DataRecord &getValue() { return value; }
	/// Return the data record associated with this node (const version)
	inline const DataRecord &getValue() const { return value; }
	/// Set the data record associated with this node
	inline void setValue(const DataRecord &val) { value = val; }
};

/**
 * \brief Generic multi-dimensional kd-tree data structure for point data
 *
 * Organizes a list of point data in a hierarchical manner. For data
 * with spatial extents, \ref GenericKDTree and \ref ShapeKDTree will be
 * more appropriate.
 *
 * \tparam KDNode Underlying node data structure. See \ref BasicKDNode as
 * an example for the required public interface
 */
template <typename KDNode> class TKDTree {
public:
	typedef KDNode                           node_type;
	typedef typename KDNode::point_type      point_type;
	typedef typename point_type::value_type  value_type;
	typedef typename point_type::vector_type vector_type;
	typedef TAABB<point_type>                aabb_type;

	/// Supported tree construction heuristics
	enum EHeuristic {
		/// Create a balanced tree by splitting along the median
		EBalanced = 0,
		
		/// Create a left-balanced tree
		ELeftBalanced,

		/**
		 * \brief Use the sliding midpoint tree construction rule. This 
		 * ensures that cells do not become overly elongated.
		 */
		ESlidingMidpoint,

		/**
		 * \brief Choose the split plane by optimizing a cost heuristic
		 * based on the ratio of voxel volumes. Note that the implementation
		 * here is not particularly optimized, and it furthermore runs in
		 * time O(n (log n)^2) instead of O(n log n)
		 */
		EVoxelVolume
	};

	/// Result data type for k-nn queries
	struct SearchResult {
		Float distSquared;
		uint32_t index;

		inline SearchResult(Float distSquared, uint32_t index)
			: distSquared(distSquared), index(index) { }

		std::string toString() const {
			std::ostringstream oss;
			oss << "SearchResult[distance=" << std::sqrt(distSquared) 
				<< ", index=" << index << "]";
			return oss.str();
		}

		inline bool operator==(const SearchResult &r) const {
			return distSquared == r.distSquared &&
				index == r.index;
		}
	};

	/// Comparison functor for nearest-neighbor search queries
	struct SearchResultComparator : public 
		std::binary_function<SearchResult, SearchResult, bool> {
	public:
		inline bool operator()(const SearchResult &a, const SearchResult &b) const {
			return a.distSquared < b.distSquared;
		}
	};

public:
	/** 
	 * \brief Create an empty KD-tree that can hold the specified
	 * number of points
	 */
	inline TKDTree(size_t nodes, EHeuristic heuristic = ESlidingMidpoint)
		: m_nodes(nodes), m_heuristic(heuristic), m_depth(0) { }

	/// Return one of the KD-tree nodes by index
	inline KDNode &operator[](size_t idx) { return m_nodes[idx]; }
	/// Return one of the KD-tree nodes by index (const version)
	inline const KDNode &operator[](size_t idx) const { return m_nodes[idx]; }

	/// Return the AABB of the underlying point data
	inline const aabb_type &getAABB() const { return m_aabb; }
	/// Return the depth of the constructed KD-tree
	inline size_t getDepth() const { return m_depth; }
	/// Return the size of the kd-tree
	inline size_t getSize() const { return m_nodes.size(); }

	/// Construct the KD-tree hierarchy
	void build() {
		m_aabb.reset();

		BOOST_FOREACH(KDNode &node, m_nodes) {
			m_aabb.expandBy(node.getPosition());
		}

		m_depth = 0;
		build(1, m_nodes.begin(), m_nodes.end());
	}

	/**
	 * \brief Run a k-nearest-neighbor search query
	 *
	 * \param p Search position
	 * \param k Maximum number of search results
	 * \param result Index list of search results
	 * \param searchRadius   Maximum search radius (this can be used to
	 *      restrict the knn query to a subset of the data)
	 * \return The number of used traversal steps
	 */
	size_t nnSearch(const point_type &p, size_t k, std::vector<SearchResult> &results, 
			Float searchRadius = std::numeric_limits<Float>::infinity()) const {
		uint32_t *stack = (uint32_t *) alloca((m_depth+1) * sizeof(uint32_t));
		size_t index = 0, stackPos = 1, traversalSteps = 0;
		bool isHeap = false;
		Float distSquared = searchRadius*searchRadius;
		stack[0] = 0;

		results.clear();
		results.reserve(k+1);
	
		while (stackPos > 0) {
			const KDNode &node = m_nodes[index];
			++traversalSteps;
			int nextIndex;
	
			/* Recurse on inner nodes */
			if (!node.isLeaf()) {
				Float distToPlane = p[node.getAxis()] 
					- node.getPosition()[node.getAxis()];
	
				uint32_t first, second;
				bool searchBoth = distToPlane*distToPlane <= distSquared;

				if (distToPlane > 0) {
					first = node.getRightIndex(index);
					second = searchBoth ? node.getLeftIndex(index) : 0;
				} else {
					first = node.getLeftIndex(index);
					second = searchBoth ? node.getRightIndex(index) : 0;
				}

				if (first != 0 && second != 0) {
					nextIndex = first;
					stack[stackPos++] = second;
				} else if (first != 0) {
					nextIndex = first;
				} else if (second != 0) {
					nextIndex = second;
				} else {
					nextIndex = stack[--stackPos];
				}
			} else {
				nextIndex = stack[--stackPos];
			}
	
			/* Check if the current point is within the query's search radius */
			const Float pointDistSquared = (node.getPosition() - p).lengthSquared();
	
			if (pointDistSquared < distSquared) {
				/* Switch to a max-heap when the available search 
				   result space is exhausted */
				if (results.size() < k) {
					/* There is still room, just add the point to
					   the search result list */
					results.push_back(SearchResult(pointDistSquared, index));
				} else {
					if (!isHeap) {
						/* Establish the max-heap property */
						std::make_heap(results.begin(), results.end(), 
								SearchResultComparator());
						isHeap = true;
					}
	
					/* Add the new point, remove the one that is farthest away */
					results.push_back(SearchResult(pointDistSquared, index));
					std::push_heap(results.begin(), results.end(), SearchResultComparator());
					std::pop_heap(results.begin(), results.end(), SearchResultComparator());
					results.pop_back();
	
					/* Reduce the search radius accordingly */
					distSquared = results[0].distSquared;
				}
			}
			index = nextIndex;
		}
		return traversalSteps;
	}

	/**
	 * \brief Execute a search query and run the specified functor on them,
	 * while potentially modifying nodes within the search radius
	 *
	 * The functor must have an operator() implementation, which accepts
	 * a \a KDNode as its argument.
	 *
	 * \param p Search position
	 * \param result Index list of search results
	 * \param searchRadius Search radius 
	 * \return The number of used traversal steps
	 */
	template <typename Functor> size_t executeModifier(const point_type &p,
			Float searchRadius, Functor &functor) {
		uint32_t *stack = (uint32_t *) alloca((m_depth+1) * sizeof(uint32_t));
		size_t index = 0, stackPos = 1, traversalSteps = 0;
		Float distSquared = searchRadius*searchRadius;
		stack[0] = 0;

		while (stackPos > 0) {
			KDNode &node = m_nodes[index];
			++traversalSteps;
			int nextIndex;
	
			/* Recurse on inner nodes */
			if (!node.isLeaf()) {
				Float distToPlane = p[node.getAxis()] 
					- node.getPosition()[node.getAxis()];
	
				uint32_t first, second;
				bool searchBoth = distToPlane*distToPlane <= distSquared;

				if (distToPlane > 0) {
					first = node.getRightIndex(index);
					second = searchBoth ? node.getLeftIndex(index) : 0;
				} else {
					first = node.getLeftIndex(index);
					second = searchBoth ? node.getRightIndex(index) : 0;
				}

				if (first != 0 && second != 0) {
					nextIndex = first;
					stack[stackPos++] = second;
				} else if (first != 0) {
					nextIndex = first;
				} else if (second != 0) {
					nextIndex = second;
				} else {
					nextIndex = stack[--stackPos];
				}
			} else {
				nextIndex = stack[--stackPos];
			}
	
			/* Check if the current point is within the query's search radius */
			const Float pointDistSquared = (node.getPosition() - p).lengthSquared();
	
			if (pointDistSquared < distSquared)
				functor(node);

			index = nextIndex;
		}
		return traversalSteps;
	}

	/**
	 * \brief Execute a search query and run the specified functor on them
	 *
	 * The functor must have an operator() implementation, which accepts
	 * a constant reference to a \a KDNode as its argument.
	 *
	 * \param p Search position
	 * \param result Index list of search results
	 * \param searchRadius  Search radius 
	 * \return The number of used traversal steps
	 */
	template <typename Functor> size_t executeQuery(const point_type &p,
			Float searchRadius, Functor &functor) const {
		uint32_t *stack = (uint32_t *) alloca((m_depth+1) * sizeof(uint32_t));
		size_t index = 0, stackPos = 1, traversalSteps = 0;
		Float distSquared = searchRadius*searchRadius;
		stack[0] = 0;

		while (stackPos > 0) {
			const KDNode &node = m_nodes[index];
			++traversalSteps;
			int nextIndex;
	
			/* Recurse on inner nodes */
			if (!node.isLeaf()) {
				Float distToPlane = p[node.getAxis()] 
					- node.getPosition()[node.getAxis()];
	
				uint32_t first, second;
				bool searchBoth = distToPlane*distToPlane <= distSquared;

				if (distToPlane > 0) {
					first = node.getRightIndex(index);
					second = searchBoth ? node.getLeftIndex(index) : 0;
				} else {
					first = node.getLeftIndex(index);
					second = searchBoth ? node.getRightIndex(index) : 0;
				}

				if (first != 0 && second != 0) {
					nextIndex = first;
					stack[stackPos++] = second;
				} else if (first != 0) {
					nextIndex = first;
				} else if (second != 0) {
					nextIndex = second;
				} else {
					nextIndex = stack[--stackPos];
				}
			} else {
				nextIndex = stack[--stackPos];
			}
	
			/* Check if the current point is within the query's search radius */
			const Float pointDistSquared = (node.getPosition() - p).lengthSquared();
	
			if (pointDistSquared < distSquared)
				functor(node);

			index = nextIndex;
		}
		return traversalSteps;
	}


	/**
	 * \brief Run a search query
	 *
	 * \param p Search position
	 * \param result Index list of search results
	 * \param searchRadius  Search radius 
	 * \return The number of used traversal steps
	 */
	size_t search(const point_type &p, Float searchRadius, std::vector<uint32_t> &results) const {
		uint32_t *stack = (uint32_t *) alloca((m_depth+1) * sizeof(uint32_t));
		size_t index = 0, stackPos = 1, traversalSteps = 0;
		Float distSquared = searchRadius*searchRadius;
		stack[0] = 0;

		results.clear();
	
		while (stackPos > 0) {
			const KDNode &node = m_nodes[index];
			++traversalSteps;
			int nextIndex;
	
			/* Recurse on inner nodes */
			if (!node.isLeaf()) {
				Float distToPlane = p[node.getAxis()] 
					- node.getPosition()[node.getAxis()];
	
				uint32_t first, second;
				bool searchBoth = distToPlane*distToPlane <= distSquared;

				if (distToPlane > 0) {
					first = node.getRightIndex(index);
					second = searchBoth ? node.getLeftIndex(index) : 0;
				} else {
					first = node.getLeftIndex(index);
					second = searchBoth ? node.getRightIndex(index) : 0;
				}

				if (first != 0 && second != 0) {
					nextIndex = first;
					stack[stackPos++] = second;
				} else if (first != 0) {
					nextIndex = first;
				} else if (second != 0) {
					nextIndex = second;
				} else {
					nextIndex = stack[--stackPos];
				}
			} else {
				nextIndex = stack[--stackPos];
			}
	
			/* Check if the current point is within the query's search radius */
			const Float pointDistSquared = (node.getPosition() - p).lengthSquared();
	
			if (pointDistSquared < distSquared) 
				results.push_back(index);

			index = nextIndex;
		}
		return traversalSteps;
	}

protected:
	struct CoordinateOrdering : public std::binary_function<KDNode, KDNode, bool> {
	public:
		inline CoordinateOrdering(int axis) : m_axis(axis) { }
		inline bool operator()(const KDNode &n1, const KDNode &n2) const {
			return n1.getPosition()[m_axis] < n2.getPosition()[m_axis];
		}
	private:
		int m_axis;
	};

	struct LessThanOrEqual : public std::unary_function<KDNode, bool> {
	public:
		inline LessThanOrEqual(int axis, value_type value) : m_axis(axis), m_value(value) { }
		inline bool operator()(const KDNode &n1) const {
			return n1.getPosition()[m_axis] <= m_value;
		}
	private:
		int m_axis;
		value_type m_value;
	};

	void build(size_t depth,
			  typename std::vector<KDNode>::iterator rangeStart, 
			  typename std::vector<KDNode>::iterator rangeEnd) {
		m_depth = std::max(depth, m_depth);
		if (rangeEnd-rangeStart <= 0) {
			SLog(EError, "Internal error!");
		} else if (rangeEnd-rangeStart == 1) {
			/* Create a leaf node */
			rangeStart->setLeaf(true);
			return;
		}

		int axis = 0;
		typename std::vector<KDNode>::iterator split;

		switch (m_heuristic) {
			case EBalanced: {
					/* Split along the median */
					split = rangeStart + (rangeEnd-rangeStart)/2;
					axis = m_aabb.getLargestAxis();
					std::nth_element(rangeStart, split, rangeEnd, CoordinateOrdering(axis));
				};
				break;

			case ELeftBalanced: {
					size_t treeSize = rangeEnd-rangeStart;
					/* Layer 0 contains one node */
					size_t p = 1;

					/* Traverse downwards until the first incompletely
					   filled tree level is encountered */
					while (2*p <= treeSize)
						p *= 2;

					/* Calculate the number of filled slots in the last level */
					size_t remaining = treeSize - p + 1;

					if (2*remaining < p) {
						/* Case 2: The last level contains too few nodes. Remove
						   overestimate from the left subtree node count and add
						   the remaining nodes */
						p = (p >> 1) + remaining;
					}

					axis = m_aabb.getLargestAxis();
					
					split = rangeStart + (p - 1);
					std::nth_element(rangeStart, split, rangeEnd,
						CoordinateOrdering(axis));
				};
				break;

			case ESlidingMidpoint: {
					/* Sliding midpoint rule: find a split that is close to the spatial median */
					axis = m_aabb.getLargestAxis();

					value_type midpoint = (value_type) 0.5f 
						* (m_aabb.max[axis]+m_aabb.min[axis]);

					size_t nLT = std::count_if(rangeStart, rangeEnd,
							LessThanOrEqual(axis, midpoint));

					/* Re-adjust the split to pass through a nearby point */
					split = rangeStart + nLT;

					if (split == rangeStart)
						++split;
					else if (split == rangeEnd)
						--split;

					std::nth_element(rangeStart, split, rangeEnd,
						CoordinateOrdering(axis));
				};
				break;
	
			case EVoxelVolume: {
					Float bestCost = std::numeric_limits<Float>::infinity();

					for (int dim=0; dim<point_type::dim; ++dim) {
						std::sort(rangeStart, rangeEnd, CoordinateOrdering(dim));

						size_t numLeft = 1, numRight = rangeEnd-rangeStart-2;
						aabb_type leftAABB(m_aabb), rightAABB(m_aabb);
						Float invVolume = 1.0f / m_aabb.getVolume();
						for (typename std::vector<KDNode>::iterator it = rangeStart+1; it != rangeEnd; ++it) {
							++numLeft; --numRight;
							leftAABB.max[dim] = it->getPosition()[dim];
							rightAABB.min[dim] = it->getPosition()[dim];

							Float cost = (numLeft * leftAABB.getVolume()
								+ numRight * rightAABB.getVolume()) * invVolume;
							if (cost < bestCost) {
								bestCost = cost;
								axis = dim;
								split = it;
							}
						}
					}
					std::nth_element(rangeStart, split, rangeEnd,
						CoordinateOrdering(axis));
				};
				break;
		}

		value_type splitPos = split->getPosition()[axis];
		split->setAxis(axis);
		if (split+1 != rangeEnd) 
			split->setRightIndex(rangeStart - m_nodes.begin(), (uint32_t) (split + 1 - m_nodes.begin()));
		else 
			split->setRightIndex(rangeStart - m_nodes.begin(), 0);

		split->setLeftIndex(rangeStart - m_nodes.begin(), rangeStart + 1 - m_nodes.begin());
		split->setLeaf(false);
		std::iter_swap(rangeStart, split);

		/* Recursively build the children */
		value_type temp = m_aabb.max[axis];
		m_aabb.max[axis] = splitPos;
		build(depth+1, rangeStart+1, split+1);
		m_aabb.max[axis] = temp;

		if (split+1 != rangeEnd) {
			temp = m_aabb.min[axis];
			m_aabb.min[axis] = splitPos;
			build(depth+1, split+1, rangeEnd);
			m_aabb.min[axis] = temp;
		}
	}
protected:
	std::vector<KDNode> m_nodes;
	aabb_type m_aabb;
	EHeuristic m_heuristic;
	size_t m_depth;
};

MTS_NAMESPACE_END

#endif /* __KDTREE_H */

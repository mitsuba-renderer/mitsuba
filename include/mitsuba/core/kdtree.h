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

MTS_NAMESPACE_BEGIN

/**
 * \brief Basic point node for use with \ref TKDTree.
 *
 * \tparam PointType Underlying point data type (e.g. \ref TPoint3<float>)
 * \tparam DataRecord Custom payload to be attached to each node
 */
template <typename PointType, typename DataRecord> struct BasicKDNode {
	PointType position;
	BasicKDNode *left, *right;
	DataRecord data;
	uint8_t axis;

	/// Initialize a KD-tree node with the given data record
	inline BasicKDNode(const DataRecord &data) : position(0), 
		left(NULL), right(NULL), data(data), axis(0) { }

	/// Return a pointer to the left child of this node
	inline BasicKDNode *getLeft() { return left; }
	/// Return a pointer to the left child of this node (const version)
	inline const BasicKDNode *getLeft() const { return left; }
	/// Set the left child of this node
	inline void setLeft(BasicKDNode *node) { left = node; }

	/// Return a pointer to the right child of this node
	inline BasicKDNode *getRight() { return right; }
	/// Return a pointer to the right child of this node (const version)
	inline const BasicKDNode *getRight() const { return right; }
	/// Set the right child of this node
	inline void setRight(BasicKDNode *node) { right = node; }

	/// Return the split axis associated with this node
	inline uint8_t getAxis() const { return axis; }
	/// Set the split axis associated with this node
	inline void setAxis(uint8_t value) { axis = value; }

	/// Return the position associated with this node
	inline const PointType &getPosition() const { return position; }
	/// Set the position associated with this node
	inline void setPosition(const PointType &value) { position = value; }

	/// Return the data record associated with this node
	inline DataRecord &getData() { return data; }
	/// Return the data record associated with this node (const version)
	inline const DataRecord &getData() const { return data; }
	/// Set the data record associated with this node
	inline void setData(const DataRecord &value) { data = value; }
};

/**
 * \brief Generic multi-dimensional kd-tree data structure for point data
 * using the sliding midpoint tree construction rule. This ensures that
 * cells do not become overly elongated.
 *
 * Organizes a list of point data in a hierarchical manner. For data
 * with spatial extents, \ref GenericKDTree and \ref ShapeKDTree will be
 * more appropriate.
 *
 * \tparam PointType Underlying point data type (e.g. \ref TPoint3<float>)
 * 
 * \tparam KDNode Underlying node data structure. See \ref BasicKDNode as
 * an example for the required public interface
 */
template <typename PointType, typename KDNode> class TKDTree {
	typedef typename PointType::value_type  value_type;
	typedef typename PointType::vector_type vector_type;
	typedef TAABB<PointType>                aabb_type;

	/** 
	 * \brief Create an empty KD-tree that can hold the specified
	 * number of points
	 */
	inline TKDTree(size_t nodes) : m_nodes(nodes) {}

	/// Return one of the KD-tree nodes by index (const version)
	inline const KDNode &operator[](size_t idx) const { return m_nodes[idx]; }
	/// Return one of the KD-tree nodes by index
	inline KDNode &operator[](size_t idx) { return m_nodes[idx]; }


	void build() {
		m_aabb.reset();

		BOOST_FOREACH(KDNode &node, m_nodes) {
			m_aabb.expandBy(node.getPosition());
		}

		build(m_nodes.begin(), m_nodes.end());
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

	KDNode *build(typename std::vector<KDNode>::iterator rangeStart, 
				  typename std::vector<KDNode>::iterator rangeEnd) {
		SAssert(rangeEnd > rangeStart);

		if (rangeEnd-rangeStart <= 1) {
			/* Create a leaf node */
			return;
		}
	
		/* Find a split that is close to the spatial median */
		int axis = m_aabb.getLargestAxis();
		value_type midpoint = (value_type) 0.5f * (m_aabb.max[axis]+m_aabb.min[axis]);

		size_t nLT = std::count_if(rangeStart, rangeEnd,
				LessThanOrEqual(axis, midpoint));

		/* Re-adjust the split to pass through a nearby photon */
		typename std::vector<KDNode>::iterator split = rangeStart + nLT;
		std::nth_element(rangeStart, split, rangeEnd,
			CoordinateOrdering(axis));
		value_type splitPos = split->getPosition()[axis];
		split->setAxis(axis);

		/* Recursively build the children */
		value_type temp = m_aabb.max[axis];

		m_aabb.max[axis] = splitPos;
		split->setLeft(build(rangeStart, split));
		m_aabb.max[axis] = temp;

		temp = m_aabb.min[axis];
		m_aabb.min[axis] = splitPos;
		split->setRight(build(split+1, rangeEnd));
		m_aabb.min[axis] = temp;

		return split;
	}
protected:
	std::vector<KDNode> m_nodes;
	aabb_type m_aabb;
};

MTS_NAMESPACE_END

#endif /* __KDTREE_H */

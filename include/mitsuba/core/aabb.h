/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#if !defined(__AABB_H)
#define __AABB_H

#include <mitsuba/core/bsphere.h>

MTS_NAMESPACE_BEGIN


/**
 * \brief Generic multi-dimensional bounding box data structure
 *
 * Maintains a component-wise minimum and maximum position and provides
 * various convenience functions to query or change them.
 *
 * \tparam T Underlying point data type (e.g. \c TPoint3<float>)
 * \ingroup libcore
 */
template <typename T> struct TAABB {
	typedef T                       point_type;
	typedef typename T::value_type  value_type;
	typedef typename T::vector_type vector_type;

	/** 
	 * \brief Create a new invalid bounding box
	 * 
	 * Initializes the components of the minimum 
	 * and maximum position to \f$\infty\f$ and \f$-\infty\f$,
	 * respectively.
	 */
	inline TAABB() {
		reset();
	}

	/// Unserialize a bounding box from a binary data stream
	inline TAABB(Stream *stream) {
		min = point_type(stream);
		max = point_type(stream);
	}

	/// Create a collapsed AABB from a single point
	inline TAABB(const point_type &p) 
		: min(p), max(p) { }

	/// Create a bounding box from two positions
	inline TAABB(const point_type &min, const point_type &max)
		: min(min), max(max) {
#if defined(MTS_DEBUG)
		for (int i=0; i<point_type::dim; ++i) 
			SAssert(min[i] <= max[i]);
#endif
	}

	/// Equality test
	inline bool operator==(const TAABB &aabb) const {
		return min == aabb.min && max == aabb.max;
	}

	/// Inequality test
	inline bool operator!=(const TAABB &aabb) const {
		return min != aabb.min || max != aabb.max;
	}

	/// Clip to another bounding box
	inline void clip(const TAABB &aabb) {
		for (int i=0; i<point_type::dim; ++i) {
			min[i] = std::max(min[i], aabb.min[i]);
			max[i] = std::min(max[i], aabb.max[i]);
		}
	}

	/** 
	 * \brief Mark the bounding box as invalid.
	 * 
	 * This operation sets the components of the minimum 
	 * and maximum position to \f$\infty\f$ and \f$-\infty\f$,
	 * respectively.
	 */
	inline void reset() {
		min = point_type( std::numeric_limits<value_type>::infinity());
		max = point_type(-std::numeric_limits<value_type>::infinity());
	}

	/// Calculate the n-dimensional volume of the bounding box
	inline value_type getVolume() const {
		vector_type diff = max-min;
		value_type result = diff[0];
		for (int i=1; i<point_type::dim; ++i)
			result *= diff[i];
		return result;
	}
	
	/// Calculate the n-1 dimensional volume of the boundary
	inline Float getSurfaceArea() const {
		vector_type d = max - min;
		Float result = 0.0f;
		for (int i=0; i<point_type::dim; ++i) {
			Float term = 1.0f;
			for (int j=0; j<point_type::dim; ++j) {
				if (i == j)
					continue;
				term *= d[j];
			}
			result += term;
		}
		return 2.0f * result;
	}

	/// Return the center point
	inline point_type getCenter() const {
		return (max + min) * (value_type) 0.5;
	}

	/// Check whether a point lies on or inside the bounding box
	inline bool contains(const point_type &vec) const {
		for (int i=0; i<point_type::dim; ++i)
			if (vec[i] < min[i] || vec[i] > max[i])
				return false;
		return true;
	}

	/// Check whether a given bounding box is contained within this one
	inline bool contains(const TAABB &aabb) const {
		if (!isValid())
			return false;
		for (int i=0; i<point_type::dim; ++i)
			if (aabb.min[i] < min[i] || aabb.max[i] > max[i])
				return false;
		return true;
	}

	/// Axis-aligned bounding box overlap test
	inline bool overlaps(const TAABB &aabb) const {
		for (int i=0; i<point_type::dim; ++i) 
			if (max[i] < aabb.min[i] || min[i] > aabb.max[i])
				return false;
		return true;
	}

	/// Expand the bounding box to contain another point
	inline void expandBy(const point_type &p) {
		for (int i=0; i<point_type::dim; ++i) {
			min[i] = std::min(min[i], p[i]);
			max[i] = std::max(max[i], p[i]);
		}
	}

	/// Expand the bounding box to contain another bounding box
	inline void expandBy(const TAABB &aabb) {
		for (int i=0; i<point_type::dim; ++i) {
			min[i] = std::min(min[i], aabb.min[i]);
			max[i] = std::max(max[i], aabb.max[i]);
		}
	}

	/// Calculate the squared point-AABB distance
	inline value_type squaredDistanceTo(const point_type &p) const {
		value_type result = 0;
		for (int i=0; i<point_type::dim; ++i) {
			value_type value = 0;
			if (p[i] < min[i])
				value = min[i] - p[i];
			else if (p[i] > max[i])
				value = p[i] - max[i];
			result += value*value;
		}
		return result;
	}

	/// Calculate the point-AABB distance
	inline value_type distanceTo(const point_type &p) const {
		return std::sqrt(squaredDistanceTo(p));
	}

	/// Calculate the minimum squared AABB-AABB distance
	inline value_type squaredDistanceTo(const TAABB &aabb) const {
		value_type result = 0;

		for (int i=0; i<point_type::dim; ++i) {
			value_type value = 0;
			if (aabb.max[i] < min[i])
				value = min[i] - aabb.max[i];
			else if (aabb.min[i] > max[i])
				value = aabb.min[i] - max[i];
			result += value*value;
		}
		return result;
	}

	/// Calculate the minimum AABB-AABB distance
	inline value_type distanceTo(const TAABB &aabb) const {
		return std::sqrt(squaredDistanceTo(aabb));
	}

	/// Return whether this bounding box is valid
	inline bool isValid() const {
		for (int i=0; i<point_type::dim; ++i) 
			if (max[i] < min[i])
				return false;
		return true;
	}

	/**
	 * \brief Return whether or not this bounding box 
	 * covers anything at all.
	 *
	 * A bounding box which only covers a single point
	 * is considered nonempty.
	 */
	inline bool isEmpty() const {
		for (int i=0; i<point_type::dim; ++i) {
			if (max[i] > min[i])
				return false;
		}
		return true;
	}

	/// Return the axis index with the largest associated side length
	inline int getLargestAxis() const {
		vector_type d = max - min;
		int largest = 0;

		for (int i=1; i<point_type::dim; ++i)
			if (d[i] > d[largest])
				largest = i;
		return largest;
	}

	/// Return the axis index with the shortest associated side length
	inline int getShortestAxis() const {
		vector_type d = max - min;
		int shortest = 0;

		for (int i=1; i<point_type::dim; ++i)
			if (d[i] < d[shortest])
				shortest = i;
		return shortest;
	}

	/**
	 * \brief Calculate the bounding box extents
	 * \return max-min
	 */
	inline vector_type getExtents() const {
		return max - min;
	}

	/// Serialize this bounding box to a binary data stream
	inline void serialize(Stream *stream) const {
		min.serialize(stream);
		max.serialize(stream);
	}

	/// Return a string representation of the bounding box
	std::string toString() const {
		std::ostringstream oss;
		oss << "AABB[";
		if (!isValid()) {
			oss << "invalid";
		} else {
			oss << "min=" << min.toString()
				<< ", max=" << max.toString();
		}
		oss	<< "]";
		return oss.str();
	}

	point_type min; ///< Component-wise minimum 
	point_type max; ///< Component-wise maximum 
};


/**
 * \brief Axis-aligned bounding box data structure in three dimensions
 * 
 * Maintains a component-wise minimum and maximum position and provides
 * various convenience functions to query or change them.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
struct MTS_EXPORT_CORE AABB : public TAABB<Point> {
public:
	/** 
	 * \brief Create a new invalid bounding box
	 * 
	 * Initializes the components of the minimum 
	 * and maximum position to \f$\infty\f$ and \f$-\infty\f$,
	 * respectively.
	 */
	inline AABB() : TAABB<Point>() { }

	/// Unserialize a bounding box from a binary data stream
	inline AABB(Stream *stream) : TAABB<Point>(stream) { }

	/// Create a collapsed AABB from a single point
	inline AABB(const Point &p) : TAABB<Point>(p) { }

	/// Create a bounding box from two positions
	inline AABB(const point_type &min, const point_type &max) 
		: TAABB<Point>(min, max) {
	}

	/// Construct from a TAABB<Point>
	inline AABB(const TAABB<Point> &aabb) 
		: TAABB<Point>(aabb) { }

	/// Calculate the surface area of the bounding box
	inline Float getSurfaceArea() const {
		Vector d = max - min;
		return (Float) 2.0 * (d.x*d.y + d.x*d.z + d.y*d.z);
	}

	/**
	 * \brief Return the position of a bounding box corner
	 * \param corner Requested corner index (0..7)
	 */
	Point getCorner(uint8_t corner) const;

	/**
	 * \brief Bounding sphere-box overlap test
	 *
	 * Implements the technique proposed by Jim Arvo in
	 * "A simple method for box-sphere intersection testing"
	 * (Graphics Gems, 1990)
	 */
	bool overlaps(const BSphere &sphere) const;

	/** \brief Calculate the near and far ray-AABB intersection
	 * points (if they exist).
	 *
	 * \remark In the Python bindings, this function returns the
	 * \c nearT and \c farT values as a tuple (or \c None, when no
	 * intersection was found)
	 */
	FINLINE bool rayIntersect(const Ray &ray, Float &nearT, Float &farT) const {
		nearT = -std::numeric_limits<Float>::infinity();
		farT  = std::numeric_limits<Float>::infinity();

		/* For each pair of AABB planes */
		for (int i=0; i<3; i++) {
			const Float direction = ray.d[i];
			const Float origin = ray.o[i];
			const Float minVal = min[i], maxVal = max[i];

			if (direction == 0) {
				/* The ray is parallel to the planes */
				if (origin < minVal || origin > maxVal)
					return false;
			} else {
				/* Calculate intersection distances */
				Float t1 = (minVal - origin) * ray.dRcp[i];
				Float t2 = (maxVal - origin) * ray.dRcp[i];

				if (t1 > t2) {
					Float tmp = t1;
					t1 = t2;
					t2 = tmp;
				}

				nearT = std::max(nearT, t1);
				farT = std::min(farT, t2);

				if (nearT > farT)
					return false;
			}
		}
		return true;
	}

#ifdef MTS_SSE
	/**
	 * \brief Intersect against a packet of four rays. 
	 * \return \c false if none of the rays intersect.
	 */
	FINLINE bool rayIntersectPacket(const RayPacket4 &ray, RayInterval4 &interval) const;
#endif

	/// Create a bounding sphere, which contains the axis-aligned box
	BSphere getBSphere() const;
};

MTS_NAMESPACE_END

#ifdef MTS_SSE
#include <mitsuba/core/aabb_sse.h>
#endif

#endif /* __AABB_H */

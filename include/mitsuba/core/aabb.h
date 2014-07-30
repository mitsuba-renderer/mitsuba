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
#if !defined(__MITSUBA_CORE_AABB_H_)
#define __MITSUBA_CORE_AABB_H_

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
	typedef T                           PointType;
	typedef typename T::Scalar          Scalar;
	typedef typename T::VectorType      VectorType;
	typedef TRay<PointType, VectorType> RayType;

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
		min = PointType(stream);
		max = PointType(stream);
	}

	/// Create a collapsed AABB from a single point
	inline TAABB(const PointType &p)
		: min(p), max(p) { }

	/// Create a bounding box from two positions
	inline TAABB(const PointType &min, const PointType &max)
		: min(min), max(max) {
#if defined(MTS_DEBUG)
		for (int i=0; i<PointType::dim; ++i)
			SAssert(min[i] <= max[i]);
#endif
	}

	/// Copy constructor
	inline TAABB(const TAABB &aabb)
		: min(aabb.min), max(aabb.max) { }

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
		for (int i=0; i<PointType::dim; ++i) {
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
		min = PointType( std::numeric_limits<Scalar>::infinity());
		max = PointType(-std::numeric_limits<Scalar>::infinity());
	}

	/// Calculate the n-dimensional volume of the bounding box
	inline Scalar getVolume() const {
		VectorType diff = max-min;
		Scalar result = diff[0];
		for (int i=1; i<PointType::dim; ++i)
			result *= diff[i];
		return result;
	}

	/// Calculate the n-1 dimensional volume of the boundary
	inline Float getSurfaceArea() const {
		VectorType d = max - min;
		Float result = 0.0f;
		for (int i=0; i<PointType::dim; ++i) {
			Float term = 1.0f;
			for (int j=0; j<PointType::dim; ++j) {
				if (i == j)
					continue;
				term *= d[j];
			}
			result += term;
		}
		return 2.0f * result;
	}

	/// Return the center point
	inline PointType getCenter() const {
		return (max + min) * (Scalar) 0.5;
	}

	/// Return the position of one of the corners (in <tt>0..2^dim-1</tt>)
	inline PointType getCorner(int index) const {
		PointType result;
		for (int d=0; d<PointType::dim; ++d) {
			if (index & (1 << d))
				result[d] = max[d];
			else
				result[d] = min[d];
		}
		return result;
	}

	/// Return a child bounding box in a interval-, quad-, octree, etc.
	inline TAABB getChild(int index) const {
		TAABB result(getCenter());

		for (int d=0; d<PointType::dim; ++d) {
			if (index & (1 << d))
				result.max[d] = max[d];
			else
				result.min[d] = min[d];
		}

		return result;
	}

	/// Check whether a point lies on or inside the bounding box
	inline bool contains(const PointType &p) const {
		for (int i=0; i<PointType::dim; ++i)
			if (p[i] < min[i] || p[i] > max[i])
				return false;
		return true;
	}

	/// Check whether a given bounding box is contained within this one
	inline bool contains(const TAABB &aabb) const {
		if (!isValid())
			return false;
		for (int i=0; i<PointType::dim; ++i)
			if (aabb.min[i] < min[i] || aabb.max[i] > max[i])
				return false;
		return true;
	}

	/// Axis-aligned bounding box overlap test
	inline bool overlaps(const TAABB &aabb) const {
		for (int i=0; i<PointType::dim; ++i)
			if (max[i] < aabb.min[i] || min[i] > aabb.max[i])
				return false;
		return true;
	}

	/// Expand the bounding box to contain another point
	inline void expandBy(const PointType &p) {
		for (int i=0; i<PointType::dim; ++i) {
			min[i] = std::min(min[i], p[i]);
			max[i] = std::max(max[i], p[i]);
		}
	}

	/// Expand the bounding box to contain another bounding box
	inline void expandBy(const TAABB &aabb) {
		for (int i=0; i<PointType::dim; ++i) {
			min[i] = std::min(min[i], aabb.min[i]);
			max[i] = std::max(max[i], aabb.max[i]);
		}
	}

	/// Calculate the squared point-AABB distance
	inline Scalar squaredDistanceTo(const PointType &p) const {
		Scalar result = 0;
		for (int i=0; i<PointType::dim; ++i) {
			Scalar value = 0;
			if (p[i] < min[i])
				value = min[i] - p[i];
			else if (p[i] > max[i])
				value = p[i] - max[i];
			result += value*value;
		}
		return result;
	}

	/// Calculate the point-AABB distance
	inline Scalar distanceTo(const PointType &p) const {
		return std::sqrt(squaredDistanceTo(p));
	}

	/// Calculate the minimum squared AABB-AABB distance
	inline Scalar squaredDistanceTo(const TAABB &aabb) const {
		Scalar result = 0;

		for (int i=0; i<PointType::dim; ++i) {
			Scalar value = 0;
			if (aabb.max[i] < min[i])
				value = min[i] - aabb.max[i];
			else if (aabb.min[i] > max[i])
				value = aabb.min[i] - max[i];
			result += value*value;
		}
		return result;
	}

	/// Calculate the minimum AABB-AABB distance
	inline Scalar distanceTo(const TAABB &aabb) const {
		return std::sqrt(squaredDistanceTo(aabb));
	}

	/// Return whether this bounding box is valid
	inline bool isValid() const {
		for (int i=0; i<PointType::dim; ++i)
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
		for (int i=0; i<PointType::dim; ++i) {
			if (max[i] > min[i])
				return false;
		}
		return true;
	}

	/// Return the axis index with the largest associated side length
	inline int getLargestAxis() const {
		VectorType d = max - min;
		int largest = 0;

		for (int i=1; i<PointType::dim; ++i)
			if (d[i] > d[largest])
				largest = i;
		return largest;
	}

	/// Return the axis index with the shortest associated side length
	inline int getShortestAxis() const {
		VectorType d = max - min;
		int shortest = 0;

		for (int i=1; i<PointType::dim; ++i)
			if (d[i] < d[shortest])
				shortest = i;
		return shortest;
	}

	/**
	 * \brief Calculate the bounding box extents
	 * \return max-min
	 */
	inline VectorType getExtents() const {
		return max - min;
	}

	/** \brief Calculate the near and far ray-AABB intersection
	 * points (if they exist).
	 *
	 * The parameters \c nearT and \c farT are used to return the
	 * ray distances to the intersections (including negative distances).
	 * Any previously contained value is overwritten, even if there was
	 * no intersection.
	 *
	 * \remark In the Python bindings, this function returns the
	 * \c nearT and \c farT values as a tuple (or \c None, when no
	 * intersection was found)
	 */
	FINLINE bool rayIntersect(const RayType &ray, Float &nearT, Float &farT) const {
		nearT = -std::numeric_limits<Float>::infinity();
		farT  = std::numeric_limits<Float>::infinity();

		/* For each pair of AABB planes */
		for (int i=0; i<PointType::dim; i++) {
			const Float origin = ray.o[i];
			const Float minVal = min[i], maxVal = max[i];

			if (ray.d[i] == 0) {
				/* The ray is parallel to the planes */
				if (origin < minVal || origin > maxVal)
					return false;
			} else {
				/* Calculate intersection distances */
				Float t1 = (minVal - origin) * ray.dRcp[i];
				Float t2 = (maxVal - origin) * ray.dRcp[i];

				if (t1 > t2)
					std::swap(t1, t2);

				nearT = std::max(t1, nearT);
				farT = std::min(t2, farT);

				if (!(nearT <= farT))
					return false;
			}
		}

		return true;
	}

	/** \brief Calculate the overlap between an axis-aligned bounding box
	 * and a ray segment
	 *
	 * This function is an extended version of the simpler \ref rayIntersect command
	 * provided above. The first change is that input values passed via
	 * the \c nearT and \c farT parameters are considered to specify a query interval.
	 *
	 * This interval is intersected against the bounding box, returning the remaining
	 * interval using the \c nearT and \c farT parameters. Furthermore, the
	 * interval endpoints are also returned as 3D positions via the \c near and
	 * \c far parameters. Special care is taken to reduce round-off errors.
	 *
	 * \remark In the Python bindings, this function has the signature
	 * <tt>(nearT, farT, near, far) = rayIntersect(ray, nearT, farT)</tt>.
	 * It returns \c None when no intersection was found.
	 */
	FINLINE bool rayIntersect(const RayType &ray, Float &nearT, Float &farT, PointType &near, PointType &far) const {
		int nearAxis = -1, farAxis = -1;

		/* For each pair of AABB planes */
		for (int i=0; i<PointType::dim; i++) {
			const Float origin = ray.o[i];
			const Float minVal = min[i], maxVal = max[i];

			if (ray.d[i] == 0) {
				/* The ray is parallel to the planes */
				if (origin < minVal || origin > maxVal)
					return false;
			} else {
				/* Calculate intersection distances */
				Float t1 = (minVal - origin) * ray.dRcp[i];
				Float t2 = (maxVal - origin) * ray.dRcp[i];

				bool flip = t1 > t2;
				if (flip)
					std::swap(t1, t2);

				if (t1 > nearT) {
					nearT = t1;
					nearAxis = flip ? (i + PointType::dim) : i;
				}

				if (t2 < farT) {
					farT = t2;
					farAxis = flip ? i : (i + PointType::dim);
				}
			}
		}

		if (!(nearT <= farT))
			return false;

		near = ray(nearT); far = ray(farT);

		/* Avoid roundoff errors on the component where the intersection took place */
		if (nearAxis >= 0)
			near[nearAxis % PointType::dim] = ((Float *) this)[nearAxis];

		if (farAxis >= 0)
			far[farAxis % PointType::dim]   = ((Float *) this)[farAxis];

		return true;
	}

	/// Serialize this bounding box to a binary data stream
	inline void serialize(Stream *stream) const {
		min.serialize(stream);
		max.serialize(stream);
	}

	/// Return a string representation of the bounding box
	std::string toString() const {
		std::ostringstream oss;
		oss << "AABB" << PointType::dim << "[";
		if (!isValid()) {
			oss << "invalid";
		} else {
			oss << "min=" << min.toString()
				<< ", max=" << max.toString();
		}
		oss	<< "]";
		return oss.str();
	}

	PointType min; ///< Component-wise minimum
	PointType max; ///< Component-wise maximum
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
struct AABB : public TAABB<Point> {
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
	inline AABB(const PointType &min, const PointType &max)
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
	MTS_EXPORT_CORE Point getCorner(uint8_t corner) const;

	/**
	 * \brief Bounding sphere-box overlap test
	 *
	 * Implements the technique proposed by Jim Arvo in
	 * "A simple method for box-sphere intersection testing"
	 * (Graphics Gems, 1990)
	 */
	MTS_EXPORT_CORE bool overlaps(const BSphere &sphere) const;

#ifdef MTS_SSE
	/**
	 * \brief Intersect against a packet of four rays.
	 * \return \c false if none of the rays intersect.
	 */
	FINLINE bool rayIntersectPacket(const RayPacket4 &ray, RayInterval4 &interval) const;
#endif

	/// Create a bounding sphere, which contains the axis-aligned box
	MTS_EXPORT_CORE BSphere getBSphere() const;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_AABB_H_ */

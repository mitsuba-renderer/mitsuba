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

#if !defined(__AABB_H)
#define __AABB_H

#include <mitsuba/core/bsphere.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Axis-aligned bounding box data structure in three dimensions
 * 
 * Maintains a component-wise minimum and maximum position and provides
 * various convenience functions to query or change them.
 */
struct MTS_EXPORT_CORE AABB {
public:
	Point min; ///< Component-wise minimum 
	Point max; ///< Component-wise maximum 

	/** 
	 * \brief Construct an invalid bounding box
	 * 
	 * The minimum and maximum positions will be
	 * initialized to \f$(\infty,\infty,\infty)\f$ 
	 * and \f$(-\infty, -\infty, -\infty)\f$, respectively.
	 */
	inline AABB() {
		reset();
	}

	/// Unserialize a bounding box from a binary data stream
	inline AABB(Stream *stream) {
		min = Point(stream);
		max = Point(stream);
	}
	
	/// Create a collapsed AABB from a single point
	inline AABB(const Point &p) 
		: min(p), max(p) { }

	/// Create a bounding box from two 3D positions
	inline AABB(const Point &min, const Point &max)
		: min(min), max(max) {
		SAssert(min.x <= max.x
			&& min.y <= max.y
			&& min.z <= max.z);
	}

	/// Copy-constructor
	inline AABB(const AABB &aabb) 
		: min(aabb.min), max(aabb.max) {
	}

	/// Assignment operator
	inline AABB &operator=(const AABB &aabb) {
		min = aabb.min;
		max = aabb.max;
		return *this;
	}

	/// Equality test
	inline bool operator==(const AABB &aabb) const {
		return min == aabb.min && max == aabb.max;
	}

	/// Inequality test
	inline bool operator!=(const AABB &aabb) const {
		return min != aabb.min || max != aabb.max;
	}

	/** 
	 * \brief Mark the bounding box as invalid.
	 * 
	 * This operation sets the
	 * minimum position to \f$(\infty,\infty,\infty)\f$ and the
	 * maximum position to \f$(-\infty, -\infty, -\infty)\f$.
	 */
	inline void reset() {
		const Float inf = std::numeric_limits<Float>::infinity();
		min = Point(inf, inf, inf);
		max = Point(-inf, -inf, -inf);
	}

	/// Clip to another bounding box
	inline void clip(const AABB &aabb) {
		min.x = std::max(min.x, aabb.min.x);
		min.y = std::max(min.y, aabb.min.y);
		min.z = std::max(min.z, aabb.min.z);
		max.x = std::min(max.x, aabb.max.x);
		max.y = std::min(max.y, aabb.max.y);
		max.z = std::min(max.z, aabb.max.z);
	}

	/// Return the center point
	inline Point getCenter() const {
		return (max + min) * (Float) 0.5;
	}

	/// Calculate the volume of the bounding box
	inline Float getVolume() const {
		Float x = max.x - min.x;
		Float y = max.y - min.y;
		Float z = max.z - min.z;

		return x*x + y*y + z*z;
	}

	/// Calculate the surface area of the bounding box
	inline Float getSurfaceArea() const {
		Vector d = max - min;
		return (Float) 2.0 * (d.x*d.y + d.x*d.z + d.y*d.z);
	}

	/**
	 * \brief Calculate the bounding box extents
	 * \return max-min
	 */
	inline Vector getExtents() const {
		return max - min;
	}

	/// Return the axis index with the largest associated side length
	inline int getLargestAxis() const {
		Vector d = max - min;
		if (d.x >= d.y && d.x >= d.z)
			return 0;
		else if (d.y >= d.z)
			return 1;
		else
			return 2;
	}

	/// Return the axis index with the smallest associated side length
	inline int getSmallestAxis() const {
		Vector d = max - min;
		if (d.x <= d.y && d.x <= d.z)
			return 0;
		else if (d.y <= d.z)
			return 1;
		else
			return 2;
	}

	/// Return whether this bounding box is valid
	inline bool isValid() const {
		return max.x >= min.x && max.y >= min.y && max.z >= min.z;
	}

	/**
	 * \brief Return whether or not this bounding box covers a 
	 * nonzero amount of space
	 */
	inline bool isEmpty() const {
		return max.x <= min.x
			&& max.y <= min.y
			&& max.z <= min.z;
	}


	/// Return the component-wise minimum point of the bounding box
	inline const Point &getMinimum() const {
		return min;
	}

	/// Return the component-wise maximum point of the bounding box
	inline const Point &getMaximum() const {
		return max;
	}

	/**
	 * \brief Return the position of a bounding box corner
	 * \param corner Requested corner index (0..7)
	 */
	Point getCorner(uint8_t corner) const;

	/// Check whether a point lies on or inside the bounding box
	bool contains(const Point &vec) const;

	/// Check whether a given bounding box is contained within this one
	bool contains(const AABB &aabb) const;

	/**
	 * \brief Bounding sphere-box overlap test
	 *
	 * Implements the technique proposed by Jim Arvo in
	 * "A simple method for box-sphere intersection testing"
	 * (Graphics Gems, 1990)
	 */
	bool overlaps(const BSphere &sphere) const;

	/// Axis-aligned bounding box overlap test
	bool overlaps(const AABB &saabb) const;

	/// Expand the bounding box to contain another point
	void expandBy(const Point &vec);

	/// Expand the bounding box to contain another bounding box
	void expandBy(const AABB &aabb);

	/// Calculate the point-AABB distance
	Float distanceTo(const Point &p) const;

	/** \brief Calculate the near and far ray-AABB intersection
	 * points (if they exist).
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
	 * \return \a false if none of the rays intersect.
	 */
	FINLINE bool rayIntersectPacket(const RayPacket4 &ray, RayInterval4 &interval) const;
#endif

	/// Create a bounding sphere, which contains the axis-aligned box
	BSphere getBSphere() const;

	/// Serialize this bounding box to a binary data stream
	inline void serialize(Stream *stream) const {
		min.serialize(stream);
		max.serialize(stream);
	}

	/// Return a string representation of the bounding box
	std::string toString() const;
};

MTS_NAMESPACE_END

#ifdef MTS_SSE
#include <mitsuba/core/aabb_sse.h>
#endif

#endif /* __AABB_H */

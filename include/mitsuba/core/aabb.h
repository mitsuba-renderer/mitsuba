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

/** \brief Axis-aligned bounding box data structure
 */
struct MTS_EXPORT_CORE AABB {
public:
	Point min;
	Point max;

	/// Construct an invalid bounding box
	inline AABB() {
		reset();
	}

	/// Unserialize an AABB from a stream
	inline AABB(Stream *stream) {
		min = Point(stream);
		max = Point(stream);
	}

	/// Create a bounding box from two 3-dimension vectors
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

	/// Equal operator
	inline bool operator==(const AABB &aabb) const {
		return min == aabb.min && max == aabb.max;
	}

	/// Not equal operator
	inline bool operator!=(const AABB &aabb) const {
		return min != aabb.min || max != aabb.max;
	}

	/// Mark the bounding box as invalid
	inline void reset() {
		const Float inf = std::numeric_limits<Float>::infinity();
		min = Point(inf, inf, inf);
		max = Point(-inf, -inf, -inf);
	}

	/// Calculate the volume of the bounding box
	inline Float getVolume() const {
		Float x = max.x - min.x;
		Float y = max.y - min.y;
		Float z = max.z - min.z;

		return x*x + y*y + z*z;
	}

	/// Clip to another AABB
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
		return max*.5f + min*.5f;
	}

	/// Calculate the surface area of the bounding box
	inline Float getSurfaceArea() const {
		Vector d = max - min;
		return 2.0f * (d.x*d.y + d.x*d.z + d.y*d.z);
	}
	
	/// Calculate the AABB extents
	inline Vector getExtents() const {
		return max - min;
	}

	/// Return the axis with the largest corresponding AABB side
	inline int getLargestAxis() const {
		Vector d = max - min;
		if (d.x >= d.y && d.x >= d.z)
			return 0;
		else if (d.y >= d.z)
			return 1;
		else
			return 2;
	}

	/// Return the axis with the smallest corresponding AABB side
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
	 * Return whether this bounding box covers a non-zero
	 * amount of space
	 */
	inline bool isEmpty() const {
		return max.x <= min.x
			&& max.y <= min.y
			&& max.z <= min.z;
	}


	/// Return the minimum vector of the bounding box
	inline const Point &getMinimum() const {
		return min;
	}

	/// Return the maximum vector of the bounding box
	inline const Point &getMaximum() const {
		return max;
	}

	/// Return the middle point
	Point getMidPoint() const;

	/** \brief Return the vector coordinates of a bounding
	 * box corner
	 * @param corner Corner index (0..7)
	 */
	Point getCorner(uint8_t corner) const;

	/// Checks whether a vector is inside the bounding box
	bool contains(const Point &vec) const;

	/// Bounding sphere-AABB overlap test
	inline bool overlaps(const BSphere &sphere) const {
		Float distance = 0;
		for (int i=0; i<3; ++i) {
			if (sphere.center[i] < min[i]) {
				Float d = sphere.center[i]-min[i];
				distance += d*d;
			} else if (sphere.center[i] > max[i]) {
				Float d = sphere.center[i]-max[i];
				distance += d*d;
			}
		}
		return distance < sphere.radius*sphere.radius;
	}

	/// Expands the bounding box to contain another vector
	void expandBy(const Point &vec);

	/// Expands the bounding box to contain another bounding box
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
	 * Intersect against a packet of four rays. Returns false if none of 
	 * the rays intersect.
	 */
	FINLINE bool rayIntersectPacket(const RayPacket4 &ray, RayInterval4 &interval) const;
#endif

	/// Serialize this AABB to a stream
	inline void serialize(Stream *stream) const {
		min.serialize(stream);
		max.serialize(stream);
	}

	/// Returns a string representation of the bounding box
	std::string toString() const;
};

MTS_NAMESPACE_END

#ifdef MTS_SSE
#include <mitsuba/core/aabb_sse.h>
#endif

#endif /* __AABB_H */

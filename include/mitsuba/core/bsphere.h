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
#if !defined(__MITSUBA_CORE_BSPHERE_H_)
#define __MITSUBA_CORE_BSPHERE_H_

#include <mitsuba/core/ray.h>

MTS_NAMESPACE_BEGIN

/** \brief Bounding sphere data structure in three dimensions
 *
 * \ingroup libcore
 * \ingroup libpython
 */
struct BSphere {
	Point center;
	Float radius;

	/// Construct a bounding sphere at the origin having radius zero
	inline BSphere() : center(0.0f), radius(0.0f) { }

	/// Unserialize a bounding sphere from a binary data stream
	inline BSphere(Stream *stream) {
		center = Point(stream);
		radius = stream->readFloat();
	}

	/// Create a bounding sphere from a given center point and radius
	inline BSphere(const Point &center, Float radius)
		: center(center), radius(radius) {
	}

	/// Copy constructor
	inline BSphere(const BSphere &boundingSphere)
		: center(boundingSphere.center), radius(boundingSphere.radius) {
	}

	/// Return whether this bounding sphere has a radius of zero or less.
	inline bool isEmpty() const {
		return radius <= 0.0f;
	}

	/// Expand the bounding sphere radius to contain another point.
	inline void expandBy(const Point p) {
		radius = std::max(radius, (p-center).length());
	}

	/// Check whether the specified point is inside or on the sphere
	inline bool contains(const Point p) const {
		return (p - center).length() <= radius;
	}

	/// Equality test
	inline bool operator==(const BSphere &boundingSphere) const {
		return center == boundingSphere.center && radius == boundingSphere.radius;
	}

	/// Inequality test
	inline bool operator!=(const BSphere &boundingSphere) const {
		return center != boundingSphere.center || radius != boundingSphere.radius;
	}

	/**
	 * \brief Calculate the intersection points with the given ray
	 * \return \c true if the ray intersects the bounding sphere
	 *
	 * \remark In the Python bindings, this function returns the
	 * \c nearT and \c farT values as a tuple (or \c None, when no
	 * intersection was found)
	 */
	inline bool rayIntersect(const Ray &ray, Float &nearHit, Float &farHit) const {
		Vector o = ray.o - center;
		Float A = ray.d.lengthSquared();
		Float B = 2 * dot(o, ray.d);
		Float C = o.lengthSquared() - radius*radius;

		return solveQuadratic(A, B, C, nearHit, farHit);
	}

	/// Serialize this bounding sphere to a binary data stream
	inline void serialize(Stream *stream) const {
		center.serialize(stream);
		stream->writeFloat(radius);
	}

	/// Return a string representation of the bounding sphere
	inline std::string toString() const {
		std::ostringstream oss;
		oss << "BSphere[center = " << center.toString()
			<< ", radius = " << radius << "]";
		return oss.str();
	}
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_BSPHERE_H_ */

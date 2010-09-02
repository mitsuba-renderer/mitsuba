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

#if !defined(__NORMAL_H)
#define __NORMAL_H

#include <mitsuba/core/vector.h>

MTS_NAMESPACE_BEGIN

/** \brief Simple three-dimensional normal class using floating point values.
 * This class is different from the Vector class in how it is treated
 * by linear transformations.
 */
class Normal {
public:
	Float x, y, z;

	inline Normal(Float _x = 0.0f, Float _y = 0.0f, Float _z = 0.0f)
		: x(_x), y(_y), z(_z) {
	}

	explicit inline Normal(const Vector &v)
		: x(v.x), y(v.y), z(v.z) {
	}

	inline Normal(Stream *stream) {
		x = stream->readFloat();
		y = stream->readFloat();
		z = stream->readFloat();
	}

	inline Normal operator+(const Normal &v) const {
		return Normal(x + v.x, y + v.y, z + v.z);
	}

	inline Normal operator-(const Normal &v) const {
		return Normal(x - v.x, y - v.y, z - v.z);
	}

	inline Normal& operator+=(const Normal &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}

	inline Normal& operator-=(const Normal &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}

	inline Normal operator*(Float f) const {
		return Normal(x*f, y*f, z*f);
	}

	inline Normal &operator*=(Float f) {
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}

	inline Normal operator-() const {
		return Normal(-x, -y, -z);
	}

	inline Normal operator/(Float f) const {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Normal: Division by zero!");
#endif
		Float r = 1.0f / f;
		return Normal(x * r, y * r, z * r);
	}

	inline Normal &operator/=(Float f) {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Normal: Division by zero!");
#endif
		Float r = 1.0f / f;
		x *= r;
		y *= r;
		z *= r;
		return *this;
	}

	inline Float operator[](int i) const {
		return (&x)[i];
	}

	inline Float &operator[](int i) {
		return (&x)[i];
	}

	inline bool isZero() const {
		return x==0 && y == 0 && z == 0;
	}

	inline Float lengthSquared() const {
		return x*x + y*y + z*z;
	}

	inline Float length() const {
		return std::sqrt(lengthSquared());
	}

	inline bool operator==(const Normal &v) const {
		return (v.x == x && v.y == y && v.z == z);
	}
	
	inline bool operator!=(const Normal &v) const {
		return !operator==(v);
	}

	inline void serialize(Stream *stream) const {
		stream->writeFloat(x);
		stream->writeFloat(y);
		stream->writeFloat(z);
	}

	inline std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << ", " << y << ", " << z << "]";
		return oss.str();
	}
};

inline Float dot(const Normal &v1, const Normal &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline Float dot(const Vector &v1, const Normal &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline Float dot(const Normal &v1, const Vector &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline Float absDot(const Normal &v1, const Normal &v2) {
	return std::abs(dot(v1, v2));
}

inline Float absDot(const Vector &v1, const Normal &v2) {
	return std::abs(dot(v1, v2));
}

inline Float absDot(const Normal &v1, const Vector &v2) {
	return std::abs(dot(v1, v2));
}

inline Normal normalize(const Normal &n) {
	Float length = n.length();
#ifdef MTS_DEBUG
	if (length == 0.0f) {
		SLog(EWarn, "Zero-length normal encountered!");
		return n;
	} else {
		return n / length;
	}
#else
	return n / length;
#endif
}

inline Vector::Vector(const Normal &n)
 : x(n.x), y(n.y), z(n.z) {
}

MTS_NAMESPACE_END

#endif /* __NORMAL_H */

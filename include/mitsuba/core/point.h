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

#if !defined(__POINT_H)
#define __POINT_H

#include <mitsuba/core/vector.h>

MTS_NAMESPACE_BEGIN

/** \brief Simple three-dimensional point class using floating point values.
 * This class is different from the Vector class in how it is treated
 * by linear transformations.
 */
class Point {
public:
	Float x, y, z;

	inline Point(Float _x = 0.0f, Float _y = 0.0f, Float _z = 0.0f)
		: x(_x), y(_y), z(_z) {
	}

	inline Point(Stream *stream) {
		x = stream->readFloat();
		y = stream->readFloat();
		z = stream->readFloat();
	}

	explicit inline Point(const Vector &v) 
		: x(v.x), y(v.y), z(v.z) {
	}

	inline Point operator+(const Vector &v) const {
		return Point(x + v.x, y + v.y, z + v.z);
	}

	inline Point operator-(const Vector &v) const {
		return Point(x - v.x, y - v.y, z - v.z);
	}

	inline Vector operator-(const Point &p) const {
		return Vector(x - p.x, y - p.y, z - p.z);
	}
	
	inline Point operator+(const Point &p) const {
			return Point(x + p.x, y + p.y, z + p.z);
		}
	
	inline Point& operator+=(const Point &p) {
		x += p.x; y += p.y; z += p.z;
		return *this;
	}

	inline Point& operator+=(const Vector &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}

	inline Point& operator-=(const Vector &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}

	inline Point operator*(Float f) const {
		return Point(x*f, y*f, z*f);
	}

	inline Point &operator*=(Float f) {
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}

	inline Point operator/(Float f) const {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Point: Division by zero!");
#endif
		Float r = 1.0f / f;
		return Point(x * r, y * r, z * r);
	}

	inline Point &operator/=(Float f) {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Point: Division by zero!");
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

	inline bool operator==(const Point &v) const {
		return (v.x == x && v.y == y && v.z == z);
	}
	
	inline bool operator!=(const Point &v) const {
		return !operator==(v);
	}
	
	inline bool isZero() const {
		return x==0 && y == 0 && z == 0;
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

class Point2 {
public:
	Float x, y;

	inline Point2(Float _x = 0.0f, Float _y = 0.0f)
		: x(_x), y(_y) {
	}

	inline Point2(Stream *stream) {
		x = stream->readFloat();
		y = stream->readFloat();
	}

	inline Point2 operator*(Float f) const {
		return Point2(x*f, y*f);
	}

	inline Point2 &operator*=(Float f) {
		x *= f;
		y *= f;
		return *this;
	}
	
	inline Point2 operator+(const Point2 &p) const {
		return Point2(x + p.x, y + p.y);
	}
	
	inline Point2 operator+(const Vector2 &p) const {
		return Point2(x + p.x, y + p.y);
	}

	inline Point2& operator+=(const Point2 &p) {
		x += p.x; y += p.y;
		return *this;
	}
	
	inline Point2& operator+=(const Vector2 &p) {
		x += p.x; y += p.y;
		return *this;
	}

	inline Vector2 operator-(const Point2 &p) const {
		return Vector2(x - p.x, y - p.y);
	}
	
	inline Point2 operator-(const Vector2 &p) const {
		return Point2(x - p.x, y - p.y);
	}

	inline Point2& operator-=(const Point2 &p) {
		x -= p.x; y -= p.y;
		return *this;
	}
	
	inline Point2& operator-=(const Vector2 &p) {
		x -= p.x; y -= p.y;
		return *this;
	}

	inline Point2 operator/(Float f) const {
		Float r = 1.0f / f;
		return Point2(x * r, y * r);
	}

	inline Point2 &operator/=(Float f) {
		Float r = 1.0f / f;
		x *= r;
		y *= r;
		return *this;
	}

	inline Float operator[](int i) const {
		return (&x)[i];
	}

	inline Float &operator[](int i) {
		return (&x)[i];
	}

	inline bool operator==(const Point2 &v) const {
		return (v.x == x && v.y == y);
	}
	
	inline bool operator!=(const Point2 &v) const {
		return !operator==(v);
	}
	
	inline bool isZero() const {
		return x==0 && y == 0;
	}

	inline void serialize(Stream *stream) const {
		stream->writeFloat(x);
		stream->writeFloat(y);
	}

	inline std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << ", " << y << "]";
		return oss.str();
	}
};

class Point3i {
public:
	int x, y, z;

	inline Point3i(int _x = 0, int _y = 0, int _z = 0)
		: x(_x), y(_y), z(_z) {
	}

	inline Point3i(Stream *stream) {
		x = stream->readInt();
		y = stream->readInt();
		z = stream->readInt();
	}

	explicit inline Point3i(const Vector3i &v) 
		: x(v.x), y(v.y), z(v.z) {
	}

	inline Point3i operator+(const Vector3i &v) const {
		return Point3i(x + v.x, y + v.y, z + v.z);
	}

	inline Point3i operator-(const Vector3i &v) const {
		return Point3i(x - v.x, y - v.y, z - v.z);
	}

	inline Vector3i operator-(const Point3i &p) const {
		return Vector3i(x - p.x, y - p.y, z - p.z);
	}
	
	inline Point3i operator+(const Point3i &p) const {
			return Point3i(x + p.x, y + p.y, z + p.z);
		}
	
	inline Point3i& operator+=(const Point3i &p) {
		x += p.x; y += p.y; z += p.z;
		return *this;
	}

	inline Point3i& operator+=(const Vector3i &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}

	inline Point3i& operator-=(const Vector3i &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}

	inline Point3i operator*(int f) const {
		return Point3i(x*f, y*f, z*f);
	}

	inline Point3i &operator*=(int f) {
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}

	inline Point3i operator/(int i) const {
#ifdef MTS_DEBUG
		if (i == 0)
			SLog(EWarn, "Point3i: Division by zero!");
#endif
		return Point3i(x / i, y / i, z / i);
	}

	inline Point3i &operator/=(int i) {
#ifdef MTS_DEBUG
		if (i == 0)
			SLog(EWarn, "Point3i: Division by zero!");
#endif
		x /= i;
		y /= i;
		z /= i;
		return *this;
	}

	inline int operator[](int i) const {
		return (&x)[i];
	}

	inline int &operator[](int i) {
		return (&x)[i];
	}

	inline bool operator==(const Point3i &v) const {
		return (v.x == x && v.y == y && v.z == z);
	}
	
	inline bool operator!=(const Point3i &v) const {
		return !operator==(v);
	}
	
	inline bool isZero() const {
		return x==0 && y == 0 && z == 0;
	}

	inline void serialize(Stream *stream) const {
		stream->writeInt(x);
		stream->writeInt(y);
		stream->writeInt(z);
	}

	inline std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << ", " << y << ", " << z << "]";
		return oss.str();
	}
};

class Point2i {
public:
	int x, y;

	inline Point2i(int _x = 0, int _y = 0)
		: x(_x), y(_y) {
	}
	
	explicit inline Point2i(const Vector2i &v) 
		: x(v.x), y(v.y) {
	}
	
	inline Point2i(Stream *stream) {
		x = stream->readInt();
		y = stream->readInt();
	}

	inline Point2i operator*(int i) const {
		return Point2i(x*i, y*i);
	}

	inline Point2i &operator*=(int i) {
		x *= i;
		y *= i;
		return *this;
	}
	
	inline Point2i operator+(const Point2i &p) const {
		return Point2i(x + p.x, y + p.y);
	}
	
	inline Point2i operator+(const Vector2i &v) const {
		return Point2i(x + v.x, y + v.y);
	}

	inline Point2i& operator+=(const Point2i &p) {
		x += p.x; y += p.y;
		return *this;
	}
	
	inline Point2i& operator+=(const Vector2i &v) {
		x += v.x; y += v.y;
		return *this;
	}

	inline Point2i operator-(const Point2i &p) const {
		return Point2i(x - p.x, y - p.y);
	}
	
	inline Point2i operator-(const Vector2i &v) const {
		return Point2i(x - v.x, y - v.y);
	}

	inline Point2i& operator-=(const Point2i &p) {
		x -= p.x; y -= p.y;
		return *this;
	}

	inline Point2i& operator-=(const Vector2i &v) {
		x -= v.x; y -= v.y;
		return *this;
	}

	inline Point2i operator/(int i) const {
		return Point2i(x / i, y / i);
	}

	inline Point2i &operator/=(int i) {
		x /= i;
		y /= i;
		return *this;
	}

	int operator[](int i) const {
		return (&x)[i];
	}

	int &operator[](int i) {
		return (&x)[i];
	}

	inline bool operator==(const Point2i &v) const {
		return (v.x == x && v.y == y);
	}
	
	inline bool operator!=(const Point2i &v) const {
		return !operator==(v);
	}

	inline bool isZero() const {
		return x==0 && y == 0;
	}

	inline void serialize(Stream *stream) const {
		stream->writeInt(x);
		stream->writeInt(y);
	}

	inline std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << ", " << y << "]";
		return oss.str();
	}
};

inline Vector::Vector(const Point3i &p)
 : x((Float) p.x), y((Float) p.y), z((Float) p.z) {
}

inline Vector3i::Vector3i(const Point3i &p)
 : x(p.x), y(p.y), z(p.z) {
}

inline Vector::Vector(const Point &p)
 : x(p.x), y(p.y), z(p.z) {
}

inline Vector2::Vector2(const Point2 &p)
 : x(p.x), y(p.y) {
}

inline Vector2i::Vector2i(const Point2i &p)
 : x(p.x), y(p.y) {
}

inline Float distance(const Point &p1, const Point &p2) {
	return (p1 - p2).length();
}

inline Float distanceSquared(const Point &p1, const Point &p2) {
	return (p1 - p2).lengthSquared();
}

MTS_NAMESPACE_END

#endif /* __POINT_H */


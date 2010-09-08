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

#if !defined(__VECTOR_H)
#define __VECTOR_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/stream.h>

MTS_NAMESPACE_BEGIN

class Normal;
class Point;
class Point2;
class Point2i;
class Point3i;
class Vector3i;

/** \brief Simple three-dimensional vector class using floating point values.
 */
class Vector {
public:
	Float x, y, z;

	inline Vector(Float _x = 0.0f, Float _y = 0.0f, Float _z = 0.0f)
		: x(_x), y(_y), z(_z) {
	}
	
	inline Vector(Stream *stream) {
		x = stream->readFloat();
		y = stream->readFloat();
		z = stream->readFloat();
	}
	
	inline Vector(const Normal &n);
	inline Vector(const Point &p);
	inline Vector(const Vector3i &v);
	inline Vector(const Point3i &v);

	inline Vector operator+(const Vector &v) const {
		return Vector(x + v.x, y + v.y, z + v.z);
	}

	inline Vector operator-(const Vector &v) const {
		return Vector(x - v.x, y - v.y, z - v.z);
	}

	inline Vector& operator+=(const Vector &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}

	inline Vector& operator-=(const Vector &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}

	inline Vector operator*(Float f) const {
		return Vector(x*f, y*f, z*f);
	}

	inline Vector &operator*=(Float f) {
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}

	inline Vector operator-() const {
		return Vector(-x, -y, -z);
	}

	inline Vector operator/(Float f) const {
#ifdef MTS_DEBUG
		if (f == 0) {
			SLog(EWarn, "Vector: Division by zero!");
			//exit(-1);
		}
#endif
		Float r = 1.0f / f;
		return Vector(x * r, y * r, z * r);
	}

	inline Vector &operator/=(Float f) {
#ifdef MTS_DEBUG
		if (f == 0) {
			SLog(EWarn, "Vector: Division by zero!");
			//exit(-1);
		}
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

	inline Float lengthSquared() const {
		return x*x + y*y + z*z;
	}

	inline Float length() const {
		return std::sqrt(lengthSquared());
	}
	
	inline bool isZero() const {
		return x==0 && y == 0 && z == 0;
	}

	inline bool operator==(const Vector &v) const {
		return (v.x == x && v.y == y && v.z == z);
	}
	
	inline bool operator!=(const Vector &v) const {
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

inline Vector operator*(Float f, const Vector &v) {
	return v*f;
}

inline Float dot(const Vector &v1, const Vector &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline Float absDot(const Vector &v1, const Vector &v2) {
	return std::abs(dot(v1, v2));
}

inline Vector cross(const Vector &v1, const Vector &v2) {
	/* Left-handed vector cross product */
	return Vector(
		(v1.y * v2.z) - (v1.z * v2.y), 
		(v1.z * v2.x) - (v1.x * v2.z),
		(v1.x * v2.y) - (v1.y * v2.x)
	);
}

inline Vector normalize(const Vector &v) {
	return v / v.length();
}

/** \brief Simple two-dimensional vector class using single precision
 */
class Vector2 {
public:
	Float x, y;

	inline Vector2(Float _x = 0.0f, Float _y = 0.0f)
		: x(_x), y(_y) {
	}
	
	inline Vector2(const Point2 &p);
	
	inline Vector2(Stream *stream) {
		x = stream->readFloat();
		y = stream->readFloat();
	}
	
	inline Vector2 operator+(const Vector2 &v) const {
		return Vector2(x + v.x, y + v.y);
	}

	inline Vector2 operator-(const Vector2 &v) const {
		return Vector2(x - v.x, y - v.y);
	}

	inline Vector2& operator+=(const Vector2 &v) {
		x += v.x; y += v.y; 
		return *this;
	}

	inline Vector2& operator-=(const Vector2 &v) {
		x -= v.x; y -= v.y;
		return *this;
	}

	inline Vector2 operator*(Float f) const {
		return Vector2(x*f, y*f);
	}

	inline Vector2 &operator*=(Float f) {
		x *= f;
		y *= f;
		return *this;
	}

	inline Vector2 operator-() const {
		return Vector2(-x, -y);
	}

	inline Vector2 operator/(Float f) const {
		Float r = 1.0f / f;
		return Vector2(x * r, y * r);
	}

	inline Vector2 &operator/=(Float f) {
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

	inline Float lengthSquared() const {
		return x*x + y*y;
	}

	inline Float length() const {
		return std::sqrt(lengthSquared());
	}

	inline bool isZero() const {
		return x==0 && y == 0;
	}

	inline bool operator==(const Vector2 &v) const {
		return (v.x == x && v.y == y);
	}
	
	inline bool operator!=(const Vector2 &v) const {
		return !operator==(v);
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

inline Vector2 normalize(const Vector2 &v) {
	return v / v.length();
}

/** \brief Simple two-dimensional vector class using integer values.
 */
class Vector2i {
public:
	int x, y;

	inline Vector2i(int _x = 0, int _y = 0)
		: x(_x), y(_y) {
	}
	
	inline Vector2i(Stream *stream) {
		x = stream->readInt();
		y = stream->readInt();
	}
	
	
	inline Vector2i(const Point2i &p);
	
	inline Vector2i operator+(const Vector2i &v) const {
		return Vector2i(x + v.x, y + v.y);
	}

	inline Vector2i operator-(const Vector2i &v) const {
		return Vector2i(x - v.x, y - v.y);
	}

	inline Vector2i& operator+=(const Vector2i &v) {
		x += v.x; y += v.y; 
		return *this;
	}

	inline Vector2i& operator-=(const Vector2i &v) {
		x -= v.x; y -= v.y;
		return *this;
	}

	inline Vector2i operator*(int i) const {
		return Vector2i(x*i, y*i);
	}

	inline Vector2i &operator*=(int i) {
		x *= i;
		y *= i;
		return *this;
	}

	inline Vector2i operator-() const {
		return Vector2i(-x, -y);
	}

	inline Vector2i operator/(int i) const {
		return Vector2i(x / i, y / i);
	}

	inline Vector2i &operator/=(int i) {
		x /= i;
		y /= i;
		return *this;
	}

	inline int operator[](int i) const {
		return (&x)[i];
	}

	inline int &operator[](int i) {
		return (&x)[i];
	}

	inline int lengthSquared() const {
		return x*x + y*y;
	}

	int length() const {
		return (int) std::sqrt((float) lengthSquared());
	}
	
	inline bool isZero() const {
		return x==0 && y == 0;
	}

	inline bool operator==(const Vector2i &v) const {
		return (v.x == x && v.y == y);
	}
	
	inline bool operator!=(const Vector2i &v) const {
		return !operator==(v);
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

/** \brief Simple three-dimensional vector class using integer values.
 */
class Vector3i {
public:
	int x, y, z;
	
	inline Vector3i(int _x = 0, int _y = 0, int _z = 0)
		: x(_x), y(_y), z(_z) {
	}
	
	inline Vector3i(const Vector3i &v) 
		: x(v.x), y(v.y), z(v.z) {
	}
	
	inline Vector3i(const Point3i &v);
	
	inline Vector3i(Stream *stream) {
		x = stream->readInt();
		y = stream->readInt();
		z = stream->readInt();
	}
	
	explicit inline Vector3i(const Vector &v) 
		: x((int) v.x), y((int) v.y), z((int) v.z) {
	}

	inline Vector3i operator+(const Vector3i &v) const {
		return Vector3i(x + v.x, y + v.y, z + v.z);
	}

	inline Vector3i operator-(const Vector3i &v) const {
		return Vector3i(x - v.x, y - v.y, z - v.z);
	}

	inline Vector3i& operator+=(const Vector3i &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}

	inline Vector3i& operator-=(const Vector3i &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}

	inline Vector3i operator*(int i) const {
		return Vector3i(x*i, y*i, z*i);
	}

	inline Vector3i &operator*=(int i) {
		x *= i;
		y *= i;
		z *= i;
		return *this;
	}

	inline Vector3i operator-() const {
		return Vector3i(-x, -y, -z);
	}

	inline Vector3i operator/(int i) const {
		return Vector3i(x / i, y / i, z / i);
	}

	inline Vector3i &operator/=(int i) {
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

	inline int lengthSquared() const {
		return x*x + y*y + z*z;
	}
	
	int length() const {
		return (int) std::sqrt((float) lengthSquared());
	}

	inline bool isZero() const {
		return x==0 && y == 0 && z == 0;
	}

	inline bool operator==(const Vector3i &v) const {
		return (v.x == x && v.y == y && v.z == z);
	}
	
	inline bool operator!=(const Vector3i &v) const {
		return !operator==(v);
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

inline Vector::Vector(const Vector3i &v) 
 : x((Float) v.x), y((Float) v.y), z((Float) v.z) {
}

class Vector4 {
public:
	Float x, y, z, w;

	inline Vector4(Float _x = 0.0f, Float _y = 0.0f, Float _z = 0.0f, Float _w = 0.0f)
		: x(_x), y(_y), z(_z), w(_w) {
	}

	inline Vector4(Stream *stream) {
		x = stream->readFloat();
		y = stream->readFloat();
		z = stream->readFloat();
		w = stream->readFloat();
	}

	inline Vector4 operator+(const Vector4 &v) const {
		return Vector4(x + v.x, y + v.y, z + v.z, w + v.w);
	}

	inline Vector4 operator-(const Vector4 &v) const {
		return Vector4(x - v.x, y - v.y, z - v.z, w - v.w);
	}

	inline Vector4& operator+=(const Vector4 &v) {
		x += v.x; y += v.y; z += v.z; w += v.w;
		return *this;
	}

	inline Vector4& operator-=(const Vector4 &v) {
		x -= v.x; y -= v.y; z -= v.z; w -= v.w;
		return *this;
	}

	inline Vector4 operator*(Float f) const {
		return Vector4(x*f, y*f, z*f, w*f);
	}

	inline Vector4 &operator*=(Float f) {
		x *= f;
		y *= f;
		z *= f;
		w *= f;
		return *this;
	}

	inline Vector4 operator-() const {
		return Vector4(-x, -y, -z, -w);
	}

	inline Vector4 operator/(Float f) const {
#ifdef MTS_DEBUG
		if (f == 0) 
			SLog(EWarn, "Vector4: Division by zero!");
#endif
		Float r = 1.0f / f;
		return Vector4(x * r, y * r, z * r, w * r);
	}

	inline Vector4 &operator/=(Float f) {
#ifdef MTS_DEBUG
		if (f == 0) 
			SLog(EWarn, "Vector4: Division by zero!");
#endif
		Float r = 1.0f / f;
		x *= r;
		y *= r;
		z *= r;
		w *= r;
		return *this;
	}

	inline Float operator[](int i) const {
		return (&x)[i];
	}

	inline Float &operator[](int i) {
		return (&x)[i];
	}

	inline Float lengthSquared() const {
		return x*x + y*y + z*z + w*w;
	}

	inline Float length() const {
		return std::sqrt(lengthSquared());
	}
	
	inline bool isZero() const {
		return x==0 && y == 0 && z == 0 && w == 0;
	}

	inline bool operator==(const Vector4 &v) const {
		return (v.x == x && v.y == y && v.z == z && v.w == w);
	}
	
	inline bool operator!=(const Vector4 &v) const {
		return !operator==(v);
	}

	inline void serialize(Stream *stream) const {
		stream->writeFloat(x);
		stream->writeFloat(y);
		stream->writeFloat(z);
		stream->writeFloat(w);
	}

	inline std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << ", " << y << ", " << z << ", " << w << "]";
		return oss.str();
	}
};

inline Vector4 operator*(Float f, const Vector4 &v) {
	return v*f;
}

inline Float dot(const Vector4 &v1, const Vector4 &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.w * v2.w;
}

inline Float absDot(const Vector4 &v1, const Vector4 &v2) {
	return std::abs(dot(v1, v2));
}

inline Vector4 normalize(const Vector4 &v) {
	return v / v.length();
}

MTS_NAMESPACE_END

#endif /* __VECTOR_H */

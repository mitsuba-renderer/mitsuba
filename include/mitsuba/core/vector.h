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
#if !defined(__MITSUBA_CORE_VECTOR_H_)
#define __MITSUBA_CORE_VECTOR_H_

#include <mitsuba/core/stream.h>

MTS_NAMESPACE_BEGIN
/**
 * \headerfile mitsuba/core/vector.h mitsuba/mitsuba.h
 * \brief Parameterizable one-dimensional vector data structure
 * \ingroup libcore
 */
template <typename T> struct TVector1 {
	typedef T          Scalar;
	typedef TPoint1<T> PointType;

	T x;

	/// Number of dimensions
	const static int dim = 1;

	/** \brief Construct a new vector without initializing it.
	 *
	 * This construtor is useful when the vector will either not
	 * be used at all (it might be part of a larger data structure)
	 * or initialized at a later point in time. Always make sure
	 * that one of the two is the case! Otherwise your program will do
	 * computations involving uninitialized memory, which will probably
	 * lead to a difficult-to-find bug.
	 */
#if !defined(MTS_DEBUG_UNINITIALIZED)
	TVector1() { }
#else
	TVector1() { x = std::numeric_limits<T>::quiet_NaN(); }
#endif

	/// Initialize the vector with the specified value
	TVector1(T x) : x(x) {  }

	/// Initialize the vector with the components of another vector data structure
	template <typename T1> explicit TVector1(const TVector1<T1> &v)
		: x((T) v.x) { }

	/// Initialize the vector with the components of a point data structure
	template <typename T1> explicit TVector1(const TPoint1<T1> &p)
		: x((T) p.x) { }

	/// Unserialize a vector from a binary data stream
	explicit TVector1(Stream *stream) {
		x = stream->readElement<T>();
	}

	/// Add two vectors and return the result
	TVector1 operator+(const TVector1 &v) const {
		return TVector1(x + v.x);
	}

	/// Subtract two vectors and return the result
	TVector1 operator-(const TVector1 &v) const {
		return TVector1(x - v.x);
	}

	/// Add another vector to the current one
	TVector1& operator+=(const TVector1 &v) {
		x += v.x;
		return *this;
	}

	/// Subtract another vector from the current one
	TVector1& operator-=(const TVector1 &v) {
		x -= v.x;
		return *this;
	}

	/// Multiply the vector by the given scalar and return the result
	TVector1 operator*(T f) const {
		return TVector1(x*f);
	}

	/// Multiply the vector by the given scalar
	TVector1 &operator*=(T f) {
		x *= f;
		return *this;
	}

	/// Return a negated version of the vector
	TVector1 operator-() const {
		return TVector1(-x);
	}

	/// Divide the vector by the given scalar and return the result
	TVector1 operator/(T f) const {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Vector1: Division by zero!");
#endif
		return TVector1(x / f);
	}

	/// Divide the vector by the given scalar
	TVector1 &operator/=(T f) {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Vector1: Division by zero!");
#endif
		x /= f;
		return *this;
	}

	/// Index into the vector's components
	T &operator[](int i) {
		return (&x)[i];
	}

	/// Index into the vector's components (const version)
	T operator[](int i) const {
		return (&x)[i];
	}

	/// Return the squared 1-norm of this vector
	T lengthSquared() const {
		return x*x;
	}

	/// Return the 1-norm of this vector
	T length() const {
		return std::abs(x);
	}

	/// Return whether or not this vector is identically zero
	bool isZero() const {
		return x == 0;
	}

	/// Equality test
	bool operator==(const TVector1 &v) const {
		return v.x == x;
	}

	/// Inequality test
	bool operator!=(const TVector1 &v) const {
		return v.x != x;
	}

	/// Serialize this vector to a binary data stream
	void serialize(Stream *stream) const {
		stream->writeElement<T>(x);
	}

	/// Implicit conversion to Scalar
	operator Scalar() const { return x; }

	/// Return a readable string representation of this vector
	std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << "]";
		return oss.str();
	}
};

template <typename T> inline TVector1<T> operator*(T f, const TVector1<T> &v) {
	return v*f;
}

template <typename T> inline T dot(const TVector1<T> &v1, const TVector1<T> &v2) {
	return v1.x * v2.x;
}

template <typename T> inline T absDot(const TVector1<T> &v1, const TVector1<T> &v2) {
	return std::abs(dot(v1, v2));
}

template <typename T> inline TVector1<T> normalize(const TVector1<T> &v) {
	return v / v.length();
}

template <typename T> inline TVector1<T> normalizeStrict(const TVector1<T> &v, const char *errMsg) {
	Float length = v.length();
	if (length == 0)
		SLog(EError, "normalizeStrict(): %s", errMsg);
	return v / length;
}

template <> inline TVector1<int> TVector1<int>::operator/(int s) const {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Vector1i: Division by zero!");
#endif
	return TVector1(x/s);
}

template <> inline TVector1<int> &TVector1<int>::operator/=(int s) {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Vector1i: Division by zero!");
#endif

	x /= s;
	return *this;
}

namespace detail {
// Helper structures to define "length" as a floating point value

template <typename T, bool IsInteger = false> struct VectorLength {
	typedef T type;
};

template <typename T> struct VectorLength<T, true> {
	typedef Float type;
};

template <> struct VectorLength<uint64_t, true> {
	typedef double type;
};

template <> struct VectorLength<int64_t, true> {
	typedef double type;
};

}

/**
 * \headerfile mitsuba/core/vector.h mitsuba/mitsuba.h
 * \brief Parameterizable two-dimensional vector data structure
 * \ingroup libcore
 */
template <typename T> struct TVector2 {
	typedef T          Scalar;
	typedef TPoint2<T> PointType;
	typedef typename detail::VectorLength<T, std::numeric_limits<T>::is_integer>::type LengthType;

	T x, y;

	/// Number of dimensions
	const static int dim = 2;

	/** \brief Construct a new vector without initializing it.
	 *
	 * This construtor is useful when the vector will either not
	 * be used at all (it might be part of a larger data structure)
	 * or initialized at a later point in time. Always make sure
	 * that one of the two is the case! Otherwise your program will do
	 * computations involving uninitialized memory, which will probably
	 * lead to a difficult-to-find bug.
	 */
#if !defined(MTS_DEBUG_UNINITIALIZED)
	TVector2() { }
#else
	TVector2() { x = y = std::numeric_limits<T>::quiet_NaN(); }
#endif

	/// Initialize the vector with the specified X and Z components
	TVector2(T x, T y) : x(x), y(y) {  }

	/// Initialize all components of the the vector with the specified value
	explicit TVector2(T val) : x(val), y(val) { }

	/// Initialize the vector with the components of another vector data structure
	template <typename T2> explicit TVector2(const TVector2<T2> &v)
		: x((T) v.x), y((T) v.y) { }

	/// Initialize the vector with the components of a point data structure
	template <typename T2> explicit TVector2(const TPoint2<T2> &p)
		: x((T) p.x), y((T) p.y) { }

	/// Unserialize a vector from a binary data stream
	explicit TVector2(Stream *stream) {
		x = stream->readElement<T>();
		y = stream->readElement<T>();
	}

	/// Add two vectors and return the result
	TVector2 operator+(const TVector2 &v) const {
		return TVector2(x + v.x, y + v.y);
	}

	/// Subtract two vectors and return the result
	TVector2 operator-(const TVector2 &v) const {
		return TVector2(x - v.x, y - v.y);
	}

	/// Add another vector to the current one
	TVector2& operator+=(const TVector2 &v) {
		x += v.x; y += v.y;
		return *this;
	}

	/// Subtract another vector from the current one
	TVector2& operator-=(const TVector2 &v) {
		x -= v.x; y -= v.y;
		return *this;
	}

	/// Multiply the vector by the given scalar and return the result
	TVector2 operator*(T f) const {
		return TVector2(x*f, y*f);
	}

	/// Multiply the vector by the given scalar
	TVector2 &operator*=(T f) {
		x *= f; y *= f;
		return *this;
	}

	/// Return a negated version of the vector
	TVector2 operator-() const {
		return TVector2(-x, -y);
	}

	/// Divide the vector by the given scalar and return the result
	TVector2 operator/(T f) const {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Vector2: Division by zero!");
#endif
		T recip = (T) 1 / f;
		return TVector2(x * recip, y * recip);
	}

	/// Divide the vector by the given scalar
	TVector2 &operator/=(T f) {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Vector2: Division by zero!");
#endif
		T recip = (T) 1 / f;
		x *= recip; y *= recip;
		return *this;
	}

	/// Index into the vector's components
	T &operator[](int i) {
		return (&x)[i];
	}

	/// Index into the vector's components (const version)
	T operator[](int i) const {
		return (&x)[i];
	}

	/// Return the squared 2-norm of this vector
	T lengthSquared() const {
		return x*x + y*y;
	}

	/// Return the 2-norm of this vector
	LengthType length() const {
		return (LengthType) std::sqrt((LengthType) lengthSquared());
	}

	/// Return whether or not this vector is identically zero
	bool isZero() const {
		return x == 0 && y == 0;
	}

	/// Equality test
	bool operator==(const TVector2 &v) const {
		return (v.x == x && v.y == y);
	}

	/// Inequality test
	bool operator!=(const TVector2 &v) const {
		return v.x != x || v.y != y;
	}

	/// Serialize this vector to a binary data stream
	void serialize(Stream *stream) const {
		stream->writeElement<T>(x);
		stream->writeElement<T>(y);
	}

	/// Return a readable string representation of this vector
	std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << ", " << y << "]";
		return oss.str();
	}
};

template <typename T> inline TVector2<T> operator*(T f, const TVector2<T> &v) {
	return v*f;
}

template <typename T> inline T dot(const TVector2<T> &v1, const TVector2<T> &v2) {
	return v1.x * v2.x + v1.y * v2.y;
}

template <typename T> inline T absDot(const TVector2<T> &v1, const TVector2<T> &v2) {
	return std::abs(dot(v1, v2));
}

template <typename T> inline T det(const TVector2<T> &v1, const TVector2<T> &v2) {
	return v1.x * v2.y - v1.y * v2.x;
}

template <typename T> inline TVector2<T> normalize(const TVector2<T> &v) {
	return v / v.length();
}

template <typename T> inline TVector2<T> normalizeStrict(const TVector2<T> &v, const char *errMsg) {
	Float length = v.length();
	if (length == 0)
		SLog(EError, "normalizeStrict(): %s", errMsg);
	return v / length;
}

template <> inline TVector2<int> TVector2<int>::operator/(int s) const {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Vector2i: Division by zero!");
#endif
	return TVector2(x/s, y/s);
}

template <> inline TVector2<int> &TVector2<int>::operator/=(int s) {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Vector2i: Division by zero!");
#endif

	x /= s;
	y /= s;
	return *this;
}

/**
 * \headerfile mitsuba/core/vector.h mitsuba/mitsuba.h
 * \brief Parameterizable three-dimensional vector data structure
 * \ingroup libcore
 */
template <typename T> struct TVector3 {
	typedef T          Scalar;
	typedef TPoint3<T> PointType;
	typedef typename detail::VectorLength<T, std::numeric_limits<T>::is_integer>::type LengthType;

	T x, y, z;

	/// Number of dimensions
	const static int dim = 3;

	/** \brief Construct a new vector without initializing it.
	 *
	 * This construtor is useful when the vector will either not
	 * be used at all (it might be part of a larger data structure)
	 * or initialized at a later point in time. Always make sure
	 * that one of the two is the case! Otherwise your program will do
	 * computations involving uninitialized memory, which will probably
	 * lead to a difficult-to-find bug.
	 */
#if !defined(MTS_DEBUG_UNINITIALIZED)
	TVector3() { }
#else
	TVector3() { x = y = z = std::numeric_limits<T>::quiet_NaN(); }
#endif

	/// Initialize the vector with the specified X, Y and Z components
	TVector3(T x, T y, T z) : x(x), y(y), z(z) {  }

	/// Initialize all components of the the vector with the specified value
	explicit TVector3(T val) : x(val), y(val), z(val) { }

	/// Initialize the vector with the components of another vector data structure
	template <typename T2> explicit TVector3(const TVector3<T2> &v)
		: x((T) v.x), y((T) v.y), z((T) v.z) { }

	/// Initialize the vector with the components of a point data structure
	template <typename T2> explicit TVector3(const TPoint3<T2> &p)
		: x((T) p.x), y((T) p.y), z((T) p.z) { }

	/// Unserialize a vector from a binary data stream
	explicit TVector3(Stream *stream) {
		x = stream->readElement<T>();
		y = stream->readElement<T>();
		z = stream->readElement<T>();
	}

	/// Add two vectors and return the result
	TVector3 operator+(const TVector3 &v) const {
		return TVector3(x + v.x, y + v.y, z + v.z);
	}

	/// Subtract two vectors and return the result
	TVector3 operator-(const TVector3 &v) const {
		return TVector3(x - v.x, y - v.y, z - v.z);
	}

	/// Add another vector to the current one
	TVector3& operator+=(const TVector3 &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}

	/// Subtract another vector from the current one
	TVector3& operator-=(const TVector3 &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}

	/// Multiply the vector by the given scalar and return the result
	TVector3 operator*(T f) const {
		return TVector3(x*f, y*f, z*f);
	}

	/// Multiply the vector by the given scalar
	TVector3 &operator*=(T f) {
		x *= f; y *= f; z *= f;
		return *this;
	}

	/// Return a negated version of the vector
	TVector3 operator-() const {
		return TVector3(-x, -y, -z);
	}

	/// Divide the vector by the given scalar and return the result
	TVector3 operator/(T f) const {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Vector3: Division by zero!");
#endif
		T recip = (T) 1 / f;
		return TVector3(x * recip, y * recip, z * recip);
	}

	/// Divide the vector by the given scalar
	TVector3 &operator/=(T f) {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Vector3: Division by zero!");
#endif
		T recip = (T) 1 / f;
		x *= recip; y *= recip; z *= recip;
		return *this;
	}

	/// Index into the vector's components
	T &operator[](int i) {
		return (&x)[i];
	}

	/// Index into the vector's components (const version)
	T operator[](int i) const {
		return (&x)[i];
	}

	/// Return the squared 2-norm of this vector
	T lengthSquared() const {
		return x*x + y*y + z*z;
	}

	/// Return the 2-norm of this vector
	LengthType length() const {
		return (LengthType) std::sqrt((LengthType) lengthSquared());
	}

	/// Return whether or not this vector is identically zero
	bool isZero() const {
		return x == 0 && y == 0 && z == 0;
	}

	/// Equality test
	bool operator==(const TVector3 &v) const {
		return (v.x == x && v.y == y && v.z == z);
	}

	/// Inequality test
	bool operator!=(const TVector3 &v) const {
		return v.x != x || v.y != y || v.z != z;
	}

	/// Serialize this vector to a binary data stream
	void serialize(Stream *stream) const {
		stream->writeElement<T>(x);
		stream->writeElement<T>(y);
		stream->writeElement<T>(z);
	}

	/// Return a readable string representation of this vector
	std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << ", " << y << ", " << z << "]";
		return oss.str();
	}
};

template <typename T> inline TVector3<T> operator*(T f, const TVector3<T> &v) {
	return v*f;
}

template <typename T> inline T dot(const TVector3<T> &v1, const TVector3<T> &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <typename T> inline T absDot(const TVector3<T> &v1, const TVector3<T> &v2) {
	return std::abs(dot(v1, v2));
}

template <typename T> inline TVector3<T> cross(const TVector3<T> &v1, const TVector3<T> &v2) {
	return TVector3<T>(
		(v1.y * v2.z) - (v1.z * v2.y),
		(v1.z * v2.x) - (v1.x * v2.z),
		(v1.x * v2.y) - (v1.y * v2.x)
	);
}

template <typename T> inline TVector3<T> normalize(const TVector3<T> &v) {
	return v / v.length();
}

template <typename T> inline TVector3<T> normalizeStrict(const TVector3<T> &v, const char *errMsg) {
	Float length = v.length();
	if (length == 0)
		SLog(EError, "normalizeStrict(): %s", errMsg);
	return v / length;
}

template <> inline TVector3<int> TVector3<int>::operator/(int s) const {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Vector3i: Division by zero!");
#endif
	return TVector3(x/s, y/s, z/s);
}

template <> inline TVector3<int> &TVector3<int>::operator/=(int s) {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Vector3i: Division by zero!");
#endif

	x /= s;
	y /= s;
	z /= s;
	return *this;
}


/**
 * \headerfile mitsuba/core/vector.h mitsuba/mitsuba.h
 * \brief Parameterizable four-dimensional vector data structure
 * \ingroup libcore
 */
template <typename T> struct TVector4 {
	typedef T          Scalar;
	typedef TPoint4<T> PointType;
	typedef typename detail::VectorLength<T, std::numeric_limits<T>::is_integer>::type LengthType;

	T x, y, z, w;

	/// Number of dimensions
	const static int dim = 4;

	/** \brief Construct a new vector without initializing it.
	 *
	 * This construtor is useful when the vector will either not
	 * be used at all (it might be part of a larger data structure)
	 * or initialized at a later point in time. Always make sure
	 * that one of the two is the case! Otherwise your program will do
	 * computations involving uninitialized memory, which will probably
	 * lead to a difficult-to-find bug.
	 */
#if !defined(MTS_DEBUG_UNINITIALIZED)
	TVector4() { }
#else
	TVector4() { x = y = z = w = std::numeric_limits<T>::quiet_NaN(); }
#endif

	/// Initialize the vector with the specified X, Y and Z components
	TVector4(T x, T y, T z, T w) : x(x), y(y), z(z), w(w) {  }

	/// Initialize all components of the the vector with the specified value
	explicit TVector4(T val) : x(val), y(val), z(val), w(val) { }

	/// Initialize the vector with the components of another vector data structure
	template <typename T2> explicit TVector4(const TVector4<T2> &v)
		: x((T) v.x), y((T) v.y), z((T) v.z), w((T) v.w) { }

	/// Initialize the vector with the components of a point data structure
	template <typename T2> explicit TVector4(const TPoint4<T2> &p)
		: x((T) p.x), y((T) p.y), z((T) p.z), w((T) p.w) { }

	/// Unserialize a vector from a binary data stream
	explicit TVector4(Stream *stream) {
		x = stream->readElement<T>();
		y = stream->readElement<T>();
		z = stream->readElement<T>();
		w = stream->readElement<T>();
	}

	/// Add two vectors and return the result
	TVector4 operator+(const TVector4 &v) const {
		return TVector4(x + v.x, y + v.y, z + v.z, w + v.w);
	}

	/// Subtract two vectors and return the result
	TVector4 operator-(const TVector4 &v) const {
		return TVector4(x - v.x, y - v.y, z - v.z, w - v.w);
	}

	/// Add another vector to the current one
	TVector4& operator+=(const TVector4 &v) {
		x += v.x; y += v.y; z += v.z; w += v.w;
		return *this;
	}

	/// Subtract another vector from the current one
	TVector4& operator-=(const TVector4 &v) {
		x -= v.x; y -= v.y; z -= v.z; w -= v.w;
		return *this;
	}

	/// Multiply the vector by the given scalar and return the result
	TVector4 operator*(T f) const {
		return TVector4(x * f, y * f, z * f, w * f);
	}

	/// Multiply the vector by the given scalar
	TVector4 &operator*=(T f) {
		x *= f; y *= f; z *= f; w *= f;
		return *this;
	}

	/// Return a negated version of the vector
	TVector4 operator-() const {
		return TVector4(-x, -y, -z, -w);
	}

	/// Divide the vector by the given scalar and return the result
	TVector4 operator/(T f) const {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Vector4: Division by zero!");
#endif
		T recip = (T) 1 / f;
		return TVector4(x * recip, y * recip, z * recip, w * recip);
	}

	/// Divide the vector by the given scalar
	TVector4 &operator/=(T f) {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Vector4: Division by zero!");
#endif
		T recip = (T) 1 / f;
		x *= recip; y *= recip; z *= recip; w *= recip;
		return *this;
	}

	/// Index into the vector's components
	T &operator[](int i) {
		return (&x)[i];
	}

	/// Index into the vector's components (const version)
	T operator[](int i) const {
		return (&x)[i];
	}

	/// Return the squared 2-norm of this vector
	T lengthSquared() const {
		return x*x + y*y + z*z + w*w;
	}

	/// Return the 2-norm of this vector
	LengthType length() const {
		return (LengthType) std::sqrt((LengthType) lengthSquared());
	}

	/// Return whether or not this vector is identically zero
	bool isZero() const {
		return x == 0 && y == 0 && z == 0 && w == 0;
	}

	/// Equality test
	bool operator==(const TVector4 &v) const {
		return (v.x == x && v.y == y && v.z == z && v.w == w);
	}

	/// Inequality test
	bool operator!=(const TVector4 &v) const {
		return v.x != x || v.y != y || v.z != z || v.w != w;
	}

	/// Serialize this vector to a binary data stream
	void serialize(Stream *stream) const {
		stream->writeElement<T>(x);
		stream->writeElement<T>(y);
		stream->writeElement<T>(z);
		stream->writeElement<T>(w);
	}

	/// Return a readable string representation of this vector
	std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << ", " << y << ", " << z << ", " << w << "]";
		return oss.str();
	}
};

template <typename T> inline TVector4<T> operator*(T f, const TVector4<T> &v) {
	return v*f;
}

template <typename T> inline T dot(const TVector4<T> &v1, const TVector4<T> &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.w * v2.w;
}

template <typename T> inline T absDot(const TVector4<T> &v1, const TVector4<T> &v2) {
	return std::abs(dot(v1, v2));
}

template <typename T> inline TVector4<T> normalize(const TVector4<T> &v) {
	return v / v.length();
}

template <typename T> inline TVector4<T> normalizeStrict(const TVector4<T> &v, const char *errMsg) {
	Float length = v.length();
	if (length == 0)
		SLog(EError, "normalizeStrict(): %s", errMsg);
	return v / length;
}

template <> inline TVector4<int> TVector4<int>::operator/(int s) const {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Vector4i: Division by zero!");
#endif
	return TVector4(x/s, y/s, z/s, w/s);
}

template <> inline TVector4<int> &TVector4<int>::operator/=(int s) {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Vector4i: Division by zero!");
#endif

	x /= s;
	y /= s;
	z /= s;
	w /= s;
	return *this;
}

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_VECTOR_H_ */

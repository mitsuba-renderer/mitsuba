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
#if !defined(__MITSUBA_CORE_POINT_H_)
#define __MITSUBA_CORE_POINT_H_

#include <mitsuba/core/vector.h>

MTS_NAMESPACE_BEGIN

/**
 * \headerfile mitsuba/core/point.h mitsuba/mitsuba.h
 * \brief Parameterizable one-dimensional point data structure
 *
 * \ingroup libcore
 */
template <typename T> struct TPoint1 {
	typedef T           Scalar;
	typedef TVector1<T> VectorType;

	T x;

	/// Number of dimensions
	const static int dim = 1;

	/** \brief Construct a new point without initializing it.
	 *
	 * This construtor is useful when the point will either not
	 * be used at all (it might be part of a larger data structure)
	 * or initialized at a later point in time. Always make sure
	 * that one of the two is the case! Otherwise your program will do
	 * computations involving uninitialized memory, which will probably
	 * lead to a difficult-to-find bug.
	 */
#if !defined(MTS_DEBUG_UNINITIALIZED)
	TPoint1() { }
#else
	TPoint1() { x = std::numeric_limits<T>::quiet_NaN(); }
#endif

	/// Initialize the point with the specified value
	TPoint1(T x) : x(x) {  }

	/// Initialize the point with the components of another point
	template <typename T1> explicit TPoint1(const TPoint1<T1> &p)
		: x((T) p.x) { }

	/// Initialize the point with the components of a vector data structure
	template <typename T1> explicit TPoint1(const TVector1<T1> &v)
		: x((T) v.x) { }

	/// Unserialize a point from a binary data stream
	explicit TPoint1(Stream *stream) {
		x = stream->readElement<T>();
	}

	/// Add a vector to a point and return the result
	TPoint1 operator+(const TVector1<T> &v) const {
		return TPoint1(x + v.x);
	}

	/// Add two points and return the result (e.g. to compute a weighted position)
	TPoint1 operator+(const TPoint1 &p) const {
		return TPoint1(x + p.x);
	}

	/// Add a vector to this one (e.g. to compute a weighted position)
	TPoint1& operator+=(const TVector1<T> &v) {
		x += v.x;
		return *this;
	}

	/// Add a point to this one (e.g. to compute a weighted position)
	TPoint1& operator+=(const TPoint1 &p) {
		x += p.x;
		return *this;
	}

	/// Subtract a vector from this point
	TPoint1 operator-(const TVector1<T> &v) const {
		return TPoint1(x - v.x);
	}

	/// Subtract two points from each other and return the difference as a vector
	TVector1<T> operator-(const TPoint1 &p) const {
		return TVector1<T>(x - p.x);
	}

	/// Subtract a vector from this point
	TPoint1& operator-=(const TVector1<T> &v) {
		x -= v.x;
		return *this;
	}

	/// Scale the point's coordinates by the given scalar and return the result
	TPoint1 operator*(T f) const {
		return TPoint1(x * f);
	}

	/// Scale the point's coordinates by the given scalar
	TPoint1 &operator*=(T f) {
		x *= f;
		return *this;
	}

	/// Return a version of the point, which has been flipped along the origin
	TPoint1 operator-() const {
		return TPoint1(-x);
	}

	/// Divide the point's coordinates by the given scalar and return the result
	TPoint1 operator/(T f) const {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Point1: Division by zero!");
#endif
		return TPoint1(x / f);
	}

	/// Divide the point's coordinates by the given scalar
	TPoint1 &operator/=(T f) {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Point1: Division by zero!");
#endif
		x /= f;
		return *this;
	}

	/// Index into the point's components
	T &operator[](int i) {
		return (&x)[i];
	}

	/// Index into the point's components (const version)
	T operator[](int i) const {
		return (&x)[i];
	}

	/// Return whether or not this point is identically zero
	bool isZero() const {
		return x == 0;
	}

	/// Equality test
	bool operator==(const TPoint1 &v) const {
		return (v.x == x);
	}

	/// Inequality test
	bool operator!=(const TPoint1 &v) const {
		return v.x != x;
	}

	/// Serialize this point to a binary data stream
	void serialize(Stream *stream) const {
		stream->writeElement<T>(x);
	}

	/// Implicit conversion to Scalar
	operator Scalar() const { return x; }

	/// Return a readable string representation of this point
	std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << "]";
		return oss.str();
	}
};

template <typename T> inline TPoint1<T> operator*(T f, const TPoint1<T> &v) {
	return v*f;
}

template <typename T> inline T distance(const TPoint1<T> &p1, const TPoint1<T> &p2) {
	return std::abs(p2.x-p2.x);
}

template <typename T> inline T distanceSquared(const TPoint1<T> &p1, const TPoint1<T> &p2) {
	return (p1-p2).lengthSquared();
}

template <> inline TPoint1<int> TPoint1<int>::operator/(int s) const {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Point1i: Division by zero!");
#endif
	return TPoint1(x/s);
}

template <> inline TPoint1<int> &TPoint1<int>::operator/=(int s) {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Point1i: Division by zero!");
#endif
	x /= s;
	return *this;
}

/**
 * \headerfile mitsuba/core/point.h mitsuba/mitsuba.h
 * \brief Parameterizable two-dimensional point data structure
 *
 * \ingroup libcore
 */
template <typename T> struct TPoint2 {
	typedef T           Scalar;
	typedef TVector2<T> VectorType;

	T x, y;

	/// Number of dimensions
	const static int dim = 2;

	/** \brief Construct a new point without initializing it.
	 *
	 * This construtor is useful when the point will either not
	 * be used at all (it might be part of a larger data structure)
	 * or initialized at a later point in time. Always make sure
	 * that one of the two is the case! Otherwise your program will do
	 * computations involving uninitialized memory, which will probably
	 * lead to a difficult-to-find bug.
	 */
#if !defined(MTS_DEBUG_UNINITIALIZED)
	TPoint2() { }
#else
	TPoint2() { x = y = std::numeric_limits<T>::quiet_NaN(); }
#endif

	/// Initialize the point with the specified X, Y and Z components
	TPoint2(T x, T y) : x(x), y(y) {  }

	/// Initialize the point with the components of another point
	template <typename T2> explicit TPoint2(const TPoint2<T2> &p)
		: x((T) p.x), y((T) p.y) { }

	/// Initialize the point with the components of a vector data structure
	template <typename T2> explicit TPoint2(const TVector2<T2> &v)
		: x((T) v.x), y((T) v.y) { }

	/// Initialize all components of the the point with the specified value
	explicit TPoint2(T val) : x(val), y(val) { }

	/// Unserialize a point from a binary data stream
	explicit TPoint2(Stream *stream) {
		x = stream->readElement<T>();
		y = stream->readElement<T>();
	}

	/// Add a vector to a point and return the result
	TPoint2 operator+(const TVector2<T> &v) const {
		return TPoint2(x + v.x, y + v.y);
	}

	/// Add two points and return the result (e.g. to compute a weighted position)
	TPoint2 operator+(const TPoint2 &p) const {
		return TPoint2(x + p.x, y + p.y);
	}

	/// Add a vector to this one (e.g. to compute a weighted position)
	TPoint2& operator+=(const TVector2<T> &v) {
		x += v.x; y += v.y;
		return *this;
	}

	/// Add a point to this one (e.g. to compute a weighted position)
	TPoint2& operator+=(const TPoint2 &p) {
		x += p.x; y += p.y;
		return *this;
	}

	/// Subtract a vector from this point
	TPoint2 operator-(const TVector2<T> &v) const {
		return TPoint2(x - v.x, y - v.y);
	}

	/// Subtract two points from each other and return the difference as a vector
	TVector2<T> operator-(const TPoint2 &p) const {
		return TVector2<T>(x - p.x, y - p.y);
	}

	/// Subtract a vector from this point
	TPoint2& operator-=(const TVector2<T> &v) {
		x -= v.x; y -= v.y;
		return *this;
	}

	/// Scale the point's coordinates by the given scalar and return the result
	TPoint2 operator*(T f) const {
		return TPoint2(x * f, y * f);
	}

	/// Scale the point's coordinates by the given scalar
	TPoint2 &operator*=(T f) {
		x *= f; y *= f;
		return *this;
	}

	/// Return a version of the point, which has been flipped along the origin
	TPoint2 operator-() const {
		return TPoint2(-x, -y);
	}

	/// Divide the point's coordinates by the given scalar and return the result
	TPoint2 operator/(T f) const {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Point2: Division by zero!");
#endif
		T recip = (T) 1 / f;
		return TPoint2(x * recip, y * recip);
	}

	/// Divide the point's coordinates by the given scalar
	TPoint2 &operator/=(T f) {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Point2: Division by zero!");
#endif
		T recip = (T) 1 / f;
		x *= recip; y *= recip;
		return *this;
	}

	/// Index into the point's components
	T &operator[](int i) {
		return (&x)[i];
	}

	/// Index into the point's components (const version)
	T operator[](int i) const {
		return (&x)[i];
	}

	/// Return whether or not this point is identically zero
	bool isZero() const {
		return x == 0 && y == 0;
	}

	/// Equality test
	bool operator==(const TPoint2 &v) const {
		return (v.x == x && v.y == y);
	}

	/// Inequality test
	bool operator!=(const TPoint2 &v) const {
		return v.x != x || v.y != y;
	}

	/// Serialize this point to a binary data stream
	void serialize(Stream *stream) const {
		stream->writeElement<T>(x);
		stream->writeElement<T>(y);
	}

	/// Return a readable string representation of this point
	std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << ", " << y << "]";
		return oss.str();
	}
};

template <typename T> inline TPoint2<T> operator*(T f, const TPoint2<T> &v) {
	return v*f;
}

template <typename T> inline T distance(const TPoint2<T> &p1, const TPoint2<T> &p2) {
	return (p1-p2).length();
}

template <typename T> inline T distanceSquared(const TPoint2<T> &p1, const TPoint2<T> &p2) {
	return (p1-p2).lengthSquared();
}

template <> inline TPoint2<int> TPoint2<int>::operator/(int s) const {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Point2i: Division by zero!");
#endif
	return TPoint2(x/s, y/s);
}

template <> inline TPoint2<int> &TPoint2<int>::operator/=(int s) {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Point2i: Division by zero!");
#endif
	x /= s;
	y /= s;
	return *this;
}

/**
 * \headerfile mitsuba/core/point.h mitsuba/mitsuba.h
 * \brief Parameterizable three-dimensional point data structure
 * \ingroup libcore
 */
template <typename T> struct TPoint3 {
	typedef T           Scalar;
	typedef TVector3<T> VectorType;

	T x, y, z;

	/// Number of dimensions
	const static int dim = 3;

	/** \brief Construct a new point without initializing it.
	 *
	 * This construtor is useful when the point will either not
	 * be used at all (it might be part of a larger data structure)
	 * or initialized at a later point in time. Always make sure
	 * that one of the two is the case! Otherwise your program will do
	 * computations involving uninitialized memory, which will probably
	 * lead to a difficult-to-find bug.
	 */
#if !defined(MTS_DEBUG_UNINITIALIZED)
	TPoint3() { }
#else
	TPoint3() { x = y = z = std::numeric_limits<T>::quiet_NaN(); }
#endif

	/// Initialize the point with the specified X, Y and Z components
	TPoint3(T x, T y, T z) : x(x), y(y), z(z) {  }

	/// Initialize the point with the components of another point
	template <typename T2> explicit TPoint3(const TPoint3<T2> &p)
		: x((T) p.x), y((T) p.y), z((T) p.z) { }

	/// Initialize the point with the components of a vector data structure
	template <typename T2> explicit TPoint3(const TVector3<T2> &v)
		: x((T) v.x), y((T) v.y), z((T) v.z) { }

	/// Initialize all components of the the point with the specified value
	explicit TPoint3(T val) : x(val), y(val), z(val) { }

	/// Unserialize a point from a binary data stream
	explicit TPoint3(Stream *stream) {
		x = stream->readElement<T>();
		y = stream->readElement<T>();
		z = stream->readElement<T>();
	}

	/// Add a vector to a point and return the result
	TPoint3 operator+(const TVector3<T> &v) const {
		return TPoint3(x + v.x, y + v.y, z + v.z);
	}

	/// Add two points and return the result (e.g. to compute a weighted position)
	TPoint3 operator+(const TPoint3 &p) const {
		return TPoint3(x + p.x, y + p.y, z + p.z);
	}

	/// Add a vector to this one (e.g. to compute a weighted position)
	TPoint3& operator+=(const TVector3<T> &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}

	/// Add a point to this one (e.g. to compute a weighted position)
	TPoint3& operator+=(const TPoint3 &p) {
		x += p.x; y += p.y; z += p.z;
		return *this;
	}

	/// Subtract a vector from this point
	TPoint3 operator-(const TVector3<T> &v) const {
		return TPoint3(x - v.x, y - v.y, z - v.z);
	}

	/// Subtract two points from each other and return the difference as a vector
	TVector3<T> operator-(const TPoint3 &p) const {
		return TVector3<T>(x - p.x, y - p.y, z - p.z);
	}

	/// Subtract a vector from this point
	TPoint3& operator-=(const TVector3<T> &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}

	/// Scale the point's coordinates by the given scalar and return the result
	TPoint3 operator*(T f) const {
		return TPoint3(x * f, y * f, z * f);
	}

	/// Scale the point's coordinates by the given scalar
	TPoint3 &operator*=(T f) {
		x *= f; y *= f; z *= f;
		return *this;
	}

	/// Return a version of the point, which has been flipped along the origin
	TPoint3 operator-() const {
		return TPoint3(-x, -y, -z);
	}

	/// Divide the point's coordinates by the given scalar and return the result
	TPoint3 operator/(T f) const {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Point3: Division by zero!");
#endif
		T recip = (T) 1 / f;
		return TPoint3(x * recip, y * recip, z * recip);
	}

	/// Divide the point's coordinates by the given scalar
	TPoint3 &operator/=(T f) {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Point3: Division by zero!");
#endif
		T recip = (T) 1 / f;
		x *= recip; y *= recip; z *= recip;
		return *this;
	}

	/// Index into the point's components
	T &operator[](int i) {
		return (&x)[i];
	}

	/// Index into the point's components (const version)
	T operator[](int i) const {
		return (&x)[i];
	}

	/// Return whether or not this point is identically zero
	bool isZero() const {
		return x == 0 && y == 0 && z == 0;
	}

	/// Equality test
	bool operator==(const TPoint3 &v) const {
		return (v.x == x && v.y == y && v.z == z);
	}

	/// Inequality test
	bool operator!=(const TPoint3 &v) const {
		return v.x != x || v.y != y || v.z != z;
	}

	/// Serialize this point to a binary data stream
	void serialize(Stream *stream) const {
		stream->writeElement<T>(x);
		stream->writeElement<T>(y);
		stream->writeElement<T>(z);
	}


	/// Return a readable string representation of this point
	std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << ", " << y << ", " << z << "]";
		return oss.str();
	}
};

template <typename T> inline TPoint3<T> operator*(T f, const TPoint3<T> &v) {
	return v*f;
}

template <typename T> inline T distance(const TPoint3<T> &p1, const TPoint3<T> &p2) {
	return (p1-p2).length();
}

template <typename T> inline T distanceSquared(const TPoint3<T> &p1, const TPoint3<T> &p2) {
	return (p1-p2).lengthSquared();
}

template <> inline TPoint3<int> TPoint3<int>::operator/(int s) const {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Point3i: Division by zero!");
#endif
	return TPoint3(x/s, y/s, z/s);
}

template <> inline TPoint3<int> &TPoint3<int>::operator/=(int s) {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Point3i: Division by zero!");
#endif
	x /= s;
	y /= s;
	z /= s;
	return *this;
}

/**
 * \headerfile mitsuba/core/point.h mitsuba/mitsuba.h
 * \brief Parameterizable four-dimensional point data structure
 * \ingroup libcore
 */
template <typename T> struct TPoint4 {
	typedef T           Scalar;
	typedef TVector4<T> VectorType;

	T x, y, z, w;

	/// Number of dimensions
	const static int dim = 4;

	/** \brief Construct a new point without initializing it.
	 *
	 * This construtor is useful when the point will either not
	 * be used at all (it might be part of a larger data structure)
	 * or initialized at a later point in time. Always make sure
	 * that one of the two is the case! Otherwise your program will do
	 * computations involving uninitialized memory, which will probably
	 * lead to a difficult-to-find bug.
	 */
#if !defined(MTS_DEBUG_UNINITIALIZED)
	TPoint4() { }
#else
	TPoint4() { x = y = z = w = std::numeric_limits<T>::quiet_NaN(); }
#endif

	/// Initialize the point with the specified X, Y and Z components
	TPoint4(T x, T y, T z, T w) : x(x), y(y), z(z), w(w) {  }

	/// Initialize the point with the components of another point
	template <typename T2> explicit TPoint4(const TPoint4<T2> &p)
		: x((T) p.x), y((T) p.y), z((T) p.z), w((T) p.w) { }

	/// Initialize the point with the components of a vector data structure
	template <typename T2> explicit TPoint4(const TVector4<T2> &v)
		: x((T) v.x), y((T) v.y), z((T) v.z), w((T) v.w) { }

	/// Initialize all components of the the point with the specified value
	explicit TPoint4(T val) : x(val), y(val), z(val), w(val) { }

	/// Unserialize a point from a binary data stream
	explicit TPoint4(Stream *stream) {
		x = stream->readElement<T>();
		y = stream->readElement<T>();
		z = stream->readElement<T>();
		w = stream->readElement<T>();
	}

	/// Add a vector to a point and return the result
	TPoint4 operator+(const TVector4<T> &v) const {
		return TPoint4(x + v.x, y + v.y, z + v.z, w + v.w);
	}

	/// Add two points and return the result (e.g. to compute a weighted position)
	TPoint4 operator+(const TPoint4 &p) const {
		return TPoint4(x + p.x, y + p.y, z + p.z, w + p.w);
	}

	/// Add a vector to this one (e.g. to compute a weighted position)
	TPoint4& operator+=(const TVector4<T> &v) {
		x += v.x; y += v.y; z += v.z; w += v.w;
		return *this;
	}

	/// Add a point to this one (e.g. to compute a weighted position)
	TPoint4& operator+=(const TPoint4 &p) {
		x += p.x; y += p.y; z += p.z; w += p.w;
		return *this;
	}

	/// Subtract a vector from this point
	TPoint4 operator-(const TVector4<T> &v) const {
		return TPoint4(x - v.x, y - v.y, z - v.z, w - v.w);
	}

	/// Subtract two points from each other and return the difference as a vector
	TVector4<T> operator-(const TPoint4 &p) const {
		return TVector4<T>(x - p.x, y - p.y, z - p.z, w - p.w);
	}

	/// Subtract a vector from this point
	TPoint4& operator-=(const TVector4<T> &v) {
		x -= v.x; y -= v.y; z -= v.z; w -= v.w;
		return *this;
	}

	/// Scale the point's coordinates by the given scalar and return the result
	TPoint4 operator*(T f) const {
		return TPoint4(x * f, y * f, z * f, w * f);
	}

	/// Scale the point's coordinates by the given scalar
	TPoint4 &operator*=(T f) {
		x *= f; y *= f; z *= f; w *= f;
		return *this;
	}

	/// Return a version of the point, which has been flipped along the origin
	TPoint4 operator-() const {
		return TPoint4(-x, -y, -z, -w);
	}

	/// Divide the point's coordinates by the given scalar and return the result
	TPoint4 operator/(T f) const {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Point4: Division by zero!");
#endif
		T recip = (T) 1 / f;
		return TPoint4(x * recip, y * recip, z * recip, w * recip);
	}

	/// Divide the point's coordinates by the given scalar
	TPoint4 &operator/=(T f) {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "Point4: Division by zero!");
#endif
		T recip = (T) 1 / f;
		x *= recip; y *= recip; z *= recip; w *= recip;
		return *this;
	}

	/// Index into the point's components
	T &operator[](int i) {
		return (&x)[i];
	}

	/// Index into the point's components (const version)
	T operator[](int i) const {
		return (&x)[i];
	}

	/// Return whether or not this point is identically zero
	bool isZero() const {
		return x == 0 && y == 0 && z == 0 && w == 0;
	}

	/// Equality test
	bool operator==(const TPoint4 &v) const {
		return (v.x == x && v.y == y && v.z == z && v.w == w);
	}

	/// Inequality test
	bool operator!=(const TPoint4 &v) const {
		return v.x != x || v.y != y || v.z != z || v.w != w;
	}

	/// Serialize this point to a binary data stream
	void serialize(Stream *stream) const {
		stream->writeElement<T>(x);
		stream->writeElement<T>(y);
		stream->writeElement<T>(z);
		stream->writeElement<T>(w);
	}

	/// Return a readable string representation of this point
	std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << ", " << y << ", " << z << ", " << w << "]";
		return oss.str();
	}
};

template <typename T> inline TPoint4<T> operator*(T f, const TPoint4<T> &v) {
	return v*f;
}

template <typename T> inline T distance(const TPoint4<T> &p1, const TPoint4<T> &p2) {
	return (p1-p2).length();
}

template <typename T> inline T distanceSquared(const TPoint4<T> &p1, const TPoint4<T> &p2) {
	return (p1-p2).lengthSquared();
}

template <> inline TPoint4<int> TPoint4<int>::operator/(int s) const {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Point4i: Division by zero!");
#endif
	return TPoint4(x/s, y/s, z/s, w/s);
}

template <> inline TPoint4<int> &TPoint4<int>::operator/=(int s) {
#ifdef MTS_DEBUG
	if (s == 0)
		SLog(EWarn, "Point4i: Division by zero!");
#endif

	x /= s;
	y /= s;
	z /= s;
	w /= s;
	return *this;
}

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_POINT_H_ */


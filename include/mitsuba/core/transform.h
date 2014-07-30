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
#if !defined(__MITSUBA_CORE_TRANSFORM_H_)
#define __MITSUBA_CORE_TRANSFORM_H_

#include <mitsuba/core/matrix.h>
#include <mitsuba/core/ray.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Encapsulates a 4x4 linear transformation and its inverse
 * \ingroup libcore
 * \ingroup libpython
 */
struct MTS_EXPORT_CORE Transform {
public:
	/// Create an identity transformation
	Transform() {
		m_transform.setIdentity();
		m_invTransform.setIdentity();
	}

	/// Unserialize a transformation from a stream
	inline Transform(Stream *stream) {
		m_transform = Matrix4x4(stream);
		m_invTransform = Matrix4x4(stream);
	}

	/** \brief Create a transform from the given matrix
	 * and calculate the inverse
	 */
	Transform(const Matrix4x4 &trafo)
		: m_transform(trafo) {
		bool success = m_transform.invert(m_invTransform);
		if (!success)
			SLog(EError, "Unable to invert singular matrix %s", trafo.toString().c_str());
	}

	/// Create a transform from the given matrices
	Transform(const Matrix4x4 &trafo, const Matrix4x4 &invTrafo)
		: m_transform(trafo), m_invTransform(invTrafo) {
	}

	/// Return the inverse transform
	Transform inverse() const {
		return Transform(m_invTransform, m_transform);
	}

	/// Matrix-matrix multiplication
	Transform operator*(const Transform &t) const;

	/// Return the determinant of the upper left 3x3 submatrix
	inline Float det3x3() const {
		return m_transform.det3x3();
	}

	/// Test for a scale component
	inline bool hasScale() const {
		for (int i=0; i<3; ++i) {
			for (int j=i; j<3; ++j) {
				Float sum = 0;
				for (int k=0; k<3; ++k)
					sum += m_transform.m[i][k] * m_transform.m[j][k];

				if (i == j && std::abs(sum-1) > 1e-3f)
					return true;
				else if (i != j && std::abs(sum) > 1e-3f)
					return true;
			}
		}
		return false;
	}

	/// Test if this is the identity transformation
	inline bool isIdentity() const {
		return m_transform.isIdentity();
	}

	/// Test if this is the zero transformation
	inline bool isZero() const {
		return m_transform.isZero();
	}

	/**
	 * \brief Matrix-vector multiplication for points in 3d space
	 *
	 * \remark In the Python bindings, this is function implemented as
	 * the multiplication operator (\c __mul__).
	 */
	inline Point operator()(const Point &p) const {
		Float x = m_transform.m[0][0] * p.x + m_transform.m[0][1] * p.y
		        + m_transform.m[0][2] * p.z + m_transform.m[0][3];
		Float y = m_transform.m[1][0] * p.x + m_transform.m[1][1] * p.y
		        + m_transform.m[1][2] * p.z + m_transform.m[1][3];
		Float z = m_transform.m[2][0] * p.x + m_transform.m[2][1] * p.y
		        + m_transform.m[2][2] * p.z + m_transform.m[2][3];
		Float w = m_transform.m[3][0] * p.x + m_transform.m[3][1] * p.y
		        + m_transform.m[3][2] * p.z + m_transform.m[3][3];
#ifdef MTS_DEBUG
		if (w == 0)
			SLog(EWarn, "w==0 in Transform::operator(Point &)");
#endif
		if (w == 1.0f)
			return Point(x, y, z);
		else
			return Point(x, y, z) / w;
	}

	/// Transform a point by an affine / non-projective matrix
	inline Point transformAffine(const Point &p) const {
		Float x = m_transform.m[0][0] * p.x + m_transform.m[0][1] * p.y
		        + m_transform.m[0][2] * p.z + m_transform.m[0][3];
		Float y = m_transform.m[1][0] * p.x + m_transform.m[1][1] * p.y
		        + m_transform.m[1][2] * p.z + m_transform.m[1][3];
		Float z = m_transform.m[2][0] * p.x + m_transform.m[2][1] * p.y
		        + m_transform.m[2][2] * p.z + m_transform.m[2][3];
		return Point(x,y,z);
	}

	/// Transform a point by a affine / non-projective matrix (no temporaries)
	inline void transformAffine(const Point &p, Point &dest) const {
		dest.x = m_transform.m[0][0] * p.x + m_transform.m[0][1] * p.y
		       + m_transform.m[0][2] * p.z + m_transform.m[0][3];
		dest.y = m_transform.m[1][0] * p.x + m_transform.m[1][1] * p.y
		       + m_transform.m[1][2] * p.z + m_transform.m[1][3];
		dest.z = m_transform.m[2][0] * p.x + m_transform.m[2][1] * p.y
		       + m_transform.m[2][2] * p.z + m_transform.m[2][3];
	}

	/**
	 * \brief Matrix-vector multiplication for points in 3d space (no temporaries)
	 * \remark This function is not available in the Python bindings
	 */
    inline void operator()(const Point &p, Point &dest) const {
		dest.x = m_transform.m[0][0] * p.x + m_transform.m[0][1] * p.y
		       + m_transform.m[0][2] * p.z + m_transform.m[0][3];
		dest.y = m_transform.m[1][0] * p.x + m_transform.m[1][1] * p.y
		       + m_transform.m[1][2] * p.z + m_transform.m[1][3];
		dest.z = m_transform.m[2][0] * p.x + m_transform.m[2][1] * p.y
		       + m_transform.m[2][2] * p.z + m_transform.m[2][3];
		Float w = m_transform.m[3][0] * p.x + m_transform.m[3][1] * p.y
		        + m_transform.m[3][2] * p.z + m_transform.m[3][3];

#ifdef MTS_DEBUG
		if (w == 0)
			SLog(EWarn, "w==0 in Transform::operator(Point &, Point &)");
#endif
		if (w != 1.0f)
			dest /= w;
	}

	/**
	 * \brief Matrix-vector multiplication for vectors in 3d space
	 * \remark In the Python bindings, this is function implemented as
	 * the multiplication operator (\c __mul__).
	 */
    inline Vector operator()(const Vector &v) const {
		Float x = m_transform.m[0][0] * v.x + m_transform.m[0][1] * v.y
		        + m_transform.m[0][2] * v.z;
		Float y = m_transform.m[1][0] * v.x + m_transform.m[1][1] * v.y
		        + m_transform.m[1][2] * v.z;
		Float z = m_transform.m[2][0] * v.x + m_transform.m[2][1] * v.y
		        + m_transform.m[2][2] * v.z;
		return Vector(x, y, z);
	}

	/**
	 * \brief Matrix-vector multiplication for vectors in 3d space (no temporaries)
	 * \remark This function is not available in the Python bindings
	 */
    inline void operator()(const Vector &v, Vector &dest) const {
		dest.x = m_transform.m[0][0] * v.x + m_transform.m[0][1] * v.y
		       + m_transform.m[0][2] * v.z;
		dest.y = m_transform.m[1][0] * v.x + m_transform.m[1][1] * v.y
		       + m_transform.m[1][2] * v.z;
		dest.z = m_transform.m[2][0] * v.x + m_transform.m[2][1] * v.y
		       + m_transform.m[2][2] * v.z;
	}

	/**
	 * \brief Matrix-normal multiplication
	 * \remark In the Python bindings, this is function implemented as
	 * the multiplication operator (\c __mul__).
	 */
    inline Normal operator()(const Normal &v) const {
		Float x = m_invTransform.m[0][0] * v.x + m_invTransform.m[1][0] * v.y
		        + m_invTransform.m[2][0] * v.z;
		Float y = m_invTransform.m[0][1] * v.x + m_invTransform.m[1][1] * v.y
		        + m_invTransform.m[2][1] * v.z;
		Float z = m_invTransform.m[0][2] * v.x + m_invTransform.m[1][2] * v.y
		        + m_invTransform.m[2][2] * v.z;
		return Normal(x, y, z);
	}

	/**
	 * \brief Matrix-normal multiplication (no temporaries)
	 * \remark This function is not available in the Python bindings
	 */
    inline void operator()(const Normal &v, Normal &dest) const {
		dest.x = m_invTransform.m[0][0] * v.x + m_invTransform.m[1][0] * v.y
		       + m_invTransform.m[2][0] * v.z;
		dest.y = m_invTransform.m[0][1] * v.x + m_invTransform.m[1][1] * v.y
		       + m_invTransform.m[2][1] * v.z;
		dest.z = m_invTransform.m[0][2] * v.x + m_invTransform.m[1][2] * v.y
		       + m_invTransform.m[2][2] * v.z;
	}

	/**
	 * \brief 4D matrix-vector multiplication
	 * \remark In the Python bindings, this is function implemented as
	 * the multiplication operator (\c __mul__).
	 */
	inline Vector4 operator()(const Vector4 &v) const {
		Float x = m_transform.m[0][0] * v.x + m_transform.m[0][1] * v.y
		        + m_transform.m[0][2] * v.z + m_transform.m[0][3] * v.w;
		Float y = m_transform.m[1][0] * v.x + m_transform.m[1][1] * v.y
		        + m_transform.m[1][2] * v.z + m_transform.m[1][3] * v.w;
		Float z = m_transform.m[2][0] * v.x + m_transform.m[2][1] * v.y
		        + m_transform.m[2][2] * v.z + m_transform.m[2][3] * v.w;
		Float w = m_transform.m[3][0] * v.x + m_transform.m[3][1] * v.y
		        + m_transform.m[3][2] * v.z + m_transform.m[3][3] * v.w;
		return Vector4(x,y,z,w);
	}

	/**
	 * \brief 4D matrix-vector multiplication (no temporaries)
	 * \remark This function is not available in the Python bindings
	 */
	inline void operator()(const Vector4 &v, Vector4 &dest) const {
		dest.x = m_transform.m[0][0] * v.x + m_transform.m[0][1] * v.y
		       + m_transform.m[0][2] * v.z + m_transform.m[0][3] * v.w;
		dest.y = m_transform.m[1][0] * v.x + m_transform.m[1][1] * v.y
		       + m_transform.m[1][2] * v.z + m_transform.m[1][3] * v.w;
		dest.z = m_transform.m[2][0] * v.x + m_transform.m[2][1] * v.y
		       + m_transform.m[2][2] * v.z + m_transform.m[2][3] * v.w;
		dest.w = m_transform.m[3][0] * v.x + m_transform.m[3][1] * v.y
		       + m_transform.m[3][2] * v.z + m_transform.m[3][3] * v.w;
	}

	/**
	 * \brief Transform a ray. Assumes that there is no scaling (no temporaries)
	 * \remark This function is not available in the Python bindings
	 */
	inline void operator()(const Ray &a, Ray &b) const {
		b.mint = a.mint;
		b.maxt = a.maxt;
		operator()(a.o, b.o);
		operator()(a.d, b.d);
#ifdef MTS_DEBUG_FP
		bool state = disableFPExceptions();
#endif
		/* Re-compute the reciprocal */
		b.dRcp.x = 1.0f / b.d.x;
		b.dRcp.y = 1.0f / b.d.y;
		b.dRcp.z = 1.0f / b.d.z;
#ifdef MTS_DEBUG_FP
		restoreFPExceptions(state);
#endif
		b.time = a.time;
	}

	/**
	 * \brief Transform a ray
	 * \remark In the Python bindings, this is function implemented as
	 * the multiplication operator (\c __mul__).
	 */
	inline Ray operator()(const Ray &ray) const {
		Ray result;
		operator()(ray, result);
		return result;
	}

	/// Transform a ray by an affine / non-projective matrix (no temporaries)
	inline void transformAffine(const Ray &a, Ray &b) const {
		b.mint = a.mint;
		b.maxt = a.maxt;
		transformAffine(a.o, b.o);
		operator()(a.d, b.d);
#ifdef MTS_DEBUG_FP
		bool state = disableFPExceptions();
#endif
		/* Re-compute the reciprocal */
		b.dRcp.x = 1.0f / b.d.x;
		b.dRcp.y = 1.0f / b.d.y;
		b.dRcp.z = 1.0f / b.d.z;
#ifdef MTS_DEBUG_FP
		restoreFPExceptions(state);
#endif
	}

	/// Transform a ray by an affine / non-projective matrix (no temporaries)
	inline Ray transformAffine(const Ray &ray) const {
		Ray result;
		transformAffine(ray, result);
		return result;
	}

	/// Return the underlying matrix
	inline const Matrix4x4 &getMatrix() const { return m_transform; }

	/// Return the underlying inverse matrix (const version)
	inline const Matrix4x4 &getInverseMatrix() const { return m_invTransform; }

	/// Create a translation transformation
	static Transform translate(const Vector &v);

	/// Create a rotation transformation around an arbitrary axis. The angle is specified in degrees
	static Transform rotate(const Vector &axis, Float angle);

	/// Create a scale transformation
	static Transform scale(const Vector &v);

	/** \brief Create a perspective transformation.
	 *   (Maps [near, far] to [0, 1])
	 * \param fov Field of view in degrees
	 * \param clipNear Near clipping plane
	 * \param clipFar Far clipping plane
	 */
	static Transform perspective(Float fov, Float clipNear, Float clipFar);

	/** \brief Create a perspective transformation for OpenGL.
	 *   (Maps [-near, -far] to [-1, 1])
	 * \param fov Field of view in degrees
	 * \param clipNear Near clipping plane distance
	 * \param clipFar Far clipping plane distance
	 */
	static Transform glPerspective(Float fov, Float clipNear, Float clipFar);

	/** \brief Create a perspective transformation for OpenGL.
	 * \param left Left clipping plane coordinate
	 * \param right Right clipping plane coordinate
	 * \param top Top clipping plane coordinate
	 * \param bottom Bottom clipping plane coordinate
	 * \param nearVal Near clipping plane distance
	 * \param farVal Far clipping plane distance
	 */
	static Transform glFrustum(Float left, Float right, Float bottom, Float top, Float nearVal, Float farVal);

	/** \brief Create an orthographic transformation, which maps Z to [0,1]
	 * and leaves the X and Y coordinates untouched.
	 * \param clipNear Near clipping plane
	 * \param clipFar Far clipping plane
	 */
	static Transform orthographic(Float clipNear, Float clipFar);

	/** \brief Create an orthographic transformation for OpenGL
	 * \param clipNear Near clipping plane
	 * \param clipFar Far clipping plane
	 */
	static Transform glOrthographic(Float clipNear, Float clipFar);

	/** \brief Create an orthographic transformation for OpenGL
	 *
	 * Slightly extended variant which also handles non-unity clipping
	 * planes on the X and Y axes and matches the 'glOrtho' spec.
	 */
	static Transform glOrthographic(Float clipLeft, Float clipRight,
		Float clipBottom, Float clipTop, Float clipNear, Float clipFar);

	/** \brief Create a look-at camera transformation
	 * \param p Camera position
	 * \param t Target vector
	 * \param u Up vector
	 */
	static Transform lookAt(const Point &p, const Point &t, const Vector &u);

	/** \brief Create an orthogonal transformation that takes
	 * the standard to the supplied frame
	 */
	static Transform fromFrame(const Frame &frame);

	/// Serialize a transformation to a stream
	inline void serialize(Stream *stream) const {
		m_transform.serialize(stream);
		m_invTransform.serialize(stream);
	}

	/// Equality comparison operator
	inline bool operator==(const Transform &trafo) const {
		return m_transform == trafo.m_transform;
	}

	/// Inequality comparison operator
	inline bool operator!=(const Transform &trafo) const {
		return m_transform != trafo.m_transform;
	}

	/// Return a string representation
	std::string toString() const;
private:
	Matrix4x4 m_transform;
	Matrix4x4 m_invTransform;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_TRANSFORM_H_ */

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

#if !defined(__TRANSFORM_H)
#define __TRANSFORM_H

#include <mitsuba/core/ray.h>
#include <mitsuba/core/normal.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief 4x4 matrix data type
 * (thin reference counted array wrapper).
 */
class MTS_EXPORT_CORE Matrix4x4 : public Object {
public:
	/// Initialize with the identity matrix
	Matrix4x4();

	/// Unserialize a matrix from a stream
	Matrix4x4(Stream *stream);

	/// Copy-constructor
	Matrix4x4(const Matrix4x4 *mat);
	
	/// Copy-constructor
	Matrix4x4(Float m[4][4]);

	/// Initialize with the given values
	Matrix4x4(
		Float a00, Float a01, Float a02, Float a03,
		Float a10, Float a11, Float a12, Float a13,
		Float a20, Float a21, Float a22, Float a23,
		Float a30, Float a31, Float a32, Float a33
	);

	/// Return the determinant of the upper left 3x3 sub-matrix
	Float det3x3() const;

	/// Perform a symmetric 4x4 eigendecomposition into Q and D.
	void symmEigenDecomp(Matrix4x4 *Q, Vector4 &d);

	/// Return the inverse of this matrix
	ref<Matrix4x4> inverse() const;

	/// Return the transpose of this matrix
	ref<Matrix4x4> transpose() const;

	/// Serialize the matrix to a stream
	void serialize(Stream *stream) const;

	/// Return a string representation
	std::string toString() const;

	Float m[4][4];
	
	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Matrix4x4() { }
};
	
struct Frame;

/**
 * \brief Transform class encapsulating a linear affine
 * transformation (using homogenous coordinates) and its inverse
 */
class MTS_EXPORT_CORE Transform {
public:
	/// Create an identity transformation
	Transform();
	
	/// Unserialize a transformation from a stream
	inline Transform(Stream *stream) { 
		m_transform = new Matrix4x4(stream);
		m_invTransform = new Matrix4x4(stream);
	}

	/** \brief Create a transform from the given matrix
	 * and calculate the inverse
	 */
	Transform(const Matrix4x4 *trafo);
	
	/// Create a transform from the given matrices
	Transform(const Matrix4x4 *trafo, const Matrix4x4 *invTrafo);

	/// Return the inverse transform
	Transform inverse() const;

	/// Matrix-matrix multiplication
	Transform operator*(const Transform &t) const;

	/// Return the determinant of the upper left 3x3 submatrix
	inline Float det3x3() const {
		return m_transform->det3x3();
	}

	/// Test for a scale component
	inline bool hasScale() const {
		for (int i=0; i<3; ++i) {
			for (int j=i; j<3; ++j) {
				Float sum = 0;
				for (int k=0; k<3; ++k)
					sum += m_transform->m[i][k] * m_transform->m[j][k];

				if (i == j && std::abs(sum-1) > Epsilon)
					return true;
				else if (i != j && std::abs(sum) > Epsilon)
					return true;
			}
		}
		return false;
	}

	/// Test if this is the identity matrix
	inline bool isIdentity() const {
		for (int i=0; i<4; ++i) {
			for (int j=0; j<4; ++j) {
				if (m_transform->m[i][j] != ((i==j) ? 1 : 0))
					return false;
			}
		}
		return true;
	}

	/// Matrix-vector multiplication for points in 3d space
	inline Point operator()(const Point &p) const {
		Float x = m_transform->m[0][0] * p.x + m_transform->m[0][1] * p.y
				+ m_transform->m[0][2] * p.z + m_transform->m[0][3];
		Float y = m_transform->m[1][0] * p.x + m_transform->m[1][1] * p.y
				+ m_transform->m[1][2] * p.z + m_transform->m[1][3];
		Float z = m_transform->m[2][0] * p.x + m_transform->m[2][1] * p.y
				+ m_transform->m[2][2] * p.z + m_transform->m[2][3];
		Float w = m_transform->m[3][0] * p.x + m_transform->m[3][1] * p.y
				+ m_transform->m[3][2] * p.z + m_transform->m[3][3];
#ifdef MTS_DEBUG
		if (w == 0)
			SLog(EWarn, "w==0 in Transform::operator(Point &)");
#endif
		if (w == 1.0f)
			return Point(x, y, z);
		else
			return Point(x, y, z) / w;
	}

	/// Transform a point by a non-projective matrix
	inline Point transformBasic(const Point &p) const {
		Float x = m_transform->m[0][0] * p.x + m_transform->m[0][1] * p.y
				+ m_transform->m[0][2] * p.z + m_transform->m[0][3];
		Float y = m_transform->m[1][0] * p.x + m_transform->m[1][1] * p.y
				+ m_transform->m[1][2] * p.z + m_transform->m[1][3];
		Float z = m_transform->m[2][0] * p.x + m_transform->m[2][1] * p.y
				+ m_transform->m[2][2] * p.z + m_transform->m[2][3];
		return Point(x,y,z);
	}

	/// Matrix-vector multiplication for points in 3d space (no temporaries)
    inline void operator()(const Point &p, Point &dest) const {
		dest.x = m_transform->m[0][0] * p.x + m_transform->m[0][1] * p.y
			   + m_transform->m[0][2] * p.z + m_transform->m[0][3];
		dest.y = m_transform->m[1][0] * p.x + m_transform->m[1][1] * p.y
			   + m_transform->m[1][2] * p.z + m_transform->m[1][3];
		dest.z = m_transform->m[2][0] * p.x + m_transform->m[2][1] * p.y
			   + m_transform->m[2][2] * p.z + m_transform->m[2][3];
		Float w = m_transform->m[3][0] * p.x + m_transform->m[3][1] * p.y
				+ m_transform->m[3][2] * p.z + m_transform->m[3][3];

#ifdef MTS_DEBUG
		if (w == 0)
			SLog(EWarn, "w==0 in Transform::operator(Point &, Point &)");
#endif
		if (w != 1.0f)
			dest /= w;
	}

	/// Matrix-vector multiplication for vectors in 3d space
    inline Vector operator()(const Vector &v) const {
		Float x = m_transform->m[0][0] * v.x + m_transform->m[0][1] * v.y
				+ m_transform->m[0][2] * v.z;
		Float y = m_transform->m[1][0] * v.x + m_transform->m[1][1] * v.y
				+ m_transform->m[1][2] * v.z;
		Float z = m_transform->m[2][0] * v.x + m_transform->m[2][1] * v.y
				+ m_transform->m[2][2] * v.z;
		return Vector(x, y, z);
	}

	/// Matrix-vector multiplication for vectors in 3d space (no temporaries)
    inline void operator()(const Vector &v, Vector &dest) const {
		dest.x = m_transform->m[0][0] * v.x + m_transform->m[0][1] * v.y
			   + m_transform->m[0][2] * v.z;
		dest.y = m_transform->m[1][0] * v.x + m_transform->m[1][1] * v.y
			   + m_transform->m[1][2] * v.z;
		dest.z = m_transform->m[2][0] * v.x + m_transform->m[2][1] * v.y
			   + m_transform->m[2][2] * v.z;
	}

	/// Matrix-normal multiplication
    inline Normal operator()(const Normal &v) const {
		Float x = m_invTransform->m[0][0] * v.x + m_invTransform->m[1][0] * v.y
				+ m_invTransform->m[2][0] * v.z;
		Float y = m_invTransform->m[0][1] * v.x + m_invTransform->m[1][1] * v.y
				+ m_invTransform->m[2][1] * v.z;
		Float z = m_invTransform->m[0][2] * v.x + m_invTransform->m[1][2] * v.y
				+ m_invTransform->m[2][2] * v.z;
		return Normal(x, y, z);
	}

	/// Matrix-normal multiplication (no temporaries)
    inline void operator()(const Normal &v, Normal &dest) const {
		dest.x = m_invTransform->m[0][0] * v.x + m_invTransform->m[1][0] * v.y
			   + m_invTransform->m[2][0] * v.z;
		dest.y = m_invTransform->m[0][1] * v.x + m_invTransform->m[1][1] * v.y
			   + m_invTransform->m[2][1] * v.z;
		dest.z = m_invTransform->m[0][2] * v.x + m_invTransform->m[1][2] * v.y
			   + m_invTransform->m[2][2] * v.z;
	}

	/// Transform a ray. Assumes that there is no scaling
	inline void operator()(const Ray &a, Ray &b) const {
		b.mint = a.mint;
		b.maxt = a.maxt;
		operator()(a.o, b.o);
		operator()(a.d, b.d);
#ifdef MTS_DEBUG_FP
		disable_fpexcept();
#endif
		/* Re-compute the reciprocal */
		b.dRcp.x = 1.0f / b.d.x;
		b.dRcp.y = 1.0f / b.d.y;
		b.dRcp.z = 1.0f / b.d.z;
#ifdef MTS_DEBUG_FP
		enable_fpexcept();
#endif
	}

	/// Return the underlying matrix
	inline const Matrix4x4 *getMatrix() const { return m_transform.get(); }

	/// Create a translation transformation
	static Transform translate(const Vector &v);

	/// Create a rotation transformation around an arbitrary axis. The angle is specified in degrees
	static Transform rotate(const Vector &axis, Float angle);

	/// Create a scale transformation
	static Transform scale(const Vector &v);
	
	/** \brief Create a perspective transformation.
	 *   (Maps [near, far] to [0, 1])
	 * @param fov Field of view in degrees
	 * @param clipNear Near clipping plane
	 * @param clipFar Far clipping plane
	 */
	static Transform perspective(Float fov, Float clipNear, Float clipFar);
	
	/** \brief Create a perspective transformation for OpenGL.
	 *   (Maps [near, far] to [-1, 1])
	 * @param fov Field of view in degrees
	 * @param clipNear Near clipping plane distance
	 * @param clipFar Far clipping plane distance
	 */
	static Transform glPerspective(Float fov, Float clipNear, Float clipFar);

	/** \brief Create a perspective transformation for OpenGL.
	 * @param left Left clipping plane coordinate
	 * @param right Right clipping plane coordinate
	 * @param top Top clipping plane coordinate
	 * @param bottom Bottom clipping plane coordinate
	 * @param clipNear Near clipping plane distance
	 * @param clipFar Far clipping plane distance
	 */
	static Transform glFrustum(Float left, Float right, Float bottom, Float top, Float nearVal, Float farVal);

	/** \brief Create an orthographic transformation, which maps Z to [0,1]
	 * and leaves the X and Y coordinates untouched.
	 * @param clipNear Near clipping plane
	 * @param clipFar Far clipping plane
	 */
	static Transform orthographic(Float clipNear, Float clipFar);
	
	/** \brief Create an orthographic transformation for OpenGL
	 * @param clipNear Near clipping plane
	 * @param clipFar Far clipping plane
	 */
	static Transform glOrthographic(Float clipNear, Float clipFar);

	/** \brief Create a look-at camera transformation
	 * @param p Camera position
	 * @param t Target vector
	 * @param u Up vector
	 */
	static Transform lookAt(const Point &p, const Point &t, const Vector &u);

	/** \brief Create an orthogonal transformation that takes
	 * the standard to the supplied frame
	 */
	static Transform fromFrame(const Frame &frame);

	/// Serialize a transformation to a stream
	inline void serialize(Stream *stream) const {
		m_transform->serialize(stream);
		m_invTransform->serialize(stream);
	}

	/// Return a string representation
	std::string toString() const;
private:
	ref<const Matrix4x4> m_transform;
	ref<const Matrix4x4> m_invTransform;
};

/// Matrix-matrix multiplication
inline ref<Matrix4x4> operator*(const ref<const Matrix4x4> &m1, const ref<const Matrix4x4> &m2) {
	ref<Matrix4x4> retval = new Matrix4x4();
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			retval->m[i][j] = m1->m[i][0] * m2->m[0][j] +
							  m1->m[i][1] * m2->m[1][j] +
							  m1->m[i][2] * m2->m[2][j] +
							  m1->m[i][3] * m2->m[3][j];
	return retval;
}

MTS_NAMESPACE_END

#endif /* __TRANSFORM_H */

/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#include <mitsuba/core/transform.h>
#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN

// -----------------------------------------------------------------------
//  Linear transformation class
// -----------------------------------------------------------------------

Transform Transform::operator*(const Transform &t) const {
	return Transform(m_transform * t.m_transform,
		t.m_invTransform * m_invTransform);
}

Transform Transform::translate(const Vector &v) {
	Matrix4x4 trafo(
		1, 0, 0, v.x,
		0, 1, 0, v.y,
		0, 0, 1, v.z,
		0, 0, 0, 1
	);
	Matrix4x4 invTrafo(
		1, 0, 0, -v.x,
		0, 1, 0, -v.y,
		0, 0, 1, -v.z,
		0, 0, 0, 1
	);
	return Transform(trafo, invTrafo);
}

Transform Transform::scale(const Vector &v) {
	Matrix4x4 trafo(
		v.x, 0,   0,   0,
		0,   v.y, 0,   0,
		0,   0,   v.z, 0,
		0,   0,   0,   1
	);
	Matrix4x4 invTrafo(
		1.0f/v.x, 0,        0,        0,
		0,        1.0f/v.y, 0,        0,
		0,        0,        1.0f/v.z, 0,
		0,        0,        0,        1
	);
	return Transform(trafo, invTrafo);
}

Transform Transform::rotate(const Vector &axis, Float angle) {
	/* Make sure that the axis is normalized */
	Vector naxis = normalize(axis);

	Float rotAngle = degToRad(angle);
	Float sinAngle = std::sin(rotAngle);
	Float cosAngle = std::cos(rotAngle);
	Matrix4x4 result;

	result.m[0][0] = naxis.x * naxis.x + (1.f - naxis.x * naxis.x) * cosAngle;
	result.m[0][1] = naxis.x * naxis.y * (1.f - cosAngle) - naxis.z * sinAngle;
	result.m[0][2] = naxis.x * naxis.z * (1.f - cosAngle) + naxis.y * sinAngle;
	result.m[0][3] = 0;

	result.m[1][0] = naxis.x * naxis.y * (1.f - cosAngle) + naxis.z * sinAngle;
	result.m[1][1] = naxis.y * naxis.y + (1.f - naxis.y * naxis.y) * cosAngle;
	result.m[1][2] = naxis.y * naxis.z * (1.f - cosAngle) - naxis.x * sinAngle;
	result.m[1][3] = 0;

	result.m[2][0] = naxis.x * naxis.z * (1.f - cosAngle) - naxis.y * sinAngle;
	result.m[2][1] = naxis.y * naxis.z * (1.f - cosAngle) + naxis.x * sinAngle;
	result.m[2][2] = naxis.z * naxis.z + (1.f - naxis.z * naxis.z) * cosAngle;
	result.m[2][3] = 0;

	result.m[3][0] = 0;
	result.m[3][1] = 0;
	result.m[3][2] = 0;
	result.m[3][3] = 1;

	/* The matrix is orthonormal */
	Matrix4x4 transp;
	result.transpose(transp);
	return Transform(result, transp);
}

Transform Transform::perspective(Float fov, Float clipNear, Float clipFar) {
	/* Project vectors in camera space onto a plane at z=1:
	 *
	 *  xProj = x / z
	 *  yProj = y / z
	 *  zProj = (far * (z - near)) / (z * (far-near))
	 *  
	 *  Camera-space depths are not mapped linearly!
	 */

	Float recip = 1.0f / (clipFar - clipNear);

	Matrix4x4 trafo(
		1,   0,   0,   0,
		0,   1,   0,   0,
		0,   0,   clipFar * recip, -clipNear * clipFar * recip,
		0,   0,   1,   0
	);

	/* Perform a scale so that the field of view is mapped
	 * to the interval [-1, 1] */
	Float cot = 1.0f / std::tan(degToRad(fov / 2.0f));

	return Transform::scale(Vector(cot, cot, 1.0f)) * Transform(trafo);
}

Transform Transform::glPerspective(Float fov, Float clipNear, Float clipFar) {
	Float recip = 1.0f / (clipNear - clipFar);
	Float cot = 1.0f / std::tan(degToRad(fov / 2.0f));

	Matrix4x4 trafo(
		cot,   0,     0,   0,
		0,     cot,   0,   0,
		0,     0,     (clipFar + clipNear) * recip,  2 * clipFar * clipNear * recip,
		0,     0,     -1,   0
	);

	return Transform(trafo);
}

Transform Transform::glFrustum(Float left, Float right, Float bottom, Float top, Float nearVal, Float farVal) {
	Float invFMN = 1 / (farVal-nearVal);
	Float invTMB = 1 / (top-bottom);
	Float invRML = 1 / (right-left);

	Matrix4x4 trafo(
		2*nearVal*invRML, 0, (right+left)*invRML, 0,
		0, 2*nearVal*invTMB, (top+bottom)*invTMB, 0,
		0, 0, -(farVal + nearVal) * invFMN, -2*farVal*nearVal*invFMN,
		0, 0, -1, 0
	);

	return Transform(trafo);
}

Transform Transform::orthographic(Float clipNear, Float clipFar) {
	return scale(Vector(1.0f, 1.0f, 1.0f / (clipFar - clipNear))) *
		   translate(Vector(0.0f, 0.0f, -clipNear));
}

Transform Transform::glOrthographic(Float clipNear, Float clipFar) {
	Float a = -2.0f / (clipFar - clipNear),
	      b = (clipFar + clipNear) / (clipFar - clipNear);

	Matrix4x4 trafo(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, a, b,
		0, 0, 0, 1
	);
	return Transform(trafo);
}


Transform Transform::lookAt(const Point &p, const Point &t, const Vector &up) {
	Matrix4x4 result;

	Vector dirct = normalize(t-p);
	Vector right = normalize(cross(dirct, up));

	/* Generate a new, orthogonalized up vector */
	Vector newUp = cross(right, dirct);

	/* Store as columns */
	result.m[0][0] = right.x; result.m[1][0] = right.y; result.m[2][0] = right.z; result.m[3][0] = 0;
	result.m[0][1] = newUp.x; result.m[1][1] = newUp.y; result.m[2][1] = newUp.z; result.m[3][1] = 0;
	result.m[0][2] = dirct.x; result.m[1][2] = dirct.y; result.m[2][2] = dirct.z; result.m[3][2] = 0;
	result.m[0][3] = p.x;     result.m[1][3] = p.y;     result.m[2][3] = p.z;     result.m[3][3] = 1;

	return Transform(result);
}

Transform Transform::fromFrame(const Frame &frame) {
	Matrix4x4 result(
		frame.s.x, frame.t.x, frame.n.x, 0,
		frame.s.y, frame.t.y, frame.n.y, 0, 
		frame.s.z, frame.t.z, frame.n.z, 0,
		0, 0, 0, 1
	);
	/* The matrix is orthonormal */
	Matrix4x4 transp;
	result.transpose(transp);
	return Transform(result, transp);
}

std::string Transform::toString() const {
	return m_transform.toString();
}

/**
 * \brief Tridiagonalizes a matrix using a Householder transformation
 * and returns the result.
 *
 * Based on "Geometric Tools" by David Eberly.
 */
static inline void tred3(Matrix3x3 &m, Float *diag, Float *subd) {
	Float m00 = m(0, 0), m01 = m(0, 1), m02 = m(0, 2),
	      m11 = m(1, 1), m12 = m(1, 2), m22 = m(2, 2);

	diag[0] = m00;
	subd[2] = 0;
	if (std::abs(m02) > std::numeric_limits<Float>::epsilon()) {
		Float length = std::sqrt(m01*m01 + m02*m02),
			  invLength = 1 / length;
		m01 *= invLength;
		m02 *= invLength;
		Float q = 2*m01*m12 + m02*(m22 - m11);
		diag[1] = m11 + m02*q;
		diag[2] = m22 - m02*q;
		subd[0] = length;
		subd[1] = m12 - m01*q;
		m = Matrix3x3(
			1, 0, 0,
			0, m01, m02,
			0, m02, -m01);
	} else {
		/* The matrix is already tridiagonal,
		   return an identity transformation matrix */
		diag[1] = m11;
		diag[2] = m22;
		subd[0] = m01;
		subd[1] = m12;
		m = Matrix3x3(
			1, 0, 0,
			0, 1, 0,
			0, 0, 1);
	}
}

/**
 * Compute the eigenvalues and eigenvectors of a 3x3 symmetric tridiagonal
 * matrix using the QL algorithm with implicit shifts
 *
 * Based on "Geometric Tools" by David Eberly and Numerical Recipes by
 * Teukolsky, et al. Originally, this code goes back to the Handbook 
 * for Automatic Computation by Wilkionson and Reinsch, as well as
 * the corresponding EISPACK routine.
 */

static inline bool ql3(Matrix3x3 &m, Float *diag, Float *subd) {
	const int maxIter = 32;

	for (int i = 0; i < 3; ++i) {
		int k;
		for (k = 0; k < maxIter; ++k) {
			int j;
			for (j = i; j < 2; ++j) {
				Float tmp = std::abs(diag[j]) + std::abs(diag[j+1]);
				if (std::abs(subd[j]) + tmp == tmp) 
					break;
			}
			if (j == i)
				break;

			Float value0 = (diag[i + 1] - diag[i])/(2*subd[i]);
			Float value1 = std::sqrt(value0*value0 + 1);
			value0 = diag[j] - diag[i] + subd[i] / 
				((value0 < 0) ? (value0 - value1) : (value0 + value1));

			Float sn = 1, cs = 1, value2 = 0;
			for (int l = j - 1; l >= i; --l) {
				Float value3 = sn*subd[l], value4 = cs*subd[l];
				if (std::abs(value3) >= std::abs(value0)) {
					cs = value0 / value3;
					value1 = std::sqrt(cs*cs + 1);
					subd[l + 1] = value3*value1;
					sn = 1.0f / value1;
					cs *= sn;
				} else {
					sn = value3 / value0;
					value1 = std::sqrt(sn*sn + 1);
					subd[l + 1] = value0*value1;
					cs = 1.0f / value1;
					sn *= cs;
				}
				value0 = diag[l + 1] - value2;
				value1 = (diag[l] - value0)*sn + 2*value4*cs;
				value2 = sn*value1;
				diag[l + 1] = value0 + value2;
				value0 = cs*value1 - value4;

				for (int t = 0; t < 3; ++t) {
					value3 = m(t, l + 1);
					m(t, l + 1) = sn*m(t, l) + cs*value3;
					m(t, l) = cs*m(t, l) - sn*value3;
				}
			}
			diag[i] -= value2;
			subd[i] = value0;
			subd[j] = 0;
		}
		if (k == maxIter)
			return false;
	}
	return true;
}

/// Fast 3x3 eigenvalue decomposition
bool eig3(Matrix3x3 &m, Float lambda[3]) {
	Float subd[3];

	/* Reduce to Hessenberg form */
	tred3(m, lambda, subd);

	/* Iteratively find the eigenvalues and eigenvectors
	   using implicitly shifted QL iterations */
	if (!ql3(m, lambda, subd))
		return false;

	/* Sort the eigenvalues in decreasing order */
	for (int i=0; i<2; ++i) {
		/* Locate the maximum eigenvalue. */
		int largest = i;
		Float maxValue = lambda[largest];
		for (int j = i+1; j<3; ++j) {
			if (lambda[j] > maxValue) {
				largest = j;
				maxValue = lambda[largest];
			}
		}

		if (largest != i) {
			// Swap the eigenvalues.
			std::swap(lambda[i], lambda[largest]);

			// Swap the eigenvectors
			for (int j=0; j<3; ++j) 
				std::swap(m(j, i), m(j, largest));
		}
	}
	
	return true;
}


MTS_NAMESPACE_END

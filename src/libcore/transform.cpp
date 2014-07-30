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
	Float sinTheta, cosTheta;

	/* Make sure that the axis is normalized */
	Vector naxis = normalize(axis);
	math::sincos(degToRad(angle), &sinTheta, &cosTheta);

	Matrix4x4 result;
	result(0, 0) = naxis.x * naxis.x + (1.0f - naxis.x * naxis.x) * cosTheta;
	result(0, 1) = naxis.x * naxis.y * (1.0f - cosTheta) - naxis.z * sinTheta;
	result(0, 2) = naxis.x * naxis.z * (1.0f - cosTheta) + naxis.y * sinTheta;
	result(0, 3) = 0;

	result(1, 0) = naxis.x * naxis.y * (1.0f - cosTheta) + naxis.z * sinTheta;
	result(1, 1) = naxis.y * naxis.y + (1.0f - naxis.y * naxis.y) * cosTheta;
	result(1, 2) = naxis.y * naxis.z * (1.0f - cosTheta) - naxis.x * sinTheta;
	result(1, 3) = 0;

	result(2, 0) = naxis.x * naxis.z * (1.0f - cosTheta) - naxis.y * sinTheta;
	result(2, 1) = naxis.y * naxis.z * (1.0f - cosTheta) + naxis.x * sinTheta;
	result(2, 2) = naxis.z * naxis.z + (1.0f - naxis.z * naxis.z) * cosTheta;
	result(2, 3) = 0;

	result(3, 0) = 0;
	result(3, 1) = 0;
	result(3, 2) = 0;
	result(3, 3) = 1;

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

	/* Perform a scale so that the field of view is mapped
	 * to the interval [-1, 1] */
	Float cot = 1.0f / std::tan(degToRad(fov / 2.0f));

	Matrix4x4 trafo(
		cot,  0,    0,   0,
		0,    cot,  0,   0,
		0,    0,    clipFar * recip, -clipNear * clipFar * recip,
		0,    0,    1,   0
	);


	return Transform(trafo);
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
	      b = -(clipFar + clipNear) / (clipFar - clipNear);

	Matrix4x4 trafo(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, a, b,
		0, 0, 0, 1
	);
	return Transform(trafo);
}

Transform Transform::glOrthographic(Float clipLeft, Float clipRight,
		Float clipBottom, Float clipTop, Float clipNear, Float clipFar) {
	Float fx =  2.0f / (clipRight - clipLeft),
	      fy =  2.0f / (clipTop - clipBottom),
	      fz = -2.0f / (clipFar - clipNear),
	      tx = -(clipRight + clipLeft) / (clipRight - clipLeft),
	      ty = -(clipTop + clipBottom) / (clipTop - clipBottom),
	      tz = -(clipFar + clipNear) / (clipFar - clipNear);

	Matrix4x4 trafo(
		fx,  0,  0,  tx,
		 0, fy,  0,  ty,
		 0,  0, fz,  tz,
		 0,  0,  0,   1
	);

	return Transform(trafo);
}

Transform Transform::lookAt(const Point &p, const Point &t, const Vector &up) {
	Vector dir = normalizeStrict(t-p, "lookAt(): 'origin' and 'target' coincide!");
	Vector left = normalizeStrict(cross(up, dir), "lookAt(): the forward and upward direction must be linearly independent!");
	Vector newUp = cross(dir, left);

	Matrix4x4 result, inverse;
	result(0, 0)  = left.x;  result(1, 0)  = left.y;  result(2, 0)  = left.z;  result(3, 0)  = 0;
	result(0, 1)  = newUp.x; result(1, 1)  = newUp.y; result(2, 1)  = newUp.z; result(3, 1)  = 0;
	result(0, 2)  = dir.x;   result(1, 2)  = dir.y;   result(2, 2)  = dir.z;   result(3, 2)  = 0;
	result(0, 3)  = p.x;     result(1, 3)  = p.y;     result(2, 3)  = p.z;     result(3, 3)  = 1;

	/* The inverse is simple to compute for this matrix, so do it directly here */
	Vector q(
		result(0, 0) * p.x + result(1, 0) * p.y + result(2, 0) * p.z,
		result(0, 1) * p.x + result(1, 1) * p.y + result(2, 1) * p.z,
		result(0, 2) * p.x + result(1, 2) * p.y + result(2, 2) * p.z);

	inverse(0, 0) = left.x; inverse(1, 0) = newUp.x; inverse(2, 0) = dir.x; inverse(3, 0) = 0;
	inverse(0, 1) = left.y; inverse(1, 1) = newUp.y; inverse(2, 1) = dir.y; inverse(3, 1) = 0;
	inverse(0, 2) = left.z; inverse(1, 2) = newUp.z; inverse(2, 2) = dir.z; inverse(3, 2) = 0;
	inverse(0, 3) = -q.x;   inverse(1, 3) = -q.y;    inverse(2, 3) = -q.z;  inverse(3, 3) = 1;

	return Transform(result, inverse);
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
 * \brief Tridiagonalizes a 3x3 matrix using a single
 * Householder transformation and returns the result.
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

/**
 * Compute the roots of the cubic polynomial.  Double-precision arithmetic
 * is used to minimize the effects due to subtractive cancellation.  The
 * roots are returned in increasing order.\
 *
 * Based on "Geometric Tools" by David Eberly.
 */
static void eig3_evals(const Matrix3x3 &A, double root[3]) {
	const double msInv3 = 1.0 / 3.0, msRoot3 = std::sqrt(3.0);

	// Convert the unique matrix entries to double precision.
	double a00 = (double) A(0, 0), a01 = (double) A(0, 1), a02 = (double) A(0, 2),
	       a11 = (double) A(1, 1), a12 = (double) A(1, 2), a22 = (double) A(2, 2);

	// The characteristic equation is x^3 - c2*x^2 + c1*x - c0 = 0.  The
	// eigenvalues are the roots to this equation, all guaranteed to be
	// real-valued, because the matrix is symmetric.
	double c0 = a00*a11*a22 + 2.0*a01*a02*a12 - a00*a12*a12 -
		a11*a02*a02 - a22*a01*a01;

	double c1 = a00*a11 - a01*a01 + a00*a22 - a02*a02 +
		a11*a22 - a12*a12;

	double c2 = a00 + a11 + a22;

	// Construct the parameters used in classifying the roots of the equation
	// and in solving the equation for the roots in closed form.
	double c2Div3 = c2*msInv3;
	double aDiv3 = (c1 - c2*c2Div3)*msInv3;
	if (aDiv3 > 0.0)
		aDiv3 = 0.0;

	double halfMB = 0.5*(c0 + c2Div3*(2.0*c2Div3*c2Div3 - c1));

	double q = halfMB*halfMB + aDiv3*aDiv3*aDiv3;
	if (q > 0.0)
		q = 0.0;

	// Compute the eigenvalues by solving for the roots of the polynomial.
	double magnitude = std::sqrt(-aDiv3);
	double angle = std::atan2(std::sqrt(-q), halfMB)*msInv3;
	double cs = std::cos(angle);
	double sn = std::sin(angle);
	double root0 = c2Div3 + 2.0*magnitude*cs;
	double root1 = c2Div3 - magnitude*(cs + msRoot3*sn);
	double root2 = c2Div3 - magnitude*(cs - msRoot3*sn);

	// Sort in increasing order.
	if (root1 <= root0) {
		root[0] = root0;
		root[1] = root1;
	} else {
		root[0] = root1;
		root[1] = root0;
	}

	if (root2 <= root[1]) {
		root[2] = root2;
	} else {
		root[2] = root[1];
		if (root2 <= root[0]) {
			root[1] = root2;
		} else {
			root[1] = root[0];
			root[0] = root2;
		}
	}
}

/**
 * Determine if M has positive rank.  The maximum-magnitude entry of M is
 * returned.  The row in which it is contained is also returned.
 *
 * Based on "Geometric Tools" by David Eberly.
 */
static bool eig3_rank(const Matrix3x3 &m, Float& maxEntry, Vector &maxRow) {
    // Locate the maximum-magnitude entry of the matrix.
    maxEntry = -1.0f;
    int maxRowIndex = -1;
    for (int row = 0; row < 3; ++row) {
        for (int col = row; col < 3; ++col) {
            Float absValue = std::abs(m(row, col));
            if (absValue > maxEntry) {
                maxEntry = absValue;
                maxRowIndex = row;
            }
        }
    }

    // Return the row containing the maximum, to be used for eigenvector
    // construction.
    maxRow = m.row(maxRowIndex);

    return maxEntry >= std::numeric_limits<Float>::epsilon();
}

/**
 * Compute the eigenvectors
 *
 * Based on "Geometric Tools" by David Eberly.
 */
static void eig3_evecs(Matrix3x3& A, const Float lambda[3], const Vector& U2, int i0, int i1, int i2) {
	Vector U0, U1;
	coordinateSystem(normalize(U2), U0, U1);

	// V[i2] = c0*U0 + c1*U1,  c0^2 + c1^2=1
	// e2*V[i2] = c0*A*U0 + c1*A*U1
	// e2*c0 = c0*U0.Dot(A*U0) + c1*U0.Dot(A*U1) = d00*c0 + d01*c1
	// e2*c1 = c0*U1.Dot(A*U0) + c1*U1.Dot(A*U1) = d01*c0 + d11*c1
	Vector tmp = A*U0, evecs[3];

	Float p00 = lambda[i2] - dot(U0, tmp),
		  p01 = dot(U1, tmp),
		  p11 = lambda[i2] - dot(U1, A*U1),
		  maxValue = std::abs(p00),
		  absValue = std::abs(p01),
		  invLength;

	if (absValue > maxValue)
		maxValue = absValue;

	int row = 0;
	absValue = std::abs(p11);
	if (absValue > maxValue) {
		maxValue = absValue;
		row = 1;
	}

	if (maxValue >= std::numeric_limits<Float>::epsilon()) {
		if (row == 0) {
			invLength = 1/std::sqrt(p00*p00 + p01*p01);
			p00 *= invLength;
			p01 *= invLength;
			evecs[i2] = p01*U0 + p00*U1;
		} else {
			invLength = 1/std::sqrt(p11*p11 + p01*p01);
			p11 *= invLength;
			p01 *= invLength;
			evecs[i2] = p11*U0 + p01*U1;
		}
	} else {
		if (row == 0)
			evecs[i2] = U1;
		else
			evecs[i2] = U0;
	}

	// V[i0] = c0*U2 + c1*Cross(U2,V[i2]) = c0*R + c1*S
	// e0*V[i0] = c0*A*R + c1*A*S
	// e0*c0 = c0*R.Dot(A*R) + c1*R.Dot(A*S) = d00*c0 + d01*c1
	// e0*c1 = c0*S.Dot(A*R) + c1*S.Dot(A*S) = d01*c0 + d11*c1
	Vector S = cross(U2, evecs[i2]);
	tmp = A*U2;
	p00 = lambda[i0] - dot(U2, tmp);
	p01 = dot(S, tmp);
	p11 = lambda[i0] - dot(S, A*S);
	maxValue = std::abs(p00);
	row = 0;
	absValue = std::abs(p01);
	if (absValue > maxValue)
		maxValue = absValue;
	absValue = std::abs(p11);
	if (absValue > maxValue) {
		maxValue = absValue;
		row = 1;
	}

	if (maxValue >= std::numeric_limits<Float>::epsilon()) {
		if (row == 0) {
			invLength = 1/std::sqrt(p00*p00 + p01*p01);
			p00 *= invLength;
			p01 *= invLength;
			evecs[i0] = p01*U2 + p00*S;
		} else {
			invLength = 1/std::sqrt(p11*p11 + p01*p01);
			p11 *= invLength;
			p01 *= invLength;
			evecs[i0] = p11*U2 + p01*S;
		}
	} else {
		if (row == 0)
			evecs[i0] = S;
		else
			evecs[i0] = U2;
	}

	evecs[i1] = cross(evecs[i2], evecs[i0]);
	A = Matrix3x3(
		evecs[0].x, evecs[1].x, evecs[2].x,
		evecs[0].y, evecs[1].y, evecs[2].y,
		evecs[0].z, evecs[1].z, evecs[2].z
	);
}

/**
 * \brief Fast non-iterative 3x3 eigenvalue decomposition
 *
 * Based on "Geometric Tools" by David Eberly.
 */
void eig3_noniter(Matrix3x3 &A, Float lambda[3]) {
	// Compute the eigenvalues using double-precision arithmetic.
	double root[3];
	eig3_evals(A, root);
	lambda[0] = (Float) root[0];
	lambda[1] = (Float) root[1];
	lambda[2] = (Float) root[2];
	Matrix3x3 eigs;

	Float maxEntry[3];
	Vector maxRow[3];
	for (int i = 0; i < 3; ++i) {
		Matrix3x3 M(A);
		M(0, 0) -= lambda[i];
		M(1, 1) -= lambda[i];
		M(2, 2) -= lambda[i];
		if (!eig3_rank(M, maxEntry[i], maxRow[i])) {
			A.setIdentity();
			return;
		}
	}

	Float totalMax = maxEntry[0];
	int i = 0;
	if (maxEntry[1] > totalMax) {
		totalMax = maxEntry[1];
		i = 1;
	} if (maxEntry[2] > totalMax) {
		i = 2;
	}

	if (i == 0) {
		maxRow[0] = normalize(maxRow[0]);
		eig3_evecs(A, lambda, maxRow[0], 1, 2, 0);
	} else if (i == 1) {
		maxRow[1] = normalize(maxRow[1]);
		eig3_evecs(A, lambda, maxRow[1], 2, 0, 1);
	} else {
		maxRow[2] = normalize(maxRow[2]);
		eig3_evecs(A, lambda, maxRow[2], 0, 1, 2);
	}
}

MTS_NAMESPACE_END

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

#include <mitsuba/core/transform.h>
#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN

// -----------------------------------------------------------------------
//  Linear transformation class
// -----------------------------------------------------------------------

Transform::Transform() {
	m_transform = new Matrix4x4();
	m_invTransform = new Matrix4x4();
}

Transform::Transform(const Matrix4x4 *trafo)
	: m_transform(trafo) {
	m_invTransform = m_transform->inverse();
}

Transform::Transform(const Matrix4x4 *trafo, const Matrix4x4 *invTrafo)
	: m_transform(trafo), m_invTransform(invTrafo) {
}

Transform Transform::inverse() const {
	return Transform(m_invTransform.get(), m_transform.get());
}

Transform Transform::operator*(const Transform &t) const {
	ref<Matrix4x4> trafo = m_transform * t.m_transform;
	ref<Matrix4x4> invTrafo = t.m_invTransform * m_invTransform;
	return Transform(trafo, invTrafo);
}

Transform Transform::translate(const Vector &v) {
	ref<Matrix4x4> trafo = new Matrix4x4(
		1, 0, 0, v.x,
		0, 1, 0, v.y,
		0, 0, 1, v.z,
		0, 0, 0, 1
	);
	ref<Matrix4x4> invTrafo = new Matrix4x4(
		1, 0, 0, -v.x,
		0, 1, 0, -v.y,
		0, 0, 1, -v.z,
		0, 0, 0, 1
	);
	return Transform(trafo, invTrafo);
}

Transform Transform::scale(const Vector &v) {
	ref<Matrix4x4> trafo = new Matrix4x4(
		v.x, 0,   0,   0,
		0,   v.y, 0,   0,
		0,   0,   v.z, 0,
		0,   0,   0,   1
	);
	ref<Matrix4x4> invTrafo = new Matrix4x4(
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

	Float m[4][4], rotAngle = degToRad(angle);
	Float sinAngle = std::sin(rotAngle);
	Float cosAngle = std::cos(rotAngle);

	m[0][0] = naxis.x * naxis.x + (1.f - naxis.x * naxis.x) * cosAngle;
	m[0][1] = naxis.x * naxis.y * (1.f - cosAngle) - naxis.z * sinAngle;
	m[0][2] = naxis.x * naxis.z * (1.f - cosAngle) + naxis.y * sinAngle;
	m[0][3] = 0;

	m[1][0] = naxis.x * naxis.y * (1.f - cosAngle) + naxis.z * sinAngle;
	m[1][1] = naxis.y * naxis.y + (1.f - naxis.y * naxis.y) * cosAngle;
	m[1][2] = naxis.y * naxis.z * (1.f - cosAngle) - naxis.x * sinAngle;
	m[1][3] = 0;

	m[2][0] = naxis.x * naxis.z * (1.f - cosAngle) - naxis.y * sinAngle;
	m[2][1] = naxis.y * naxis.z * (1.f - cosAngle) + naxis.x * sinAngle;
	m[2][2] = naxis.z * naxis.z + (1.f - naxis.z * naxis.z) * cosAngle;
	m[2][3] = 0;

	m[3][0] = 0;
	m[3][1] = 0;
	m[3][2] = 0;
	m[3][3] = 1;

	/* The matrix is orthonormal */
	ref<Matrix4x4> matrix = new Matrix4x4(m);
	return Transform(matrix, matrix->transpose());
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

	ref<Matrix4x4> trafo = new Matrix4x4(
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

	ref<Matrix4x4> trafo = new Matrix4x4(
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

	ref<Matrix4x4> trafo = new Matrix4x4(
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

	ref<Matrix4x4> trafo = new Matrix4x4(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, a, b,
		0, 0, 0, 1
	);
	return Transform(trafo);
}


Transform Transform::lookAt(const Point &p, const Point &t, const Vector &up) {
	Float m[4][4];

	Vector dirct = normalize(t-p);
	Vector right = normalize(cross(dirct, up));

	/* Generate a new, orthogonalized up vector */
	Vector newUp = cross(right, dirct);

	/* Store as columns */
	m[0][0] = right.x; m[1][0] = right.y; m[2][0] = right.z; m[3][0] = 0;
	m[0][1] = newUp.x; m[1][1] = newUp.y; m[2][1] = newUp.z; m[3][1] = 0;
	m[0][2] = dirct.x; m[1][2] = dirct.y; m[2][2] = dirct.z; m[3][2] = 0;
	m[0][3] = p.x;     m[1][3] = p.y;     m[2][3] = p.z;     m[3][3] = 1;

	ref<Matrix4x4> mtx = new Matrix4x4(m);

	return Transform(mtx, mtx->inverse());
}

Transform Transform::fromFrame(const Frame &frame) {
	ref<Matrix4x4> mtx = new Matrix4x4(
		frame.s.x, frame.t.x, frame.n.x, 0,
		frame.s.y, frame.t.y, frame.n.y, 0, 
		frame.s.z, frame.t.z, frame.n.z, 0,
		0, 0, 0, 1
	);
	return Transform(mtx, mtx->transpose());
}

std::string Transform::toString() const {
	return m_transform->toString();
}

// -----------------------------------------------------------------------
//  Matrix4x4 reference counted array wrapper
// -----------------------------------------------------------------------

Matrix4x4::Matrix4x4() {
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			m[i][j] = (i == j) ? 1.0f : 0.0f;
}

Matrix4x4::Matrix4x4(Float value) {
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			m[i][j] = value;
}

Matrix4x4::Matrix4x4(Float _m[4][4]) {
	memcpy(m, _m, sizeof(Float) * 16);
}

Matrix4x4::Matrix4x4(const Matrix4x4 *mat) {
	memcpy(m, mat->m, sizeof(Float) * 16);
}

Matrix4x4::Matrix4x4(Stream *stream) {
	stream->readFloatArray(&m[0][0], 16);
}

Matrix4x4::Matrix4x4(
	Float a00, Float a01, Float a02, Float a03,
	Float a10, Float a11, Float a12, Float a13,
	Float a20, Float a21, Float a22, Float a23,
	Float a30, Float a31, Float a32, Float a33)
{
	m[0][0] = a00; m[0][1] = a01; m[0][2] = a02; m[0][3] = a03;
	m[1][0] = a10; m[1][1] = a11; m[1][2] = a12; m[1][3] = a13;
	m[2][0] = a20; m[2][1] = a21; m[2][2] = a22; m[2][3] = a23;
	m[3][0] = a30; m[3][1] = a31; m[3][2] = a32; m[3][3] = a33;
}

std::string Matrix4x4::toString() const {
	std::ostringstream oss;
	oss << "Matrix4x4[" << std::endl;
	for (int i=0; i<4; i++) {
		oss << "  [";
		for (int j=0; j<4; j++) {
			oss << m[i][j];
			if (j != 3)
				oss << ", ";
		}
		oss << "]";

		if (i != 3)
			oss << ",";
		oss << std::endl;
	}
	oss << "]";
	return oss.str();
}

ref<Matrix4x4> Matrix4x4::transpose() const {
	ref<Matrix4x4> matrix = new Matrix4x4();
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			matrix->m[i][j] = m[j][i];
	return matrix;
}

ref<Matrix4x4> Matrix4x4::inverse() const {
	int indxc[4], indxr[4];
	int ipiv[4] = { 0, 0, 0, 0 };
	Float minv[4][4];
	memcpy(minv, m, 4*4*sizeof(Float));

	for (int i = 0; i < 4; i++) {
		int irow = -1, icol = -1;
		Float big = 0.;
		for (int j = 0; j < 4; j++) {
			if (ipiv[j] != 1) {
				for (int k = 0; k < 4; k++) {
					if (ipiv[k] == 0) {
						if (std::abs(minv[j][k]) >= big) {
							big = std::abs(minv[j][k]);
							irow = j;
							icol = k;
						}
					}
					else if (ipiv[k] > 1)
						Log(EError, "Singular matrix in Matrix4x4::inverse");
				}
			}
		}
		++ipiv[icol];
		if (irow != icol) {
			for (int k = 0; k < 4; ++k)
				std::swap(minv[irow][k], minv[icol][k]);
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (minv[icol][icol] == 0.)
			Log(EError, "Singular matrix in Matrix4x4::inverse");
		Float pivinv = 1.f / minv[icol][icol];
		minv[icol][icol] = 1.f;
		for (int j = 0; j < 4; j++)
			minv[icol][j] *= pivinv;
		for (int j = 0; j < 4; j++) {
			if (j != icol) {
				Float save = minv[j][icol];
				minv[j][icol] = 0;
				for (int k = 0; k < 4; k++)
					minv[j][k] -= minv[icol][k]*save;
			}
		}
	}
	for (int j = 3; j >= 0; j--) {
		if (indxr[j] != indxc[j]) {
			for (int k = 0; k < 4; k++)
				std::swap(minv[k][indxr[j]], minv[k][indxc[j]]);
		}
	}
	return new Matrix4x4(minv);
}

Float Matrix4x4::det3x3() const {
	return ((m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]))
		  - (m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]))
		  + (m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])));
}

void Matrix4x4::serialize(Stream *stream) const {
	stream->writeFloatArray(&m[0][0], 16);
}

static void tred2(Float V[4][4], Float d[4], Float e[4]);
static void tql2(Float V[4][4], Float d[4], Float e[4]);

void Matrix4x4::symmEigenDecomp(Matrix4x4 *Q, Vector4 &d) {
	Float e[4];
	*Q = *this;
	tred2(Q->m, (Float *) &d, e);
	tql2(Q->m, (Float *) &d, e);
}

/* Eigendecomposition code for symmetric 3x3 matrices
   Copied from the public domain Java Matrix library JAMA. */

inline Float hypot2(Float x, Float y){
	return std::sqrt(x * x + y * y);
}

// Symmetric Householder reduction to tridiagonal form.
static void tred2(Float V[4][4], Float d[4], Float e[4]) {
	//  This is derived from the Algol procedures tred2 by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	const int n = 4;
	for (int j = 0; j < n; j++)
		d[j] = V[n - 1][j];
	
	// Householder reduction to tridiagonal form.
	for (int i = n - 1; i > 0; i--) {
		// Scale to avoid under/overflow.
		Float scale = 0.0f, h = 0.0f;

		for (int k = 0; k < i; k++)
			scale = scale + std::abs(d[k]);

		if (scale == 0.0f) {
			e[i] = d[i - 1];
			for (int j = 0; j < i; j++) {
				d[j] = V[i - 1][j];
				V[i][j] = 0.0f;
				V[j][i] = 0.0f;
			}
		} else {
			// Generate Householder vector.

			for (int k = 0; k < i; k++) {
				d[k] /= scale;
				h += d[k] * d[k];
			}
			Float f = d[i - 1], 
				  g = std::sqrt(h);

			if (f > 0)
				g = -g;
			e[i] = scale * g;
			h = h - f * g;
			d[i - 1] = f - g;

			for (int j = 0; j < i; j++)
				e[j] = 0.0f;

			 // Apply similarity transformation to remaining columns.
			for (int j = 0; j < i; j++) {
				f = d[j];
				V[j][i] = f;
				g = e[j] + V[j][j] * f;
				for (int k = j + 1; k <= i - 1; k++) {
					g += V[k][j] * d[k];
					e[k] += V[k][j] * f;
				}
				e[j] = g;
			}

			f = 0.0f;
			for (int j = 0; j < i; j++) {
				e[j] /= h;
				f += e[j] * d[j];
			}
			Float hh = f / (h + h);

			for (int j = 0; j < i; j++)
				e[j] -= hh * d[j];
			for (int j = 0; j < i; j++) {
				f = d[j];
				g = e[j];
				for (int k = j; k <= i - 1; k++)
					V[k][j] -= (f * e[k] + g * d[k]);
				d[j] = V[i - 1][j];
				V[i][j] = 0.0f;
			}
		}
		d[i] = h;
	}

	// Accumulate transformations.
	for (int i = 0; i < n - 1; i++) {
		V[n - 1][i] = V[i][i];
		V[i][i] = 1.0f;
		Float h = d[i + 1];

		if (h != 0.0f) {
			for (int k = 0; k <= i; k++)
				d[k] = V[k][i + 1] / h;
			for (int j = 0; j <= i; j++) {
				Float g = 0.0f;

				for (int k = 0; k <= i; k++)
					g += V[k][i + 1] * V[k][j];
				for (int k = 0; k <= i; k++)
					V[k][j] -= g * d[k];
			}
		}
		for (int k = 0; k <= i; k++)
			V[k][i + 1] = 0.0f;
	}
	for (int j = 0; j < n; j++) {
		d[j] = V[n - 1][j];
		V[n - 1][j] = 0.0f;
	}
	V[n - 1][n - 1] = 1.0f;
	e[0] = 0.0;
}

// Symmetric tridiagonal QL algorithm.
static void tql2(Float V[4][4], Float d[4], Float e[4]) {
	//  This is derived from the Algol procedures tql2, by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	const int n = 4;

	for (int i = 1; i < n; i++)
		e[i - 1] = e[i];

	e[n - 1] = 0.0f;
	Float f = 0.0f, tst1 = 0.0f;

#if defined(SINGLE_PRECISION)
	Float eps = pow(2.0f, -23.0f);
#else
	Float eps = pow(2.0, -52.0);
#endif

	for (int l = 0; l < n; l++) {
		// Find small subdiagonal element
		tst1 = std::max(tst1, std::abs(d[l]) + std::abs(e[l]));
		int m = l;

		while (m < n) {
			if (std::abs(e[m]) <= eps * tst1)
				break;
			m++;
		}

		// If m == l, d[l] is an eigenvalue,
		// otherwise, iterate.
		if (m > l) {
			int iter = 0;

			do {
				iter = iter + 1;    // (Could check iteration count here.)

				// Compute implicit shift
				Float g = d[l];
				Float p = (d[l + 1] - g) / (2.0f * e[l]);
				Float r = hypot2(p, 1.0f);

				if (p < 0)
					r = -r;
				d[l] = e[l] / (p + r);
				d[l + 1] = e[l] * (p + r);
				Float dl1 = d[l + 1];

				Float h = g - d[l];

				for (int i = l + 2; i < n; i++)
					d[i] -= h;
				f = f + h;

				// Implicit QL transformation.
				p = d[m];
				Float c = 1.0f;
				Float c2 = c, c3 = c;
				Float el1 = e[l + 1];
				Float s = 0.0f, s2 = 0.0f;

				for (int i = m - 1; i >= l; i--) {
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c * e[i];
					h = c * p;
					r = hypot2(p, e[i]);
					e[i + 1] = s * r;
					s = e[i] / r;
					c = p / r;
					p = c * d[i] - s * g;
					d[i + 1] = h + s * (c * g + s * d[i]);

					// Accumulate transformation.
					for (int k = 0; k < n; k++) {
						h = V[k][i + 1];
						V[k][i + 1] =
							s * V[k][i] + c * h;
						V[k][i] = c * V[k][i] - s * h;
					}
				}
				p = -s * s2 * c3 * el1 * e[l] / dl1;
				e[l] = s * p;
				d[l] = c * p;
				// Check for convergence.
			} while (std::abs(e[l]) > eps * tst1);
		}
		d[l] = d[l] + f;
		e[l] = 0.0f;
	}

	// Sort eigenvalues and corresponding vectors.
	for (int i = 0; i < n - 1; i++) {
		int k = i;

		Float p = d[i];

		for (int j = i + 1; j < n; j++) {
			if (d[j] < p) {
				k = j;
				p = d[j];
			}
		}

		if (k != i) {
			d[k] = d[i];
			d[i] = p;
			for (int j = 0; j < n; j++) {
				p = V[j][i];
				V[j][i] = V[j][k];
				V[j][k] = p;
			}
		}
	}
}

MTS_IMPLEMENT_CLASS(Matrix4x4, false, Object)
MTS_NAMESPACE_END

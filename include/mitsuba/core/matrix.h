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

#if !defined(_MATRIX_H)
#define _MATRIX_H

#include <mitsuba/mitsuba.h>
#include <boost/static_assert.hpp>

MTS_NAMESPACE_BEGIN

/**
 * Generic fixed-size dense matrix class using a row-major storage format
 */
template <int M, int N, typename T> struct Matrix {
public:
	T m[M][N];

	/** 
	 * \brief Construct a new MxN matrix without initializing it.
	 * 
	 * This construtor is useful when the matrix will either not
	 * be used at all (it might be part of a larger data structure)
	 * or initialized at a later point in time. Always make sure
	 * that one of the two is the case! Otherwise your program will do
	 * computations involving uninitialized memory, which will probably
	 * lead to a difficult-to-find bug.
	 */
#if !defined(MTS_DEBUG_UNINITIALIZED)
	Matrix() { }
#else
	Matrix() {
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				m[i][j] = std::numeric_limits<double>::quiet_NaN();
	}
#endif

	/// Initialize the matrix with constant entries
	explicit inline Matrix(T value) {
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				m[i][j] = value;
	}

	/// Initialize the matrix from a given MxN array
	explicit inline Matrix(T _m[M][N]) {
		memcpy(m, _m, sizeof(T) * M * N);
	}

	/// Unserialize a matrix from a stream
	explicit inline Matrix(Stream *stream) {
		stream->readArray(&m[0][0], M * N);
	}

	/// Initialize with the identity matrix
	void setIdentity() {
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				m[i][j] = (i == j) ? 1.0f : 0.0f;
	}

	/// Indexing operator
	inline T &operator()(int i, int j) { return m[i][j]; }

	/// Indexing operator (const verions)
	inline const T & operator()(int i, int j) const { return m[i][j]; }

	/// Equality operator
	inline bool operator==(const Matrix &mat) const {
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				if (m[i][j] != mat.m[i][j])
					return false;
		return true;
	}

	/// Inequality operator
	inline bool operator!=(const Matrix &mat) const {
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				if (m[i][j] != mat.m[i][j])
					return true;
		return false;
	}

	/// Matrix addition (returns a temporary)
	inline Matrix operator+(const Matrix &mat) const {
		Matrix result;
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				result.m[i][j] = m[i][j] + mat.m[i][j];
		return result;
	}

	/// Matrix-scalar addition (returns a temporary)
	inline Matrix operator+(T value) const {
		Matrix result;
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				result.m[i][j] = m[i][j] + value;
		return result;
	}

	/// Matrix addition
	inline const Matrix &operator+=(const Matrix &mat) {
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				m[i][j] += mat.m[i][j];
		return *this;
	}

	/// Matrix-scalar addition
	inline const Matrix &operator+=(T value) {
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				m[i][j] += value;
		return *this;
	}

	/// Matrix subtraction (returns a temporary)
	inline Matrix operator-(const Matrix &mat) const {
		Matrix result;
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				result.m[i][j] = m[i][j] - mat.m[i][j];
		return result;
	}

	/// Matrix-scalar subtraction (returns a temporary)
	inline Matrix operator-(T value) const {
		Matrix result;
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				result.m[i][j] = m[i][j] - value;
		return result;
	}

	/// Matrix subtraction
	inline const Matrix &operator-=(const Matrix &mat) {
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				m[i][j] -= mat.m[i][j];
		return *this;
	}

	/// Matrix-scalar subtraction
	inline const Matrix &operator-(T value) {
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				m[i][j] -= value;
		return *this;
	}

	/// Component-wise negation
	inline Matrix operator-() const {
		Matrix result;
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				result.m[i][j] = -m[i][j];
		return result;
	}

	/// Scalar multiplication (creates a temporary)
	inline Matrix operator*(T value) const {
		Matrix result;
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				result.m[i][j] = m[i][j]*value;
		return result;
	}

	/// Scalar multiplication
	inline const Matrix& operator*=(T value) {
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				m[i][j] *= value;
		return *this;
	}

	/// Matrix multiplication (creates a temporary)
	template <int M2, int N2> inline Matrix<M, N2, T> operator*(const Matrix<M2, N2, T> &mat) const {
		BOOST_STATIC_ASSERT(N == M2);
		Matrix<M, N2, T> result;
		for (int i=0; i<M; ++i) {
			for (int j=0; j<N2; ++j) {
				T sum = 0;
				for (int k=0; k<N; ++k)
					sum += m[i][k] * mat.m[k][j];
				result.m[i][j] = sum;
			}
		}
		return result;
	}

	/// Matrix multiplication (for square matrices)
	inline const Matrix &operator*=(const Matrix &mat) {
		BOOST_STATIC_ASSERT(M == N);
		Matrix temp = *this * mat;
		*this = temp;
		return *this;
	}

	/// Compute the inverse of this matrix using the Gauss-Jordan algorithm
	bool invert(Matrix &target) const;

	/// Perform a symmetric eigendecomposition into Q and D.
	void symmEigenDecomp(Matrix &Q, T *d) const {
		BOOST_STATIC_ASSERT(M == N);
		T e[M];
		Q = *this;
		tred2(Q.m, d, e);
		tql2(Q.m, d, e);
	}

	/// Compute the transpose of this matrix
	void transpose(Matrix<N, M, T> &target) const {
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				target.m[i][j] = m[j][i];
	}

	/// Serialize the matrix to a stream
	void serialize(Stream *stream) const {
		stream->writeArray(&m[0][0], M*N);
	}

	/// Return a string representation
	std::string toString() const {
		std::ostringstream oss;
		oss << "Matrix[" << M << "x" << N << ","<< std::endl;
		for (int i=0; i<4; ++i) {
			oss << "  [";
			for (int j=0; j<4; ++j) {
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
protected:
	/// Symmetric Householder reduction to tridiagonal form.
	static void tred2(T V[M][N], T d[M], T e[M]);

	/// Symmetric tridiagonal QL algorithm.
	static void tql2(T V[M][N], T d[N], T e[N]);
};

/**
 * \brief Basic 4x4 matrix data type
 */
struct MTS_EXPORT_CORE Matrix4x4 : public Matrix<4, 4, Float> {
	Matrix4x4() { }

	/// Initialize the matrix with constant entries
	explicit inline Matrix4x4(Float value) : Matrix<4, 4, Float>(value) { }

	/// Initialize the matrix from a given MxN array
	explicit inline Matrix4x4(Float _m[4][4]) : Matrix<4, 4, Float>(_m) { }

	/// Unserialize a matrix from a stream
	explicit inline Matrix4x4(Stream *stream) : Matrix<4, 4, Float>(stream) { }

	/// Initialize with the given values
	inline Matrix4x4(
		Float a00, Float a01, Float a02, Float a03,
		Float a10, Float a11, Float a12, Float a13,
		Float a20, Float a21, Float a22, Float a23,
		Float a30, Float a31, Float a32, Float a33) {
		m[0][0] = a00; m[0][1] = a01; m[0][2] = a02; m[0][3] = a03;
		m[1][0] = a10; m[1][1] = a11; m[1][2] = a12; m[1][3] = a13;
		m[2][0] = a20; m[2][1] = a21; m[2][2] = a22; m[2][3] = a23;
		m[3][0] = a30; m[3][1] = a31; m[3][2] = a32; m[3][3] = a33;
	}

	/// Return the determinant of the upper left 3x3 sub-matrix
	inline Float det3x3() const {
		return ((m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]))
			  - (m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]))
			  + (m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])));
	}
};

MTS_NAMESPACE_END

#include "matrix.inl"

#endif /* _MATRIX_H */

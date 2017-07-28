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
#if !defined(__MITSUBA_CORE_MATRIX_H_)
#define __MITSUBA_CORE_MATRIX_H_

#include <mitsuba/mitsuba.h>
#include <boost/static_assert.hpp>

MTS_NAMESPACE_BEGIN

/**
 * \brief Generic fixed-size dense matrix class using a row-major storage format
 *
 * \ingroup libcore
 */
template <int M, int N, typename T> struct Matrix {
public:
    T m[M][N];

    /**
     * \brief Construct a new MxN matrix without initializing it.
     *
     * This constructor is useful when the matrix will either not
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
    explicit inline Matrix(const T _m[M][N]) {
        memcpy(m, _m, sizeof(T) * M * N);
    }

    /// Initialize the matrix from a given (flat) MxN array in row-major order
    explicit inline Matrix(const T _m[M*N]) {
        memcpy(m, _m, sizeof(T) * M * N);
    }

    /// Unserialize a matrix from a stream
    explicit inline Matrix(Stream *stream) {
        stream->readArray(&m[0][0], M * N);
    }

    /// Copy constructor
    inline Matrix(const Matrix &mtx) {
        memcpy(m, mtx.m, sizeof(T) * M * N);
    }

    /// Initialize with the identity matrix
    void setIdentity() {
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                m[i][j] = (i == j) ? 1.0f : 0.0f;
    }

    /// Initialize with zeroes
    void setZero() {
        memset(m, 0, sizeof(T) * M * N);
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

    /// Assignment operator
    inline Matrix &operator=(const Matrix &mat) {
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                m[i][j] = mat.m[i][j];
        return *this;
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

    /// Matrix-scalar addition
    inline const Matrix &operator-=(T value) {
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

    /// Scalar division (creates a temporary)
    inline Matrix operator/(T value) const {
        Matrix result;
#ifdef MTS_DEBUG
        if (value == 0)
            SLog(EWarn, "Matrix: Division by zero!");
#endif
        Float recip = 1/value;
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                result.m[i][j] = m[i][j]*recip;
        return result;
    }

    /// Scalar division
    inline const Matrix& operator/=(T value) {
#ifdef MTS_DEBUG
        if (value == 0)
            SLog(EWarn, "Matrix: Division by zero!");
#endif
        Float recip = 1/value;
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                m[i][j] *= recip;
        return *this;
    }

    /// Matrix multiplication (for square matrices)
    inline const Matrix &operator*=(const Matrix &mat) {
        BOOST_STATIC_ASSERT(M == N);
        Matrix temp = *this * mat;
        *this = temp;
        return *this;
    }

    /// Compute the trace of a square matrix
    inline Float trace() const {
        BOOST_STATIC_ASSERT(M == N);
        Float sum = 0;
        for (int i=0; i<M; ++i)
            sum += m[i][i];
        return sum;
    }

    /// Compute the Frobenius norm
    inline Float frob() const {
        Float val = 0;
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                val += m[i][j] * m[i][j];
        return std::sqrt(val);
    }

    /**
     * \brief Compute the LU decomposition of a matrix
     *
     * For an m-by-n matrix A with m >= n, the LU decomposition is an
     * m-by-n unit lower triangular matrix L, an n-by-n upper triangular
     * matrix U,
     *
     * and a permutation vector piv of length m so that A(piv,:) = L*U.
     * If m < n, then L is m-by-m and U is m-by-n.
     *
     * The LU decomposition with pivoting always exists, even if the matrix is
     * singular, so the constructor will never fail.
     * The primary use of the
     *
     * LU decomposition is in the solution of square systems of simultaneous
     * linear equations.
     *
     * \param Target matrix (the L and U parts will be stored
     *        together in a packed format)
     * \param piv Storage for the permutation created by the pivoting
     * \param pivsign Sign of the permutation
     * \return \c true if the matrix was nonsingular.
     *
     * Based on the implementation in JAMA.
     */
    bool lu(Matrix &LU, int piv[M], int &pivsign) const;

    /**
     * Compute the Cholesky decomposition of a symmetric
     * positive definite matrix
     *
     * \param L Target matrix (a lower triangular matrix such
     *    that A=L*L')
     * \return \c false If the matrix is not symmetric positive
     *    definite.
     * Based on the implementation in JAMA.
     */
    bool chol(Matrix &L) const;

    /**
     * Solve A*X==B, where \a this is a Cholesky decomposition of \a A created by \ref chol()
     *
     * \param B A matrix with as many rows as \a A and any number of columns
     * \param X A matrix such that L*L'*X == B
     *
     * Based on the implementation in JAMA.
     */
    template <int K> void cholSolve(const Matrix<M, K, T> &B,
            Matrix<M, K, T> &X) const;

    /**
     * Solve A*X==B, where \a this is a LU decomposition of \a A created by \ref lu()
     *
     * \param B A matrix with as many rows as \a A and any number of  columns
     * \param X A matrix such that L*U*X == B(piv, :)
     * \param piv Pivot vector returned by \ref lu()
     *
     * Based on the implementation in JAMA.
     */
    template <int K> void luSolve(const Matrix<M, K, T> &B,
            Matrix<M, K, T> &X, int piv[M]) const;

    /**
     * \brief Compute the determinant of a decomposed matrix
     * created by \ref lu()
     *
     * \param pivsign The sign of the pivoting permutation returned
     *        by \ref lu()
     *
     * Based on the implementation in JAMA.
     */
    T luDet(int pivsign) const;

    /**
     * \brief Compute the determinant of a decomposed matrix
     * created by \ref chol()
     */

    T cholDet() const;

    /// Check if the matrix is identically zero
    inline bool isZero() const {
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                if (m[i][j] != 0)
                    return false;
        return true;
    }

    /// Test if this is the identity matrix
    inline bool isIdentity() const {
        for (int i=0; i<M; ++i) {
            for (int j=0; j<N; ++j) {
                if (m[i][j] != ((i==j) ? 1 : 0))
                    return false;
            }
        }
        return true;
    }

    /**
     * \brief Compute the determinant of a square matrix (internally
     * creates a LU decomposition)
     */
    inline T det() const {
        Matrix LU;
        int piv[M], pivsign;
        if (!lu(LU, piv, pivsign))
            return 0.0f;
        return LU.luDet(pivsign);
    }

    /// Compute the inverse of a square matrix using the Gauss-Jordan algorithm
    bool invert(Matrix &target) const;

    /**
     * \brief Perform a symmetric eigendecomposition of a square matrix
     * into Q and D.
     *
     * Based on the implementation in JAMA.
     */
    inline void symEig(Matrix &Q, T d[M]) const {
        BOOST_STATIC_ASSERT(M == N);
        T e[M];
        Q = *this;
        tred2(Q.m, d, e);
        tql2(Q.m, d, e);
    }

    /// Compute the transpose of this matrix
    inline void transpose(Matrix<N, M, T> &target) const {
        for (int i=0; i<M; ++i)
            for (int j=0; j<N; ++j)
                target.m[i][j] = m[j][i];
    }

    /// Serialize the matrix to a stream
    inline void serialize(Stream *stream) const {
        stream->writeArray(&m[0][0], M*N);
    }

    /// Return a string representation
    std::string toString() const {
        std::ostringstream oss;
        oss << "Matrix" << M << "x" << N << "["<< std::endl;
        for (int i=0; i<M; ++i) {
            oss << "  ";
            for (int j=0; j<N; ++j) {
                oss << m[i][j];
                if (j != N-1)
                    oss << ", ";
            }
            if (i != M-1)
                oss << ";";
            oss << std::endl;
        }
        oss << "]";
        return oss.str();
    }
protected:
    /// Symmetric Householder reduction to tridiagonal form.
    static void tred2(T V[M][N], T d[N], T e[N]);

    /// Symmetric tridiagonal QL algorithm.
    static void tql2(T V[M][N], T d[N], T e[N]);
};

/**
 * \brief Basic 2x2 matrix data type
 * \ingroup libcore
 */
struct MTS_EXPORT_CORE Matrix2x2 : public Matrix<2, 2, Float> {
public:
    inline Matrix2x2() { }

    /// Initialize the matrix with constant entries
    explicit inline Matrix2x2(Float value) : Matrix<2, 2, Float>(value) { }

    /// Initialize the matrix from a given 2x2 array
    explicit inline Matrix2x2(const Float _m[2][2]) : Matrix<2, 2, Float>(_m) { }

    /// Initialize the matrix from a given (float) 2x2 array in row-major order
    explicit inline Matrix2x2(const Float _m[4]) : Matrix<2, 2, Float>(_m) { }

    /// Initialize the matrix from two 2D column vectors
    explicit inline Matrix2x2(const Vector2 &v1, const Vector2 &v2) {
        m[0][0] = v1.x; m[0][1] = v2.x;
        m[1][0] = v1.y; m[1][1] = v2.y;
    }

    /// Unserialize a matrix from a stream
    explicit inline Matrix2x2(Stream *stream) : Matrix<2, 2, Float>(stream) { }

    /// Copy constructor
    inline Matrix2x2(const Matrix<2, 2, Float> &mtx) : Matrix<2, 2, Float>(mtx) { }

    /// Initialize with the given values
    inline Matrix2x2(Float a00, Float a01, Float a10, Float a11) {
        m[0][0] = a00; m[0][1] = a01;
        m[1][0] = a10; m[1][1] = a11;
    }

    /// Return the determinant (Faster than Matrix::det)
    inline Float det() const {
        return m[0][0]*m[1][1] - m[0][1]*m[1][0];
    }

    /// Compute the inverse (Faster than Matrix::invert)
    FINLINE bool invert2x2(Matrix2x2 &target) const {
        Float det = m[0][0]*m[1][1] - m[0][1]*m[1][0];
        if (std::abs(det) <= RCPOVERFLOW)
            return false;
        Float invDet = 1/det;
        target.m[0][0] =  m[1][1] * invDet;
        target.m[0][1] = -m[0][1] * invDet;
        target.m[1][1] =  m[0][0] * invDet;
        target.m[1][0] = -m[1][0] * invDet;
        return true;
    }
    FINLINE bool invert(Matrix2x2 &target) const {
        return invert2x2(target);
    }

     /// Compute the inverse with det (Faster than Matrix::invert)
     FINLINE bool invert2x2(Matrix2x2 &target, Float& det) const {
        det = m[0][0]*m[1][1] - m[0][1]*m[1][0];
        if (std::abs(det) <= RCPOVERFLOW)
            return false;
        Float invDet = 1/det;
        target.m[0][0] = m[1][1] * invDet;
        target.m[0][1] = -m[0][1] * invDet;
        target.m[1][1] = m[0][0] * invDet;
        target.m[1][0] = -m[1][0] * invDet;
        return true;
     }

    /// Matrix-vector multiplication
    inline Vector2 operator*(const Vector2 &v) const {
        return Vector2(
            m[0][0] * v.x + m[0][1] * v.y,
            m[1][0] * v.x + m[1][1] * v.y
        );
    }

    /// Scalar multiplication (creates a temporary)
    inline Matrix2x2 operator*(Float value) const {
        Matrix2x2 result;
        for (int i=0; i<2; ++i)
            for (int j=0; j<2; ++j)
                result.m[i][j] = m[i][j]*value;
        return result;
    }

    /// Assignment operator
    inline Matrix2x2 &operator=(const Matrix<2, 2, Float> &mat) {
        for (int i=0; i<2; ++i)
            for (int j=0; j<2; ++j)
                m[i][j] = mat.m[i][j];
        return *this;
    }

    /// Return a row by index
    inline Vector2 row(int i) const {
        return Vector2(m[i][0], m[i][1]);
    }

    /// Return a column by index
    inline Vector2 col(int i) const {
        return Vector2(m[0][i], m[1][i]);
    }
};

/**
 * \brief Basic 3x3 matrix data type
 * \ingroup libcore
 */
struct MTS_EXPORT_CORE Matrix3x3 : public Matrix<3, 3, Float> {
public:
    inline Matrix3x3() { }

    /// Initialize the matrix with constant entries
    explicit inline Matrix3x3(Float value) : Matrix<3, 3, Float>(value) { }

    /// Initialize the matrix from a given 3x3 array
    explicit inline Matrix3x3(const Float _m[3][3]) : Matrix<3, 3, Float>(_m) { }

    /// Initialize the matrix from a given (float) 3x3 array in row-major order
    explicit inline Matrix3x3(const Float _m[9]) : Matrix<3, 3, Float>(_m) { }

    /// Initialize the matrix from three 3D column vectors
    explicit inline Matrix3x3(const Vector &v1, const Vector &v2, const Vector &v3) {
        m[0][0] = v1.x; m[0][1] = v2.x; m[0][2] = v3.x;
        m[1][0] = v1.y; m[1][1] = v2.y; m[1][2] = v3.y;
        m[2][0] = v1.z; m[2][1] = v2.z; m[2][2] = v3.z;
    }

    /// Unserialize a matrix from a stream
    explicit inline Matrix3x3(Stream *stream) : Matrix<3, 3, Float>(stream) { }

    /// Copy constructor
    inline Matrix3x3(const Matrix<3, 3, Float> &mtx) : Matrix<3, 3, Float>(mtx) { }

    /// Initialize with the given values
    inline Matrix3x3(Float a00, Float a01, Float a02,
            Float a10, Float a11, Float a12,
            Float a20, Float a21, Float a22) {
        m[0][0] = a00; m[0][1] = a01; m[0][2] = a02;
        m[1][0] = a10; m[1][1] = a11; m[1][2] = a12;
        m[2][0] = a20; m[2][1] = a21; m[2][2] = a22;
    }

    /// Return the determinant (Faster than Matrix::det())
    inline Float det() const {
        return ((m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]))
              - (m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]))
              + (m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])));
    }

    /// Matrix-vector multiplication
    inline Vector operator*(const Vector &v) const {
        return Vector(
            m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z);
    }

    /// multiply the matrix by transpose of v, returns the transpose of the line vector.
    inline Vector preMult(const Vector &v) const {
        return Vector(v.x * m[0][0] + v.y * m[1][0] + v.z * m[2][0],
                      v.x * m[0][1] + v.y * m[1][1] + v.z * m[2][1],
                      v.x * m[0][2] + v.y * m[1][2] + v.z * m[2][2]);
    }

    /// Scalar multiplication (creates a temporary)
    inline Matrix3x3 operator*(Float value) const {
        Matrix3x3 result;
        for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
                result.m[i][j] = m[i][j]*value;
        return result;
    }


    /// Assignment operator
    inline Matrix3x3 &operator=(const Matrix<3, 3, Float> &mat) {
        for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
                m[i][j] = mat.m[i][j];
        return *this;
    }

    /// Return a row by index
    inline Vector3 row(int i) const {
        return Vector3(
            m[i][0], m[i][1], m[i][2]
        );
    }

    /// Return a column by index
    inline Vector3 col(int i) const {
        return Vector3(
            m[0][i], m[1][i], m[2][i]
        );
    }
};


/**
 * \brief Basic 4x4 matrix data type
 * \ingroup libcore
 * \ingroup libpython
 */
struct MTS_EXPORT_CORE Matrix4x4 : public Matrix<4, 4, Float> {
    inline Matrix4x4() { }

    /// Initialize the matrix with constant entries
    explicit inline Matrix4x4(Float value) : Matrix<4, 4, Float>(value) { }

    /// Initialize the matrix from a given 4x4 array
    explicit inline Matrix4x4(const Float _m[4][4]) : Matrix<4, 4, Float>(_m) { }

    /// Initialize the matrix from a given (float) 4x4 array in row-major order
    explicit inline Matrix4x4(const Float _m[16]) : Matrix<4, 4, Float>(_m) { }

    /// Initialize the matrix from four 4D column vectors
    explicit inline Matrix4x4(const Vector4 &v1, const Vector4 &v2, const Vector4 &v3, const Vector4 &v4) {
        m[0][0] = v1.x; m[0][1] = v2.x; m[0][2] = v3.x; m[0][3] = v4.x;
        m[1][0] = v1.y; m[1][1] = v2.y; m[1][2] = v3.y; m[1][3] = v4.y;
        m[2][0] = v1.z; m[2][1] = v2.z; m[2][2] = v3.z; m[2][3] = v4.z;
        m[3][0] = v1.w; m[3][1] = v2.w; m[3][2] = v3.w; m[3][3] = v4.w;
    }

    /// Unserialize a matrix from a stream
    explicit inline Matrix4x4(Stream *stream) : Matrix<4, 4, Float>(stream) { }

    /// Copy constructor
    inline Matrix4x4(const Matrix<4, 4, Float> &mtx) : Matrix<4, 4, Float>(mtx) { }

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

    /// Matrix-vector multiplication
    inline Vector4 operator*(const Vector4 &v) const {
        return Vector4(
            m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3] * v.w,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3] * v.w,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3] * v.w,
            m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3] * v.w
        );
    }

    /// Scalar multiplication (creates a temporary)
    inline Matrix4x4 operator*(Float value) const {
        Matrix4x4 result;
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                result.m[i][j] = m[i][j]*value;
        return result;
    }

    /// Assignment operator
    inline Matrix4x4 &operator=(const Matrix<4, 4, Float> &mat) {
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                m[i][j] = mat.m[i][j];
        return *this;
    }

    /// Return a row by index
    inline Vector4 row(int i) const {
        return Vector4(
            m[i][0], m[i][1], m[i][2], m[i][3]
        );
    }

    /// Return a column by index
    inline Vector4 col(int i) const {
        return Vector4(
            m[0][i], m[1][i], m[2][i], m[3][i]
        );
    }
};

/// Matrix multiplication (creates a temporary)
template <typename T, int M1, int N1, int M2, int N2> inline Matrix<M1, N2, T>
        operator*(const Matrix<M1, N1, T> &mat1, const Matrix<M2, N2, T> &mat2) {
    BOOST_STATIC_ASSERT(N1 == M2);
    Matrix<M1, N2, T> result;
    for (int i=0; i<M1; ++i) {
        for (int j=0; j<N2; ++j) {
            T sum = 0;
            for (int k=0; k<N1; ++k)
                sum += mat1.m[i][k] * mat2.m[k][j];
            result.m[i][j] = sum;
        }
    }
    return result;
}

template <typename T, int M, int N> inline Matrix<M, N, T> operator*(T f, const Matrix<M, N, T> &m) {
    return m*f;
}

/**
 * \brief Fast 3x3 eigenvalue decomposition
 *
 * \param m
 *    Matrix in question -- will be replaced with the eigenvectors
 * \param lambda
 *    Parameter used to returns the eigenvalues
 * \return
 *    \c true upon success.
 */
extern MTS_EXPORT_CORE bool eig3(Matrix3x3 &m, Float lambda[3]);

/**
 * \brief Fast non-iterative 3x3 eigenvalue decomposition
 *
 * \param m
 *    Matrix in question -- will be replaced with the eigenvectors
 * \param lambda
 *    Parameter used to returns the eigenvalues
 */
extern MTS_EXPORT_CORE void eig3_noniter(Matrix3x3 &m, Float lambda[3]);

MTS_NAMESPACE_END

#include "matrix.inl"

#endif /* __MITSUBA_CORE_MATRIX_H_ */

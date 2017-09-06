MTS_NAMESPACE_BEGIN

/* Implementations are based on the public domain JAMA library */

template <int M, int N, typename T> bool Matrix<M, N, T>::chol(Matrix &L) const {
    BOOST_STATIC_ASSERT(M == N);

    for (int j = 0; j < N; j++) {
        T *Lrowj = L.m[j], d = 0;
        for (int k = 0; k < j; k++) {
            T *Lrowk = L.m[k];
            T s = 0;
            for (int i = 0; i < k; i++)
                s += Lrowk[i]*Lrowj[i];
            Lrowj[k] = s = (m[j][k] - s)/L.m[k][k];
            d = d + s*s;
            if (m[k][j] != m[j][k])
                return false;
        }
        d = m[j][j] - d;
        if (d <= 0)
            return false;
        L.m[j][j] = std::sqrt(std::max(d, (T) 0));
        for (int k = j+1; k < N; k++)
            L.m[j][k] = 0;
    }

    return true;
}

template <int M, int N, typename T> template <int K> void Matrix<M, N, T>::cholSolve(const Matrix<M, K, T> &B,
            Matrix<M, K, T> &X) const {
    BOOST_STATIC_ASSERT(M == N);

    memcpy(X.m, B.m, sizeof(T)*M*K);

    // Solve L*Y = B;
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < K; j++) {
            for (int i = 0; i < k ; i++)
                X.m[k][j] -= X.m[i][j]*m[k][i];
            X.m[k][j] /= m[k][k];
        }
    }

    // Solve L'*X = Y;
    for (int k = N-1; k >= 0; k--) {
        for (int j = 0; j < K; j++) {
            for (int i = k+1; i < N ; i++)
                X.m[k][j] -= X.m[i][j]*m[i][k];
            X.m[k][j] /= m[k][k];
        }
    }
}

template <int M, int N, typename T> bool Matrix<M, N, T>::lu(Matrix &LU,
        int piv[M], int &pivsign) const {
    LU = *this;

    // Computes L and U with the "daxpy"-based elimination algorithm
    for (int i = 0; i < M; i++)
        piv[i] = i;
    pivsign = 1;

    // Main loop.
    for (int k = 0; k < N; k++) {
        // Find pivot.
        int p = k;
        for (int i = k+1; i < M; i++)
            if (std::abs(LU.m[i][k]) > std::abs(LU.m[p][k]))
               p = i;

         // Exchange if necessary.
        if (p != k) {
            for (int j = 0; j < N; j++)
                std::swap(LU.m[p][j], LU.m[k][j]);
            std::swap(piv[p], piv[k]);
            pivsign = -pivsign;
        }

        // Compute multipliers and eliminate k-th column.
        if (LU.m[k][k] != 0) {
            for (int i = k+1; i < M; i++) {
                LU.m[i][k] /= LU.m[k][k];
                for (int j = k+1; j < N; j++)
                    LU.m[i][j] -= LU.m[i][k]*LU.m[k][j];
            }
        }
    }
    for (int j = 0; j < N; j++)
        if (LU.m[j][j] == 0)
            return false;
    return true;
}

template <int M, int N, typename T> template <int K> void Matrix<M, N, T>::luSolve(const Matrix<M, K, T> &B,
            Matrix<M, K, T> &X, int piv[M]) const {
    BOOST_STATIC_ASSERT(M == N);

    // Copy right hand side with pivoting
    for (int i=0; i<M; ++i)
        for (int j=0; j<K; ++j)
            X.m[i][j] = B.m[piv[i]][j];

    // Solve L*Y = B(piv,:)
    for (int k = 0; k < N; k++)
        for (int i = k+1; i < N; i++)
            for (int j = 0; j < K; j++)
                X.m[i][j] -= X.m[k][j]*m[i][k];

    // Solve U*X = Y;
    for (int k = N-1; k >= 0; k--) {
        for (int j = 0; j < K; j++)
            X.m[k][j] /= m[k][k];

        for (int i = 0; i < k; i++)
            for (int j = 0; j < K; j++)
                X.m[i][j] -= X.m[k][j]*m[i][k];
    }
}

template <int M, int N, typename T> T Matrix<M, N, T>::luDet(int pivsign) const {
    BOOST_STATIC_ASSERT(M == N);
    T d = (T) pivsign;
    for (int j = 0; j < N; j++)
        d *= m[j][j];
    return d;
}

template <int M, int N, typename T> T Matrix<M, N, T>::cholDet() const {
    BOOST_STATIC_ASSERT(M == N);
    T d = m[0][0];
    for (int j = 1; j < N; j++)
        d *= m[j][j];
    return d*d;
}

template <int M, int N, typename T> bool Matrix<M, N, T>::invert(Matrix &target) const {
    BOOST_STATIC_ASSERT(M == N);

    int indxc[N], indxr[N];
    int ipiv[N];
    memset(ipiv, 0, sizeof(int)*N);
    memcpy(target.m, m, M*N*sizeof(T));

    for (int i = 0; i < N; i++) {
        int irow = -1, icol = -1;
        T big = 0;
        for (int j = 0; j < N; j++) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < N; k++) {
                    if (ipiv[k] == 0) {
                        if (std::abs(target.m[j][k]) >= big) {
                            big = std::abs(target.m[j][k]);
                            irow = j;
                            icol = k;
                        }
                    } else if (ipiv[k] > 1) {
                        return false;
                    }
                }
            }
        }
        ++ipiv[icol];
        if (irow != icol) {
            for (int k = 0; k < N; ++k)
                std::swap(target.m[irow][k], target.m[icol][k]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (target.m[icol][icol] == 0)
            return false;
        T pivinv = 1.f / target.m[icol][icol];
        target.m[icol][icol] = 1.f;
        for (int j = 0; j < N; j++)
            target.m[icol][j] *= pivinv;
        for (int j = 0; j < N; j++) {
            if (j != icol) {
                T save = target.m[j][icol];
                target.m[j][icol] = 0;
                for (int k = 0; k < N; k++)
                    target.m[j][k] -= target.m[icol][k]*save;
            }
        }
    }
    for (int j = N-1; j >= 0; j--) {
        if (indxr[j] != indxc[j]) {
            for (int k = 0; k < N; k++)
                std::swap(target.m[k][indxr[j]], target.m[k][indxc[j]]);
        }
    }
    return true;
}

// Symmetric Householder reduction to tridiagonal form.
template <int M, int N, typename T> void
    Matrix<M, N, T>::tred2(T V[M][N], T d[N], T e[N]) {
    //  This is derived from the Algol procedures tred2 by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.

    for (int j = 0; j < N; j++)
        d[j] = V[N - 1][j];

    // Householder reduction to tridiagonal form.
    for (int i = N - 1; i > 0; i--) {
        // Scale to avoid under/overflow.
        T scale = 0.0f, h = 0.0f;

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
            T f = d[i - 1],
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
            T hh = f / (h + h);

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
    for (int i = 0; i < N - 1; i++) {
        V[N - 1][i] = V[i][i];
        V[i][i] = 1.0f;
        T h = d[i + 1];

        if (h != 0.0f) {
            for (int k = 0; k <= i; k++)
                d[k] = V[k][i + 1] / h;
            for (int j = 0; j <= i; j++) {
                T g = 0.0f;

                for (int k = 0; k <= i; k++)
                    g += V[k][i + 1] * V[k][j];
                for (int k = 0; k <= i; k++)
                    V[k][j] -= g * d[k];
            }
        }
        for (int k = 0; k <= i; k++)
            V[k][i + 1] = 0.0f;
    }
    for (int j = 0; j < N; j++) {
        d[j] = V[N - 1][j];
        V[N - 1][j] = 0.0f;
    }
    V[N - 1][N - 1] = 1.0f;
    e[0] = 0.0;
}

// Symmetric tridiagonal QL algorithm.
template <int M, int N, typename T> void
    Matrix<M, N, T>::tql2(T V[M][N], T d[N], T e[N]) {
    //  This is derived from the Algol procedures tql2, by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.

    for (int i = 1; i < N; i++)
        e[i - 1] = e[i];

    e[N - 1] = 0.0f;
    T f = 0.0f, tst1 = 0.0f;

#if defined(SINGLE_PRECISION)
    T eps = pow(2.0f, -23.0f);
#else
    T eps = pow(2.0, -52.0);
#endif

    for (int l = 0; l < N; l++) {
        // Find small subdiagonal element
        tst1 = std::max(tst1, std::abs(d[l]) + std::abs(e[l]));
        int m = l;

        while (m < N) {
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
                T g = d[l];
                T p = (d[l + 1] - g) / (2.0f * e[l]);
                T r = math::hypot2((T) 1, p);

                if (p < 0)
                    r = -r;
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                T dl1 = d[l + 1];

                T h = g - d[l];

                for (int i = l + 2; i < N; i++)
                    d[i] -= h;
                f = f + h;

                // Implicit QL transformation.
                p = d[m];
                T c = 1.0f;
                T c2 = c, c3 = c;
                T el1 = e[l + 1];
                T s = 0.0f, s2 = 0.0f;

                for (int i = m - 1; i >= l; i--) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = math::hypot2(p, e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);

                    // Accumulate transformation.
                    for (int k = 0; k < N; k++) {
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
    for (int i = 0; i < N - 1; i++) {
        int k = i;

        T p = d[i];

        for (int j = i + 1; j < N; j++) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }

        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (int j = 0; j < N; j++) {
                p = V[j][i];
                V[j][i] = V[j][k];
                V[j][k] = p;
            }
        }
    }
}

MTS_NAMESPACE_END

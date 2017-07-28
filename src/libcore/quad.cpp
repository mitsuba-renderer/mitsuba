#include <mitsuba/core/quad.h>
#include <boost/bind.hpp>

MTS_NAMESPACE_BEGIN

float legendreP(int l, float x_) {
    SAssert(l >= 0);

    if (l == 0) {
        return (float) 1.0f;
    } else if (l == 1) {
        return x_;
    } else {
        /* Evaluate the recurrence in double precision */
        double x = (double) x_;
        double Lppred = 1.0, Lpred = x, Lcur = 0.0;

        for (int k = 2; k <= l; ++k) {
            Lcur = ((2*k-1) * x * Lpred - (k - 1) * Lppred) / k;
            Lppred = Lpred; Lpred = Lcur;
        }

        return (float) Lcur;
    }
}

double legendreP(int l, double x) {
    SAssert(l >= 0);

    if (l == 0) {
        return (double) 1.0f;
    } else if (l == 1) {
        return x;
    } else {
        double Lppred = 1.0, Lpred = x, Lcur = 0.0;

        for (int k = 2; k <= l; ++k) {
            Lcur = ((2*k-1) * x * Lpred - (k - 1) * Lppred) / k;
            Lppred = Lpred; Lpred = Lcur;
        }

        return Lcur;
    }
}

std::pair<float, float> legendrePD(int l, float x_) {
    SAssert(l >= 0);

    if (l == 0) {
        return std::make_pair((float) 1.0f, (float) 0.0f);
    } else if (l == 1) {
        return std::make_pair(x_, (float) 1.0f);
    } else {
        /* Evaluate the recurrence in double precision */
        double x = (double) x_;
        double Lppred = 1.0, Lpred = x, Lcur = 0.0,
               Dppred = 0.0, Dpred = 1.0, Dcur = 0.0;

        for (int k = 2; k <= l; ++k) {
            Lcur = ((2*k-1) * x * Lpred - (k - 1) * Lppred) / k;
            Dcur = Dppred + (2*k-1) * Lpred;
            Lppred = Lpred; Lpred = Lcur;
            Dppred = Dpred; Dpred = Dcur;
        }

        return std::make_pair((float) Lcur, (float) Dcur);
    }
}

std::pair<double, double> legendrePD(int l, double x) {
    SAssert(l >= 0);

    if (l == 0) {
        return std::make_pair(1.0, 0.0);
    } else if (l == 1) {
        return std::make_pair(x, 1.0);
    } else {
        double Lppred = 1.0, Lpred = x, Lcur = 0.0,
               Dppred = 0.0, Dpred = 1.0, Dcur = 0.0;

        for (int k = 2; k <= l; ++k) {
            Lcur = ((2*k-1) * x * Lpred - (k - 1) * Lppred) / k;
            Dcur = Dppred + (2*k-1) * Lpred;
            Lppred = Lpred; Lpred = Lcur;
            Dppred = Dpred; Dpred = Dcur;
        }

        return std::make_pair(Lcur, Dcur);
    }
}

/// Evaluate the function legendrePD(l+1, x) - legendrePD(l-1, x)
static std::pair<double, double> legendreQ(int l, double x) {
    SAssert(l >= 1);

    if (l == 1) {
        return std::make_pair(0.5 * (3*x*x-1) - 1, 3*x);
    } else {
        /* Evaluate the recurrence in double precision */
        double Lppred = 1.0, Lpred = x, Lcur = 0.0,
               Dppred = 0.0, Dpred = 1.0, Dcur = 0.0;

        for (int k = 2; k <= l; ++k) {
            Lcur = ((2*k-1) * x * Lpred - (k-1) * Lppred) / k;
            Dcur = Dppred + (2*k-1) * Lpred;
            Lppred = Lpred; Lpred = Lcur;
            Dppred = Dpred; Dpred = Dcur;
        }

        double Lnext = ((2*l+1) * x * Lpred - l * Lppred) / (l+1);
        double Dnext = Dppred + (2*l+1) * Lpred;

        return std::make_pair(Lnext - Lppred, Dnext - Dppred);
    }
}

double legendreP(int l, int m, double x) {
    double p_mm = 1;

    if (m > 0) {
        double somx2 = std::sqrt((1 - x) * (1 + x));
        double fact = 1;
        for (int i=1; i<=m; i++) {
            p_mm *= (-fact) * somx2;
            fact += 2;
        }
    }

    if (l == m)
        return p_mm;

    double p_mmp1 = x * (2*m + 1) * p_mm;
    if (l == m+1)
        return p_mmp1;

    double p_ll = 0;
    for (int ll=m+2; ll <= l; ++ll) {
        p_ll = ((2*ll-1)*x*p_mmp1 - (ll+m-1) * p_mm) / (ll-m);
        p_mm = p_mmp1;
        p_mmp1 = p_ll;
    }

    return p_ll;
}

float legendreP(int l, int m, float x) {
    /* Evaluate the recurrence in double precision */
    double p_mm = 1;

    if (m > 0) {
        double somx2 = std::sqrt((1 - x) * (1 + x));
        double fact = 1;
        for (int i=1; i<=m; i++) {
            p_mm *= (-fact) * somx2;
            fact += 2;
        }
    }

    if (l == m)
        return (float) p_mm;

    double p_mmp1 = x * (2*m + 1) * p_mm;
    if (l == m+1)
        return (float) p_mmp1;

    double p_ll = 0;
    for (int ll=m+2; ll <= l; ++ll) {
        p_ll = ((2*ll-1)*x*p_mmp1 - (ll+m-1) * p_mm) / (ll-m);
        p_mm = p_mmp1;
        p_mmp1 = p_ll;
    }

    return (float) p_ll;
}

void gaussLegendre(int n, Float *nodes, Float *weights) {
    if (n-- < 1)
        SLog(EError, "gaussLegendre(): n must be >= 1");

    if (n == 0) {
        nodes[0] = 0;
        weights[0] = 2;
    } else if (n == 1) {
        nodes[0] = (Float) -std::sqrt(1.0/3.0);
        nodes[1] = -nodes[0];
        weights[0] = weights[1] = 1;
    }

    int m = (n+1)/2;
    for (int i=0; i<m; ++i) {
        /* Initial guess for this root using that of a Chebyshev polynomial */

        double x = -std::cos((double) (2*i + 1) / (double) (2*n + 2) * M_PI);
        int it = 0;

        while (true) {
            if (++it > 20)
                SLog(EError, "gaussLegendre(%i): did not converge after 20 iterations!", n);

            /* Search for the interior roots of P_{n+1}(x) using Newton's method. */
            std::pair<double, double> L = legendrePD(n+1, x);
            double step = L.first / L.second;
            x -= step;

            if (std::abs(step) <= 4 * std::abs(x) * std::numeric_limits<double>::epsilon())
                break;
        }

        std::pair<double, double> L = legendrePD(n+1, x);
        weights[i] = weights[n-i] = (Float) (2.0 / ((1-x*x) * (L.second*L.second)));
        nodes[i] = (Float) x; nodes[n-i] = (Float) -x;
        SAssert(i == 0 || x > nodes[i-1]);
    }

    if ((n % 2) == 0) {
        std::pair<double, double> L = legendrePD(n+1, 0.0);
        weights[n/2] = (Float) (2.0 / (L.second*L.second));
        nodes[n/2] = 0;
    }
}

void gaussLobatto(int n, Float *nodes, Float *weights) {
    if (n-- < 2)
        SLog(EError, "gaussLobatto(): n must be >= 2");

    nodes[0] = -1;
    nodes[n] =  1;
    weights[0] = weights[n] = (Float) 2 / (Float) (n * (n+1));

    int m = (n+1)/2;
    for (int i=1; i<m; ++i) {
        /* Initial guess for this root -- see "On the Legendre-Gauss-Lobatto Points
           and Weights" by Seymor V. Parter, Journal of Sci. Comp., Vol. 14, 4, 1999 */

        double x = -std::cos((i + 0.25) * M_PI / n - 3/(8*n*M_PI * (i + 0.25)));
        int it = 0;

        while (true) {
            if (++it > 20)
                SLog(EError, "gaussLobatto(%i): did not converge after 20 iterations!", n);

            /* Search for the interior roots of P_n'(x) using Newton's method. The same
               roots are also shared by P_{n+1}-P_{n-1}, which is nicer to evaluate. */

            std::pair<double, double> Q = legendreQ(n, x);
            double step = Q.first / Q.second;
            x -= step;

            if (std::abs(step) <= 4 * std::abs(x) * std::numeric_limits<double>::epsilon())
                break;
        }

        double Ln = legendreP(n, x);
        weights[i] = weights[n-i] = (Float) (2.0 / ((n * (n+1)) * Ln * Ln));
        nodes[i] = (Float) x; nodes[n-i] = (Float) -x;
        SAssert(x > nodes[i-1]);
    }

    if ((n % 2) == 0) {
        double Ln = legendreP(n, 0.0);
        weights[n/2] = (Float) (2.0 / ((n * (n+1)) * Ln * Ln));
        nodes[n/2] = 0.0;
    }
}


/*!
 \brief integral of a one-dimensional function using an adaptive
 Gauss-Lobatto integral

 Copyright (C) 2008 Klaus Spanderen

 This code is based on code in QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE. See the license for more details.
*/

const Float GaussLobattoIntegrator::m_alpha = (Float) std::sqrt(2.0/3.0);
const Float GaussLobattoIntegrator::m_beta  = (Float) (1.0/std::sqrt(5.0));
const Float GaussLobattoIntegrator::m_x1    = (Float) 0.94288241569547971906;
const Float GaussLobattoIntegrator::m_x2    = (Float) 0.64185334234578130578;
const Float GaussLobattoIntegrator::m_x3    = (Float) 0.23638319966214988028;

GaussLobattoIntegrator::GaussLobattoIntegrator(size_t maxEvals,
    Float absError, Float relError, bool useConvergenceEstimate, bool warn)
    : m_absError(absError),
      m_relError(relError),
      m_maxEvals(maxEvals),
      m_useConvergenceEstimate(useConvergenceEstimate),
      m_warn(warn) {
    if (m_absError == 0 && m_relError == 0)
        SLog(EError, "GaussLobattoIntegrator:: Absolute and relative "
            "error requirements can't both be zero!");
}

Float GaussLobattoIntegrator::integrate(
        const boost::function<Float (Float)>& f, Float a, Float b, size_t *_evals) const {
    Float factor = 1;
    size_t evals = 0;
    if (a == b) {
        return 0;
    } else if (b < a) {
        std::swap(a, b);
        factor = -1;
    }
    const Float absTolerance = calculateAbsTolerance(f, a, b, evals);
    evals += 2;
    Float result = factor * adaptiveGaussLobattoStep(f, a, b, f(a), f(b), absTolerance, evals);
    if (evals >= m_maxEvals && m_warn)
        SLog(EWarn, "GaussLobattoIntegrator: Maximum number of evaluations reached!");
    if (_evals)
        *_evals = evals;
    return result;
}

Float GaussLobattoIntegrator::calculateAbsTolerance(
        const boost::function<Float (Float)>& f, Float a, Float b, size_t &evals) const {
    const Float m = (a+b)/2;
    const Float h = (b-a)/2;
    const Float y1 = f(a);
    const Float y3 = f(m-m_alpha*h);
    const Float y5 = f(m-m_beta*h);
    const Float y7 = f(m);
    const Float y9 = f(m+m_beta*h);
    const Float y11= f(m+m_alpha*h);
    const Float y13= f(b);

    Float acc = h*((Float) 0.0158271919734801831*(y1+y13)
                 + (Float) 0.0942738402188500455*(f(m-m_x1*h)+f(m+m_x1*h))
                 + (Float) 0.1550719873365853963*(y3+y11)
                 + (Float) 0.1888215739601824544*(f(m-m_x2*h)+ f(m+m_x2*h))
                 + (Float) 0.1997734052268585268*(y5+y9)
                 + (Float) 0.2249264653333395270*(f(m-m_x3*h)+f(m+m_x3*h))
                 + (Float) 0.2426110719014077338*y7);
    evals += 13;

    Float r = 1.0;
    if (m_useConvergenceEstimate) {
        const Float integral2 = (h/6)*(y1+y13+5*(y5+y9));
        const Float integral1 = (h/1470)*
            (77*(y1+y13) + 432*(y3+y11) + 625*(y5+y9) + 672*y7);

        if (std::abs(integral2-acc) != 0.0)
            r = std::abs(integral1-acc)/std::abs(integral2-acc);
        if (r == 0.0 || r > 1.0)
            r = 1.0;
    }
    Float result = std::numeric_limits<Float>::infinity();

    if (m_relError != 0 && acc != 0)
        result = acc * std::max(m_relError,
            std::numeric_limits<Float>::epsilon())
            / (r*std::numeric_limits<Float>::epsilon());

    if (m_absError != 0)
        result = std::min(result, m_absError
            / (r*std::numeric_limits<Float>::epsilon()));

    return result;
}

Float GaussLobattoIntegrator::adaptiveGaussLobattoStep(
                                 const boost::function<Float (Float)>& f,
                                 Float a, Float b, Float fa, Float fb,
                                 Float acc, size_t &evals) const {
    const Float h=(b-a)/2;
    const Float m=(a+b)/2;

    const Float mll=m-m_alpha*h;
    const Float ml =m-m_beta*h;
    const Float mr =m+m_beta*h;
    const Float mrr=m+m_alpha*h;

    const Float fmll= f(mll);
    const Float fml = f(ml);
    const Float fm  = f(m);
    const Float fmr = f(mr);
    const Float fmrr= f(mrr);

    const Float integral2=(h/6)*(fa+fb+5*(fml+fmr));
    const Float integral1=(h/1470)*(77*(fa+fb)
        + 432*(fmll+fmrr) + 625*(fml+fmr) + 672*fm);

    evals += 5;

    if (evals >= m_maxEvals)
        return integral1;

    Float dist = acc + (integral1-integral2);
    if (dist==acc || mll<=a || b<=mrr) {
        return integral1;
    } else {
        return  adaptiveGaussLobattoStep(f,a,mll,fa,fmll,acc,evals)
              + adaptiveGaussLobattoStep(f,mll,ml,fmll,fml,acc,evals)
              + adaptiveGaussLobattoStep(f,ml,m,fml,fm,acc,evals)
              + adaptiveGaussLobattoStep(f,m,mr,fm,fmr,acc,evals)
              + adaptiveGaussLobattoStep(f,mr,mrr,fmr,fmrr,acc,evals)
              + adaptiveGaussLobattoStep(f,mrr,b,fmrr,fb,acc,evals);
    }
}

/* Adaptive multidimensional integration of a vector of const Integrand &s.
 *
 * Copyright (c) 2005-2010 Steven G. Johnson
 *
 * Portions (see comments) based on HIntLib (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 2002-2005 Rudolf Schuerer.
 *     (http://www.cosy.sbg.ac.at/~rschuer/hintlib/)
 *
 * Portions (see comments) based on GNU GSL (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 1996-2000 Brian Gough.
 *     (http://www.gnu.org/software/gsl/)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* Adaptive multidimensional integration on hypercubes (or, really,
   hyper-rectangles) using cubature rules.

   A cubature rule takes a function and a hypercube and evaluates
   the function at a small number of points, returning an estimate
   of the integral as well as an estimate of the error, and also
   a suggested dimension of the hypercube to subdivide.

   Given such a rule, the adaptive integration is simple:

   1) Evaluate the cubature rule on the hypercube(s).
      Stop if converged.

   2) Pick the hypercube with the largest estimated error,
      and divide it in two along the suggested dimension.

   3) Goto (1).

 The basic algorithm is based on the adaptive cubature described in

     A. C. Genz and A. A. Malik, "An adaptive algorithm for numeric
     integration over an N-dimensional rectangular region,"
     J. Comput. Appl. Math. 6 (4), 295-302 (1980).

 and subsequently extended to integrating a vector of const Integrand &s in

     J. Berntsen, T. O. Espelid, and A. Genz, "An adaptive algorithm
     for the approximate calculation of multiple integrals,"
     ACM Trans. Math. Soft. 17 (4), 437-451 (1991).

 Note, however, that we do not use any of code from the above authors
 (in part because their code is Fortran 77, but mostly because it is
 under the restrictive ACM copyright license).  I did make use of some
 GPL code from Rudolf Schuerer's HIntLib and from the GNU Scientific
 Library as listed in the copyright notice above, on the other hand.

 I am also grateful to Dmitry Turbiner <dturbiner@alum.mit.edu>, who
 implemented an initial prototype of the "vectorized" functionality
 for evaluating multiple points in a single call (as opposed to
 multiple functions in a single call).  (Although Dmitry implemented
 a working version, I ended up re-implementing this feature from
 scratch as part of a larger code-cleanup, and in order to have
 a single code path for the vectorized and non-vectorized APIs.  I
 subsequently implemented the algorithm by Gladwell to extract
 even more parallelism by evalutating many hypercubes at once.)
*/

/***************************************************************************/
/* Basic datatypes */

typedef NDIntegrator::VectorizedIntegrand VectorizedIntegrand;

typedef struct {
    Float val, err;
} esterr;

static Float relError(esterr ee) {
    return (ee.val == 0 ? std::numeric_limits<Float>::infinity() :
        std::abs(ee.err / ee.val));
}

static Float errMax(unsigned int fdim, const esterr *ee) {
    Float errmax = 0;
    unsigned int k;
    for (k = 0; k < fdim; ++k)
        if (ee[k].err > errmax) errmax = ee[k].err;
    return errmax;
}

typedef struct {
    unsigned int dim;
    Float *data;    /* length 2*dim = center followed by half-widths */
    Float vol;  /* cache volume = product of widths */
} hypercube;

static Float compute_vol(const hypercube *h) {
    unsigned int i;
    Float vol = 1;
    for (i = 0; i < h->dim; ++i)
        vol *= 2 * h->data[i + h->dim];
    return vol;
}

static hypercube make_hypercube(unsigned int dim, const Float *center, const Float *halfwidth) {
    unsigned int i;
    hypercube h;
    h.dim = dim;
    h.data = (Float *) malloc(sizeof(Float) * dim * 2);
    h.vol = 0;
    if (h.data) {
        for (i = 0; i < dim; ++i) {
            h.data[i] = center[i];
            h.data[i + dim] = halfwidth[i];
        }
        h.vol = compute_vol(&h);
    }
    return h;
}

static hypercube make_hypercube_range(unsigned int dim, const Float *xmin, const Float *xmax) {
    hypercube h = make_hypercube(dim, xmin, xmax);
    unsigned int i;
    if (h.data) {
        for (i = 0; i < dim; ++i) {
            h.data[i] = 0.5f * (xmin[i] + xmax[i]);
            h.data[i + dim] = 0.5f * (xmax[i] - xmin[i]);
        }
        h.vol = compute_vol(&h);
    }
    return h;
}

static void destroy_hypercube(hypercube *h) {
    free(h->data);
    h->dim = 0;
}

typedef struct {
    hypercube h;
    unsigned int splitDim;
    unsigned int fdim; /* dimensionality of vector const Integrand & */
    esterr *ee; /* array of length fdim */
    Float errmax; /* max ee[k].err */
} region;

static region make_region(const hypercube *h, unsigned int fdim) {
    region R;
    R.h = make_hypercube(h->dim, h->data, h->data + h->dim);
    R.splitDim = 0;
    R.fdim = fdim;
    R.ee = R.h.data ? (esterr *) malloc(sizeof(esterr) * fdim) : NULL;
    return R;
}

static void destroy_region(region *R) {
    destroy_hypercube(&R->h);
    free(R->ee);
    R->ee = 0;
}

static bool cut_region(region *R, region *R2) {
    unsigned int d = R->splitDim, dim = R->h.dim;
    *R2 = *R;
    R->h.data[d + dim] *= 0.5f;
    R->h.vol *= 0.5f;
    R2->h = make_hypercube(dim, R->h.data, R->h.data + dim);
    if (!R2->h.data)
        return NDIntegrator::EFailure;
    R->h.data[d] -= R->h.data[d + dim];
    R2->h.data[d] += R->h.data[d + dim];
    R2->ee = (esterr *) malloc(sizeof(esterr) * R2->fdim);
    return R2->ee == NULL;
}

struct rule_s; /* forward declaration */

typedef NDIntegrator::EResult (*evalError_func)(struct rule_s *r,
                  unsigned int fdim, const VectorizedIntegrand &f,
                  unsigned int nR, region *R);
typedef void (*destroy_func)(struct rule_s *r);

typedef struct rule_s {
    unsigned int dim, fdim;         /* the dimensionality & number of functions */
    unsigned int num_points;       /* number of evaluation points */
    unsigned int num_regions; /* max number of regions evaluated at once */
    Float *pts; /* points to eval: num_regions * num_points * dim */
    Float *vals; /* num_regions * num_points * fdim */
    evalError_func evalError;
    destroy_func destroy;
} rule;

static void destroy_rule(rule *r) {
    if (r) {
        if (r->destroy)
            r->destroy(r);
        free(r->pts);
        free(r);
    }
}

static NDIntegrator::EResult alloc_rule_pts(rule *r, unsigned int num_regions) {
    if (num_regions > r->num_regions) {
        free(r->pts);
        r->pts = r->vals = NULL;
        r->num_regions = 0;
        /* allocate extra so that repeatedly calling alloc_rule_pts with
           growing num_regions only needs a logarithmic number of allocations */
        num_regions *= 2;
        r->pts = (Float *) malloc(sizeof(Float) *
                 (num_regions * r->num_points * (r->dim + r->fdim)));
        if (r->fdim + r->dim > 0 && !r->pts)
            return NDIntegrator::EFailure;
        r->vals = r->pts + num_regions * r->num_points * r->dim;
        r->num_regions = num_regions;
    }
    return NDIntegrator::ESuccess;
}

static rule *make_rule(size_t sz, /* >= sizeof(rule) */
               unsigned int dim, unsigned int fdim, unsigned int num_points,
               evalError_func evalError, destroy_func destroy) {
    rule *r;

    if (sz < sizeof(rule))
        return NULL;
    r = (rule *) malloc(sz);
    if (!r)
        return NULL;
    r->pts = r->vals = NULL;
    r->num_regions = 0;
    r->dim = dim; r->fdim = fdim;
    r->num_points = num_points;
    r->evalError = evalError;
    r->destroy = destroy;
    return r;
}

/* note: all regions must have same fdim */
static int eval_regions(unsigned int nR, region *R,
            const VectorizedIntegrand & f, rule *r)
{
    unsigned int iR;
    if (nR == 0)
        return NDIntegrator::ESuccess; /* nothing to evaluate */
    if (r->evalError(r, R->fdim, f, nR, R))
        return NDIntegrator::EFailure;
    for (iR = 0; iR < nR; ++iR)
        R[iR].errmax = errMax(R->fdim, R[iR].ee);
    return NDIntegrator::ESuccess;
}

/***************************************************************************/
/* Functions to loop over points in a hypercube. */

/* Based on orbitrule.cpp in HIntLib-0.0.10 */

/* ls0 returns the least-significant 0 bit of n (e.g. it returns
   0 if the LSB is 0, it returns 1 if the 2 LSBs are 01, etcetera). */
static unsigned int ls0(unsigned int n)
{
#if defined(__GNUC__) && \
    ((__GNUC__ == 3 && __GNUC_MINOR__ >= 4) || __GNUC__ > 3)
    return __builtin_ctz(~n); /* gcc builtin for version >= 3.4 */
#else
    const unsigned int bits[256] = {
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6,
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7,
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6,
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
        0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 8,
    };
    unsigned int bit = 0;
    while ((n & 0xff) == 0xff) {
        n >>= 8;
        bit += 8;
    }
    return bit + bits[n & 0xff];
#endif
}

/**
 *  Evaluate the integration points for all 2^n points (+/-r,...+/-r)
 *
 *  A Gray-code ordering is used to minimize the number of coordinate updates
 *  in p, although this doesn't matter as much now that we are saving all pts.
 */
static void evalR_Rfs(Float *pts, unsigned int dim, Float *p, const Float *c, const Float *r) {
    unsigned int signs = 0; /* 0/1 bit = +/- for corresponding element of r[] */

    /* We start with the point where r is ADDed in every coordinate
       (this implies signs=0). */
    for (unsigned int i = 0; i < dim; ++i)
        p[i] = c[i] + r[i];

    /* Loop through the points in Gray-code ordering */
    for (unsigned i = 0;; ++i) {
        unsigned int mask, d;
        memcpy(pts, p, sizeof(Float) * dim); pts += dim;
        d = ls0(i); /* which coordinate to flip */
        if (d >= dim)
            break;

        /* flip the d-th bit and add/subtract r[d] */
        mask = 1U << d;
        signs ^= mask;
        p[d] = (signs & mask) ? c[d] - r[d] : c[d] + r[d];
    }
}

static void evalRR0_0fs(Float *pts, unsigned int dim, Float *p, const Float *c, const Float *r) {
    for (unsigned i = 0; i < dim - 1; ++i) {
        p[i] = c[i] - r[i];
        for (unsigned j = i + 1; j < dim; ++j) {
            p[j] = c[j] - r[j];
            memcpy(pts, p, sizeof(Float) * dim); pts += dim;
            p[i] = c[i] + r[i];
            memcpy(pts, p, sizeof(Float) * dim); pts += dim;
            p[j] = c[j] + r[j];
            memcpy(pts, p, sizeof(Float) * dim); pts += dim;
            p[i] = c[i] - r[i];
            memcpy(pts, p, sizeof(Float) * dim); pts += dim;
            p[j] = c[j];    /* Done with j -> Restore p[j] */
        }
        p[i] = c[i];        /* Done with i -> Restore p[i] */
    }
}

static void evalR0_0fs4d(Float *pts, unsigned int dim, Float *p, const Float *c,
             const Float *r1, const Float *r2) {
    memcpy(pts, p, sizeof(Float) * dim); pts += dim;
    for (unsigned i = 0; i < dim; i++) {
        p[i] = c[i] - r1[i];
        memcpy(pts, p, sizeof(Float) * dim); pts += dim;
        p[i] = c[i] + r1[i];
        memcpy(pts, p, sizeof(Float) * dim); pts += dim;
        p[i] = c[i] - r2[i];
        memcpy(pts, p, sizeof(Float) * dim); pts += dim;
        p[i] = c[i] + r2[i];
        memcpy(pts, p, sizeof(Float) * dim); pts += dim;
        p[i] = c[i];
    }
}

#define num0_0(dim) (1U)
#define numR0_0fs(dim) (2 * (dim))
#define numRR0_0fs(dim) (2 * (dim) * (dim-1))
#define numR_Rfs(dim) (1U << (dim))

/***************************************************************************/
/* Based on rule75genzmalik.cpp in HIntLib-0.0.10: An embedded
   cubature rule of degree 7 (embedded rule degree 5) due to A. C. Genz
   and A. A. Malik.  See:

         A. C. Genz and A. A. Malik, "An imbedded [sic] family of fully
         symmetric numerical integration rules," SIAM
         J. Numer. Anal. 20 (3), 580-588 (1983).
*/

typedef struct {
     rule parent;

     /* temporary arrays of length dim */
     Float *widthLambda, *widthLambda2, *p;

     /* dimension-dependent constants */
     Float weight1, weight3, weight5;
     Float weightE1, weightE3;
} rule75genzmalik;

#define real(x) ((Float)(x))
#define to_int(n) ((int)(n))

static int isqr(int x)
{
     return x * x;
}

static void destroy_rule75genzmalik(rule *r_)
{
     rule75genzmalik *r = (rule75genzmalik *) r_;
     free(r->p);
}

static NDIntegrator::EResult rule75genzmalik_evalError(rule *r_, unsigned int fdim, const VectorizedIntegrand &f, unsigned int nR, region *R) {
    /* lambda2 = sqrt(9/70), lambda4 = sqrt(9/10), lambda5 = sqrt(9/19) */
    const Float lambda2 = (Float) 0.3585685828003180919906451539079374954541;
    const Float lambda4 = (Float) 0.9486832980505137995996680633298155601160;
    const Float lambda5 = (Float) 0.6882472016116852977216287342936235251269;
    const Float weight2 = (Float) (980.0 / 6561.0);
    const Float weight4 = (Float) (200.0 / 19683.0);
    const Float weightE2 = (Float) (245.0 / 486.0);
    const Float weightE4 = (Float) (25.0 / 729.0);
    const Float ratio = (lambda2 * lambda2) / (lambda4 * lambda4);

    rule75genzmalik *r = (rule75genzmalik *) r_;
    unsigned int i, j, dim = r_->dim, npts = 0;
    Float *diff, *pts, *vals;

    if (alloc_rule_pts(r_, nR))
        return NDIntegrator::EFailure;
    pts = r_->pts; vals = r_->vals;

    for (unsigned int iR = 0; iR < nR; ++iR) {
        const Float *center = R[iR].h.data;
        const Float *halfwidth = R[iR].h.data + dim;

        for (i = 0; i < dim; ++i)
            r->p[i] = center[i];

        for (i = 0; i < dim; ++i)
            r->widthLambda2[i] = halfwidth[i] * lambda2;
        for (i = 0; i < dim; ++i)
            r->widthLambda[i] = halfwidth[i] * lambda4;

        /* Evaluate points in the center, in (lambda2,0,...,0) and
            (lambda3=lambda4, 0,...,0).  */
        evalR0_0fs4d(pts + npts*dim, dim, r->p, center,
            r->widthLambda2, r->widthLambda);
        npts += num0_0(dim) + 2 * numR0_0fs(dim);

        /* Calculate points for (lambda4, lambda4, 0, ...,0) */
        evalRR0_0fs(pts + npts*dim, dim, r->p, center, r->widthLambda);
        npts += numRR0_0fs(dim);

        /* Calculate points for (lambda5, lambda5, ..., lambda5) */
        for (i = 0; i < dim; ++i)
            r->widthLambda[i] = halfwidth[i] * lambda5;
        evalR_Rfs(pts + npts*dim, dim, r->p, center, r->widthLambda);
        npts += numR_Rfs(dim);
    }

    /* Evaluate the const Integrand & function(s) at all the points */
    f((size_t) npts, pts, vals);

    /* we are done with the points, and so we can re-use the pts
       array to store the maximum difference diff[i] in each dimension
       for each hypercube */
    diff = pts;
    for (i = 0; i < dim * nR; ++i)
        diff[i] = 0;

    for (j = 0; j < fdim; ++j) {
        for (unsigned int iR = 0; iR < nR; ++iR) {
            Float result, res5th;
            Float val0, sum2=0, sum3=0, sum4=0, sum5=0;
            unsigned int k, k0 = 0;

            /* accumulate j-th function values into j-th integrals
               NOTE: this relies on the ordering of the eval functions
               above, as well as on the internal structure of
               the evalR0_0fs4d function */

            val0 = vals[0]; /* central point */
            k0 += 1;

            for (k = 0; k < dim; ++k) {
                Float v0 = vals[k0 + 4*k];
                Float v1 = vals[(k0 + 4*k) + 1];
                Float v2 = vals[(k0 + 4*k) + 2];
                Float v3 = vals[(k0 + 4*k) + 3];

                sum2 += v0 + v1;
                sum3 += v2 + v3;

                diff[iR * dim + k] +=
                    std::abs(v0 + v1 - 2*val0 - ratio * (v2 + v3 - 2*val0));
            }
            k0 += 4*k;

            for (k = 0; k < numRR0_0fs(dim); ++k)
                sum4 += vals[k0 + k];
            k0 += k;

            for (k = 0; k < numR_Rfs(dim); ++k)
                sum5 += vals[k0 + k];

            /* Calculate fifth and seventh order results */
            result = R[iR].h.vol * (r->weight1 * val0 + weight2 * sum2 + r->weight3 * sum3 + weight4 * sum4 + r->weight5 * sum5);
            res5th = R[iR].h.vol * (r->weightE1 * val0 + weightE2 * sum2 + r->weightE3 * sum3 + weightE4 * sum4);

            R[iR].ee[j].val = result;
            R[iR].ee[j].err = std::abs(res5th - result);

            vals += r_->num_points;
        }
    }

    /* figure out dimension to split: */
    for (unsigned int iR = 0; iR < nR; ++iR) {
        Float maxdiff = 0;
        unsigned int dimDiffMax = 0;

        for (i = 0; i < dim; ++i) {
            if (diff[iR*dim + i] > maxdiff) {
                maxdiff = diff[iR*dim + i];
                dimDiffMax = i;
            }
        }
        R[iR].splitDim = dimDiffMax;
    }
    return NDIntegrator::ESuccess;
}

static rule *make_rule75genzmalik(unsigned int dim, unsigned int fdim) {
    rule75genzmalik *r;

    if (dim < 2) return NULL; /* this rule does not support 1d integrals */

    /* Because of the use of a bit-field in evalR_Rfs, we are limited
       to be < 32 dimensions (or however many bits are in unsigned).
       This is not a practical limitation...long before you reach
       32 dimensions, the Genz-Malik cubature becomes excruciatingly
       slow and is superseded by other methods (e.g. Monte-Carlo). */
    if (dim >= sizeof(unsigned) * 8)
        return NULL;

    r = (rule75genzmalik *) make_rule(sizeof(rule75genzmalik),
            dim, fdim, num0_0(dim) + 2 * numR0_0fs(dim)
            + numRR0_0fs(dim) + numR_Rfs(dim),
            rule75genzmalik_evalError,
            destroy_rule75genzmalik);
     if (!r)
         return NULL;

    r->weight1 = (real(12824 - 9120 * to_int(dim) + 400 * isqr(to_int(dim))) / real(19683));
    r->weight3 = real(1820 - 400 * to_int(dim)) / real(19683);
    r->weight5 = real(6859) / real(19683) / real(1U << dim);
    r->weightE1 = (real(729 - 950 * to_int(dim) + 50 * isqr(to_int(dim))) / real(729));
    r->weightE3 = real(265 - 100 * to_int(dim)) / real(1458);
    r->p = (Float *) malloc(sizeof(Float) * dim * 3);
    if (!r->p) {
        destroy_rule((rule *) r);
        return NULL;
    }
    r->widthLambda = r->p + dim;
    r->widthLambda2 = r->p + 2 * dim;
    return (rule *) r;
}

/***************************************************************************/
/* 1d 15-point Gaussian quadrature rule, based on qk15.c and qk.c in
   GNU GSL (which in turn is based on QUADPACK). */

static NDIntegrator::EResult rule15gauss_evalError(rule *r,
                 unsigned int fdim, const VectorizedIntegrand & f,
                 unsigned int nR, region *R) {
     /* Gauss quadrature weights and kronrod quadrature abscissae and
        weights as evaluated with 80 decimal digit arithmetic by
        L. W. Fullerton, Bell Labs, Nov. 1981. */
    const unsigned int n = 8;
    const Float xgk[8] = {  /* abscissae of the 15-point kronrod rule */
        (Float) 0.991455371120812639206854697526329,
        (Float) 0.949107912342758524526189684047851,
        (Float) 0.864864423359769072789712788640926,
        (Float) 0.741531185599394439863864773280788,
        (Float) 0.586087235467691130294144838258730,
        (Float) 0.405845151377397166906606412076961,
        (Float) 0.207784955007898467600689403773245,
        (Float) 0.000000000000000000000000000000000
        /* xgk[1], xgk[3], ... abscissae of the 7-point gauss rule.
           xgk[0], xgk[2], ... to optimally extend the 7-point gauss rule */
    };
    static const Float wg[4] = {  /* weights of the 7-point gauss rule */
        (Float) 0.129484966168869693270611432679082,
        (Float) 0.279705391489276667901467771423780,
        (Float) 0.381830050505118944950369775488975,
        (Float) 0.417959183673469387755102040816327
    };
    static const Float wgk[8] = { /* weights of the 15-point kronrod rule */
        (Float) 0.022935322010529224963732008058970,
        (Float) 0.063092092629978553290700663189204,
        (Float) 0.104790010322250183839876322541518,
        (Float) 0.140653259715525918745189590510238,
        (Float) 0.169004726639267902826583426598550,
        (Float) 0.190350578064785409913256402421014,
        (Float) 0.204432940075298892414161999234649,
        (Float) 0.209482141084727828012999174891714
    };
    unsigned int j, npts = 0;
    Float *pts, *vals;

    if (alloc_rule_pts(r, nR))
        return NDIntegrator::EFailure;

    pts = r->pts; vals = r->vals;

    for (unsigned int iR = 0; iR < nR; ++iR) {
        const Float center = R[iR].h.data[0];
        const Float halfwidth = R[iR].h.data[1];

        pts[npts++] = center;

        for (j = 0; j < (n - 1) / 2; ++j) {
            int j2 = 2*j + 1;
            Float w = halfwidth * xgk[j2];
            pts[npts++] = center - w;
            pts[npts++] = center + w;
        }
        for (j = 0; j < n/2; ++j) {
            int j2 = 2*j;
            Float w = halfwidth * xgk[j2];
            pts[npts++] = center - w;
            pts[npts++] = center + w;
        }

        R[iR].splitDim = 0; /* no choice but to divide 0th dimension */
    }

    f((size_t) npts, pts, vals);

    for (unsigned int k = 0; k < fdim; ++k) {
        for (unsigned int iR = 0; iR < nR; ++iR) {
            const Float halfwidth = R[iR].h.data[1];
            Float result_gauss = vals[0] * wg[n/2 - 1];
            Float result_kronrod = vals[0] * wgk[n - 1];
            Float result_abs = std::abs(result_kronrod);
            Float result_asc, mean, err;

            /* accumulate integrals */
            npts = 1;
            for (j = 0; j < (n - 1) / 2; ++j) {
                int j2 = 2*j + 1;
                Float v = vals[npts] + vals[npts+1];
                result_gauss += wg[j] * v;
                result_kronrod += wgk[j2] * v;
                result_abs += wgk[j2] * (std::abs(vals[npts]) + std::abs(vals[npts+1]));
                npts += 2;
            }
            for (j = 0; j < n/2; ++j) {
                int j2 = 2*j;
                result_kronrod += wgk[j2] * (vals[npts] + vals[npts+1]);
                result_abs += wgk[j2] * (std::abs(vals[npts]) + std::abs(vals[npts+1]));
                npts += 2;
            }

            /* integration result */
            R[iR].ee[k].val = result_kronrod * halfwidth;

            /* error estimate (from GSL, probably dates back to QUADPACK
            ... not completely clear to me why we don't just use
            std::abs(result_kronrod - result_gauss) * halfwidth */
            mean = result_kronrod * 0.5f;
            result_asc = wgk[n - 1] * std::abs(vals[0] - mean);
            npts = 1;
            for (j = 0; j < (n - 1) / 2; ++j) {
                int j2 = 2*j + 1;
                result_asc += wgk[j2] * (std::abs(vals[npts]-mean)
                         + std::abs(vals[npts+1]-mean));
                npts += 2;
            }
            for (j = 0; j < n/2; ++j) {
                int j2 = 2*j;
                result_asc += wgk[j2] * (std::abs(vals[npts]-mean)
                         + std::abs(vals[npts+1]-mean));
                npts += 2;
            }
            err = std::abs(result_kronrod - result_gauss) * halfwidth;
            result_abs *= halfwidth;
            result_asc *= halfwidth;
            if (result_asc != 0 && err != 0) {
                /* Recommended error estimate for the 7-15 G-K rule */
                Float scale = std::pow((200 * err / result_asc), (Float) 1.5);
                err = (scale < 1) ? result_asc * scale : result_asc;
            }
            #if 0
                /* This seems a bit excessive (and creates problems for single
                   precision code) */
                if (result_abs > std::numeric_limits<Float>::min() / (50 * std::numeric_limits<Float>::epsilon())) {
                    Float min_err = 50 * std::numeric_limits<Float>::epsilon() * result_abs;
                    if (min_err > err)
                        err = min_err;
                }
            #endif
            R[iR].ee[k].err = err;

            /* increment vals to point to next batch of results */
            vals += 15;
        }
    }
    return NDIntegrator::ESuccess;
}

static rule *make_rule15gauss(unsigned int dim, unsigned int fdim) {
     if (dim != 1) return NULL; /* this rule is only for 1d integrals */

     return make_rule(sizeof(rule), dim, fdim, 15, rule15gauss_evalError, 0);
}

/***************************************************************************/
/* binary heap implementation (ala _Introduction to Algorithms_ by
   Cormen, Leiserson, and Rivest), for use as a priority queue of
   regions to integrate. */

typedef region heap_item;
#define KEY(hi) ((hi).errmax)

typedef struct {
    unsigned int n, nalloc;
    heap_item *items;
    unsigned int fdim;
    esterr *ee; /* array of length fdim of the total const Integrand & & error */
} heap;

static void heap_resize(heap *h, unsigned int nalloc) {
    h->nalloc = nalloc;
    h->items = (heap_item *) realloc(h->items, sizeof(heap_item) * nalloc);
}

static heap heap_alloc(unsigned int nalloc, unsigned int fdim) {
    heap h;
    unsigned int i;
    h.n = 0;
    h.nalloc = 0;
    h.items = 0;
    h.fdim = fdim;
    h.ee = (esterr *) malloc(sizeof(esterr) * fdim);
    if (h.ee) {
        for (i = 0; i < fdim; ++i)
            h.ee[i].val = h.ee[i].err = 0;
        heap_resize(&h, nalloc);
    }
    return h;
}

/* note that heap_free does not deallocate anything referenced by the items */
static void heap_free(heap *h) {
    h->n = 0;
    heap_resize(h, 0);
    h->fdim = 0;
    free(h->ee);
}

static NDIntegrator::EResult heap_push(heap *h, heap_item hi) {
    int insert;
    unsigned int fdim = h->fdim;

    for (unsigned int i = 0; i < fdim; ++i) {
        h->ee[i].val += hi.ee[i].val;
        h->ee[i].err += hi.ee[i].err;
    }
    insert = h->n;
    if (++(h->n) > h->nalloc) {
        heap_resize(h, h->n * 2);
        if (!h->items)
            return NDIntegrator::EFailure;
    }
    while (insert) {
        int parent = (insert - 1) / 2;
        if (KEY(hi) <= KEY(h->items[parent]))
            break;
        h->items[insert] = h->items[parent];
        insert = parent;
    }
    h->items[insert] = hi;
    return NDIntegrator::ESuccess;
}

static NDIntegrator::EResult heap_push_many(heap *h, unsigned int ni, heap_item *hi) {
     unsigned int i;
     for (i = 0; i < ni; ++i)
      if (heap_push(h, hi[i])) return NDIntegrator::EFailure;
     return NDIntegrator::ESuccess;
}

static heap_item heap_pop(heap *h) {
    heap_item ret;
    int i, n, child;
    if (!(h->n))
        SLog(EError, "attempted to pop an empty heap\n");

    ret = h->items[0];
    h->items[i = 0] = h->items[n = --(h->n)];
    while ((child = i * 2 + 1) < n) {
        int largest;
        heap_item swap;

        if (KEY(h->items[child]) <= KEY(h->items[i]))
            largest = i;
        else
            largest = child;
        if (++child < n && KEY(h->items[largest]) < KEY(h->items[child]))
            largest = child;
        if (largest == i)
            break;
        swap = h->items[i];
        h->items[i] = h->items[largest];
        h->items[i = largest] = swap;
    }
    unsigned int fdim = h->fdim;
    for (unsigned int j = 0; j < fdim; ++j) {
        h->ee[j].val -= ret.ee[j].val;
        h->ee[j].err -= ret.ee[j].err;
    }
    return ret;
}

/***************************************************************************/

/* adaptive integration, analogous to adaptintegrator.cpp in HIntLib */

static NDIntegrator::EResult ruleadapt_integrate(rule *r, unsigned int fdim,
        const VectorizedIntegrand & f, const hypercube *h, size_t maxEval,
        Float reqAbsError, Float reqRelError, Float *val, Float *err, size_t &numEval, int parallel) {
    heap regions;
    unsigned int i, j;
    region *R = NULL; /* array of regions to evaluate */
    unsigned int nR_alloc = 0;
    esterr *ee = NULL;

    regions = heap_alloc(1, fdim);
    if (!regions.ee || !regions.items)
        goto bad;

    ee = (esterr *) malloc(sizeof(esterr) * fdim);
    if (!ee)
        goto bad;

    nR_alloc = 2;
    R = (region *) malloc(sizeof(region) * nR_alloc);
    if (!R)
        goto bad;
    R[0] = make_region(h, fdim);
    if (!R[0].ee || eval_regions(1, R, f, r) || heap_push(&regions, R[0]))
        goto bad;
    numEval += r->num_points;

    while (numEval < maxEval || !maxEval) {
        for (j = 0; j < fdim && (regions.ee[j].err <= reqAbsError ||
            relError(regions.ee[j]) <= reqRelError); ++j)
            ;
        if (j == fdim)
            break; /* convergence */

        if (parallel) {
            /* Maximize potential parallelism

               adapted from I. Gladwell, "Vectorization of one dimensional
               quadrature codes," pp. 230--238 in _Numerical Integration. Recent
               Developments, Software and Applications_, G. Fairweather and
               P. M. Keast, eds., NATO ASI Series C203, Dordrecht (1987), as
               described in J. M. Bull and T. L. Freeman, "Parallel Globally
               Adaptive Algorithms for Multi-dimensional Integration,"
               http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.42.6638

               Basically, this evaluates in one shot all regions
               that *must* be evaluated in order to reduce the
               error to the requested bound: the minimum set of
               largest-error regions whose errors push the total
               error over the bound.

               [Note: Bull and Freeman claim that the Gladwell
               approach is intrinsically inefficent because it
               "requires sorting", and propose an alternative
               algorithm that "only" requires three passes over the
               entire set of regions.  Apparently, they didn't
               realize that one could use a heap data structure, in
               which case the time to pop K biggest-error regions
               out of N is only O(K log N), much better than the
               O(N) cost of the Bull and Freeman algorithm if
               K << N, and it is also much simpler.] */
            unsigned int nR = 0;
            for (j = 0; j < fdim; ++j)
                ee[j] = regions.ee[j];
            do {
                if (nR + 2 > nR_alloc) {
                    nR_alloc = (nR + 2) * 2;
                    R = (region *) realloc(R, nR_alloc * sizeof(region));
                    if (!R)
                        goto bad;
                }
                R[nR] = heap_pop(&regions);
                for (j = 0; j < fdim; ++j)
                    ee[j].err -= R[nR].ee[j].err;
                if (cut_region(R+nR, R+nR+1))
                    goto bad;
                numEval += r->num_points * 2;
                nR += 2;
                for (j = 0; j < fdim && (ee[j].err <= reqAbsError
                    || relError(ee[j]) <= reqRelError); ++j)
                    ;
                if (j == fdim)
                    break; /* other regions have small errs */
            } while (regions.n > 0 && (numEval < maxEval || !maxEval));
            if (eval_regions(nR, R, f, r) || heap_push_many(&regions, nR, R))
                goto bad;
        } else { /* minimize number of function evaluations */
            R[0] = heap_pop(&regions); /* get worst region */
            if (cut_region(R, R+1) || eval_regions(2, R, f, r)
                || heap_push_many(&regions, 2, R))
                goto bad;
            numEval += r->num_points * 2;
        }
    }

     /* re-sum integral and errors */
    for (j = 0; j < fdim; ++j)
        val[j] = err[j] = 0;
    for (i = 0; i < regions.n; ++i) {
        for (j = 0; j < fdim; ++j) {
            val[j] += regions.items[i].ee[j].val;
            err[j] += regions.items[i].ee[j].err;
        }
        destroy_region(&regions.items[i]);
    }

    /* printf("regions.nalloc = %d\n", regions.nalloc); */
    free(ee);
    heap_free(&regions);
    free(R);
    return NDIntegrator::ESuccess;

bad:
    free(ee);
    heap_free(&regions);
    free(R);
    return NDIntegrator::EFailure;
}

static NDIntegrator::EResult integrate(unsigned fdim, const VectorizedIntegrand & f,
             unsigned dim, const Float *xmin, const Float *xmax,
             size_t maxEval, Float reqAbsError, Float reqRelError,
             Float *val, Float *err, size_t &numEval, int parallel) {
    NDIntegrator::EResult status;

    numEval = 0;
    if (fdim == 0) /* nothing to do */
        return NDIntegrator::ESuccess;
    if (dim == 0) { /* trivial integration */
        f(1, xmin, val);
        for (unsigned int i = 0; i < fdim; ++i)
            err[i] = 0;
        return NDIntegrator::ESuccess;
    }
    rule *r = dim == 1 ? make_rule15gauss(dim, fdim)
        : make_rule75genzmalik(dim, fdim);
    if (!r) {
        for (unsigned int i = 0; i < fdim; ++i) {
            val[i] = 0;
            err[i] = std::numeric_limits<Float>::infinity();
        }
        return NDIntegrator::EFailure;
    }
    hypercube h = make_hypercube_range(dim, xmin, xmax);
    status = !h.data ? NDIntegrator::EFailure
        : ruleadapt_integrate(r, fdim, f, &h,
            maxEval, reqAbsError, reqRelError,
            val, err, numEval, parallel);
    destroy_hypercube(&h);
    destroy_rule(r);
    return status;
}

class VectorizationAdapter {
public:
    VectorizationAdapter(const NDIntegrator::Integrand &integrand, size_t fdim,
            size_t dim) : m_integrand(integrand), m_fdim(fdim), m_dim(dim) {
        m_temp = new Float[m_fdim];
    }

    ~VectorizationAdapter() {
        delete[] m_temp;
    }

    void f(size_t nPt, const Float *in, Float *out) {
        for (size_t i = 0; i < nPt; ++i) {
            m_integrand(in + i*m_dim, m_temp);
            for (size_t k = 0; k < m_fdim; ++k)
                out[k*nPt + i] = m_temp[k];
        }
    }
private:
    const NDIntegrator::Integrand &m_integrand;
    size_t m_fdim, m_dim;
    Float *m_temp;
};

NDIntegrator::NDIntegrator(size_t fDim, size_t dim,
            size_t maxEvals, Float absError, Float relError)
 : m_fdim(fDim), m_dim(dim), m_maxEvals(maxEvals), m_absError(absError),
  m_relError(relError) { }

NDIntegrator::EResult NDIntegrator::integrate(const Integrand &f, const Float *min,
        const Float *max, Float *result, Float *error, size_t *_evals) const {
    VectorizationAdapter adapter(f, m_fdim, m_dim);
    size_t evals = 0;
    EResult retval = mitsuba::integrate((unsigned int) m_fdim, boost::bind(
        &VectorizationAdapter::f, &adapter, _1, _2, _3), (unsigned int) m_dim,
        min, max, m_maxEvals, m_absError, m_relError, result, error, evals, false);
    if (_evals)
        *_evals = evals;
    return retval;
}

NDIntegrator::EResult NDIntegrator::integrateVectorized(const VectorizedIntegrand &f, const Float *min,
        const Float *max, Float *result, Float *error, size_t *_evals) const {
    size_t evals = 0;
    EResult retval = mitsuba::integrate((unsigned int) m_fdim, f, (unsigned int) m_dim,
        min, max, m_maxEvals, m_absError, m_relError, result, error, evals, true);
    if (_evals)
        *_evals = evals;
    return retval;
}

MTS_NAMESPACE_END

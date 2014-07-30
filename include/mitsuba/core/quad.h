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
#if !defined(__MITSUBA_CORE_QUAD_H_)
#define __MITSUBA_CORE_QUAD_H_

#include <mitsuba/mitsuba.h>
#include <boost/function.hpp>

MTS_NAMESPACE_BEGIN

// -----------------------------------------------------------------------
//! \ingroup libcore
//! \ingroup libpython
//! @{ \name Basic tools for non-adaptive quadrature
// -----------------------------------------------------------------------

/// Evaluate the l-th Legendre polynomial using recurrence (single precision)
extern MTS_EXPORT_CORE float legendreP(int l, float x);

/// Evaluate the l-th Legendre polynomial using recurrence (double precision)
extern MTS_EXPORT_CORE double legendreP(int l, double x);

/// Evaluate an associated Legendre polynomial using recurrence (single precision)
extern MTS_EXPORT_CORE float legendreP(int l, int m, float x);

/// Evaluate an associated Legendre polynomial using recurrence (double precision)
extern MTS_EXPORT_CORE double legendreP(int l, int m, double x);

/// Evaluate the l-th Legendre polynomial and its derivative using recurrence (single precision)
extern MTS_EXPORT_CORE std::pair<float, float> legendrePD(int l, float x);

/// Evaluate the l-th Legendre polynomial and its derivative using recurrence (double precision)
extern MTS_EXPORT_CORE std::pair<double, double> legendrePD(int l, double x);

/**
 * \brief Computes the nodes and weights of a Gauss-Legendre quadrature
 * (aka "Gaussian quadrature") rule with the given number of evaluations.
 *
 * Integration is over the interval \f$[-1, 1]\f$. Gauss-Legendre quadrature
 * maximizes the order of exactly integrable polynomials achieves this up to
 * degree \f$2n-1\f$ (where \f$n\f$ is the number of function evaluations).
 *
 * This method is numerically well-behaved until about \f$n=200\f$
 * and then becomes progressively less accurate. It is generally not a
 * good idea to go much higher---in any case, a composite or
 * adaptive integration scheme will be superior for large \f$n\f$.
 *
 * \param n
 *     Desired number of evalution points
 * \param nodes
 *     Length-\c n array used to store the nodes of the quadrature rule
 * \param nodes
 *     Length-\c n array used to store the weights of the quadrature rule
 */
extern MTS_EXPORT_CORE void gaussLegendre(int n, Float *nodes, Float *weights);


/**
 * \brief Computes the nodes and weights of a Gauss-Lobatto quadrature
 * rule with the given number of evaluations.
 *
 * Integration is over the interval \f$[-1, 1]\f$. Gauss-Lobatto quadrature
 * is preferable to Gauss-Legendre quadrature whenever the endpoints of the
 * integration domain should explicitly be included. It maximizes the order
 * of exactly integrable polynomials subject to this constraint and achieves
 * this up to degree \f$2n-3\f$ (where \f$n\f$ is the number of function
 * evaluations).
 *
 * This method is numerically well-behaved until about \f$n=200\f$
 * and then becomes progressively less accurate. It is generally not a
 * good idea to go much higher---in any case, a composite or
 * adaptive integration scheme will be superior for large \f$n\f$.
 *
 * \param n
 *     Desired number of evalution points
 * \param nodes
 *     Length-\c n array used to store the nodes of the quadrature rule
 * \param nodes
 *     Length-\c n array used to store the weights of the quadrature rule
 */
extern MTS_EXPORT_CORE void gaussLobatto(int n, Float *nodes, Float *weights);

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Adaptive quadrature schemes
// -----------------------------------------------------------------------

/**
 * \brief Computes the integral of a one-dimensional function
 * using adaptive Gauss-Lobatto quadrature.
 *
 * Given a target error \f$ \epsilon \f$, the integral of
 * a function \f$ f \f$ between \f$ a \f$ and \f$ b \f$ is
 * calculated by means of the Gauss-Lobatto formula.
 *
 * References:
 * This algorithm is a C++ implementation of the algorithm outlined in
 *
 * W. Gander and W. Gautschi, Adaptive Quadrature - Revisited.
 * BIT, 40(1):84-101, March 2000. CS technical report:
 * ftp.inf.ethz.ch/pub/publications/tech-reports/3xx/306.ps.gz
 *
 * The original MATLAB version can be downloaded here
 * http://www.inf.ethz.ch/personal/gander/adaptlob.m
 *
 * This particular implementation is based on code in QuantLib,
 * a free-software/open-source library for financial quantitative
 * analysts and developers - http://quantlib.org/
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE GaussLobattoIntegrator {
public:
	typedef boost::function<Float (Float)> Integrand;

	/**
	 * Initialize a Gauss-Lobatto integration scheme
	 *
	 * \param maxEvals Maximum number of function evaluations. The
	 *    integrator will print a warning when this limit is
	 *    exceeded. It will then stop the recursion, but a few
	 *    further evaluations may still take place. Hence the limit
	 *    is not a strict one.
	 *
	 * \param absError Absolute error requirement (0 to disable)
	 * \param relError Relative error requirement (0 to disable)
	 *
	 * \param useConvergenceEstimate Estimate the convergence behavior
	 *     of the GL-quadrature by comparing the 4, 7 and 13-point
	 *     variants and increase the absolute tolerance accordingly.
	 *
	 * \param warn Should the integrator warn when the number of
	 *     function evaluations is exceeded?
	 */
	GaussLobattoIntegrator(size_t maxEvals,
						 Float absError = 0,
						 Float relError = 0,
						 bool useConvergenceEstimate = true,
						 bool warn = true);

	/**
	 * \brief Integrate the function \c f from \c a to \c b.
	 *
	 * Also returns the total number of evaluations if requested
	 */
	Float integrate(const Integrand &f, Float a, Float b,
		size_t *evals = NULL) const;
protected:
	/**
	 * \brief Perform one step of the 4-point Gauss-Lobatto rule, then
	 * compute the same integral using a 7-point Kronrod extension and
	 * compare. If the accuracy is deemed too low, recurse.
	 *
	 * \param f Function to integrate
	 * \param a Lower integration limit
	 * \param b Upper integration limit
	 * \param fa Function evaluated at the lower limit
	 * \param fb Function evaluated at the upper limit
	 * \param is Absolute tolerance in epsilons
	 */
	Float adaptiveGaussLobattoStep(const boost::function<Float (Float)>& f,
		Float a, Float b, Float fa, Float fb, Float is, size_t &evals) const;

	/**
	 * Compute the absolute error tolerance using a 13-point
	 * Gauss-Lobatto rule.
	 */
	Float calculateAbsTolerance(const boost::function<Float (Float)>& f,
		Float a, Float b, size_t &evals) const;
protected:
	Float m_absError, m_relError;
	size_t m_maxEvals;
	bool m_useConvergenceEstimate;
	bool m_warn;
	static const Float m_alpha;
	static const Float m_beta;
	static const Float m_x1;
	static const Float m_x2;
	static const Float m_x3;
};

/**
 * \brief Adaptively computes the integral of a multidimensional function using
 * either a Gauss-Kronod (1D) or a Genz-Malik (>1D) cubature rule.
 *
 * This class is a C++ wrapper around the \c cubature code by Steven G. Johnson
 * (http://ab-initio.mit.edu/wiki/index.php/Cubature)
 *
 * The original implementation is based on algorithms proposed in
 *
 * A. C. Genz and A. A. Malik, "An adaptive algorithm for numeric integration
 * over an N-dimensional rectangular region," J. Comput. Appl. Math. 6 (4),
 * 295–302 (1980).
 *
 * and
 *
 * J. Berntsen, T. O. Espelid, and A. Genz, "An adaptive algorithm for the
 * approximate calculation of multiple integrals," ACM Trans. Math. Soft. 17
 * (4), 437–451 (1991).
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE NDIntegrator {
public:
	typedef boost::function<void (const Float *, Float *)>         Integrand;
	typedef boost::function<void (size_t, const Float *, Float *)> VectorizedIntegrand;

	enum EResult {
		ESuccess = 0,
		EFailure = 1
	};

	/**
	 * Initialize the Cubature integration scheme
	 *
	 * \param fDim Number of integrands (i.e. dimensions of the image space)
	 * \param nDim Number of integration dimensions (i.e. dimensions of the
	 *      function domain)
	 * \param maxEvals Maximum number of function evaluations (0 means no
	 *      limit). The error bounds will likely be exceeded when the
	 *      integration is forced to stop prematurely. Note: the actual
	 *      number of evaluations may somewhat exceed this value.
	 * \param absError Absolute error requirement (0 to disable)
	 * \param relError Relative error requirement (0 to disable)
	 */
	NDIntegrator(size_t fDim, size_t dim,
			size_t maxEvals, Float absError = 0, Float relError = 0);

	/**
	 * \brief Integrate the function \c f over the rectangular domain
	 * bounded by \c min and \c max.
	 *
	 * The supplied function should have the interface
	 *
	 * <code>
	 * void integrand(const Float *in, Float *out);
	 * </code>
	 *
	 * The input array \c in consists of one set of input parameters
	 * having \c dim entries. The function is expected to store the
	 * results of the evaluation into the \c out array using \c fDim entries.
	 */
	EResult integrate(const Integrand &f, const Float *min, const Float *max,
			Float *result, Float *error, size_t *evals = NULL) const;

	/**
	 * \brief Integrate the function \c f over the rectangular domain
	 * bounded by \c min and \c max.
	 *
	 * This function implements a vectorized version of the above
	 * integration function, which is more efficient by evaluating
	 * the integrant in `batches'. The supplied function should
	 * have the interface
	 *
	 * <code>
	 * void integrand(int numPoints, const Float *in, Float *out);
	 * </code>
	 *
	 * Note that \c in in is not a single point, but an array of \c numPoints points
	 * (length \c numPoints x \c dim), and upon return the values of all \c fDim
	 * integrands at all \c numPoints points should be stored in \c out
	 * (length \c fDim x \c numPoints). In particular, out[i*dim + j] is the j-th
	 * coordinate of the i-th point, and the k-th function evaluation (k<fDim)
	 * for the i-th point is returned in out[k*npt + i].
	 * The size of \c numPoints will vary with the dimensionality of the problem;
	 * higher-dimensional problems will have (exponentially) larger numbers,
	 * allowing for the possibility of more parallelism. Currently, \c numPoints
	 * starts at 15 in 1d, 17 in 2d, and 33 in 3d, but as the integrator
	 * calls your integrand more and more times the value will grow. e.g. if you end
	 * up requiring several thousand points in total, \c numPoints may grow to
	 * several hundred.
	 */
	EResult integrateVectorized(const VectorizedIntegrand &f, const Float *min,
		const Float *max, Float *result, Float *error, size_t *evals = NULL) const;
protected:
	size_t m_fdim, m_dim, m_maxEvals;
	Float m_absError, m_relError;
};

//! @}
// -----------------------------------------------------------------------

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_QUAD_H_ */

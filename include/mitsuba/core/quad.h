/*
 Copyright (C) 2008 Klaus Spanderen

 This file is based on code in QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file gausslobattointegral.hpp
	\brief integral of a one-dimensional function using the adaptive
	Gauss-Lobatto integral
*/

#if !defined(__QUADRATURE_H)
#define __QUADRATURE_H

#include <mitsuba/mitsuba.h>
#include <boost/function.hpp>

MTS_NAMESPACE_BEGIN

//! Integral of a one-dimensional function
/*! Given a target accuracy \f$ \epsilon \f$, the integral of
	a function \f$ f \f$ between \f$ a \f$ and \f$ b \f$ is
	calculated by means of the Gauss-Lobatto formula
*/

/*! References:
   This algorithm is a C++ implementation of the algorithm outlined in

   W. Gander and W. Gautschi, Adaptive Quadrature - Revisited.
   BIT, 40(1):84-101, March 2000. CS technical report:
   ftp.inf.ethz.ch/pub/publications/tech-reports/3xx/306.ps.gz

   The original MATLAB version can be downloaded here
   http://www.inf.ethz.ch/personal/gander/adaptlob.m
*/

class MTS_EXPORT_CORE GaussLobattoIntegrator {
public:
	/**
	 * Initialize a Gauss-Lobatto integration scheme
	 *
	 * \param maxEvals Maximum number of function evaluations. The
	 *    integrator will throw an exception when this limit is
	 *    exceeded.
	 *
	 * \param absAccuracy Absolute accuracy requirement (0 to disable)
	 * \param relAccuracy Relative accuracy requirement (0 to disable)
	 *
	 * \param useConvergenceEstimate Estimate the convergence behavior
	 *     of the GL-quadrature by comparing the 4, 7 and 13-point
	 *     variants and increase the absolute tolerance accordingly.
	 */
	GaussLobattoIntegrator(size_t maxEvals,
						 Float absAccuracy = 0,
						 Float relAccuracy = 0,
						 bool useConvergenceEstimate = true);

	Float integrate(const boost::function<Float (Float)>& f,
					Float a, Float b);

	inline size_t getEvaluations() const { return m_evals; }
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
		Float a, Float b, Float fa, Float fb,
		Float is);

	/**
	 * Compute the absolute error tolerance using a 13-point 
	 * Gauss-Lobatto rule.
	 */
	Float calculateAbsTolerance(const boost::function<Float (Float)>& f,
		Float a, Float b);

	Float m_absAccuracy, m_relAccuracy;
	size_t m_evals, m_maxEvals;
	const bool m_useConvergenceEstimate;
	static const Float m_alpha;
	static const Float m_beta;
	static const Float m_x1;
	static const Float m_x2;
	static const Float m_x3;
};

MTS_NAMESPACE_END

#endif /* __QUADRATURE_H */

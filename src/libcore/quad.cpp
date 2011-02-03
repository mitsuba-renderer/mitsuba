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

/*! \file gausslabottointegral.cpp
    \brief integral of a one-dimensional function using the adaptive 
           Gauss-Lobatto integral
*/

#include <mitsuba/core/quad.h>

MTS_NAMESPACE_BEGIN

const Float GaussLobattoIntegrator::m_alpha = (Float) std::sqrt(2.0/3.0); 
const Float GaussLobattoIntegrator::m_beta  = (Float) (1.0/std::sqrt(5.0));
const Float GaussLobattoIntegrator::m_x1	= (Float) 0.94288241569547971906; 
const Float GaussLobattoIntegrator::m_x2	= (Float) 0.64185334234578130578;
const Float GaussLobattoIntegrator::m_x3	= (Float) 0.23638319966214988028;

GaussLobattoIntegrator::GaussLobattoIntegrator(size_t maxEvals,
	Float absAccuracy, Float relAccuracy, bool useConvergenceEstimate)
	: m_absAccuracy(absAccuracy), 
	  m_relAccuracy(relAccuracy),
	  m_maxEvals(maxEvals),
	  m_useConvergenceEstimate(useConvergenceEstimate) {
	if (m_absAccuracy == 0 && m_relAccuracy == 0)
		SLog(EError, "GaussLobattoIntegrator:: Absolute and relative "
			"accuracy requirements can't both be zero!");
}

Float GaussLobattoIntegrator::integrate(
		const boost::function<Float (Float)>& f, Float a, Float b) {
	Float factor = 1;
	if (a == b) {
		return 0;
	} else if (b < a) {
		std::swap(a, b);
		factor = -1;
	}
	m_evals = 0;
	const Float absTolerance = calculateAbsTolerance(f, a, b);
	m_evals += 2;
	return factor * adaptiveGaussLobattoStep(f, a, b, f(a), f(b), absTolerance);
}

Float GaussLobattoIntegrator::calculateAbsTolerance(
		const boost::function<Float (Float)>& f, Float a, Float b) {
	const Float m = (a+b)/2; 
	const Float h = (b-a)/2;
	const Float y1 = f(a);
	const Float y3 = f(m-m_alpha*h);
	const Float y5 = f(m-m_beta*h);
	const Float y7 = f(m);
	const Float y9 = f(m+m_beta*h);
	const Float y11= f(m+m_alpha*h);
	const Float y13= f(b);

	Float acc=h*((Float) 0.0158271919734801831*(y1+y13)
			  +(Float) 0.0942738402188500455*(f(m-m_x1*h)+f(m+m_x1*h))
			  +(Float) 0.1550719873365853963*(y3+y11)
			  +(Float) 0.1888215739601824544*(f(m-m_x2*h)+ f(m+m_x2*h))
			  +(Float) 0.1997734052268585268*(y5+y9) 
			  +(Float) 0.2249264653333395270*(f(m-m_x3*h)+f(m+m_x3*h))
			  +(Float) 0.2426110719014077338*y7);  
	m_evals += 13;

	if (acc == 0)
		SLog(EError, "GaussLobattoIntegrator: Cannot calculate absolute accuracy from relative accuracy");

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

	if (m_relAccuracy != 0)	
		result = acc * std::max(m_relAccuracy,
			std::numeric_limits<Float>::epsilon())
			/ (r*std::numeric_limits<Float>::epsilon());

	if (m_absAccuracy != 0)
		result = std::min(result, m_absAccuracy 
			/ (r*std::numeric_limits<Float>::epsilon()));

	return result;
}

Float GaussLobattoIntegrator::adaptiveGaussLobattoStep(
								 const boost::function<Float (Float)>& f,
								 Float a, Float b, Float fa, Float fb,
								 Float acc) {
	if (m_evals >= m_maxEvals)
		SLog(EError, "GaussLobattoIntegrator: Maximum number of evaluations reached!");

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

	m_evals += 5;

	Float dist = acc + (integral1-integral2);
	if(dist==acc || mll<=a || b<=mrr) {
		if (m<=a || b<=m)
			SLog(EWarn, "GaussLobattoIntegrator: Interval contains no more machine numbers!");
		return integral1;
	} else {
		return  adaptiveGaussLobattoStep(f,a,mll,fa,fmll,acc)  
			  + adaptiveGaussLobattoStep(f,mll,ml,fmll,fml,acc)
			  + adaptiveGaussLobattoStep(f,ml,m,fml,fm,acc)
			  + adaptiveGaussLobattoStep(f,m,mr,fm,fmr,acc)
			  + adaptiveGaussLobattoStep(f,mr,mrr,fmr,fmrr,acc)
			  + adaptiveGaussLobattoStep(f,mrr,b,fmrr,fb,acc);
	}
}

MTS_NAMESPACE_END

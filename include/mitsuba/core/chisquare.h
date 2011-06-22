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

#if !defined(__CHI_SQUARE_TEST_H)
#define __CHI_SQUARE_TEST_H

#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/quad.h>
#include <mitsuba/core/timer.h>
#include <boost/bind.hpp>
#include <boost/math/distributions/chi_squared.hpp>

MTS_NAMESPACE_BEGIN

/**
 * \brief Chi-square goodness-of-fit test on the sphere
 *
 * This class performs a chi-square goodness-of-fit test of the null hypothesis
 * that a specified sampling procedure produces samples that are distributed
 * according to a supplied density function. This is very useful to verify BRDF
 * and phase function sampling codes for their correctness.
 *
 * This implementation works by generating a large batch of samples, which are
 * then accumulated into rectangular bins in spherical coordinates. To obtain
 * reference bin counts, the provided density function is numerically
 * integrated over the area of each bin. Comparing the actual and reference
 * bin counts yields the desired test statistic.
 *
 * Given a probability distribution with the following interface
 *
 * <code>
 * class MyDistribution {
 *     // Sample a (optionally weighted) direction. A non-unity weight
 *     // in the return value is needed when the sampling distribution
 *     // doesn't exactly match the implementation in pdf()
 *     std::pair<Vector, Float> generateSample() const;
 *
 *     /// Compute the probability density for the specified direction
 *     Float pdf(const Vector &direction) const;
 * };
 * </code>
 *
 * the code in this class might be used as follows
 * 
 * <code>
 * MyDistribution myDistrInstance;
 * ChiSquareTest chiSqr;
 *
 * // Initialize the tables used by the chi-square test
 * chiSqr.fill(
 *    boost::bind(&MyDistribution::generateSample, myDistrInstance),
 *    boost::bind(&MyDistribution::pdf, myDistrInstance, _1)
 * );
 *
 * // Optional: dump the tables to a MATLAB file for external analysis
 * chiSqr.dumpTables("debug.m");
 *
 * // (the following assumes that the distribution has 1 parameter, e.g. exponent value)
 * if (!chiSqr.runTest(1))
 *    Log(EError, "Uh oh -- test failed, the implementation is probably incorrect!");
 * </code>
 */
class MTS_EXPORT_CORE ChiSquareTest : public Object {
public:
	/**
	 * \brief Create a new Chi-square test instance with the given
	 * resolution and sample count
	 *
	 * \param thetaBins
	 *    Number of bins wrt. latitude. The default is 10
	 *
	 * \param phiBins
	 *    Number of bins wrt. azimuth. The default is to use
	 *    twice the number of \c thetaBins
	 *
	 * \param sampleCount
	 *    Number of samples to be used when computing the bin
	 *    values. The default is \c thetaBins*phiBins*5000
	 */
	ChiSquareTest(int thetaBins = 10, int phiBins = 0, size_t sampleCount = 0);

	/// Set the log level
	inline void setLogLevel(ELogLevel logLevel) { m_logLevel = EInfo; }

	/**
	 * \brief Fill the actual and reference bin counts
	 *
	 * Please see the class documentation for a description
	 * on how to invoke this function
	 */
	void fill(
		const boost::function<std::pair<Vector, Float>()> &sampleFn,
		const boost::function<Float (const Vector &)> &pdfFn);

	/**
	 * \brief Dump the bin counts to a file using MATLAB format
	 */
	void dumpTables(const fs::path &filename);

	/**
	 * \brief Perform the actual chi-square test
	 *
	 * \param distParams
	 *     Number of parameters of the distribution in question.
	 *     Anything such as lobe width, indices of refraction, etc.
	 *     should be counted.
	 *
	 * \param pvalThresh
	 *     The implementation will reject the null hypothesis
	 *     when the computed p-value lies below this parameter
	 *     (default: 0.01f)
	 *
	 * \return \c false if the null hypothesis was rejected
	 *     and \c true otherwise.
	 */
	bool runTest(int distParams, Float pvalThresh = 0.01f);

	MTS_DECLARE_CLASS()
protected:
	/// Release all memory
	virtual ~ChiSquareTest();

	/// Functor to evaluate the pdf values in parallel using OpenMP
	static void integrand(
		const boost::function<Float (const Vector &)> &pdfFn,
			size_t nPts, const Float *in, Float *out) {
		#pragma omp parallel for
		for (int i=0; i<(int) nPts; ++i)
			out[i] = pdfFn(sphericalDirection(in[2*i], in[2*i+1])) * std::sin(in[2*i]);
	}
private:
	ELogLevel m_logLevel;
	int m_thetaBins, m_phiBins;
	size_t m_sampleCount;
	Float *m_table;
	Float *m_refTable;
};

MTS_NAMESPACE_END

#endif /* __CHI_SQUARE_TEST_H */

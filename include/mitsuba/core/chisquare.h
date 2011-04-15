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
 * // (the folowing assumes that the distribution has 1 parameter, e.g. exponent value)
 * if (!chiSqr.runTest(1))
 *    Log(EError, "Uh oh -- test failed, the implementation is probably incorrect!");
 * </code>
 */
class ChiSquareTest {
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
	ChiSquareTest(int thetaBins = 10, int phiBins = 0, size_t sampleCount = 0)
	    	: m_thetaBins(thetaBins), m_phiBins(phiBins), m_sampleCount(sampleCount) {
		if (m_phiBins == 0)
			m_phiBins = 2*m_thetaBins;
		if (m_sampleCount == 0)
			m_sampleCount = m_thetaBins * m_phiBins * 1000;
		m_table = new Float[m_thetaBins*m_phiBins];
		m_refTable = new Float[m_thetaBins*m_phiBins];
	}

	/// Release all memory
	~ChiSquareTest() {
		delete[] m_table;
		delete[] m_refTable;
	}

	/**
	 * \brief Fill the actual and reference bin counts
	 *
	 * Please see the class documentation for a description
	 * on how to invoke this function
	 */
	void fill(
		const boost::function<std::pair<Vector, Float>()> &sampleFn,
		const boost::function<Float (const Vector &)> &pdfFn) {
		memset(m_table, 0, m_thetaBins*m_phiBins*sizeof(Float));

		SLog(EInfo, "Accumulating " SIZE_T_FMT " samples into a %ix%i"
				" contingency table", m_sampleCount, m_thetaBins, m_phiBins);
		Point2 factor(m_thetaBins / M_PI, m_phiBins / (2*M_PI));

		ref<Timer> timer = new Timer();
		for (size_t i=0; i<m_sampleCount; ++i) {
			std::pair<Vector, Float> sample = sampleFn();
			Point2 sphCoords = toSphericalCoordinates(sample.first);

			int thetaBin = std::min(std::max(0,
				floorToInt(sphCoords.x * factor.x)), m_thetaBins-1);
			int phiBin = std::min(std::max(0,
				floorToInt(sphCoords.y * factor.y)), m_phiBins-1);

			m_table[thetaBin * m_phiBins + phiBin] += sample.second;
		}
		SLog(EInfo, "Done, took %i ms.", timer->getMilliseconds());
		factor = Point2(M_PI / m_thetaBins, (2*M_PI) / m_phiBins);

		SLog(EInfo, "Integrating reference contingency table");
		timer->reset();
		Float min[2], max[2];
		size_t idx = 0;

		NDIntegrator integrator(1, 2, 100000, 0, 1e-6f);
		Float maxError = 0, integral = 0;
		for (int i=0; i<m_thetaBins; ++i) {
			min[0] = i * factor.x;
			max[0] = (i+1) * factor.x;
			for (int j=0; j<m_phiBins; ++j) {
				min[1] = j * factor.y;
				max[1] = (j+1) * factor.y;
				Float result, error;
				size_t evals;

				integrator.integrateVectorized(
					boost::bind(&ChiSquareTest::integrand, pdfFn, _1, _2, _3),
					min, max, &result, &error, evals
				);

				integral += result;
				m_refTable[idx++] = result * m_sampleCount;
				maxError = std::max(maxError, error);
			}
		}
		SLog(EInfo, "Done, took %i ms (max error = %f, integral=%f).", 
				timer->getMilliseconds(), maxError, integral);
	}

	/**
	 * \brief Dump the bin counts to a file using MATLAB format
	 */
	void dumpTables(const fs::path &filename) {
		fs::ofstream out(filename);
		out << "tbl_counts = [ ";
		for (int i=0; i<m_thetaBins; ++i) {
			for (int j=0; j<m_phiBins; ++j) {
				out << m_table[i*m_phiBins+j];
				if (j+1 < m_phiBins)
					out << ", ";
			}
			if (i+1 < m_thetaBins)
				out << "; ";
		}
		out << " ];" << endl
			<< "tbl_ref = [ ";
		for (int i=0; i<m_thetaBins; ++i) {
			for (int j=0; j<m_phiBins; ++j) {
				out << m_refTable[i*m_phiBins+j];
				if (j+1 < m_phiBins)
					out << ", ";
			}
			if (i+1 < m_thetaBins)
				out << "; ";
		}
		out << " ];" << endl;
		out.close();
	}

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
	bool runTest(int distParams, Float pvalThresh = 0.01f) {
		/* Compute the chi-square statistic */
		Float chsq = 0.0f;
		Float pooledCounts = 0, pooledRef = 0;
		int pooledCells = 0;
		int df = 0;

		for (int i=0; i<m_thetaBins*m_phiBins; ++i) {
			if (m_refTable[i] < 5) {
				pooledCounts += m_table[i];
				pooledRef += m_refTable[i];
				++pooledCells;
			} else {
				Float diff = m_table[i]-m_refTable[i];
				chsq += (diff*diff) / m_refTable[i];
				++df;
			}
		}

		if (pooledCells > 0) {
			SLog(EInfo, "Pooled %i cells with an expected "
				"number of < 5 entries!", pooledCells);
			if (pooledRef < 5) {
				SLog(EWarn, "Even after pooling, the expected "
					"number of entries is < 5 (%f), expect badness!",
					pooledRef);
			}
			Float diff = pooledCounts - pooledRef;
			chsq += (diff*diff) / pooledRef;
			++df;
		}

		df -= distParams + 1;
		SLog(EInfo, "Chi-square statistic = %e (df=%i)", chsq, df);
		boost::math::chi_squared chSqDist(df);
		/* Probability of obtaining a test statistic at least
		   as extreme as the one observed under the assumption
		   that the distributions match */
		Float pval = 1 - boost::math::cdf(chSqDist, chsq);
		SLog(EInfo, "P-value = %e", pval);

		if (pval < pvalThresh) {
			SLog(EWarn, "Rejecting the null hypothesis");
			return false;
		}
		return true;

	}
protected:
	static void integrand(
		const boost::function<Float (const Vector &)> &pdfFn,
			size_t nPts, const Float *in, Float *out) {
		#pragma omp parallel for
		for (int i=0; i<(int) nPts; ++i)
			out[i] = pdfFn(sphericalDirection(in[2*i], in[2*i+1])) * std::sin(in[2*i]);
	}
private:
	int m_thetaBins, m_phiBins;
	size_t m_sampleCount;
	Float *m_table;
	Float *m_refTable;
};

MTS_NAMESPACE_END

#endif /* __CHI_SQUARE_TEST_H */

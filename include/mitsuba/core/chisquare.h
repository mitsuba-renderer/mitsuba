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
#if !defined(__MITSUBA_CORE_CHISQUARE_H_)
#define __MITSUBA_CORE_CHISQUARE_H_

#include <mitsuba/render/common.h>
#include <boost/tuple/tuple.hpp>
#include <boost/function.hpp>

MTS_NAMESPACE_BEGIN

/// Minimum expected cell frequency. Cells below this value will be pooled
#define CHISQR_MIN_EXP_FREQUENCY 5

/**
 * \brief Chi-square goodness-of-fit test on the sphere
 *
 * This class performs a chi-square goodness-of-fit test of the null hypothesis
 * that a specified sampling procedure produces samples that are distributed
 * according to a supplied density function. This is very useful to verify BRDF
 * and phase function sampling codes for their correctness. Currently, it
 * supports both 2D and discrete sampling methods and mixtures thereof.
 *
 * This implementation works by generating a large batch of samples, which are
 * then accumulated into rectangular bins in spherical coordinates. To obtain
 * reference bin counts, the provided density function is numerically
 * integrated over the area of each bin. Comparing the actual and reference
 * bin counts yields the desired test statistic.
 *
 * Given a probability distribution with the following interface
 *
 * \code
 * class MyDistribution {
 *     // Sample a (optionally weighted) direction. A non-unity weight
 *     // in the return value is needed when the sampling distribution
 *     // doesn't exactly match the implementation in pdf()
 *     boost::tuple<Vector, Float, EMeasure> generateSample() const;
 *
 *     /// Compute the probability density for the specified direction and measure
 *     Float pdf(const Vector &direction, EMeasure) const;
 * };
 * \endcode
 *
 * the code in this class might be used as follows
 *
 * \code
 * MyDistribution myDistrInstance;
 * ChiSquare chiSqr;
 *
 * // Initialize the tables used by the chi-square test
 * chiSqr.fill(
 *    boost::bind(&MyDistribution::generateSample, myDistrInstance),
 *    boost::bind(&MyDistribution::pdf, myDistrInstance, _1, _2)
 * );
 *
 * // Optional: dump the tables to a MATLAB file for external analysis
 * chiSqr.dumpTables("debug.m");
 *
 * if (!chiSqr.runTest())
 *    Log(EError, "Uh oh -- test failed, the implementation is probably incorrect!");
 * \endcode
 * \ingroup libcore
 */
class MTS_EXPORT_CORE ChiSquare : public Object {
public:
    /// Possible outcomes in \ref runTest()
    enum ETestResult {
        /// The null hypothesis was rejected
        EReject = 0,
        /// The null hypothesis was accepted
        EAccept = 1,
        /// The degrees of freedom were too low
        ELowDoF = 2
    };

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
     * \param numTests
     *    Number of independent tests that will be performed. This
     *    is used to compute the Sidak-correction factor.
     *
     * \param sampleCount
     *    Number of samples to be used when computing the bin
     *    values. The default is \c thetaBins*phiBins*5000
     */
    ChiSquare(int thetaBins = 10, int phiBins = 0,
            int numTests = 1, size_t sampleCount = 0);

    /// Get the log level
    inline ELogLevel getLogLevel() const { return m_logLevel; }

    /// Set the log level
    inline void setLogLevel(ELogLevel logLevel) { m_logLevel = logLevel; }

    /**
     * \brief Set the tolerance threshold for bins with very low
     * aggregate probabilities
     *
     * When the Chi-square test integrates the supplied probability
     * density function over the support of a bin and determines that
     * the aggregate bin probability is zero, the test would ordinarily
     * fail if as much as one sample is placed in that bin in the
     * subsequent sampling step. However, due to various numerical
     * errors in a system based on finite-precision arithmetic, it
     * may be a good idea to tolerate at least a few samples without
     * immediately rejecting the null hypothesis. This parameter
     * sets this threshold. The default value is \c number-of-samples*1e-4f
     */
    inline void setTolerance(Float tolerance) { m_tolerance = tolerance; }

    /**
     * \brief Fill the actual and reference bin counts
     *
     * Please see the class documentation for a description
     * on how to invoke this function
     */
    void fill(
        const boost::function<boost::tuple<Vector, Float, EMeasure>()> &sampleFn,
        const boost::function<Float (const Vector &, EMeasure)> &pdfFn);

    /**
     * \brief Dump the bin counts to a file using MATLAB format
     */
    void dumpTables(const fs::path &filename);

    /**
     * \brief Perform the actual chi-square test
     *
     * \param pvalThresh
     *     The implementation will reject the null hypothesis
     *     when the computed p-value lies below this parameter
     *     (default: 0.01f)
     *
     * \return A status value of type \ref ETestResult
     */
    ETestResult runTest(Float pvalThresh = 0.01f);

    MTS_DECLARE_CLASS()
protected:
    /// Release all memory
    virtual ~ChiSquare();

    /// Functor to evaluate the pdf values in parallel using OpenMP
    static void integrand(
        const boost::function<Float (const Vector &, EMeasure)> &pdfFn,
            size_t nPts, const Float *in, Float *out) {
        #if defined(MTS_OPENMP)
        #pragma omp parallel for
        #endif
        for (int i=0; i<(int) nPts; ++i)
            out[i] = pdfFn(sphericalDirection(in[2*i], in[2*i+1]), ESolidAngle)
                * std::sin(in[2*i]);
    }
private:
    ELogLevel m_logLevel;
    Float m_tolerance;
    int m_thetaBins, m_phiBins;
    int m_numTests;
    size_t m_sampleCount;
    Float *m_table;
    Float *m_refTable;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_CHISQUARE_H_ */

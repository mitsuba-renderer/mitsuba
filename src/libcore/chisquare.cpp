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

#include <mitsuba/core/chisquare.h>
#include <mitsuba/core/quad.h>
#include <mitsuba/core/timer.h>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem/fstream.hpp>
#include <set>

MTS_NAMESPACE_BEGIN

/* Simple ordering for storing vectors in a set */
struct VectorOrder {
	inline int compare(const Vector &v1, const Vector &v2) const {
		if (v1.x < v2.x) return -1;
		else if (v1.x > v2.x) return 1;
		if (v1.y < v2.y) return -1;
		else if (v1.y > v2.y) return 1;
		if (v1.z < v2.z) return -1;
		else if (v1.z > v2.z) return 1;
		return 0;
	}

	bool operator()(const Vector &v1, const Vector &v2) const {
		return compare(v1, v2) < 0;
	}
};

ChiSquare::ChiSquare(int thetaBins, int phiBins, int numTests,
		size_t sampleCount) : m_logLevel(EInfo), m_thetaBins(thetaBins),
		  m_phiBins(phiBins), m_numTests(numTests), m_sampleCount(sampleCount) {
	if (m_phiBins == 0)
		m_phiBins = 2*m_thetaBins;
	if (m_sampleCount == 0)
		m_sampleCount = m_thetaBins * m_phiBins * 1000;
	m_table = new Float[m_thetaBins*m_phiBins];
	m_refTable = new Float[m_thetaBins*m_phiBins];
	m_tolerance = m_sampleCount * 1e-4f;
}

ChiSquare::~ChiSquare() {
	delete[] m_table;
	delete[] m_refTable;
}

void ChiSquare::dumpTables(const fs::path &filename) {
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

void ChiSquare::fill(
	const boost::function<boost::tuple<Vector, Float, EMeasure>()> &sampleFn,
	const boost::function<Float (const Vector &, EMeasure measure)> &pdfFn) {
	memset(m_table, 0, m_thetaBins*m_phiBins*sizeof(Float));
	memset(m_refTable, 0, m_thetaBins*m_phiBins*sizeof(Float));

	Log(m_logLevel, "Accumulating " SIZE_T_FMT " samples into a %ix%i"
			" contingency table", m_sampleCount, m_thetaBins, m_phiBins);
	Point2 factor(m_thetaBins / M_PI, m_phiBins / (2*M_PI));

	std::set<Vector, VectorOrder> discreteDirections;

	ref<Timer> timer = new Timer();
	for (size_t i=0; i<m_sampleCount; ++i) {
		boost::tuple<Vector, Float, EMeasure> sample = sampleFn();
		Point2 sphCoords = toSphericalCoordinates(boost::get<0>(sample));

		int thetaBin = std::min(std::max(0,
			math::floorToInt(sphCoords.x * factor.x)), m_thetaBins-1);
		int phiBin = std::min(std::max(0,
			math::floorToInt(sphCoords.y * factor.y)), m_phiBins-1);
		m_table[thetaBin * m_phiBins + phiBin] += boost::get<1>(sample);
		if (boost::get<1>(sample) > 0 && boost::get<2>(sample) == EDiscrete)
			discreteDirections.insert(boost::get<0>(sample));
	}

	if (discreteDirections.size() > 0) {
		Log(EDebug, "Incorporating the disrete density over "
			SIZE_T_FMT " direction(s) into the contingency table", discreteDirections.size());
		for (std::set<Vector, VectorOrder>::const_iterator it = discreteDirections.begin();
			it != discreteDirections.end(); ++it) {
			const Vector &direction = *it;
			Point2 sphCoords = toSphericalCoordinates(direction);
			Float pdf = pdfFn(direction, EDiscrete);

			int thetaBin = std::min(std::max(0,
				math::floorToInt(sphCoords.x * factor.x)), m_thetaBins-1);
			int phiBin = std::min(std::max(0,
				math::floorToInt(sphCoords.y * factor.y)), m_phiBins-1);

			m_refTable[thetaBin * m_phiBins + phiBin] += pdf * m_sampleCount;
		}
	}

	factor = Point2(M_PI / m_thetaBins, (2*M_PI) / m_phiBins);

	Log(m_logLevel, "Done, took %i ms. Integrating reference "
		"contingency table ..", timer->getMilliseconds());
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

			integrator.integrateVectorized(
				boost::bind(&ChiSquare::integrand, pdfFn, _1, _2, _3),
				min, max, &result, &error
			);

			integral += result;
			m_refTable[idx++] += result * m_sampleCount;
			maxError = std::max(maxError, error);
		}
	}

	Log(m_logLevel, "Done, took %i ms (max error = %f, integral=%f).",
			timer->getMilliseconds(), maxError, integral);
}

struct SortedCell {
	Float expCount;
	int idx;
};

struct SortedCellFunctor {
	inline bool operator()(const SortedCell &c1, const SortedCell &c2) const {
		return c1.expCount < c2.expCount;
	}
};

ChiSquare::ETestResult ChiSquare::runTest(Float pvalThresh) {
	/* Compute the chi-square statistic */
	Float pooledCounts = 0, pooledRef = 0, chsq = 0.0f;
	int pooledCells = 0, df = 0;

	/* Process cells in order sorted by their expected counts */
	std::vector<SortedCell> cells(m_thetaBins*m_phiBins);
	for (int i=0; i<m_thetaBins*m_phiBins; ++i) {
		cells[i].expCount = m_refTable[i];
		cells[i].idx = i;
	}

	std::sort(cells.begin(), cells.end(), SortedCellFunctor());

	std::vector<SortedCell>::iterator it = cells.begin();

	while (it != cells.end()) {
		int idx = it->idx;

		if (m_refTable[idx] == 0) {
			if (m_table[idx] > m_tolerance) {
				/* Special handler for cells with an expected frequency of zero */
				Log(EWarn, "Encountered a cell (%i) with an expected frequency of zero, "
					"where the actual number of observations is %f! Rejecting the "
					"null hypothesis.", idx, m_table[idx]);
				return EReject;
			}
		} else if (m_refTable[idx] < CHISQR_MIN_EXP_FREQUENCY) {
			/* Pool cells with low expected frequencies */
			pooledCounts += m_table[idx];
			pooledRef += m_refTable[idx];
			++pooledCells;
		} else if (pooledRef > 0 && pooledRef < CHISQR_MIN_EXP_FREQUENCY) {
			/* Pool more cells until the merged cell
			   has a sufficiently high frequency */
			pooledCounts += m_table[idx];
			pooledRef += m_refTable[idx];
			++pooledCells;
		} else {
			Float diff = m_table[idx]-m_refTable[idx];
			chsq += (diff*diff) / m_refTable[idx];
			++df;
		}

		++it;
	}

	if (pooledCells > 0) {
		Log(m_logLevel, "Pooled %i cells to ensure sufficiently "
			"high expected frequencies (> %f).", pooledCells,
			(Float) CHISQR_MIN_EXP_FREQUENCY);
		Float diff = pooledCounts - pooledRef;
		chsq += (diff*diff) / pooledRef;
		++df;
	}

	/* All parameters are assumed to be known, so there is no
	   DF reduction due to model parameters */
	df -= 1;

	Log(m_logLevel, "Chi-square statistic = %e (df=%i)", chsq, df);

	if (df <= 0) {
		Log(m_logLevel, "The number of degrees of freedom (%i) is too low!", df);
		return ELowDoF;
	}

	/* Probability of obtaining a test statistic at least
	   as extreme as the one observed under the assumption
	   that the distributions match */
	boost::math::chi_squared chSqDist(df);
	Float pval = 1 - (Float) boost::math::cdf(chSqDist, chsq);

	/* Apply the Sidak correction for multiple independent hypothesis tests */
	Float alpha = 1 - std::pow(1 - pvalThresh, 1 / (Float) m_numTests);

	if (pval < alpha) {
		Log(EWarn, "Rejected the null hypothesis (P-value = %e, "
			"significance level = %e)", pval, alpha);
		return EReject;
	} else {
		Log(m_logLevel, "Accepted the null hypothesis (P-value = %e, "
			"significance level = %e)", pval, alpha);
		return EAccept;
	}
}

MTS_IMPLEMENT_CLASS(ChiSquare, false, Object)
MTS_NAMESPACE_END

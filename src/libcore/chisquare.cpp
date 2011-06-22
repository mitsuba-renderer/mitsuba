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

#include <mitsuba/core/chisquare.h>

MTS_NAMESPACE_BEGIN

ChiSquareTest::ChiSquareTest(int thetaBins, int phiBins, size_t sampleCount) 
    	: m_logLevel(EInfo), m_thetaBins(thetaBins), m_phiBins(phiBins), 
		  m_sampleCount(sampleCount) {
	if (m_phiBins == 0)
		m_phiBins = 2*m_thetaBins;
	if (m_sampleCount == 0)
		m_sampleCount = m_thetaBins * m_phiBins * 1000;
	m_table = new Float[m_thetaBins*m_phiBins];
	m_refTable = new Float[m_thetaBins*m_phiBins];
}

ChiSquareTest::~ChiSquareTest() {
	delete[] m_table;
	delete[] m_refTable;
}

bool ChiSquareTest::runTest(int distParams, Float pvalThresh) {
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

	if (pooledCells > 0 && pooledRef != 0) {
		Log(m_logLevel, "Pooled %i cells with an expected "
			"number of < 5 entries!", pooledCells);
		if (pooledRef < 5) {
			Log(EWarn, "Even after pooling %i cells, the expected "
				"number of entries is < 5 (%f), expect badness!",
				pooledCells, pooledRef);
		}
		Float diff = pooledCounts - pooledRef;
		chsq += (diff*diff) / pooledRef;
		++df;
	}

	df -= distParams + 1;
	Log(m_logLevel, "Chi-square statistic = %e (df=%i)", chsq, df);
	boost::math::chi_squared chSqDist(df);
	/* Probability of obtaining a test statistic at least
	   as extreme as the one observed under the assumption
	   that the distributions match */
	Float pval = 1 - boost::math::cdf(chSqDist, chsq);
	Log(m_logLevel, pval > 0.01 ? "P-value = %.4f" : "P-value = %e", pval);

	if (pval < pvalThresh) {
		Log(EWarn, "Rejecting the null hypothesis (P-value=%e)", pval);
		return false;
	}
	return true;
}

void ChiSquareTest::dumpTables(const fs::path &filename) {
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

void ChiSquareTest::fill(
	const boost::function<std::pair<Vector, Float>()> &sampleFn,
	const boost::function<Float (const Vector &)> &pdfFn) {
	memset(m_table, 0, m_thetaBins*m_phiBins*sizeof(Float));

	Log(m_logLevel, "Accumulating " SIZE_T_FMT " samples into a %ix%i"
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
	factor = Point2(M_PI / m_thetaBins, (2*M_PI) / m_phiBins);

	Log(m_logLevel, "Done, took %i ms. Integrating reference contingency table ..", timer->getMilliseconds());
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
	Log(m_logLevel, "Done, took %i ms (max error = %f, integral=%f).", 
			timer->getMilliseconds(), maxError, integral);
}

MTS_IMPLEMENT_CLASS(ChiSquareTest, false, Object)
MTS_NAMESPACE_END

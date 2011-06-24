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

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/chisquare.h>
#include <mitsuba/render/testcase.h>
#include <boost/bind.hpp>

/* Statistical significance level of the test. Set to
   1/2 percent by default -- we want there to be notable
   evidence before failing a test case */
#define SIGNIFICANCE_LEVEL 0.005f

MTS_NAMESPACE_BEGIN

/**
 * This testcase checks if the sampling methods of various BSDF & phase 
 * function implementations really do what they promise in their pdf()
 * methods
 */
class TestChiSquare : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test01_BSDF)
	MTS_DECLARE_TEST(test02_PhaseFunction)
	MTS_END_TESTCASE()

	/// Adapter to use BSDFs in the chi-square test
	class BSDFAdapter {
	public:
		BSDFAdapter(const BSDF *bsdf, Sampler *sampler, const Vector &wi, int component = -1)
			: m_bsdf(bsdf), m_sampler(sampler), m_wi(wi), m_component(component),
			  m_largestWeight(0) { }

		std::pair<Vector, Float> generateSample() {
			Point2 sample(m_sampler->next2D());
			Intersection its;
			BSDFQueryRecord bRec(its);
			bRec.component = m_component;
			bRec.wi = m_wi;
	
			/* Check the various sampling routines for agreement amongst each other */
			Float pdfVal;
			Spectrum f = m_bsdf->sample(bRec, pdfVal, sample);
			Spectrum sampled = m_bsdf->sample(bRec, sample);
			
			if (f.isZero() || pdfVal == 0) {
				if (!sampled.isZero()) 
					Log(EWarn, "Inconsistency: f=%s, pdf=%f, sampled f/pdf=%s",
						f.toString().c_str(), pdfVal, sampled.toString().c_str());
				return std::make_pair(bRec.wo, 0.0f);
			} else if (sampled.isZero()) {
				if (!f.isZero() && pdfVal != 0)
					Log(EWarn, "Inconsistency: f=%s, pdf=%f, sampled f/pdf=%s",
						f.toString().c_str(), pdfVal, sampled.toString().c_str());
				return std::make_pair(bRec.wo, 0.0f);
			}

			Spectrum sampled2 = f/pdfVal;
			bool mismatch = false;

			for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
				Float a = sampled[i], b = sampled2[i];
				SAssert(a >= 0 && b >= 0);
				Float min = std::min(a, b);
				Float err = std::abs(a - b);
				m_largestWeight = std::max(m_largestWeight, a * std::abs(Frame::cosTheta(bRec.wo)));

				if (min < Epsilon && err > Epsilon) // absolute error threshold
					mismatch = true;
				else if (min > Epsilon && err/min > Epsilon) // relative error threshold
					mismatch = true;
			}

			if (mismatch)
				Log(EWarn, "Inconsistency: f=%s, pdf=%f, sampled f/pdf=%s",
					f.toString().c_str(), pdfVal, sampled.toString().c_str());

			return std::make_pair(bRec.wo, 1.0f);
		}
 
		Float pdf(const Vector &wo) const {
			Intersection its;
			BSDFQueryRecord bRec(its);
			bRec.component = m_component;
			bRec.wi = m_wi;
			bRec.wo = wo;
			if (m_bsdf->f(bRec).isZero())
				return 0.0f;
			return m_bsdf->pdf(bRec);
		}

		inline Float getLargestWeight() const { return m_largestWeight; }
	private:
		ref<const BSDF> m_bsdf;
		ref<Sampler> m_sampler;
		Vector m_wi;
		int m_component;
		Float m_largestWeight;
	};

	/// Adapter to use Phase functions in the chi-square test
	class PhaseFunctionAdapter {
	public:
		PhaseFunctionAdapter(const MediumSamplingRecord &mRec,
				const PhaseFunction *phase, Sampler *sampler, const Vector &wi)
			: m_mRec(mRec), m_phase(phase), m_sampler(sampler), m_wi(wi), 
			  m_largestWeight(0) { }

		std::pair<Vector, Float> generateSample() {
			Point2 sample(m_sampler->next2D());
			PhaseFunctionQueryRecord pRec(m_mRec, m_wi);
	
			/* Check the various sampling routines for agreement amongst each other */
			Float pdfVal;
			Float f = m_phase->sample(pRec, pdfVal, m_sampler);
			Float sampled = m_phase->sample(pRec, m_sampler);

			if (f == 0 || pdfVal == 0) {
				if (sampled != 0)
					Log(EWarn, "Inconsistency: f=%f, pdf=%f, sampled f/pdf=%f",
						f, pdfVal, sampled);
				return std::make_pair(pRec.wo, 0.0f);
			} else if (sampled == 0) {
				if (f != 0 && pdfVal != 0)
					Log(EWarn, "Inconsistency: f=%f, pdf=%f, sampled f/pdf=%f",
						f, pdfVal, sampled);
				return std::make_pair(pRec.wo, 0.0f);
			}

			Float sampled2 = f/pdfVal;
			bool mismatch = false;

			SAssert(sampled >= 0 && sampled2 >= 0);
			Float min = std::min(sampled, sampled2);
			Float err = std::abs(sampled - sampled2);
			m_largestWeight = std::max(m_largestWeight, sampled);

			if (min < Epsilon && err > Epsilon) // absolute error threshold
				mismatch = true;
			else if (min > Epsilon && err/min > Epsilon) // relative error threshold
				mismatch = true;

			if (mismatch)
				Log(EWarn, "Inconsistency: f=%f, pdf=%f, sampled f/pdf=%f",
					f, pdfVal, sampled);

			return std::make_pair(pRec.wo, 1.0f);
		}
 
		Float pdf(const Vector &wo) const {
			PhaseFunctionQueryRecord pRec(m_mRec, m_wi, wo);
			if (m_phase->f(pRec) == 0)
				return 0.0f;
			return m_phase->pdf(pRec);
		}

		inline Float getLargestWeight() const { return m_largestWeight; }
	private:
		const MediumSamplingRecord &m_mRec;
		ref<const PhaseFunction> m_phase;
		ref<Sampler> m_sampler;
		Vector m_wi;
		Float m_largestWeight;
	};

	void test01_BSDF() {
		/* Load a set of BSDF instances to be tested from the following XML file */
		ref<Scene> scene = loadScene("data/tests/test_bsdf.xml");
	
		const std::vector<ConfigurableObject *> objects = scene->getReferencedObjects();
		size_t thetaBins = 10, wiSamples = 20, failureCount = 0, testCount = 0;
		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Sampler), Properties("independent")));
		ProgressReporter *progress = new ProgressReporter("Checking", wiSamples, NULL);

		Log(EInfo, "Verifying BSDF sampling routines ..");
		for (size_t i=0; i<objects.size(); ++i) {
			if (!objects[i]->getClass()->derivesFrom(MTS_CLASS(BSDF)))
				continue;

			const BSDF *bsdf = static_cast<const BSDF *>(objects[i]);
			Float largestWeight = 0;

			Log(EInfo, "Processing BSDF model %s", bsdf->toString().c_str());
			Log(EInfo, "Checking the model for %i incident directions", wiSamples);
			progress->reset();

			/* Test for a number of different incident directions */
			for (size_t j=0; j<wiSamples; ++j) {
				Vector wi;
	
				if (bsdf->getType() & BSDF::EBackSide)
					wi = squareToSphere(sampler->next2D());
				else
					wi = squareToHemispherePSA(sampler->next2D());

				BSDFAdapter adapter(bsdf, sampler, wi);
				ref<ChiSquare> chiSqr = new ChiSquare(thetaBins, 2*thetaBins, wiSamples);
				chiSqr->setLogLevel(EDebug);

				// Initialize the tables used by the chi-square test
				chiSqr->fill(
					boost::bind(&BSDFAdapter::generateSample, &adapter),
					boost::bind(&BSDFAdapter::pdf, &adapter, _1)
				);

				// (the following assumes that the distribution has 1 parameter, e.g. exponent value)
				ChiSquare::ETestResult result = chiSqr->runTest(1, SIGNIFICANCE_LEVEL);
				if (result == ChiSquare::EReject) {
					std::string filename = formatString("failure_%i.m", failureCount++);
					chiSqr->dumpTables(filename);
					failAndContinue(formatString("Uh oh, the chi-square test indicates a potential "
						"issue for wi=%s. Dumped the contingency tables to '%s' for user analysis", 
						wi.toString().c_str(), filename.c_str()));
				} else {
					succeed();
				}
				largestWeight = std::max(largestWeight, adapter.getLargestWeight());
				++testCount;
				progress->update(j+1);
			}

			if (bsdf->getComponentCount() > 1) {
				for (int comp=0; comp<bsdf->getComponentCount(); ++comp) {
					progress->reset();
					Log(EInfo, "Individually checking BSDF component %i", comp);

					/* Test for a number of different incident directions */
					for (size_t j=0; j<wiSamples; ++j) {
						Vector wi;
			
						if (bsdf->getType(comp) & BSDF::EBackSide)
							wi = squareToSphere(sampler->next2D());
						else
							wi = squareToHemispherePSA(sampler->next2D());

						BSDFAdapter adapter(bsdf, sampler, wi, comp);

						ref<ChiSquare> chiSqr = new ChiSquare(thetaBins, 2*thetaBins, wiSamples);
						chiSqr->setLogLevel(EDebug);

						// Initialize the tables used by the chi-square test
						chiSqr->fill(
							boost::bind(&BSDFAdapter::generateSample, &adapter),
							boost::bind(&BSDFAdapter::pdf, &adapter, _1)
						);

						// (the following assumes that the distribution has 1 parameter, e.g. exponent value)
						ChiSquare::ETestResult result = chiSqr->runTest(1, SIGNIFICANCE_LEVEL);
						if (result == ChiSquare::EReject) {
							std::string filename = formatString("failure_%i.m", failureCount++);
							chiSqr->dumpTables(filename);
							failAndContinue(formatString("Uh oh, the chi-square test indicates a potential "
								"issue for wi=%s. Dumped the contingency tables to '%s' for user analysis", 
								wi.toString().c_str(), filename.c_str()));
						} else {
							succeed();
						}
						largestWeight = std::max(largestWeight, adapter.getLargestWeight());
						++testCount;
						progress->update(j+1);
					}
				}
				Log(EInfo, "Done with this BSDF. The largest encountered "
						"importance weight was = %.2f", largestWeight);
				largestWeight = 0;
			}
		}
		Log(EInfo, "%i/%i BSDF checks succeeded", testCount-failureCount, testCount);
		delete progress;
	}

	void test02_PhaseFunction() {
		/* Load a set of BSDF instances to be tested from the following XML file */
		ref<Scene> scene = loadScene("data/tests/test_phase.xml");
	
		const std::vector<ConfigurableObject *> objects = scene->getReferencedObjects();
		size_t thetaBins = 10, wiSamples = 20, failureCount = 0, testCount = 0;
		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Sampler), Properties("independent")));

		ProgressReporter *progress = new ProgressReporter("Checking", wiSamples, NULL);

		Log(EInfo, "Verifying phase function sampling routines ..");
		for (size_t i=0; i<objects.size(); ++i) {
			if (!objects[i]->getClass()->derivesFrom(MTS_CLASS(PhaseFunction)))
				continue;

			const PhaseFunction *phase = static_cast<const PhaseFunction *>(objects[i]);
			Float largestWeight = 0;

			Log(EInfo, "Processing phase function model %s", phase->toString().c_str());
			Log(EInfo, "Checking the model for %i incident directions", wiSamples);
			progress->reset();
			MediumSamplingRecord mRec;

			/* Sampler fiber/particle orientation */
			mRec.orientation = squareToSphere(sampler->next2D());

			/* Test for a number of different incident directions */
			for (size_t j=0; j<wiSamples; ++j) {
				Vector wi = squareToSphere(sampler->next2D());

				PhaseFunctionAdapter adapter(mRec, phase, sampler, wi);
				ref<ChiSquare> chiSqr = new ChiSquare(thetaBins, 2*thetaBins, wiSamples);
				chiSqr->setLogLevel(EDebug);

				// Initialize the tables used by the chi-square test
				chiSqr->fill(
					boost::bind(&PhaseFunctionAdapter::generateSample, &adapter),
					boost::bind(&PhaseFunctionAdapter::pdf, &adapter, _1)
				);

				// (the following assumes that the distribution has 1 parameter, e.g. exponent value)
				ChiSquare::ETestResult result = chiSqr->runTest(1, SIGNIFICANCE_LEVEL);
				if (result == ChiSquare::EReject) {
					std::string filename = formatString("failure_%i.m", failureCount++);
					chiSqr->dumpTables(filename);
					failAndContinue(formatString("Uh oh, the chi-square test indicates a potential "
						"issue for wi=%s. Dumped the contingency tables to '%s' for user analysis", 
						wi.toString().c_str(), filename.c_str()));
				} else {
					succeed();
				}
				largestWeight = std::max(largestWeight, adapter.getLargestWeight());
				++testCount;
				progress->update(j+1);
			}

			Log(EInfo, "Done with this phase function. The largest encountered "
					"importance weight was = %.2f", largestWeight);
		}
		Log(EInfo, "%i/%i phase function checks succeeded", testCount-failureCount, testCount);
		delete progress;
	}
};

MTS_EXPORT_TESTCASE(TestChiSquare, "Chi-square test for various sampling functions")
MTS_NAMESPACE_END

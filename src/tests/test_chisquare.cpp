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
   1/4 percent by default -- we want there to be strong 
   evidence of an implementaiton error before failing 
   a test case */
#define SIGNIFICANCE_LEVEL 0.0025f

/* Relative bound on what is still accepted as roundoff 
   error -- be quite tolerant */
#if defined(SINGLE_PRECISION)
	#define ERROR_REQ 1e-2f
#else
	#define ERROR_REQ 1e-5
#endif

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

	/**
	 * Replayable fake sampler
	 */
	class FakeSampler : public Sampler {
	public:
		FakeSampler(Sampler *sampler)
			: Sampler(Properties()), m_sampler(sampler) { }

		Float next1D() {
			while (m_sampleIndex >= m_values.size())
				m_values.push_back(m_sampler->next1D());
			return m_values[m_sampleIndex++];
		}

		Point2 next2D() {
			return Point2(next1D(), next1D());
		}

		void clear() {
			m_values.clear();
			m_sampleIndex = 0;
		}

		void rewind() {
			m_sampleIndex = 0;
		}
		
		Float independent1D() { SLog(EError, "Not supported!"); return 0; }
		Point2 independent2D() { SLog(EError, "Not supported!"); return Point2(0.0f); }

		ref<Sampler> clone() {
			SLog(EError, "Not supported!");
			return NULL;
		}

		std::string toString() const { return "FakeSampler[]"; }
	private:
		ref<Sampler> m_sampler;
		std::vector<Float> m_values;
	};

	/// Adapter to use BSDFs in the chi-square test
	class BSDFAdapter {
	public:
		BSDFAdapter(const BSDF *bsdf, Sampler *sampler, const Vector &wi, 
				int component, bool passSamplerToBSDF)
			: m_bsdf(bsdf), m_sampler(sampler), m_wi(wi), m_component(component),
			  m_largestWeight(0), m_passSamplerToBSDF(passSamplerToBSDF) {
			m_fakeSampler = new FakeSampler(m_sampler);
		}

		boost::tuple<Vector, Float, EMeasure> generateSample() {
			Point2 sample(m_sampler->next2D());
			Intersection its;
			BSDFQueryRecord bRec(its);
			bRec.component = m_component;
			bRec.wi = m_wi;
			
			#if defined(MTS_DEBUG_FP)
				enableFPExceptions();
			#endif
	
			Float pdfVal, pdfVal2;
	
			/* Only make the sampler available to the BSDF when requested
			   by the testcase. This allows testing both sampling variants
			   where applicable: those that can improve by having access to 
			   an arbitrary random number stream vs. those that only use
			   a single uniform 2D sample */
	
			if (m_passSamplerToBSDF)
				bRec.sampler = m_fakeSampler;

			/* Check the various sampling routines for agreement 
			   amongst each other */
			m_fakeSampler->clear();
			Spectrum f = m_bsdf->sample(bRec, pdfVal, sample);
			m_fakeSampler->rewind();
			Spectrum sampled = m_bsdf->sample(bRec, sample);
			EMeasure measure = ESolidAngle;
			if (!f.isZero())
				measure = BSDF::getMeasure(bRec.sampledType);
			Spectrum f2 = m_bsdf->eval(bRec, measure);
			pdfVal2 = m_bsdf->pdf(bRec, measure);

			if (f.isZero() || pdfVal == 0 || pdfVal2 == 0) {
				if (!sampled.isZero())
					Log(EWarn, "Inconsistency (1): f=%s, f2=%s, pdf=%f, pdf2=%f, sampled f/pdf=%s, bRec=%s",
						f.toString().c_str(), f2.toString().c_str(), pdfVal, pdfVal2, sampled.toString().c_str(), bRec.toString().c_str());
				#if defined(MTS_DEBUG_FP)
					disableFPExceptions();
				#endif
				return boost::make_tuple(bRec.wo, 0.0f, ESolidAngle);
			} else if (sampled.isZero()) {
				if ((!f.isZero() && pdfVal != 0) || (!f2.isZero() && pdfVal2 != 0))
					Log(EWarn, "Inconsistency (2): f=%s, f2=%s, pdf=%f, pdf2=%f, sampled f/pdf=%s, bRec=%s",
						f.toString().c_str(), f2.toString().c_str(), pdfVal, pdfVal2, sampled.toString().c_str(), bRec.toString().c_str());
				#if defined(MTS_DEBUG_FP)
					disableFPExceptions();
				#endif
				return boost::make_tuple(bRec.wo, 0.0f, ESolidAngle);
			}

			Spectrum sampled2 = f/pdfVal, evaluated = f2/pdfVal2;
			if (!sampled.isValid() || !sampled2.isValid() || !evaluated.isValid()) {
				Log(EWarn, "Ooops: f=%s, f2=%s, pdf=%f, pdf2=%f, sampled f/pdf=%s, bRec=%s",
					f.toString().c_str(), f2.toString().c_str(), pdfVal, pdfVal2, sampled.toString().c_str(), bRec.toString().c_str());
				return boost::make_tuple(bRec.wo, 0.0f, ESolidAngle);
			}

			bool mismatch = false;
			for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
				Float a = sampled[i], b = sampled2[i], c = evaluated[i];
				Float min = std::min(std::min(a, b), c);
				Float err = std::max(std::max(std::abs(a - b), std::abs(a - c)), std::abs(b - c));
				m_largestWeight = std::max(m_largestWeight, a);

				if (min < ERROR_REQ && err > ERROR_REQ) // absolute error threshold
					mismatch = true;
				else if (min > ERROR_REQ && err/min > ERROR_REQ) // relative error threshold
					mismatch = true;
			}

			if (mismatch)
				Log(EWarn, "Potential inconsistency (3): f/pdf=%s (method 1), f/pdf=%s (methdod 2), sampled f/pdf=%s",
					sampled2.toString().c_str(), evaluated.toString().c_str(), sampled.toString().c_str());

			#if defined(MTS_DEBUG_FP)
				disableFPExceptions();
			#endif

			return boost::make_tuple(bRec.wo, 1.0f, measure);
		}
 
		Float pdf(const Vector &wo, EMeasure measure) {
			Intersection its;
			BSDFQueryRecord bRec(its);
			bRec.component = m_component;
			bRec.wi = m_wi;
			bRec.wo = wo;
			if (m_passSamplerToBSDF)
				bRec.sampler = m_sampler;

			#if defined(MTS_DEBUG_FP)
				enableFPExceptions();
			#endif

			if (m_bsdf->eval(bRec, measure).isZero())
				return 0.0f;

			Float result = m_bsdf->pdf(bRec, measure);

			#if defined(MTS_DEBUG_FP)
				disableFPExceptions();
			#endif
			return result;
		}

		inline Float getLargestWeight() const { return m_largestWeight; }
	private:
		ref<const BSDF> m_bsdf;
		ref<Sampler> m_sampler;
		ref<FakeSampler> m_fakeSampler;
		Vector m_wi;
		int m_component;
		Float m_largestWeight;
		bool m_passSamplerToBSDF;
	};

	/// Adapter to use Phase functions in the chi-square test
	class PhaseFunctionAdapter {
	public:
		PhaseFunctionAdapter(const MediumSamplingRecord &mRec,
				const PhaseFunction *phase, Sampler *sampler, const Vector &wi)
			: m_mRec(mRec), m_phase(phase), m_sampler(sampler), m_wi(wi), 
			  m_largestWeight(0) { }

		boost::tuple<Vector, Float, EMeasure> generateSample() {
			Point2 sample(m_sampler->next2D());
			PhaseFunctionQueryRecord pRec(m_mRec, m_wi);
			
			#if defined(MTS_DEBUG_FP)
				enableFPExceptions();
			#endif

			/* Check the various sampling routines for agreement amongst each other */
			Float pdfVal;
			Float f = m_phase->sample(pRec, pdfVal, m_sampler);
			Float sampled = m_phase->sample(pRec, m_sampler);

			if (f == 0 || pdfVal == 0) {
				if (sampled != 0)
					Log(EWarn, "Inconsistency: f=%f, pdf=%f, sampled f/pdf=%f",
						f, pdfVal, sampled);
				#if defined(MTS_DEBUG_FP)
					disableFPExceptions();
				#endif
				return boost::make_tuple(pRec.wo, 0.0f, ESolidAngle);
			} else if (sampled == 0) {
				if (f != 0 && pdfVal != 0)
					Log(EWarn, "Inconsistency: f=%f, pdf=%f, sampled f/pdf=%f",
						f, pdfVal, sampled);
				#if defined(MTS_DEBUG_FP)
					disableFPExceptions();
				#endif
				return boost::make_tuple(pRec.wo, 0.0f, ESolidAngle);
			}

			Float sampled2 = f/pdfVal;
			bool mismatch = false;

			SAssert(sampled >= 0 && sampled2 >= 0);
			Float min = std::min(sampled, sampled2);
			Float err = std::abs(sampled - sampled2);
			m_largestWeight = std::max(m_largestWeight, sampled);

			if (min < ERROR_REQ && err > ERROR_REQ) // absolute error threshold
				mismatch = true;
			else if (min > ERROR_REQ && err/min > ERROR_REQ) // relative error threshold
				mismatch = true;

			if (mismatch)
				Log(EWarn, "Inconsistency: f=%f, pdf=%f, sampled f/pdf=%f",
					f, pdfVal, sampled);

			#if defined(MTS_DEBUG_FP)
				disableFPExceptions();
			#endif
			return boost::make_tuple(pRec.wo, 1.0f, ESolidAngle);
		}
 
		Float pdf(const Vector &wo, EMeasure measure) const {
			if (measure != ESolidAngle)
				return 0.0f;

			PhaseFunctionQueryRecord pRec(m_mRec, m_wi, wo);
			#if defined(MTS_DEBUG_FP)
				enableFPExceptions();
			#endif
			if (m_phase->eval(pRec) == 0)
				return 0.0f;
			Float result = m_phase->pdf(pRec);
			#if defined(MTS_DEBUG_FP)
				disableFPExceptions();
			#endif
			return result;
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

			for (int pass=0; pass<2; ++pass) {
				Log(EInfo, "Checking the model for %i incident directions and 2D sampling", wiSamples);
				progress->reset();

				/**
				 * Do a second pass when the BSDF can optionally make use of 
				 * an infinite random number stream. The sampling method is likely
				 * quite different in this case, so everything has to be checked
				 * again.
				 */

				if (pass == 1) {
					if (!(bsdf->getType() & BSDF::ECanUseSampler))
						break;
					else
						Log(EInfo, "Checking the model again, but now with access to an infinite random number stream.");
				}

				/* Test for a number of different incident directions */
				for (size_t j=0; j<wiSamples; ++j) {
					Vector wi;
		
					if (bsdf->getType() & BSDF::EBackSide)
						wi = squareToSphere(sampler->next2D());
					else
						wi = squareToHemispherePSA(sampler->next2D());

					BSDFAdapter adapter(bsdf, sampler, wi, -1, pass == 1);
					ref<ChiSquare> chiSqr = new ChiSquare(thetaBins, 2*thetaBins, wiSamples);
					chiSqr->setLogLevel(EDebug);

					// Initialize the tables used by the chi-square test
					chiSqr->fill(
						boost::bind(&BSDFAdapter::generateSample, &adapter),
						boost::bind(&BSDFAdapter::pdf, &adapter, _1, _2)
					);

					// (the following assumes that the distribution has 1 parameter, e.g. exponent value)
					ChiSquare::ETestResult result = chiSqr->runTest(SIGNIFICANCE_LEVEL);
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
				Log(EInfo, "The largest encountered importance weight was = %.2f", largestWeight);
				largestWeight = 0;

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

							BSDFAdapter adapter(bsdf, sampler, wi, comp, pass == 1);

							ref<ChiSquare> chiSqr = new ChiSquare(thetaBins, 2*thetaBins, wiSamples);
							chiSqr->setLogLevel(EDebug);

							// Initialize the tables used by the chi-square test
							chiSqr->fill(
								boost::bind(&BSDFAdapter::generateSample, &adapter),
								boost::bind(&BSDFAdapter::pdf, &adapter, _1, _2)
							);

							// (the following assumes that the distribution has 1 parameter, e.g. exponent value)
							ChiSquare::ETestResult result = chiSqr->runTest(SIGNIFICANCE_LEVEL);
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
						Log(EInfo, "The largest encountered importance weight was = %.2f", largestWeight);
						largestWeight = 0;
					}
				}
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
					boost::bind(&PhaseFunctionAdapter::pdf, &adapter, _1, _2)
				);

				// (the following assumes that the distribution has 1 parameter, e.g. exponent value)
				ChiSquare::ETestResult result = chiSqr->runTest(SIGNIFICANCE_LEVEL);
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

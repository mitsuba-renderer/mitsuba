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

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/chisquare.h>
#include <mitsuba/core/fresolver.h>
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
 * function & emitter implementations really do what they promise in
 * their pdf() methods
 */
class TestChiSquare : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test01_BSDF)
	MTS_DECLARE_TEST(test02_PhaseFunction)
	MTS_DECLARE_TEST(test03_EmitterDirect)
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
		BSDFAdapter(const BSDF *bsdf, Sampler *sampler, const Vector &wi, int component)
			: m_bsdf(bsdf), m_sampler(sampler), m_wi(wi), m_component(component),
			  m_largestWeight(0) {
			m_fakeSampler = new FakeSampler(m_sampler);
			m_its.uv = Point2(0.0f);
			m_its.dpdu = Vector(1, 0, 0);
			m_its.dpdv = Vector(0, 1, 0);
			m_its.dudx = m_its.dvdy = 0.01f;
			m_its.dudy = m_its.dvdx = 0.00f;
			m_its.shFrame = Frame(Normal(0, 0, 1));
//			m_isSymmetric = true;
		}

		boost::tuple<Vector, Float, EMeasure> generateSample() {
			Point2 sample(m_sampler->next2D());
			BSDFSamplingRecord bRec(m_its, m_fakeSampler);
			bRec.mode = EImportance;
			bRec.component = m_component;
			bRec.wi = m_wi;

			#if defined(MTS_DEBUG_FP)
				enableFPExceptions();
			#endif

			Float pdfVal, sampledPDF;

			/* Check the various sampling routines for agreement
			   amongst each other */
			m_fakeSampler->clear();
			Spectrum sampled = m_bsdf->sample(bRec, sampledPDF, sample);
			m_fakeSampler->rewind();
			Spectrum sampled2 = m_bsdf->sample(bRec, sample);
			EMeasure measure = ESolidAngle;
			if (!sampled.isZero())
				measure = BSDF::getMeasure(bRec.sampledType);

			if (sampled.isZero() && sampled2.isZero())
				return boost::make_tuple(Vector(0.0f), 0.0f, measure);

			Spectrum f = m_bsdf->eval(bRec, measure);
			pdfVal = m_bsdf->pdf(bRec, measure);
			Spectrum manual = f/pdfVal;

#if 0
			if (m_isSymmetric) {
				/* Check for non-symmetry */
				BSDFSamplingRecord bRecRev(bRec);
				bRecRev.reverse();
				bRec.mode = EImportance;
				Spectrum fFwd = f;
				Spectrum fRev = m_bsdf->eval(bRecRev, measure);
				if (measure == ESolidAngle) {
					fFwd /= std::abs(Frame::cosTheta(bRec.wo));
					fRev /= std::abs(Frame::cosTheta(bRecRev.wo));
				}
				Float max = std::max(fFwd.max(), fRev.max());
				if (max > 0) {
					Float err = (fFwd-fRev).max() / max;
					if (err > Epsilon) {
						Log(EWarn, "Non-symmetry in %s: %s vs %s, %s", m_bsdf->toString().c_str(),
							fFwd.toString().c_str(), fRev.toString().c_str(), bRec.toString().c_str());
						m_isSymmetric = false;
					}
				}
			}
#endif

			if (!sampled.isValid() || !sampled2.isValid() || !manual.isValid()) {
				Log(EWarn, "Oops: sampled=%s, sampled2=%s, manual=%s, sampledPDF=%f, "
					"pdf=%f, f=%s, bRec=%s, measure=%i", sampled.toString().c_str(),
					sampled2.toString().c_str(), manual.toString().c_str(),
					sampledPDF, pdfVal, f.toString().c_str(), bRec.toString().c_str(),
					measure);
				return boost::make_tuple(bRec.wo, 0.0f, ESolidAngle);
			}

			bool mismatch = false;
			for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
				Float a = sampled[i], b = sampled2[i], c = manual[i];
				Float min = std::min(std::min(a, b), c);
				Float err = std::max(std::max(std::abs(a - b), std::abs(a - c)), std::abs(b - c));
				m_largestWeight = std::max(m_largestWeight, a);

				if (min < ERROR_REQ && err > ERROR_REQ) // absolute error threshold
					mismatch = true;
				else if (min > ERROR_REQ && err/min > ERROR_REQ) // relative error threshold
					mismatch = true;
			}

			if (mismatch)
				Log(EWarn, "Potential inconsistency: sampled=%s, sampled2=%s, manual=%s, sampledPDF=%f, "
					"pdf=%f, f=%s, bRec=%s, measure=%i", sampled.toString().c_str(),
					sampled2.toString().c_str(), manual.toString().c_str(),
					sampledPDF, pdfVal, f.toString().c_str(), bRec.toString().c_str(),
					measure);

			mismatch = false;
			Float min = std::min(pdfVal, sampledPDF);
			Float err = std::abs(pdfVal - sampledPDF);

			if (min < ERROR_REQ && err > ERROR_REQ) // absolute error threshold
				mismatch = true;
			else if (min > ERROR_REQ && err/min > ERROR_REQ) // relative error threshold
				mismatch = true;

			if (mismatch)
				Log(EWarn, "Potential inconsistency: pdfVal=%f, sampledPDF=%f",
					pdfVal, sampledPDF);

			#if defined(MTS_DEBUG_FP)
				disableFPExceptions();
			#endif

			return boost::make_tuple(bRec.wo, 1.0f, measure);
		}

		Float pdf(const Vector &wo, EMeasure measure) {
			BSDFSamplingRecord bRec(m_its, m_wi, wo);
			bRec.mode = EImportance;
			bRec.component = m_component;

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
//		inline bool isSymmetric() const { return m_isSymmetric; }
	private:
		Intersection m_its;
		ref<const BSDF> m_bsdf;
		ref<Sampler> m_sampler;
		ref<FakeSampler> m_fakeSampler;
		Vector m_wi;
		int m_component;
		Float m_largestWeight;
//		bool m_isSymmetric;
	};

	/// Adapter to use Phase functions in the chi-square test
	class PhaseFunctionAdapter {
	public:
		PhaseFunctionAdapter(const MediumSamplingRecord &mRec,
				const PhaseFunction *phase, Sampler *sampler, const Vector &wi)
			: m_mRec(mRec), m_phase(phase), m_sampler(sampler), m_wi(wi),
			  m_largestWeight(0) {
			m_fakeSampler = new FakeSampler(m_sampler);
		}

		boost::tuple<Vector, Float, EMeasure> generateSample() {
			PhaseFunctionSamplingRecord pRec(m_mRec, m_wi);

			#if defined(MTS_DEBUG_FP)
				enableFPExceptions();
			#endif
			Float sampledPDF, pdfVal;

			/* Check the various sampling routines for agreement amongst each other */
			m_fakeSampler->clear();
			Float sampled = m_phase->sample(pRec, sampledPDF, m_fakeSampler);
			m_fakeSampler->rewind();
			Float sampled2 = m_phase->sample(pRec, m_fakeSampler);
			Float f = m_phase->eval(pRec);
			pdfVal = m_phase->pdf(pRec);
			Float manual = f/pdfVal;

			if (std::isnan(sampled) || std::isnan(sampled2) || std::isnan(manual) ||
				sampled < 0 || sampled2 < 0 || manual < 0) {
				Log(EWarn, "Oops: sampled=%f, sampled2=%f, manual=%f, sampledPDF=%f, "
					"pdf=%f, f=%f, pRec=%s", sampled, sampled2, manual,
					sampledPDF, pdfVal, f, pRec.toString().c_str());
				return boost::make_tuple(pRec.wo, 0.0f, ESolidAngle);
			}

			bool mismatch = false;
			Float min = std::min(std::min(sampled, sampled2), manual);
			Float err = std::max(std::max(std::abs(sampled - sampled2),
					std::abs(sampled - manual)), std::abs(sampled2 - manual));
			m_largestWeight = std::max(m_largestWeight, sampled);

			if (min < ERROR_REQ && err > ERROR_REQ) // absolute error threshold
				mismatch = true;
			else if (min > ERROR_REQ && err/min > ERROR_REQ) // relative error threshold
				mismatch = true;

			if (mismatch)
				Log(EWarn, "Potential inconsistency: sampled=%f, sampled2=%f, manual=%s, "
					"sampledPDF=%f, pdf=%f, f=%f, pRec=%s", sampled, sampled2, manual,
					sampledPDF, pdfVal, f, pRec.toString().c_str());

			mismatch = false;
			min = std::min(pdfVal, sampledPDF);
			err = std::abs(pdfVal - sampledPDF);

			if (min < ERROR_REQ && err > ERROR_REQ) // absolute error threshold
				mismatch = true;
			else if (min > ERROR_REQ && err/min > ERROR_REQ) // relative error threshold
				mismatch = true;

			if (mismatch)
				Log(EWarn, "Potential inconsistency: pdfVal=%f, sampledPDF=%f",
					pdfVal, sampledPDF);


			return boost::make_tuple(pRec.wo,
				sampled == 0 ? 0.0f : 1.0f, ESolidAngle);
		}

		Float pdf(const Vector &wo, EMeasure measure) const {
			if (measure != ESolidAngle)
				return 0.0f;

			PhaseFunctionSamplingRecord pRec(m_mRec, m_wi, wo);
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
		ref<FakeSampler> m_fakeSampler;
		ref<const PhaseFunction> m_phase;
		ref<Sampler> m_sampler;
		Vector m_wi;
		Float m_largestWeight;
	};

	/// Adapter to use direct illumination sampling in the chi-square test
	class EmitterAdapter {
	public:
		EmitterAdapter(const Emitter *emitter, Sampler *sampler)
			: m_emitter(emitter), m_sampler(sampler), m_pRec(0.0f) {
			emitter->samplePosition(m_pRec, m_sampler->next2D());
		}

		boost::tuple<Vector, Float, EMeasure> generateSample() {
			#if defined(MTS_DEBUG_FP)
				enableFPExceptions();
			#endif

			DirectSamplingRecord dRec(Point(0.0f), 0);
			m_emitter->sampleDirect(dRec, m_sampler->next2D());

			#if defined(MTS_DEBUG_FP)
				disableFPExceptions();
			#endif

			return boost::make_tuple(dRec.d, 1.0f, dRec.measure);
		}

		Float pdf(const Vector &d, EMeasure measure) const {
			if (measure != ESolidAngle)
				return 0.0f;

			DirectSamplingRecord dRec(Point(0.0f), 0);
			dRec.d = d;
			dRec.measure = ESolidAngle;

			#if defined(MTS_DEBUG_FP)
				enableFPExceptions();
			#endif

			Float result = m_emitter->pdfDirect(dRec);

			#if defined(MTS_DEBUG_FP)
				disableFPExceptions();
			#endif

			return result;
		}

	private:
		ref<const Emitter> m_emitter;
		ref<Sampler> m_sampler;
		PositionSamplingRecord m_pRec;
	};

	void test01_BSDF() {
		/* Load a set of BSDF instances to be tested from the following XML file */
		FileResolver *resolver = Thread::getThread()->getFileResolver();
		const fs::path scenePath =
			resolver->resolveAbsolute("data/tests/test_bsdf.xml");
		ref<Scene> scene = loadScene(scenePath);

		const ref_vector<ConfigurableObject> &objects = scene->getReferencedObjects();
		int thetaBins = 10, wiSamples = 20, failureCount = 0, testCount = 0;
		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Sampler), Properties("independent")));
		ProgressReporter *progress = new ProgressReporter("Checking", wiSamples, NULL);

		Log(EInfo, "Verifying BSDF sampling routines ..");
		for (size_t i=0; i<objects.size(); ++i) {
			if (!objects[i]->getClass()->derivesFrom(MTS_CLASS(BSDF)))
				continue;

			const BSDF *bsdf = static_cast<const BSDF *>(objects[i].get());
			Float largestWeight = 0;

			Log(EInfo, "Processing BSDF model %s", bsdf->toString().c_str());

			Log(EInfo, "Checking the model for %i incident directions and 2D sampling", wiSamples);
			progress->reset();
#if 1
			/* Test for a number of different incident directions */
			for (int j=0; j<wiSamples; ++j) {
				Vector wi;

				if (bsdf->getType() & BSDF::EBackSide)
					wi = warp::squareToUniformSphere(sampler->next2D());
				else
					wi = warp::squareToCosineHemisphere(sampler->next2D());

				BSDFAdapter adapter(bsdf, sampler, wi, -1);
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

//				if (!adapter.isSymmetric())
//					Log(EWarn, "****** BSDF is non-symmetric! ******");
			}
			Log(EInfo, "The largest encountered importance weight was = %.2f", largestWeight);

#endif
			largestWeight = 0;

			if (bsdf->getComponentCount() > 1) {
				for (int comp=0; comp<bsdf->getComponentCount(); ++comp) {
					progress->reset();
					Log(EInfo, "Individually checking BSDF component %i", comp);

					/* Test for a number of different incident directions */
					for (int j=0; j<wiSamples; ++j) {
						Vector wi;

						if (bsdf->getType(comp) & BSDF::EBackSide)
							wi = warp::squareToUniformSphere(sampler->next2D());
						else
							wi = warp::squareToCosineHemisphere(sampler->next2D());

						BSDFAdapter adapter(bsdf, sampler, wi, comp);

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
		Log(EInfo, "%i/%i BSDF checks succeeded", testCount-failureCount, testCount);
		delete progress;
	}

	void test02_PhaseFunction() {
		/* Load a set of phase function instances to be tested from the following XML file */
		FileResolver *resolver = Thread::getThread()->getFileResolver();
		const fs::path scenePath =
			resolver->resolveAbsolute("data/tests/test_phase.xml");
		ref<Scene> scene = loadScene(scenePath);

		const ref_vector<ConfigurableObject> &objects = scene->getReferencedObjects();
		int thetaBins = 10, wiSamples = 20, failureCount = 0, testCount = 0;
		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Sampler), Properties("independent")));

		ProgressReporter *progress = new ProgressReporter("Checking", wiSamples, NULL);

		Log(EInfo, "Verifying phase function sampling routines ..");
		for (size_t i=0; i<objects.size(); ++i) {
			if (!objects[i]->getClass()->derivesFrom(MTS_CLASS(PhaseFunction)))
				continue;

			const PhaseFunction *phase = static_cast<const PhaseFunction *>(objects[i].get());
			Float largestWeight = 0;

			Log(EInfo, "Processing phase function model %s", phase->toString().c_str());
			Log(EInfo, "Checking the model for %i incident directions", wiSamples);
			progress->reset();
			MediumSamplingRecord mRec;

			/* Sampler fiber/particle orientation */
			mRec.orientation = warp::squareToUniformSphere(sampler->next2D());

			/* Test for a number of different incident directions */
			for (int j=0; j<wiSamples; ++j) {
				Vector wi = warp::squareToUniformSphere(sampler->next2D());

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

	void test03_EmitterDirect() {
		/* Load a set of emitter instances to be tested from the following XML file */
		FileResolver *resolver = Thread::getThread()->getFileResolver();
		const fs::path scenePath =
			resolver->resolveAbsolute("data/tests/test_emitter.xml");
		ref<Scene> scene = loadScene(scenePath);
		scene->initialize();

		const ref_vector<Emitter> &emitters = scene->getEmitters();
		int thetaBins = 10, failureCount = 0, testCount = 0;
		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Sampler), Properties("independent")));

		Log(EInfo, "Verifying emitter sampling routines ..");
		for (size_t i=0; i<emitters.size(); ++i) {
			const Emitter *emitter = emitters[i].get();

			Log(EInfo, "Processing emitter function model %s", emitter->toString().c_str());

			EmitterAdapter adapter(emitter, sampler);
			ref<ChiSquare> chiSqr = new ChiSquare(thetaBins, 2*thetaBins, 1);
			chiSqr->setLogLevel(EDebug);

			// Initialize the tables used by the chi-square test
			chiSqr->fill(
				boost::bind(&EmitterAdapter::generateSample, &adapter),
				boost::bind(&EmitterAdapter::pdf, &adapter, _1, _2)
			);
			chiSqr->dumpTables("test.m");

			// (the following assumes that the distribution has 1 parameter, e.g. exponent value)
			ChiSquare::ETestResult result = chiSqr->runTest(SIGNIFICANCE_LEVEL);
			if (result == ChiSquare::EReject) {
				std::string filename = formatString("failure_%i.m", failureCount++);
				chiSqr->dumpTables(filename);
				failAndContinue(formatString("Uh oh, the chi-square test indicates a potential "
					"issue. Dumped the contingency tables to '%s' for user analysis",
					filename.c_str()));
			} else {
				succeed();
			}
			++testCount;
		}
		Log(EInfo, "%i/%i emitter checks succeeded", testCount-failureCount, testCount);
	}
};

MTS_EXPORT_TESTCASE(TestChiSquare, "Chi-square test for various sampling functions")
MTS_NAMESPACE_END

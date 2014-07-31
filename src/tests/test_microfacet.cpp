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
#include "../bsdfs/microfacet.h"

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

class TestChiSquare : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test01_Microfacet)
	MTS_DECLARE_TEST(test02_MicrofacetVisible)
	MTS_END_TESTCASE()

	class MicrofacetAdapter {
	public:
		MicrofacetAdapter(Sampler *sampler, const MicrofacetDistribution &distr, const Vector &wi = Vector(0.0f)) : m_sampler(sampler), m_distr(distr), m_wi(wi) { }

		boost::tuple<Vector, Float, EMeasure> generateSample() {
			Float pdf;

			if (m_wi.lengthSquared() == 0) {
				Normal m = m_distr.sampleAll(m_sampler->next2D(), pdf);
				Float pdf_ref = m_distr.pdfAll(m);

				SAssert(std::isfinite(pdf) && pdf > 0);
				SAssert(std::isfinite(pdf_ref) && pdf_ref > 0);
				SAssert(std::isfinite(m.x) && std::isfinite(m.y) && std::isfinite(m.z));
				SAssert(std::abs(m.length() - 1) < 1e-4f);
				SAssert(std::abs((pdf-pdf_ref)/pdf_ref) < 1e-4f);
				return boost::make_tuple(m, 1.0f, ESolidAngle);
			} else {
				Normal m = m_distr.sampleVisible(m_wi, m_sampler->next2D());
				SAssert(std::isfinite(m.x) && std::isfinite(m.y) && std::isfinite(m.z));
				SAssert(std::abs(m.length() - 1) < 1e-4f);
				return boost::make_tuple(m, 1.0f, ESolidAngle);
			}
		}

		Float pdf(const Vector &d, EMeasure measure) const {
			if (measure != ESolidAngle)
				return 0.0f;

			Float pdf = m_wi.lengthSquared() == 0 ? m_distr.pdfAll(d)
				: m_distr.pdfVisible(m_wi, d);
			SAssert(std::isfinite(pdf) && pdf >= 0);

			return pdf;
		}

	private:
		ref<Sampler> m_sampler;
		MicrofacetDistribution m_distr;
		Vector m_wi;
	};

	void test01_Microfacet() {
		int thetaBins = 20;
		std::vector<MicrofacetDistribution> distrs;

		distrs.push_back(MicrofacetDistribution(MicrofacetDistribution::EBeckmann, 0.5f));
		distrs.push_back(MicrofacetDistribution(MicrofacetDistribution::EBeckmann, 0.5f, 0.3f));
		distrs.push_back(MicrofacetDistribution(MicrofacetDistribution::EGGX, 0.5f));
		distrs.push_back(MicrofacetDistribution(MicrofacetDistribution::EGGX, 0.5f, 0.3f));
		distrs.push_back(MicrofacetDistribution(MicrofacetDistribution::EPhong, 0.5f));
		distrs.push_back(MicrofacetDistribution(MicrofacetDistribution::EPhong, 0.5f, 0.3f));


		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Sampler), Properties("independent")));
		ref<ChiSquare> chiSqr = new ChiSquare(thetaBins, 2*thetaBins, (int) distrs.size());
		chiSqr->setLogLevel(EDebug);

		for (size_t i=0; i<distrs.size(); ++i) {
			Log(EInfo, "Testing %s", distrs[i].toString().c_str());
			// Initialize the tables used by the chi-square test
			MicrofacetAdapter adapter(sampler, distrs[i]);
			chiSqr->fill(
				boost::bind(&MicrofacetAdapter::generateSample, &adapter),
				boost::bind(&MicrofacetAdapter::pdf, &adapter, _1, _2)
			);

			// (the following assumes that the distribution has 1 parameter, e.g. exponent value)
			ChiSquare::ETestResult result = chiSqr->runTest(SIGNIFICANCE_LEVEL);
			if (result == ChiSquare::EReject) {
				std::string filename = formatString("failure_%i.m", (int) i);
				chiSqr->dumpTables(filename);
				failAndContinue(formatString("Uh oh, the chi-square test indicates a potential "
					"issue. Dumped the contingency tables to '%s' for user analysis",
					filename.c_str()));
			} else {
				//chiSqr->dumpTables(formatString("success_%i.m", (int) i));
				succeed();
			}
		}
	}

	void test02_MicrofacetVisible() {
		int thetaBins = 10;
		std::vector<std::pair<MicrofacetDistribution, Vector> > distrs;

		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Sampler), Properties("independent")));
		for (int i=0; i<10; ++i) {
			Vector wi = warp::squareToUniformHemisphere(sampler->next2D());
			distrs.push_back(std::make_pair(MicrofacetDistribution(MicrofacetDistribution::EBeckmann, 0.3f), wi));
			distrs.push_back(std::make_pair(MicrofacetDistribution(MicrofacetDistribution::EBeckmann, 0.5f, 0.3f), wi));
			distrs.push_back(std::make_pair(MicrofacetDistribution(MicrofacetDistribution::EGGX, 0.1f), wi));
			distrs.push_back(std::make_pair(MicrofacetDistribution(MicrofacetDistribution::EGGX, 0.2f, 0.3f), wi));
		}

		ref<ChiSquare> chiSqr = new ChiSquare(thetaBins, 2*thetaBins, (int) distrs.size());
		chiSqr->setLogLevel(EDebug);

		for (size_t i=0; i<distrs.size(); ++i) {
			Log(EInfo, "Testing %s (wi=%s)", distrs[i].first.toString().c_str(), distrs[i].second.toString().c_str());
			// Initialize the tables used by the chi-square test
			MicrofacetAdapter adapter(sampler, distrs[i].first, distrs[i].second);
			chiSqr->fill(
				boost::bind(&MicrofacetAdapter::generateSample, &adapter),
				boost::bind(&MicrofacetAdapter::pdf, &adapter, _1, _2)
			);

			// (the following assumes that the distribution has 1 parameter, e.g. exponent value)
			ChiSquare::ETestResult result = chiSqr->runTest(SIGNIFICANCE_LEVEL);
			if (result == ChiSquare::EReject) {
				std::string filename = formatString("failure_%i.m", (int) i);
				chiSqr->dumpTables(filename);
				failAndContinue(formatString("Uh oh, the chi-square test indicates a potential "
					"issue. Dumped the contingency tables to '%s' for user analysis",
					filename.c_str()));
			} else {
				//chiSqr->dumpTables(formatString("success_%i.m", (int) i));
				succeed();
			}
		}
	}
};

MTS_EXPORT_TESTCASE(TestChiSquare, "Chi-square test for microfacet sampling")
MTS_NAMESPACE_END

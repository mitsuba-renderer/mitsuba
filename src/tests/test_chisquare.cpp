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
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/chisquare.h>
#include <mitsuba/render/testcase.h>

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
	MTS_END_TESTCASE()

	class BSDFFunctor {
	public:
		BSDFFunctor(const BSDF *bsdf, Random *random, const Vector &wi)
				: m_bsdf(bsdf), m_random(random), m_wi(wi) { }

		std::pair<Vector, Float> generateSample() {
			Point2 sample(m_random->nextFloat(), m_random->nextFloat());
			Intersection its;
			BSDFQueryRecord bRec(its);
			bRec.wi = m_wi;
			m_bsdf->sample(bRec, sample);

			return std::make_pair(bRec.wo, 1.0f);
		}
 
		Float pdf(const Vector &wo) const {
			Intersection its;
			BSDFQueryRecord bRec(its);
			bRec.wi = m_wi;
			bRec.wo = wo;
			return m_bsdf->pdf(bRec);
		}
	private:
		ref<const BSDF> m_bsdf;
		ref<Random> m_random;
		Vector m_wi;
	};
	
	void test01_BSDF() {
		/* Load a set of BSDF instances to be tested from the following XML file */
		ref<Scene> scene = loadScene("data/tests/test_bsdf.xml");
	
		const std::vector<ConfigurableObject *> objects = scene->getReferencedObjects();
		size_t thetaBins = 10, wiSamples = 5;
		ref<Random> random = new Random();

		for (size_t i=0; i<objects.size(); ++i) {
			if (!objects[i]->getClass()->derivesFrom(MTS_CLASS(BSDF)))
				continue;

			const BSDF *bsdf = static_cast<const BSDF *>(objects[i]);

			Log(EInfo, "Checking sampling with %i different incident directions, BSDF:\n%s", 
					wiSamples, bsdf->toString().c_str());

			/* Test for a number of different incident directions */
			for (size_t j=0; j<wiSamples; ++j) {
				Vector wi;
	
				if (bsdf->getType() & (BSDF::EDiffuseTransmission | BSDF::EGlossyTransmission))
					wi = squareToSphere(Point2(random->nextFloat(), random->nextFloat()));
				else
					wi = squareToHemispherePSA(Point2(random->nextFloat(), random->nextFloat()));

				BSDFFunctor functor(bsdf, random, wi);
				ref<ChiSquareTest> chiSqr = new ChiSquareTest(thetaBins);

				// Initialize the tables used by the chi-square test
				chiSqr->fill(
					boost::bind(&BSDFFunctor::generateSample, functor),
					boost::bind(&BSDFFunctor::pdf, functor, _1)
				);

				// (the folowing assumes that the distribution has 1 parameter, e.g. exponent value)
				if (!chiSqr->runTest(1)) {
					// Optional: dump the tables to a MATLAB file for external analysis
					chiSqr->dumpTables("failure.m");
					failAndContinue("Oh oh, the chi-square test failed! Dumped the contingency tables to 'failure.m'");
				} else {
					succeed();
				}
			}
		}
	}
};

MTS_EXPORT_TESTCASE(TestChiSquare, "Chi-square test for various sampling functions")
MTS_NAMESPACE_END

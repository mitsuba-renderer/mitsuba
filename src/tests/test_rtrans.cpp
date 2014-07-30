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

#include <mitsuba/render/testcase.h>
#include <mitsuba/core/quad.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/plugin.h>
#include <boost/bind.hpp>
#include "../bsdfs/rtrans.h"

MTS_NAMESPACE_BEGIN

void transmittanceIntegrand(const BSDF *bsdf, const Vector &wi, size_t nPts, const Float *in, Float *out) {
	Intersection its;

	for (size_t i=0; i<nPts; ++i) {
		BSDFSamplingRecord bRec(its, wi, Vector(), EImportance);
		bRec.typeMask = BSDF::ETransmission;
		out[i] = bsdf->sample(bRec, Point2(in[2*i], in[2*i+1]))[0];
	}
}

void diffTransmittanceIntegrand(Float *data, size_t resolution, size_t nPts, const Float *in, Float *out) {
	for (size_t i=0; i<nPts; ++i)
		out[i] = 2 * in[i] * evalCubicInterp1D(in[i], data, resolution, 0, 1);
}

class TestRoughTransmittance : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test01_smoothTransmittance)
	MTS_DECLARE_TEST(test02_roughTransmittance)
	MTS_DECLARE_TEST(test03_roughTransmittanceFixedEta)
	MTS_DECLARE_TEST(test04_roughTransmittanceFixedEtaFixedAlpha)
	MTS_END_TESTCASE()

	Float computeDiffuseTransmittance(const char *name, Float eta, Float alpha, size_t resolution = 100) {
		Properties bsdfProps("roughdielectric");
		if (eta < 1) {
			bsdfProps.setFloat("intIOR", 1.0f);
			bsdfProps.setFloat("extIOR", 1.0f / eta);
		} else {
			bsdfProps.setFloat("extIOR", 1.0f);
			bsdfProps.setFloat("intIOR", eta);
		}

		bsdfProps.setFloat("alpha", alpha);
		bsdfProps.setString("distribution", name);
		ref<BSDF> bsdf = static_cast<BSDF *>(
				PluginManager::getInstance()->createObject(bsdfProps));

		Float *transmittances = new Float[resolution];
		Float stepSize = 1.0f / (resolution-1);
		Float error;

		NDIntegrator intTransmittance(1, 2, 50000, 0, 1e-6f);
		NDIntegrator intDiffTransmittance(1, 1, 50000, 0, 1e-6f);

		for (size_t i=0; i<resolution; ++i) {
			Float cosTheta = stepSize * i;
			Vector wi(math::safe_sqrt(1-cosTheta*cosTheta), 0, cosTheta);

			Float min[2] = {0, 0}, max[2] = {1, 1};
			intTransmittance.integrateVectorized(
				boost::bind(&transmittanceIntegrand, bsdf, wi, _1, _2, _3),
				min, max, &transmittances[i], &error, NULL);
		}

		Float Fdr;
		Float min[1] = { 0 }, max[1] = { 1 };
		intDiffTransmittance.integrateVectorized(
			boost::bind(&diffTransmittanceIntegrand, transmittances, resolution, _1, _2, _3),
			min, max, &Fdr, &error, NULL);

		delete[] transmittances;
		return Fdr;
	}

	Float computeTransmittance(const char *name, Float eta, Float alpha, Float cosTheta) {
		Properties bsdfProps("roughdielectric");
		if (cosTheta < 0) {
			cosTheta = -cosTheta;
			eta = 1.0f / eta;
		}
		if (eta < 1) {
			bsdfProps.setFloat("intIOR", 1.0f);
			bsdfProps.setFloat("extIOR", 1.0f / eta);
		} else {
			bsdfProps.setFloat("extIOR", 1.0f);
			bsdfProps.setFloat("intIOR", eta);
		}
		bsdfProps.setFloat("alpha", alpha);
		bsdfProps.setString("distribution", name);

		ref<BSDF> bsdf = static_cast<BSDF *>(
				PluginManager::getInstance()->createObject(bsdfProps));

		NDIntegrator intTransmittance(1, 2, 50000, 0, 1e-6f);

		Vector wi(math::safe_sqrt(1-cosTheta*cosTheta), 0, cosTheta);
		Float transmittance, error;

		Float min[2] = {0, 0}, max[2] = {1, 1};
		intTransmittance.integrateVectorized(
			boost::bind(&transmittanceIntegrand, bsdf, wi, _1, _2, _3),
			min, max, &transmittance, &error, NULL);

		return transmittance;
	}

	void test01_smoothTransmittance() {
		/* Smooth diffuse transmittance - compare polynomial approximations to ground truth */
		for (int i=0; i<=10; ++i) {
			Float eta = 1 + i/10.0f;

			Float f1 = fresnelDiffuseReflectance(eta, false);
			Float f2 = fresnelDiffuseReflectance(eta, true);
			Float f3 = fresnelDiffuseReflectance(1/eta, false);
			Float f4 = fresnelDiffuseReflectance(1/eta, true);

			assertEqualsEpsilon(std::abs(f1-f2), (Float) 0, 1e-3f);
			assertEqualsEpsilon(std::abs(f3-f4), (Float) 0, 1e-3f);
		}
	}

	void test02_roughTransmittance() {
		RoughTransmittance rtr(MicrofacetDistribution::EBeckmann);
		ref<Random> random = new Random();

		for (int i=0; i<50; ++i) {
			Float alpha = std::pow(random->nextFloat(), (Float) 4.0f)*4;
			Float eta = 1 + std::pow(random->nextFloat(), (Float) 4.0f)*3;
			if (alpha < 1e-5)
				alpha = 1e-5f;
			if (eta < 1+1e-5)
				eta = 1+1e-5f;
			//eta = 1/eta;

			Float refD = computeDiffuseTransmittance("beckmann", eta, alpha);
			Float datD = rtr.evalDiffuse(alpha, eta);

			cout << "Testing " << i << "/50" << endl;
			if (std::abs(refD-datD) > 1e-3f) {
				cout << endl;
				cout << "eta = " << eta << endl;
				cout << "alpha = " << alpha << endl;
				cout << "diff=" << datD-refD << " (datD=" << datD << ", ref=" << refD << ")" << endl;
			}
		}

		Float avgErr = 0.0f;
		for (int i=0; i<1000; ++i) {
			Float cosTheta = random->nextFloat();
			Float alpha = std::pow(random->nextFloat(), (Float) 4.0f)*4;
			Float eta = 1 + std::pow(random->nextFloat(), (Float) 4.0f)*3;
			if (cosTheta < 1e-5)
				cosTheta = 1e-5f;
			if (alpha < 1e-5)
				alpha = 1e-5f;
			if (eta < 1+1e-5)
				eta = 1+1e-5f;
			//eta = 1/eta;

			Float ref = computeTransmittance("beckmann", eta, alpha, cosTheta);
			Float dat = rtr.eval(cosTheta, alpha, eta);

			if (i % 20 == 0)
				cout << "Testing " << i << "/1000" << endl;
			if (std::abs(ref-dat) > 1e-3f) {
				cout << endl;
				cout << "eta = " << eta << endl;
				cout << "alpha = " << alpha << endl;
				cout << "cosTheta = " << cosTheta << endl;
				cout << "diff=" << dat-ref << " (dat=" << dat << ", ref=" << ref << ")" << endl;
			}

			avgErr += ref-dat;
		}
		avgErr /= 1000;
		cout << "Avg error = " << avgErr << endl;
	}

	void test03_roughTransmittanceFixedEta() {
		RoughTransmittance rtr(MicrofacetDistribution::EBeckmann);

		Float eta = 1.5f;
		rtr.setEta(eta);

		ref<Random> random = new Random();

		for (int i=0; i<50; ++i) {
			Float alpha = std::pow(random->nextFloat(), (Float) 4.0f)*4;
			if (alpha < 1e-5)
				alpha = 1e-5f;

			Float refD = computeDiffuseTransmittance("beckmann", eta, alpha);
			Float datD = rtr.evalDiffuse(alpha, eta);

			cout << "Testing " << i << "/50" << endl;
			if (std::abs(refD-datD) > 1e-3f) {
				cout << endl;
				cout << "alpha = " << alpha << endl;
				cout << "diff=" << datD-refD << " (datD=" << datD << ", ref=" << refD << ")" << endl;
			}
		}

		Float avgErr = 0.0f;
		for (int i=0; i<1000; ++i) {
			Float cosTheta = random->nextFloat();
			Float alpha = std::pow(random->nextFloat(), (Float) 4.0f)*4;
			if (cosTheta < 1e-5)
				cosTheta = 1e-5f;
			if (alpha < 1e-5)
				alpha = 1e-5f;

			Float ref = computeTransmittance("beckmann", eta, alpha, cosTheta);
			Float dat = rtr.eval(cosTheta, alpha, eta);

			if (i % 20 == 0)
				cout << "Testing " << i << "/1000" << endl;
			if (std::abs(ref-dat) > 1e-3f) {
				cout << endl;
				cout << "eta = " << eta << endl;
				cout << "alpha = " << alpha << endl;
				cout << "cosTheta = " << cosTheta << endl;
				cout << "diff=" << dat-ref << " (dat=" << dat << ", ref=" << ref << ")" << endl;
			}

			avgErr += ref-dat;
		}
		avgErr /= 1000;
		cout << "Avg error = " << avgErr << endl;
	}

	void test04_roughTransmittanceFixedEtaFixedAlpha() {
		ref<Timer> timer = new Timer();
		RoughTransmittance rtr(MicrofacetDistribution::EBeckmann);
		Float eta = 1.5f;
		Float alpha = 0.2f;
		rtr.setEta(eta);
		rtr.setAlpha(alpha);
		cout << "Loading and projecting took " << timer->getMilliseconds() << " ms" << endl;

		ref<Random> random = new Random();

		Float refD = computeDiffuseTransmittance("beckmann", eta, alpha);
		Float datD = rtr.evalDiffuse(alpha, eta);

		if (std::abs(refD-datD) > 1e-3f) {
			cout << endl;
			cout << "alpha = " << alpha << endl;
			cout << "diff=" << datD-refD << " (datD=" << datD << ", ref=" << refD << ")" << endl;
		}

		Float avgErr = 0.0f;
		for (int i=0; i<1000; ++i) {
			Float cosTheta = random->nextFloat();
			if (cosTheta < 1e-5)
				cosTheta = 1e-5f;

			Float ref = computeTransmittance("beckmann", eta, alpha, cosTheta);
			Float dat = rtr.eval(cosTheta, alpha, eta);

			if (i % 20 == 0)
				cout << "Testing " << i << "/1000" << endl;
			if (std::abs(ref-dat) > 1e-3f) {
				cout << endl;
				cout << "eta = " << eta << endl;
				cout << "alpha = " << alpha << endl;
				cout << "cosTheta = " << cosTheta << endl;
				cout << "diff=" << dat-ref << " (dat=" << dat << ", ref=" << ref << ")" << endl;
			}

			avgErr += ref-dat;
		}
		avgErr /= 1000;
		cout << "Avg error = " << avgErr << endl;
	}
};

MTS_EXPORT_TESTCASE(TestRoughTransmittance, "Testcase for rough transmittance computations")
MTS_NAMESPACE_END

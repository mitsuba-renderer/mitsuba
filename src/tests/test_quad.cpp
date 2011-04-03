/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/testcase.h>
#include <mitsuba/core/quad.h>
#include <boost/bind.hpp>

MTS_NAMESPACE_BEGIN

class TestQuadrature : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test01_quad)
	MTS_DECLARE_TEST(test02_simpson)
	MTS_END_TESTCASE()

	Float testF(Float t) const {
		return std::sin(t);
	}

	void test01_quad() {
		GaussLobattoIntegrator quad(1024, 0, 1e-6f);
		Float result = quad.integrate(boost::bind(&TestQuadrature::testF, this, _1), 0, 10);
		Float ref = 2 * std::pow(std::sin(5), 2);
		assertEqualsEpsilon(result, ref, 1e-6f);
	}

	/**
	 * Attempts to solve the following 1D integral equation for 't'
	 * \int_0^t density(ray.o + x * ray.d) * dx == desiredDensity.
	 * When no solution can be found in [0, maxDist] the function returns
	 * false. For convenience, the function returns the current values of sigmaT 
	 * and the albedo, as well as the 3D position 'ray.o+t*ray.d' upon
	 * success.
	 */
	bool invertDensityIntegral(const Ray &ray, Float desiredDensity) {
		/* Determine the ray segment, along which the
		   density integration should take place */
		Float mint, maxt;
		if (!m_aabb.rayIntersect(ray, mint, maxt))
			return false;
		Float t = std::max(mint, ray.mint);
		maxt = std::min(maxt, ray.maxt);
		Float length = maxt-mint;
		Point p = ray(t);

		if (length <= 0)
			return false;

		/* Compute a suitable step size */
		int nSteps = (int) std::ceil(length/m_stepSize);
		Float stepSize = maxDist/nSteps;
		Vector fullStep = ray.d * stepSize,
			   halfStep = fullStep * .5f;

		Float node1 = m_densities->lookupFloat(p),
			  accumulatedDensity = 0.0f;

		for (int i=0; i<nSteps; ++i) {
			Float node2 = m_densities->lookupFloat(p + halfStep),
				  node3 = m_densities->lookupFloat(p + fullStep);

			Float newDensity = accumulatedDensity 
				+ (node1+node2*4+node3) * (stepSize/6) * m_densityMultiplier;

			if (newDensity >= desiredDensity) {
				/* Use Newton-Bisection to find the root */
				Float a = 0, b = stepSize, x = a,
					  fa = accumulatedDensity - desiredDensity,
					  fb = newDensity - desiredDensity,
					  fx = fa;

				int it = 1;
				while (true) {
					Float dfx = df(x);

					if (std::abs(fx) < 1e-10) {
						cout << "Found root at " << x << ", fx=" << fx << endl;
						break;
					}

					x = x - fx/dfx;

					if (x <= a || x >= b || dfx == 0) {
						cout << "Doing a bisection step!" << endl;
						x = 0.5f * (b + a);
					}

					fx = f(x);

					if (fx * fa < 0)
						b = x;
					else
						a = x;
				}
			}

			Point next = p + fullStep;
			if (p == next) {
				Log(EWarn, "invertDensityIntegral(): unable to make forward progress -- "
						"round-off error issues? The step size was %f", stepSize);
				break;
			}
			accumulatedDensity + newDensity;
			node1 = node3;
			t += stepSize;
			p = next;
		}

		return false;
	}


	void test02_simpson() {
		Float mint = 0;
		Float maxt = .921;
		cout << integrateDensity(Ray(Point(0, 0, 0), Vector(0, 10, 0), mint, maxt, 0)) << endl;
	}
};

MTS_EXPORT_TESTCASE(TestQuadrature, "Testcase for quadrature routines")
MTS_NAMESPACE_END

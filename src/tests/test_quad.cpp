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
};

MTS_EXPORT_TESTCASE(TestQuadrature, "Testcase for quadrature routines")
MTS_NAMESPACE_END

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

MTS_NAMESPACE_BEGIN

class TestLinearAlgebra : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test01_eigenDecomp)
	MTS_END_TESTCASE()

	void test01_eigenDecomp() {
		ref<Matrix4x4> A = new Matrix4x4(
			1.4541, 1.1233, 1.2407, 1.2548,
			1.1233, 0.2597, 0.3819, 1.3917,
			1.2407, 0.3819, 1.1552, 1.1048,
			1.2548, 1.3917, 1.1048, 1.4712
		);
		ref<Matrix4x4> Q = new Matrix4x4();
		Vector4 d;
		A->symmEigenDecomp(Q, d);

		Vector4 refD(-0.823889076095475, 0.130902702868822,
			0.557486242256414, 4.47570013097024);

		ref<Matrix4x4> refQ = new Matrix4x4(
			0.294383137629217, 0.746207030711883, 0.191065628242818, 0.565692108217809,
			-0.789565156329591, 0.139585328270248, -0.459560990857043, 0.381977087957391,
			-0.286639718342894, -0.415237779451033, 0.738328701906588, 0.447533223711729,
			0.455810381691897, -0.501266984714151, -0.45515749947419,  0.577754235510108
		);

		assertEqualsEpsilon(d, refD, 1e-6);
		assertEqualsEpsilon(Q, refQ, 1e-6);
	}
};

MTS_EXPORT_TESTCASE(TestLinearAlgebra, "Testcase for Linear Algebra routines")
MTS_NAMESPACE_END

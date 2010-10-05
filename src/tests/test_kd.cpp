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
#include <mitsuba/render/gkdtree.h>

MTS_NAMESPACE_BEGIN

class TestKDTree : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test01_sutherlandHodgman)
	MTS_DECLARE_TEST(test02_buildSimple)
	MTS_END_TESTCASE()

	void test01_sutherlandHodgman() {
		/* Test the triangle clipping algorithm */
		Vertex vertices[3];
		vertices[0].p = Point(0, 0, 0);
		vertices[1].p = Point(1, 0, 0);
		vertices[2].p = Point(1, 1, 0);

		AABB clipAABB(
			Point(0, .5, -1),
			Point(1, 1, 1)
		);
		
		AABB expectedResult(
			Point(.5, .5, 0),
			Point(1, 1, 0)
		);

		Triangle t;
		t.idx[0] = 0; t.idx[1] = 1; t.idx[2] = 2;

		AABB clippedAABB = t.getClippedAABB(vertices, clipAABB);

		assertEquals(clippedAABB.min, expectedResult.min);
		assertEquals(clippedAABB.max, expectedResult.max);
	}

	class TriKDTree : GenericKDTree<Triangle> {
	};

	void test02_buildSimple() {
		TriKDTree tree;
		std::vector<Triangle> tris;
	}
};

MTS_EXPORT_TESTCASE(TestKDTree, "Testcase for kd-tree related code")
MTS_NAMESPACE_END

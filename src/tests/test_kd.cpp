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

#include <mitsuba/core/plugin.h>
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
		/* Test the triangle clipping algorithm on the unit triangle */
		Vertex vertices[3];
		vertices[0].p = Point(0, 0, 0);
		vertices[1].p = Point(1, 0, 0);
		vertices[2].p = Point(1, 1, 0);
		Triangle t;
		t.idx[0] = 0; t.idx[1] = 1; t.idx[2] = 2;

		/* Split the triangle in half and verify the clipped AABB */
		AABB clippedAABB = t.getClippedAABB(vertices, AABB(
			Point(0, .5, -1),
			Point(1, 1, 1)
		));

		assertEquals(Point(.5, .5, 0), clippedAABB.min);
		assertEquals(Point(1, 1, 0), clippedAABB.max);

		/* Verify that a triangle can be completely clipped away */
		clippedAABB = t.getClippedAABB(vertices, AABB(
			Point(2, 2, 2),
			Point(3, 3, 3)
		));
		assertFalse(clippedAABB.isValid());

		/* Verify that a no clipping whatsoever happens when 
		   the AABB fully contains a triangle */
		clippedAABB = t.getClippedAABB(vertices, AABB(
			Point(-1, -1, -1),
			Point(1, 1, 1)
		));
		assertEquals(Point(0, 0, 0), clippedAABB.min);
		assertEquals(Point(1, 1, 0), clippedAABB.max);

		/* Verify that a triangle within a flat cell won't be clipped away */
		clippedAABB = t.getClippedAABB(vertices, AABB(
			Point(-100,-100, 0),
			Point( 100, 100, 0)
		));
		assertEquals(Point(0, 0, 0), clippedAABB.min);
		assertEquals(Point(1, 1, 0), clippedAABB.max);

		/* Verify that a triangle just touching the clip AABB leads to a
		   collapsed point AABB */
		clippedAABB = t.getClippedAABB(vertices, AABB(
			Point(0,1, 0),
			Point(1,2, 0)
		));
		assertEquals(Point(1, 1, 0), clippedAABB.min);
		assertEquals(Point(1, 1, 0), clippedAABB.max);
	}

	void test02_buildSimple() {
		Properties bunnyProps("ply");
		bunnyProps.setString("filename", "tools/tests/bunny.ply");

		ref<TriMesh> mesh = static_cast<TriMesh *> (PluginManager::getInstance()->
				createObject(TriMesh::m_theClass, bunnyProps));
		mesh->configure();
		ref<KDTree> tree = new KDTree();
		tree->addShape(mesh);
		tree->build();
		BSphere bsphere(mesh->getBSphere());

		//Float intersectionCost, traversalCost;
		//tree.findCosts(intersectionCost, traversalCost);

		bsphere = BSphere(Point(-0.016840, 0.110154, -0.001537), .2f);

		for (int j=0; j<3; ++j) {
			ref<Random> random = new Random();
			ref<Timer> timer = new Timer();
			size_t nRays = 10000000;
			size_t nIntersections = 0;

			Log(EInfo, "Bounding sphere: %s", bsphere.toString().c_str());
			Log(EInfo, "Shooting " SIZE_T_FMT " rays ..", nRays);

			for (size_t i=0; i<nRays; ++i) {
				Point2 sample1(random->nextFloat(), random->nextFloat()),
					sample2(random->nextFloat(), random->nextFloat());
				Point p1 = bsphere.center + squareToSphere(sample1) * bsphere.radius;
				Point p2 = bsphere.center + squareToSphere(sample2) * bsphere.radius;
				Ray r(p1, normalize(p2-p1));
				Intersection its;

				if (tree->rayIntersect(r, its))
					nIntersections++;
			}

			Log(EInfo, "Found " SIZE_T_FMT " intersections in %i ms",
				nIntersections, timer->getMilliseconds());
			Log(EInfo, "%.3f MRays/s", 
				nRays / (timer->getMilliseconds() * (Float) 1000));
		}
	}
};

MTS_EXPORT_TESTCASE(TestKDTree, "Testcase for kd-tree related code")
MTS_NAMESPACE_END

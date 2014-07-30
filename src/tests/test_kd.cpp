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
#include <mitsuba/core/kdtree.h>
#include <mitsuba/render/testcase.h>
#include <mitsuba/render/skdtree.h>

MTS_NAMESPACE_BEGIN

class TestKDTree : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test01_sutherlandHodgman)
	MTS_DECLARE_TEST(test02_bunnyBenchmark)
	MTS_DECLARE_TEST(test03_pointKDTree)
	MTS_END_TESTCASE()

	void test01_sutherlandHodgman() {
		/* Test the triangle clipping algorithm on the unit triangle */
		Point vertices[3];
		vertices[0] = Point(0, 0, 0);
		vertices[1] = Point(1, 0, 0);
		vertices[2] = Point(1, 1, 0);
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

	void test02_bunnyBenchmark() {
		Properties bunnyProps("ply");
		bunnyProps.setString("filename", "data/tests/bunny.ply");

		PluginManager *pmgr = PluginManager::getInstance();
		ref<TriMesh> mesh = static_cast<TriMesh *> (
				pmgr->createObject(MTS_CLASS(TriMesh), bunnyProps));
		mesh->addChild(pmgr->createObject(Properties("diffuse")));
		mesh->configure();
		ref<ShapeKDTree> tree = new ShapeKDTree();
		tree->addShape(mesh);
		tree->build();
		BSphere bsphere = tree->getAABB().getBSphere();

		bsphere = BSphere(Point(-0.016840, 0.110154, -0.001537), .2f);

		Log(EInfo, "Bunny benchmark (http://homepages.paradise.net.nz/nickamy/benchmark.html):");
		ref<Timer> timer = new Timer();
		for (int j=0; j<3; ++j) {
			ref<Random> random = new Random();
			size_t nRays = 10000000;
			size_t nIntersections = 0;

			Log(EInfo, "  Iteration %i: shooting " SIZE_T_FMT " rays ..", j, nRays);
			timer->reset();

			for (size_t i=0; i<nRays; ++i) {
				Point2 sample1(random->nextFloat(), random->nextFloat()),
					sample2(random->nextFloat(), random->nextFloat());
				Point p1 = bsphere.center + warp::squareToUniformSphere(sample1) * bsphere.radius;
				Point p2 = bsphere.center + warp::squareToUniformSphere(sample2) * bsphere.radius;
				Ray r(p1, normalize(p2-p1), 0.0f);
				Intersection its;

				if (tree->rayIntersect(r))
					nIntersections++;
			}

			Float perc = nIntersections/(Float) nRays;
			Log(EInfo, "  Found " SIZE_T_FMT " intersections (%.3f%%) in %i ms",
				nIntersections, perc, timer->getMilliseconds());
			Log(EInfo, "  -> %.3f MRays/s",
				nRays / (timer->getMilliseconds() * (Float) 1000));
			Log(EInfo, "");
		}
	}

	void test03_pointKDTree() {
		typedef PointKDTree< SimpleKDNode<Point2, Float> > KDTree2;
		typedef PointKDTree< LeftBalancedKDNode<Point2, Float> > KDTree2Left;

		size_t nPoints = 50000, nTries = 20;
		ref<Random> random = new Random();

		for (int heuristic=0; heuristic<4; ++heuristic) {
			KDTree2 kdtree(nPoints, (KDTree2::EHeuristic) heuristic);

			for (size_t i=0; i<nPoints; ++i) {
				kdtree[i].setPosition(Point2(random->nextFloat(), random->nextFloat()));
				kdtree[i].setData(random->nextFloat());
			}

			KDTree2::SearchResult results[11];
			std::vector<KDTree2::SearchResult> resultsBF;

			if (heuristic == 0) {
				Log(EInfo, "Testing the balanced kd-tree construction heuristic");
			} else if (heuristic == 1) {
				Log(EInfo, "Testing the left-balanced kd-tree construction heuristic");
			} else if (heuristic == 2) {
				Log(EInfo, "Testing the sliding midpoint kd-tree construction heuristic");
			} else if (heuristic == 3) {
				Log(EInfo, "Testing the voxel volume kd-tree construction heuristic");
			}

			ref<Timer> timer = new Timer();
			kdtree.build(true);
			Log(EInfo, "Construction time = %i ms, depth = %i", timer->getMilliseconds(), kdtree.getDepth());

			for (int k=1; k<=10; ++k) {
				size_t nTraversals = 0;
				for (size_t it = 0; it < nTries; ++it) {
					Point2 p(random->nextFloat(), random->nextFloat());
					Float searchRadius = std::numeric_limits<Float>::infinity();
					kdtree.nnSearchCollectStatistics(p, searchRadius, k, results, nTraversals);
					resultsBF.clear();
					for (size_t j=0; j<nPoints; ++j)
						resultsBF.push_back(KDTree2::SearchResult((kdtree[j].getPosition()-p).lengthSquared(), (uint32_t) j));
					std::sort(results, results + k, KDTree2::SearchResultComparator());
					std::sort(resultsBF.begin(), resultsBF.end(), KDTree2::SearchResultComparator());
					for (int j=0; j<k; ++j)
						assertTrue(results[j] == resultsBF[j]);
				}
				Log(EInfo, "Average number of traversals for a %i-nn query = " SIZE_T_FMT, k, nTraversals / nTries);
			}
		}
		Log(EInfo, "Testing the left-balanced kd-tree construction heuristic with left-balanced nodes");
		KDTree2Left kdtree(nPoints, KDTree2Left::ELeftBalanced);

		for (size_t i=0; i<nPoints; ++i) {
			kdtree[i].setPosition(Point2(random->nextFloat(), random->nextFloat()));
			kdtree[i].setData(random->nextFloat());
		}

		std::vector<KDTree2Left::SearchResult> resultsBF;
		KDTree2Left::SearchResult results[11];

		ref<Timer> timer = new Timer();
		kdtree.build(true);
		Log(EInfo, "Construction time = %i ms, depth = %i", timer->getMilliseconds(), kdtree.getDepth());

		for (int k=1; k<=10; ++k) {
			size_t nTraversals = 0;
			for (size_t it = 0; it < nTries; ++it) {
				Point2 p(random->nextFloat(), random->nextFloat());
				Float searchRadius = std::numeric_limits<Float>::infinity();
				kdtree.nnSearchCollectStatistics(p, searchRadius, k, results, nTraversals);
				resultsBF.clear();
				for (size_t j=0; j<nPoints; ++j)
					resultsBF.push_back(KDTree2Left::SearchResult((kdtree[j].getPosition()-p).lengthSquared(), (uint32_t) j));
				std::sort(results, results + k, KDTree2Left::SearchResultComparator());
				std::sort(resultsBF.begin(), resultsBF.end(), KDTree2Left::SearchResultComparator());
				for (int j=0; j<k; ++j)
					assertTrue(results[j] == resultsBF[j]);
			}
			Log(EInfo, "Average number of traversals for a %i-nn query = " SIZE_T_FMT, k, nTraversals / nTries);
		}

		Log(EInfo, "Normal node size = " SIZE_T_FMT " bytes", sizeof(KDTree2::NodeType));
		Log(EInfo, "Left-balanced node size = " SIZE_T_FMT " bytes", sizeof(KDTree2Left::NodeType));
	}
};

MTS_EXPORT_TESTCASE(TestKDTree, "Testcase for kd-tree related code")
MTS_NAMESPACE_END

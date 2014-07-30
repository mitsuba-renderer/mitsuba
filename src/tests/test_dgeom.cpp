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
#include <mitsuba/render/testcase.h>
#include <mitsuba/render/skdtree.h>

MTS_NAMESPACE_BEGIN

class TestDGeom : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test01_trimesh_1);
	MTS_DECLARE_TEST(test02_trimesh_2);
	MTS_DECLARE_TEST(test03_trimesh_3);
	MTS_DECLARE_TEST(test04_sphere);
	MTS_DECLARE_TEST(test05_cylinder);
	MTS_END_TESTCASE()

	void test01_trimesh_1() {
		/* No normals, no UV coordinates */
		ref<TriMesh> trimesh = new TriMesh("", 1, 3);
		Triangle &tri = trimesh->getTriangles()[0];
		tri.idx[0] = 0; tri.idx[1] = 1; tri.idx[2] = 2;
		Point *vertices = trimesh->getVertexPositions();
		vertices[0] = Point(0, 0, 0);
		vertices[1] = Point(1, 0, 0);
		vertices[2] = Point(0, 1, 0);
		trimesh->configure();

		ref<ShapeKDTree> kdtree = new ShapeKDTree();
		kdtree->addShape(trimesh);
		kdtree->build();

		Intersection its;
		Ray ray(Point(0.1f, 0.2f, -1.0f), Vector(0, 0, 1), 123.0f);

		assertTrue(kdtree->rayIntersect(ray, its));
		assertEquals(its.p, Point(0.1f, 0.2f, 0.0f));
		assertEquals(its.uv, Point2(0.1f, 0.2f));
		assertEquals(its.shFrame.n, Normal(0, 0, 1));
		assertEquals(its.geoFrame.n, Normal(0, 0, 1));
		assertEquals(its.dpdu, Vector(1.0f, 0.0f, 0.0f));
		assertEquals(its.dpdv, Vector(0.0f, 1.0f, 0.0f));
		assertEquals(its.time, 123.0f);

		Vector dndu, dndv;
		its.shape->getNormalDerivative(its, dndu, dndv, false);
		assertEquals(dndu, Vector(0.0f)); assertEquals(dndv, Vector(0.0f));
		its.shape->getNormalDerivative(its, dndu, dndv, true);
		assertEquals(dndu, Vector(0.0f)); assertEquals(dndv, Vector(0.0f));
	}

	void test02_trimesh_2() {
		/* Shading normals & UV coords, but no explicit parameterization */
		ref<TriMesh> trimesh = new TriMesh("", 1, 3, true, true);
		Triangle &tri = trimesh->getTriangles()[0];
		tri.idx[0] = 0; tri.idx[1] = 1; tri.idx[2] = 2;
		Point *vertices = trimesh->getVertexPositions();
		Normal *normals = trimesh->getVertexNormals();
		Point2 *uv = trimesh->getVertexTexcoords();
		vertices[0] = Point(0, 0, 0);
		vertices[1] = Point(1, 0, 0);
		vertices[2] = Point(0, 1, 0);

		normals[0] = Normal(-0.3f, 0, 1);
		normals[1] = Normal(0.3f, 0, 1);
		normals[2] = Normal(0, 0.3f, 1);

		uv[0] = Point2(0.1f, 0.1f);
		uv[1] = Point2(1.1f, 0.1f);
		uv[2] = Point2(0.1f, 0.9f);

		trimesh->configure();

		ref<ShapeKDTree> kdtree = new ShapeKDTree();
		kdtree->addShape(trimesh);
		kdtree->build();

		Intersection its;
		Ray ray(Point(0.1f, 0.2f, -1.0f), Vector(0, 0, 1), 123.0f);

		assertTrue(kdtree->rayIntersect(ray, its));
		assertEquals(its.p, Point(0.1f, 0.2f, 0.0f));
		assertEqualsEpsilon(its.uv, Point2(0.2f, 0.26f), Epsilon);
		assertEquals(its.geoFrame.n, Normal(0, 0, 1));

		assertEqualsEpsilon(its.shFrame.n,
			normalize(normals[0]*.7f + normals[1]*.1f + normals[2]*.2f), Epsilon);

		its.shFrame.s = normalize(its.dpdu - its.shFrame.n
						* dot(its.shFrame.n, its.dpdu));

		assertEqualsEpsilon(its.dpdu, vertices[1]-vertices[0], Epsilon);
		assertEqualsEpsilon(its.dpdv, vertices[2]-vertices[0], Epsilon);
		assertEquals(its.time, 123.0f);

		Vector dndu, dndv;
		its.shape->getNormalDerivative(its, dndu, dndv, false);
		assertEquals(dndu, Vector(0.0f)); assertEquals(dndv, Vector(0.0f));
		its.shape->getNormalDerivative(its, dndu, dndv, true);
		// From mathematica
		assertEqualsEpsilon(dndu, Vector(0.571048f, 0.00614519f, 0.10242f), 1e-4f);
		assertEqualsEpsilon(dndv, Vector(0.288596f, 0.29679f, 0.0341399f), 1e-4f);
	}

	void test03_trimesh_3() {
		/* Shading normals & UV coords, but with an explicit parameterization */
		ref<TriMesh> trimesh = new TriMesh("", 1, 3, true, true);
		Triangle &tri = trimesh->getTriangles()[0];
		tri.idx[0] = 0; tri.idx[1] = 1; tri.idx[2] = 2;
		Point *vertices = trimesh->getVertexPositions();
		Normal *normals = trimesh->getVertexNormals();
		Point2 *uv = trimesh->getVertexTexcoords();
		vertices[0] = Point(0, 0, 0);
		vertices[1] = Point(1, 0, 0);
		vertices[2] = Point(0, 1, 0);

		normals[0] = Normal(-0.3f, 0, 1);
		normals[1] = Normal(0.3f, 0, 1);
		normals[2] = Normal(0, 0.3f, 1);

		uv[0] = Point2(0.0f, 0.0f);
		uv[1] = Point2(0.0f, 1.0f);
		uv[2] = Point2(1.0f, 0.0f);

		trimesh->configure();
		trimesh->computeUVTangents();

		ref<ShapeKDTree> kdtree = new ShapeKDTree();
		kdtree->addShape(trimesh);
		kdtree->build();

		Intersection its;
		Ray ray(Point(0.1f, 0.2f, -1.0f), Vector(0, 0, 1), 123.0f);

		assertTrue(kdtree->rayIntersect(ray, its));
		assertEquals(its.p, Point(0.1f, 0.2f, 0.0f));
		assertEqualsEpsilon(its.uv, Point2(0.2f, 0.1f), Epsilon);
		assertEquals(its.geoFrame.n, Normal(0, 0, 1));

		assertEqualsEpsilon(its.shFrame.n,
			normalize(normals[0]*.7f + normals[1]*.1f + normals[2]*.2f), Epsilon);

		assertEqualsEpsilon(its.shFrame.s,
			normalize(its.dpdu - its.shFrame.n * dot(its.dpdu, its.shFrame.n)), Epsilon);

		its.shFrame.s = normalize(its.dpdu - its.shFrame.n
						* dot(its.shFrame.n, its.dpdu));

		assertEquals(its.dpdu, vertices[2]-vertices[0]);
		assertEquals(its.dpdv, vertices[1]-vertices[0]);
		assertEquals(its.time, 123.0f);

		Vector dndu, dndv;
		its.shape->getNormalDerivative(its, dndu, dndv, false);
		assertEquals(dndu, Vector(0.0f)); assertEquals(dndv, Vector(0.0f));
		its.shape->getNormalDerivative(its, dndu, dndv, true);
		// From mathematica
		assertEqualsEpsilon(dndu, Vector(0.288596f, 0.29679f, 0.0341399f), 1e-4f);
		assertEqualsEpsilon(dndv, Vector(0.571048f, 0.00614519f, 0.10242f), 1e-4f);
	}

	void test04_sphere() {
		Properties props("sphere");
		Float radius = 2.0f;
		props.setFloat("radius", radius);
		props.setPoint("center", Point(0.0f));
		ref<Shape> sphere = static_cast<Shape *>(PluginManager::getInstance()->createObject(props));
		sphere->configure();

		ref<ShapeKDTree> kdtree = new ShapeKDTree();
		kdtree->addShape(sphere);
		kdtree->build();

		Ray ray(Point(3, 3, 3), normalize(Vector(-1, -1, -1)), 123.0f);
		Intersection its;

		assertTrue(kdtree->rayIntersect(ray, its));

		assertEquals(its.time, 123.0f);
		assertEqualsEpsilon(its.p, Point(ray.d*-2.0f), Epsilon);
		assertEqualsEpsilon(its.wi, Vector(0, 0, 1), Epsilon);
		assertEqualsEpsilon(its.shFrame.n, -ray.d, Epsilon);
		assertTrue(its.shFrame == its.geoFrame);

		Point2 sc = toSphericalCoordinates(Vector(its.p));
		std::swap(sc.x, sc.y);
		sc.x *= INV_TWOPI;
		sc.y *= INV_PI;

		assertEqualsEpsilon(its.uv, sc, Epsilon);

		// from mathematica
		assertEqualsEpsilon(its.dpdu, Vector(-3.6276f, 3.6276f, 0.0f)*radius, 1e-4f);
		assertEqualsEpsilon(its.dpdv, Vector(1.28255f, 1.28255f, -2.5651f)*radius, 1e-4f);

		Vector dndu, dndv;
		its.shape->getNormalDerivative(its, dndu, dndv, false);
		assertEqualsEpsilon(dndu, Vector(-3.6276f, 3.6276f, 0.0f), 1e-4f);
		assertEqualsEpsilon(dndv, Vector(1.28255f, 1.28255f, -2.5651f), 1e-4f);

		Float H, K;
		its.shape->getCurvature(its, H, K);

		assertEquals(K, 1.0f / (radius*radius));
		assertEquals(H, -1.0f / radius);
	}

	void test05_cylinder() {
		Properties props("cylinder");
		Float radius = 2.0f;
		props.setFloat("radius", radius);
		props.setPoint("p0", Point(0.0f, 0.0f, -3.0f));
		props.setPoint("p1", Point(0.0f, 0.0f,  3.0f));
		ref<Shape> sphere = static_cast<Shape *>(
				PluginManager::getInstance()->createObject(props));
		sphere->configure();

		ref<ShapeKDTree> kdtree = new ShapeKDTree();
		kdtree->addShape(sphere);
		kdtree->build();

		Ray ray(Point(3.0f, 3.0f, 3.0f), normalize(Vector(-1, -1, -1)), 123.0f);
		Intersection its;

		assertTrue(kdtree->rayIntersect(ray, its));
		Float t = 3*std::sqrt(3.0f)-std::sqrt(6.0f);
		assertEquals(its.time, 123.0f);

		assertEqualsEpsilon(its.p, ray(t), 1e-4f);
		assertEqualsEpsilon(its.t, t, 1e-4f);

		assertEqualsEpsilon(its.shFrame.n, normalize(Normal(1, 1, 0)), Epsilon);
		assertTrue(its.shFrame == its.geoFrame);
		assertEqualsEpsilon(its.uv.x, 1.0f / 8.0f, Epsilon);
		assertEqualsEpsilon(its.uv.y, (its.p.z+3)/6, Epsilon);

		// from mathematica
		assertEqualsEpsilon(its.dpdu, Vector(-8.88577, 8.88577, 0), 1e-4f);
		assertEqualsEpsilon(its.dpdv, Vector(0, 0, 6), 1e-4f);

		Vector dndu, dndv;
		its.shape->getNormalDerivative(its, dndu, dndv, false);
		assertEqualsEpsilon(dndu, Vector(-4.44288f, 4.44288f, 0.0f), 1e-4f);
		assertEqualsEpsilon(dndv, Vector(0.0f), 1e-4f);

		Float H, K;
		its.shape->getCurvature(its, H, K);
		assertEqualsEpsilon(K, 0.0f, Epsilon);
		assertEqualsEpsilon(H, -1.0f / (2*radius), Epsilon);
	}
};

MTS_EXPORT_TESTCASE(TestDGeom, "Differential geometry testcase")
MTS_NAMESPACE_END

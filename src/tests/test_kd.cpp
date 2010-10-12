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

class TriKDTree : public GenericKDTree<TriKDTree> {
public:
	TriKDTree(const Triangle *triangles,
				const Vertex *vertexBuffer,
				size_type triangleCount) 
		: m_triangles(triangles),
		  m_vertexBuffer(vertexBuffer),
		  m_triangleCount(triangleCount) {
	}
	
	bool rayIntersect(const Ray &ray, Intersection &its) const {
		uint32_t temp[MTS_KD_INTERSECTION_TEMP];
		its.t = std::numeric_limits<Float>::infinity(); 
		Float mint, maxt;
		TriAccelHandler handler(m_indices, m_triAccel);

		if (m_aabb.rayIntersect(ray, mint, maxt)) {
			if (ray.mint > mint) mint = ray.mint;
			if (ray.maxt < maxt) maxt = ray.maxt;

			if (EXPECT_TAKEN(maxt > mint)) {
				if (rayIntersectHavranCustom(handler,
							ray, mint, maxt, its.t, temp)) {
					fillIntersectionDetails(ray, its.t, temp, its);
					return true;
				}
			}
		}
		return false;
	}

	FINLINE AABB getAABB(index_type idx) const {
		return m_triangles[idx].getAABB(m_vertexBuffer);
	}

	FINLINE AABB getClippedAABB(index_type idx, const AABB &aabb) const {
		return m_triangles[idx].getClippedAABB(m_vertexBuffer, aabb);
	}

	FINLINE size_type getPrimitiveCount() const {
		return m_triangleCount;
	}
#if 0
	FINLINE EIntersectionResult intersect(const Ray &ray, index_type idx,
			Float mint, Float maxt, Float &t, void *tmp) const {
		Float tempT, tempU, tempV;
#if 0
		if (m_triangles[idx].rayIntersect(m_vertexBuffer, ray, tempU, tempV, tempT)) {
			if (tempT < mint && tempT > maxt)
				return ENo;
#else
		if (m_triAccel[idx].rayIntersect(ray, mint, maxt, tempU, tempV, tempT)) {
#endif
			index_type *indexPtr = reinterpret_cast<index_type *>(tmp);
			Float *floatPtr = reinterpret_cast<Float *>(indexPtr + 1);

			t = tempT;
			*indexPtr = idx;
			*floatPtr++ = tempU;
			*floatPtr++ = tempV;

			return EYes;
		}
		return ENo;
	}
#endif

	void build() {
		buildInternal();
		Log(EInfo, "Precomputing triangle intersection information (%s)",
				memString(sizeof(TriAccel)*m_triangleCount).c_str());
		m_triAccel = static_cast<TriAccel *>(allocAligned(m_triangleCount * sizeof(TriAccel)));
		for (size_t i=0; i<m_triangleCount; ++i) {
			const Triangle &tri = m_triangles[i];
			const Vertex &v0 = m_vertexBuffer[tri.idx[0]];
			const Vertex &v1 = m_vertexBuffer[tri.idx[1]];
			const Vertex &v2 = m_vertexBuffer[tri.idx[2]];

			m_triAccel[i].load(v0.p, v1.p, v2.p);
		}
	}

	FINLINE void fillIntersectionDetails(const Ray &ray, 
			Float t, const void *tmp, Intersection &its) const {
		its.p = ray(t);

		//const index_type *indexPtr = reinterpret_cast<const index_type *>(tmp);
		//const Float *floatPtr = reinterpret_cast<const Float *>(indexPtr + 1);

		//const Triangle &tri = m_triangles[*indexPtr];
		//const Vertex &v0 = m_vertexBuffer[tri.idx[0]];
		//const Vertex &v1 = m_vertexBuffer[tri.idx[1]];
		//const Vertex &v2 = m_vertexBuffer[tri.idx[2]];

		//const Float u = *floatPtr++, v = *floatPtr++;
		//const Vector b(1 - u - v, u, v);
/*
		its.uv = v0.uv * b.x + v1.uv * b.y + v2.uv * b.z;
		its.dpdu = v0.dpdu * b.x + v1.dpdu * b.y + v2.dpdu * b.z;
		its.dpdv = v0.dpdv * b.x + v1.dpdv * b.y + v2.dpdv * b.z;

		its.geoFrame = Frame(normalize(cross(v1.p-v0.p, v2.p-v0.p)));
		its.shFrame.n = normalize(v0.n * b.x + v1.n * b.y + v2.n * b.z);
		its.shFrame.s = normalize(its.dpdu - its.shFrame.n 
			* dot(its.shFrame.n, its.dpdu));
		its.shFrame.t = cross(its.shFrame.n, its.shFrame.s);
		its.wi = its.toLocal(-ray.d);
*/
		its.hasUVPartials = false;
	}

protected:
	struct TriAccelHandler {
		FINLINE TriAccelHandler(const index_type *indices, const TriAccel *triAccel)
			: m_indices(indices), m_triAccel(triAccel) { }

		FINLINE bool operator()(const KDNode *node, const Ray &ray, Float searchStart,
			Float searchEnd, Float &t, void *temp) const {
			bool foundIntersection = false;
			for (unsigned int entry=node->getPrimStart(),
					last = node->getPrimEnd(); entry != last; entry++) {
				const index_type primIdx = m_indices[entry];
				Float tempT, tempU, tempV;

				if (m_triAccel[primIdx].rayIntersect(ray, searchStart, searchEnd, tempU, tempV, tempT)) {
					index_type *indexPtr = reinterpret_cast<index_type *>(temp);
					Float *floatPtr = reinterpret_cast<Float *>(indexPtr + 1);
					t = searchEnd = tempT;
					*indexPtr = primIdx;
					*floatPtr++ = tempU;
					*floatPtr++ = tempV;
					foundIntersection = true;
				}
			}
			return foundIntersection;
		}
	private:
		const index_type *m_indices;
		const TriAccel *m_triAccel;
	};

private:
	const Triangle *m_triangles;
	const Vertex *m_vertexBuffer;
	TriAccel *m_triAccel;
	size_type m_triangleCount;
};


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
		TriKDTree tree(mesh->getTriangles(), 
			mesh->getVertexBuffer(), mesh->getTriangleCount());
		tree.build();
		BSphere bsphere(mesh->getBSphere());

		//Float intersectionCost, traversalCost;
		//tree.findCosts(intersectionCost, traversalCost);

		ref<KDTree> oldTree = new KDTree();
		oldTree->addShape(mesh);
		oldTree->build(); 

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

				if (tree.rayIntersect(r, its))
					nIntersections++;
			}

			Log(EInfo, "New: Found " SIZE_T_FMT " intersections in %i ms",
				nIntersections, timer->getMilliseconds());
			Log(EInfo, "New: %.3f MRays/s", 
				nRays / (timer->getMilliseconds() * (Float) 1000));

			random = new Random();
			timer->reset();
			nIntersections=0;

			for (size_t i=0; i<nRays; ++i) {
				Point2 sample1(random->nextFloat(), random->nextFloat()),
					sample2(random->nextFloat(), random->nextFloat());
				Point p1 = bsphere.center + squareToSphere(sample1) * bsphere.radius;
				Point p2 = bsphere.center + squareToSphere(sample2) * bsphere.radius;
				Ray r(p1, normalize(p2-p1));
				Intersection its;

				if (oldTree->rayIntersect(r, its))
					nIntersections++;
			}

			Log(EInfo, "Old: Found " SIZE_T_FMT " intersections in %i ms",
				nIntersections, timer->getMilliseconds());
			Log(EInfo, "Old: %.3f MRays/s", 
				nRays / (timer->getMilliseconds() * (Float) 1000));
		}
	}
};

MTS_EXPORT_TESTCASE(TestKDTree, "Testcase for kd-tree related code")
MTS_NAMESPACE_END

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

	class TriKDTree : GenericKDTree<Triangle> {
	public:
		struct AABBFunctor {
			inline AABBFunctor(const Triangle *triangles,
					           const Vertex *vertexBuffer,
							   uint32_t primitiveCount)
				: m_triangles(triangles),
				  m_vertexBuffer(vertexBuffer),
				  m_primitiveCount(primitiveCount) { }

			inline uint32_t getPrimitiveCount() const {
				return m_primitiveCount;
			}

			inline AABB operator()(uint32_t idx) const {
				return m_triangles[idx].getAABB(m_vertexBuffer);
			}

			const Triangle *m_triangles;
			const Vertex *m_vertexBuffer;
			const uint32_t m_primitiveCount;
		};

		void build(const Triangle *triangles,
				   const Vertex *vertexBuffer,
				   uint32_t primitiveCount) {
			GenericKDTree<Triangle>::build(
				AABBFunctor(triangles, vertexBuffer, primitiveCount)
			);
		}
	private:
	};


	void test02_buildSimple() {
		Properties bunnyProps("ply");
		bunnyProps.setString("filename", "tools/tests/bunny.ply");

		ref<TriMesh> mesh = static_cast<TriMesh *> (PluginManager::getInstance()->
				createObject(TriMesh::m_theClass, bunnyProps));
		mesh->configure();
		TriKDTree tree;
		tree.build(mesh->getTriangles(), mesh->getVertexBuffer(), mesh->getTriangleCount());
	}
};

MTS_EXPORT_TESTCASE(TestKDTree, "Testcase for kd-tree related code")
MTS_NAMESPACE_END

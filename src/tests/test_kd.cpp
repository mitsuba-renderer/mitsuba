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
	friend class GenericKDTree<TriKDTree>;
public:
	TriKDTree() {
		m_shapeMap.push_back(0);
	}

	void addShape(const Shape *shape) {
		Assert(!isBuilt());
		if (shape->isCompound())
			Log(EError, "Cannot add compound shapes to a kd-tree - expand them first!");
		if (shape->getClass()->derivesFrom(TriMesh::m_theClass)) {
			// Triangle meshes are expanded into individual primitives,
			// which are visible to the tree construction code. Generic
			// primitives are only handled by their AABBs
			m_shapeMap.push_back((size_type) 
				static_cast<const TriMesh *>(shape)->getTriangleCount());
			m_triangleFlag.push_back(true);
		} else {
			m_shapeMap.push_back(1);
			m_triangleFlag.push_back(false);
		}
		m_shapes.push_back(shape);
	}

	void build() {
		for (size_t i=1; i<m_shapeMap.size(); ++i)
			m_shapeMap[i] += m_shapeMap[i-1];

		buildInternal();

		ref<Timer> timer = new Timer();
		size_type primCount = getPrimitiveCount();
		Log(EDebug, "Precomputing triangle intersection information (%s)",
				memString(sizeof(TriAccel)*primCount).c_str());
		m_triAccel = static_cast<TriAccel *>(allocAligned(primCount * sizeof(TriAccel)));

		index_type idx = 0;
		for (index_type i=0; i<m_shapes.size(); ++i) {
			const Shape *shape = m_shapes[i];
			if (m_triangleFlag[i]) {
				const TriMesh *mesh = static_cast<const TriMesh *>(shape);
				const Triangle *triangles = mesh->getTriangles();
				const Vertex *vertexBuffer = mesh->getVertexBuffer();
				for (index_type j=0; j<mesh->getTriangleCount(); ++j) {
					const Triangle &tri = triangles[j];
					const Vertex &v0 = vertexBuffer[tri.idx[0]];
					const Vertex &v1 = vertexBuffer[tri.idx[1]];
					const Vertex &v2 = vertexBuffer[tri.idx[2]];
					m_triAccel[idx].load(v0.p, v1.p, v2.p);
					m_triAccel[idx].shapeIndex = i;
					m_triAccel[idx].index = j;
					++idx;
				}
			} else {
				m_triAccel[idx].shapeIndex = i;
				m_triAccel[idx].k = KNoTriangleFlag;
				++idx;
			}
		}
		Log(EDebug, "Done (%i ms)", timer->getMilliseconds());
		Log(EDebug, "");
		Assert(idx == primCount);
	}
	
protected:
	/**
	 * \brief Return the shape index corresponding to a primitive index
	 * seen by the generic kd-tree implementation. When this is a triangle
	 * mesh, the \a idx parameter is updated to the triangle index within
	 * the mesh.
	 */
	index_type findShape(index_type &idx) const {
		std::vector<index_type>::const_iterator it = std::lower_bound(
				m_shapeMap.begin(), m_shapeMap.end(), idx+1) - 1;
		idx -= *it;
		return (index_type) (it - m_shapeMap.begin());
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

	FINLINE AABB getAABB(index_type idx) const {
		index_type shapeIdx = findShape(idx);
		const Shape *shape = m_shapes[shapeIdx];
		if (m_triangleFlag[shapeIdx]) {
			const TriMesh *mesh = static_cast<const TriMesh *>(shape);
			return mesh->getTriangles()[idx].getAABB(mesh->getVertexBuffer());
		} else {
			return shape->getAABB();
		}
	}

	FINLINE AABB getClippedAABB(index_type idx, const AABB &aabb) const {
		index_type shapeIdx = findShape(idx);
		const Shape *shape = m_shapes[shapeIdx];
		if (m_triangleFlag[shapeIdx]) {
			const TriMesh *mesh = static_cast<const TriMesh *>(shape);
			return mesh->getTriangles()[idx].getClippedAABB(mesh->getVertexBuffer(), aabb);
		} else {
			return shape->getAABB();
		}
	}

	FINLINE size_type getPrimitiveCount() const {
		return m_shapeMap[m_shapeMap.size()-1];
	}
 
	FINLINE EIntersectionResult intersect(const Ray &ray, index_type idx, Float mint, 
		Float maxt, Float &t, void *temp) const {
		Float tempU, tempV, tempT;
		if (EXPECT_TAKEN(m_triAccel[idx].k != KNoTriangleFlag)) {
			const TriAccel &ta = m_triAccel[idx];
			if (ta.rayIntersect(ray, mint, maxt, tempU, tempV, tempT)) {
				index_type *indexPtr = reinterpret_cast<index_type *>(temp);
				Float *floatPtr = reinterpret_cast<Float *>((uint8_t *) indexPtr + 8);
				t = tempT;
				*indexPtr++ = ta.shapeIndex;
				*indexPtr++ = ta.index;
				*floatPtr++ = tempU;
				*floatPtr++ = tempV;
				return EYes;
			}
		} else {
			int shape = m_triAccel[idx].shapeIndex;
		}
		return ENo;
	}
private:
	std::vector<const Shape *> m_shapes;
	std::vector<bool> m_triangleFlag;
	std::vector<index_type> m_shapeMap;
	TriAccel *m_triAccel;
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
		TriKDTree tree;
		tree.addShape(mesh);
		tree.build();
		BSphere bsphere(mesh->getBSphere());

		//Float intersectionCost, traversalCost;
		//tree.findCosts(intersectionCost, traversalCost);

		ref<KDTree> oldTree = new KDTree();
		oldTree->addShape(mesh);
		oldTree->build(); 

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

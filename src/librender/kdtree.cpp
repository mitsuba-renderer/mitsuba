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

#include <mitsuba/render/kdtree.h>

MTS_NAMESPACE_BEGIN

KDTree::KDTree() : m_triAccel(NULL) {
	m_shapeMap.push_back(0);
}

KDTree::~KDTree() {
	if (m_triAccel)
		freeAligned(m_triAccel);
}

void KDTree::addShape(const Shape *shape) {
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

void KDTree::build() {
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
	Log(EDebug, "Finished -- took %i ms.", timer->getMilliseconds());
	Log(EDebug, "");
	Assert(idx == primCount);
}

bool KDTree::rayIntersect(const Ray &ray, Intersection &its) const {
	uint8_t temp[MTS_KD_INTERSECTION_TEMP];
	its.t = std::numeric_limits<Float>::infinity(); 
	Float mint, maxt;

	if (m_aabb.rayIntersect(ray, mint, maxt)) {
		if (ray.mint > mint) mint = ray.mint;
		if (ray.maxt < maxt) maxt = ray.maxt;

		if (EXPECT_TAKEN(maxt > mint)) {
			if (rayIntersectHavran<false>(ray, mint, maxt, its.t, temp)) {
				/* After having found a unique intersection, fill a proper record
				   using the temporary information collected in \ref intersect() */
				const IntersectionCache *cache = reinterpret_cast<const IntersectionCache *>(temp);
				const TriMesh *shape = static_cast<const TriMesh *>(m_shapes[cache->shapeIndex]);

				const Triangle &tri = shape->getTriangles()[cache->index];
				const Vertex *vertexBuffer = shape->getVertexBuffer();
				const Vertex &v0 = vertexBuffer[tri.idx[0]];
				const Vertex &v1 = vertexBuffer[tri.idx[1]];
				const Vertex &v2 = vertexBuffer[tri.idx[2]];
				const Vector b(1 - cache->u - cache->v, cache->u, cache->v);

				its.p = ray(its.t);
				its.geoFrame = Frame(normalize(cross(v1.p-v0.p, v2.p-v0.p)));
				its.shFrame.n = normalize(v0.n * b.x + v1.n * b.y + v2.n * b.z);
				const Vector dpdu = v0.dpdu * b.x + v1.dpdu * b.y + v2.dpdu * b.z;
				its.shFrame.s = normalize(dpdu - its.shFrame.n 
					* dot(its.shFrame.n, dpdu));
				its.shFrame.t = cross(its.shFrame.n, its.shFrame.s);
				its.dpdu = dpdu;
				its.dpdv = v0.dpdv * b.x + v1.dpdv * b.y + v2.dpdv * b.z;
				its.shFrame = its.geoFrame;
				its.uv = v0.uv * b.x + v1.uv * b.y + v2.uv * b.z;
				its.wi = its.toLocal(-ray.d);
				its.shape = shape;
				its.hasUVPartials = false;
				return true;
			}
		}
	}
	return false;
}

bool KDTree::rayIntersect(const Ray &ray) const {
	uint8_t temp[MTS_KD_INTERSECTION_TEMP];
	Float mint, maxt, t;

	if (m_aabb.rayIntersect(ray, mint, maxt)) {
		if (ray.mint > mint) mint = ray.mint;
		if (ray.maxt < maxt) maxt = ray.maxt;

		if (EXPECT_TAKEN(maxt > mint)) 
			if (rayIntersectHavran<true>(ray, mint, maxt, t, temp)) 
				return true;
	}
	return false;
}

MTS_IMPLEMENT_CLASS(KDTree, false, GenericKDTree)
MTS_NAMESPACE_END

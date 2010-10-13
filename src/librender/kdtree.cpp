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
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

KDTree::KDTree() {
#if !defined(MTS_KD_CONSERVE_MEMORY)
	m_triAccel = NULL;
#endif
	m_shapeMap.push_back(0);
}

KDTree::~KDTree() {
#if !defined(MTS_KD_CONSERVE_MEMORY)
	if (m_triAccel)
		freeAligned(m_triAccel);
#endif
}

static StatsCounter raysTraced("General", "Normal rays traced");
static StatsCounter shadowRaysTraced("General", "Shadow rays traced");

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

#if !defined(MTS_KD_CONSERVE_MEMORY)
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
	KDAssert(idx == primCount);
#endif
}

bool KDTree::rayIntersect(const Ray &ray, Intersection &its) const {
	uint8_t temp[MTS_KD_INTERSECTION_TEMP];
	its.t = std::numeric_limits<Float>::infinity(); 
	Float mint, maxt;

	++raysTraced;
	if (m_aabb.rayIntersect(ray, mint, maxt)) {
		/* Use an adaptive ray epsilon */
		Float rayMinT = ray.mint;
		if (rayMinT == Epsilon)
			rayMinT *= std::max(std::max(std::abs(ray.o.x), 
				std::abs(ray.o.y)), std::abs(ray.o.z));

		if (rayMinT > mint) mint = rayMinT;
		if (ray.maxt < maxt) maxt = ray.maxt;

		if (EXPECT_TAKEN(maxt > mint)) {
			if (rayIntersectHavran<false>(ray, mint, maxt, its.t, temp)) {
				/* After having found a unique intersection, fill a proper record
				   using the temporary information collected in \ref intersect() */
#if 1
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
#endif
				return true;
			}
		}
	}
	return false;
}

bool KDTree::rayIntersect(const Ray &ray) const {
	uint8_t temp[MTS_KD_INTERSECTION_TEMP];
	Float mint, maxt, t = std::numeric_limits<Float>::infinity();

	++shadowRaysTraced;
	if (m_aabb.rayIntersect(ray, mint, maxt)) {
		/* Use an adaptive ray epsilon */
		Float rayMinT = ray.mint;
		if (rayMinT == Epsilon)
			rayMinT *= std::max(std::max(std::abs(ray.o.x), 
				std::abs(ray.o.y)), std::abs(ray.o.z));

		if (rayMinT > mint) mint = rayMinT;
		if (ray.maxt < maxt) maxt = ray.maxt;

		if (EXPECT_TAKEN(maxt > mint)) 
			if (rayIntersectHavran<true>(ray, mint, maxt, t, temp)) 
				return true;
	}
	return false;
}


#if defined(MTS_HAS_COHERENT_RT)
static StatsCounter coherentPackets("General", "Coherent ray packets");
static StatsCounter incoherentPackets("General", "Incoherent ray packets");

void KDTree::rayIntersectPacket(const RayPacket4 &packet, 
		const RayInterval4 &rayInterval, Intersection4 &its) const {
	CoherentKDStackEntry MM_ALIGN16 stack[MTS_KD_MAXDEPTH];
	RayInterval4 MM_ALIGN16 interval;

	const KDNode * __restrict currNode = m_nodes;
	int stackIndex = 0;

	++coherentPackets;

	/* First, intersect with the kd-tree AABB to determine
	   the intersection search intervals */
	if (!m_aabb.rayIntersectPacket(packet, interval))
		return;

	interval.mint.ps = _mm_max_ps(interval.mint.ps, rayInterval.mint.ps);
	interval.maxt.ps = _mm_min_ps(interval.maxt.ps, rayInterval.maxt.ps);

	__m128 itsFound = _mm_cmpgt_ps(interval.mint.ps, interval.maxt.ps),
		   masked = itsFound;
	if (_mm_movemask_ps(itsFound) == 0xF)
		return;

	while (currNode != NULL) {
		while (EXPECT_TAKEN(!currNode->isLeaf())) {
			const uint8_t axis = currNode->getAxis();

			/* Calculate the plane intersection */
			const __m128
				splitVal = _mm_set1_ps(currNode->getSplit()),
				t = _mm_mul_ps(_mm_sub_ps(splitVal, packet.o[axis].ps),
					packet.dRcp[axis].ps);

			const __m128
				startsAfterSplit = _mm_or_ps(masked, 
					_mm_cmplt_ps(t, interval.mint.ps)),
				endsBeforeSplit = _mm_or_ps(masked,
					_mm_cmpgt_ps(t, interval.maxt.ps));

			currNode = currNode->getLeft() + packet.signs[axis][0];

			/* The interval completely completely lies on one side
			   of the split plane */
			if (EXPECT_TAKEN(_mm_movemask_ps(startsAfterSplit) == 15)) {
				currNode = currNode->getSibling();
				continue;
			}

			if (EXPECT_TAKEN(_mm_movemask_ps(endsBeforeSplit) == 15)) 
				continue;

			stack[stackIndex].node = currNode->getSibling();
			stack[stackIndex].interval.maxt =    interval.maxt;
			stack[stackIndex].interval.mint.ps = _mm_max_ps(t, interval.mint.ps);
			interval.maxt.ps =                   _mm_min_ps(t, interval.maxt.ps);
			masked = _mm_or_ps(masked, 	
					_mm_cmpgt_ps(interval.mint.ps, interval.maxt.ps));
			stackIndex++;
		}

		/* Arrived at a leaf node - intersect against primitives */
		const index_type primStart = currNode->getPrimStart();
		const index_type primEnd = currNode->getPrimEnd();

		if (EXPECT_NOT_TAKEN(primStart != primEnd)) {
#ifdef MTS_USE_TRIACCEL4
			const int count = m_packedTriangles[primStart].indirectionCount;
			primStart = m_packedTriangles[primStart].indirectionIndex;
			primEnd = primStart + count;
#endif
			__m128 
				searchStart = _mm_max_ps(rayInterval.mint.ps, 
					_mm_mul_ps(interval.mint.ps, SSEConstants::om_eps.ps)),
				searchEnd   = _mm_min_ps(rayInterval.maxt.ps, 
					_mm_mul_ps(interval.maxt.ps, SSEConstants::op_eps.ps));

			for (index_type entry=primStart; entry != primEnd; entry++) {
				const TriAccel &kdTri = m_triAccel[m_indices[entry]];
				if (EXPECT_TAKEN(kdTri.k != KNoTriangleFlag)) {
					itsFound = _mm_or_ps(itsFound, 
						kdTri.rayIntersectPacket(packet, searchStart, searchEnd, masked, its));
				} else {
					/* Not a triangle - invoke the shape's intersection routine */
					__m128 hasIts = m_shapes[kdTri.shapeIndex]->rayIntersectPacket(packet, 
							searchStart, searchEnd, masked, its);
					itsFound = _mm_or_ps(itsFound, hasIts);
					its.primIndex.pi  = mux_epi32(pstoepi32(hasIts), 
						load1_epi32(kdTri.index), its.primIndex.pi);
					its.shapeIndex.pi = mux_epi32(pstoepi32(hasIts),
						load1_epi32(kdTri.shapeIndex), its.shapeIndex.pi);
				}
				searchEnd = _mm_min_ps(searchEnd, its.t.ps);
			}
		}

		/* Abort if the tree has been traversed or if
		   intersections have been found for all four rays */
		if (_mm_movemask_ps(itsFound) == 0xF || --stackIndex < 0)
			break;

		/* Pop from the stack */
		currNode = stack[stackIndex].node;
		interval = stack[stackIndex].interval;
		masked = _mm_or_ps(itsFound, 
			_mm_cmpgt_ps(interval.mint.ps, interval.maxt.ps));
	}
}
	
void KDTree::rayIntersectPacket(const Ray *rays, Intersection *its) const {
	RayPacket4 MM_ALIGN16 packet;
	RayInterval4 MM_ALIGN16 interval(rays);
	Intersection4 MM_ALIGN16 its4;

	if (packet.load(rays)) {
		rayIntersectPacket(packet, interval, its4);

		for (int i=0; i<4; i++) {
			Intersection &it = its[i];

			it.t = its4.t.f[i];

			if (it.t != std::numeric_limits<float>::infinity()) {
				const uint32_t shapeIndex = its4.shapeIndex.i[i];
				const uint32_t primIndex = its4.primIndex.i[i];
				const Shape *shape = m_shapes[shapeIndex];

				if (EXPECT_TAKEN(primIndex != KNoTriangleFlag)) {
					const TriMesh *triMesh = static_cast<const TriMesh *>(m_shapes[shapeIndex]);
					const Triangle &t = triMesh->getTriangles()[primIndex];
					const Vertex &v0 = triMesh->getVertexBuffer()[t.idx[0]];
					const Vertex &v1 = triMesh->getVertexBuffer()[t.idx[1]];
					const Vertex &v2 = triMesh->getVertexBuffer()[t.idx[2]];
					const Vector b(1 - its4.u.f[i] - its4.v.f[i],
						its4.u.f[i], its4.v.f[i]);
					const Vector &rayD = rays[i].d;
					it.p = rays[i].o + it.t * rayD;

					Normal faceNormal(normalize(cross(v1.p-v0.p, v2.p-v0.p)));
					it.geoFrame.n = faceNormal;
					coordinateSystem(it.geoFrame.n, it.geoFrame.s, it.geoFrame.t);

					it.uv = v0.uv * b.x + v1.uv * b.y + v2.uv * b.z;
					it.dpdu = v0.dpdu * b.x + v1.dpdu * b.y + v2.dpdu * b.z;
					it.dpdv = v0.dpdv * b.x + v1.dpdv * b.y + v2.dpdv * b.z;
					it.shFrame.n = normalize(v0.n * b.x + v1.n * b.y + v2.n * b.z);

					it.shFrame.s = normalize(it.dpdu - it.shFrame.n 
						* dot(it.shFrame.n, it.dpdu));
					it.geoFrame.t = cross(it.shFrame.n, it.shFrame.s);
					it.wi = it.toLocal(-rayD);
					it.hasUVPartials = false;
					it.shape = shape;
				} else {
					/* Non-triangle shape: intersect again to fill in details */
					shape->rayIntersect(rays[i], it);
				}
			}
		}
	} else {
		rayIntersectPacketIncoherent(rays, its);
	}
}

void KDTree::rayIntersectPacketIncoherent(const RayPacket4 &packet, 
		const RayInterval4 &rayInterval, Intersection4 &its4) const {
	uint8_t temp[MTS_KD_INTERSECTION_TEMP];

	++incoherentPackets;
	for (int i=0; i<4; i++) {
		Ray ray;
		for (int axis=0; axis<3; axis++) {
			ray.o[axis] = packet.o[axis].f[i];
			ray.d[axis] = packet.d[axis].f[i];
			ray.dRcp[axis] = packet.dRcp[axis].f[i];
		}
		ray.mint = rayInterval.mint.f[i];
		ray.maxt = rayInterval.maxt.f[i];
		if (ray.mint < ray.maxt && rayIntersectHavran<false>(ray, ray.mint, ray.maxt, its4.t.f[i], temp)) {
			const IntersectionCache *cache = reinterpret_cast<const IntersectionCache *>(temp);
			its4.u.f[i] = cache->u;
			its4.v.f[i] = cache->u;
			its4.shapeIndex.i[i] = cache->shapeIndex;
			its4.primIndex.i[i] = cache->index;
		}
	}
}

#else // !MTS_HAS_COHERENT_RT
void KDTree::rayIntersectPacket(const Ray *rays, Intersection *its) const {
	rayIntersectPacketIncoherent(rays, its);
}
#endif

MTS_IMPLEMENT_CLASS(KDTree, false, GenericKDTree)
MTS_NAMESPACE_END

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

#include <mitsuba/render/skdtree.h>
#include <mitsuba/core/statistics.h>

#if defined(MTS_SSE)
#include <mitsuba/core/sse.h>
#include <mitsuba/core/aabb_sse.h>
#include <mitsuba/render/triaccel_sse.h>
#endif

MTS_NAMESPACE_BEGIN

ShapeKDTree::ShapeKDTree() {
#if !defined(MTS_KD_CONSERVE_MEMORY)
	m_triAccel = NULL;
#endif
	m_shapeMap.push_back(0);
}

ShapeKDTree::~ShapeKDTree() {
#if !defined(MTS_KD_CONSERVE_MEMORY)
	if (m_triAccel)
		freeAligned(m_triAccel);
#endif
	for (size_t i=0; i<m_shapes.size(); ++i)
		m_shapes[i]->decRef();
}

static StatsCounter raysTraced("General", "Normal rays traced");
static StatsCounter shadowRaysTraced("General", "Shadow rays traced");

void ShapeKDTree::addShape(const Shape *shape) {
	Assert(!isBuilt());
	if (shape->isCompound())
		Log(EError, "Cannot add compound shapes to a kd-tree - expand them first!");
	if (shape->getClass()->derivesFrom(MTS_CLASS(TriMesh))) {
		// Triangle meshes are expanded into individual primitives,
		// which are visible to the tree construction code. Generic
		// primitives are only handled by their AABBs
		m_shapeMap.push_back((SizeType)
			static_cast<const TriMesh *>(shape)->getTriangleCount());
		m_triangleFlag.push_back(true);
	} else {
		m_shapeMap.push_back(1);
		m_triangleFlag.push_back(false);
	}
	shape->incRef();
	m_shapes.push_back(shape);
}

void ShapeKDTree::build() {
	for (size_t i=1; i<m_shapeMap.size(); ++i)
		m_shapeMap[i] += m_shapeMap[i-1];

	SAHKDTree3D<ShapeKDTree>::buildInternal();

#if !defined(MTS_KD_CONSERVE_MEMORY)
	ref<Timer> timer = new Timer();
	SizeType primCount = getPrimitiveCount();
	Log(EDebug, "Precomputing triangle intersection information (%s)",
			memString(sizeof(TriAccel)*primCount).c_str());
	m_triAccel = static_cast<TriAccel *>(allocAligned(primCount * sizeof(TriAccel)));

	IndexType idx = 0;
	for (IndexType i=0; i<m_shapes.size(); ++i) {
		const Shape *shape = m_shapes[i];
		if (m_triangleFlag[i]) {
			const TriMesh *mesh = static_cast<const TriMesh *>(shape);
			const Triangle *triangles = mesh->getTriangles();
			const Point *positions = mesh->getVertexPositions();
			for (IndexType j=0; j<mesh->getTriangleCount(); ++j) {
				const Triangle &tri = triangles[j];
				const Point &v0 = positions[tri.idx[0]];
				const Point &v1 = positions[tri.idx[1]];
				const Point &v2 = positions[tri.idx[2]];
				m_triAccel[idx].load(v0, v1, v2);
				m_triAccel[idx].shapeIndex = i;
				m_triAccel[idx].primIndex = j;
				++idx;
			}
		} else {
			/* Create a 'fake' triangle, which redirects to a Shape */
			memset(&m_triAccel[idx], 0, sizeof(TriAccel));
			m_triAccel[idx].shapeIndex = i;
			m_triAccel[idx].k = KNoTriangleFlag;
			++idx;
		}
	}
	Log(EDebug, "Finished -- took %i ms.", timer->getMilliseconds());
	Log(m_logLevel, "");
	KDAssert(idx == primCount);
#endif
}

bool ShapeKDTree::rayIntersect(const Ray &ray, Intersection &its) const {
	uint8_t temp[MTS_KD_INTERSECTION_TEMP];
	its.t = std::numeric_limits<Float>::infinity();
	Float mint, maxt;

	#if defined(MTS_FP_DEBUG_STRICT)
		Assert(
			std::isfinite(ray.o.x) && std::isfinite(ray.o.y) && std::isfinite(ray.o.z) &&
			std::isfinite(ray.d.x) && std::isfinite(ray.d.y) && std::isfinite(ray.d.z));
	#endif

	++raysTraced;
	if (m_aabb.rayIntersect(ray, mint, maxt)) {
		/* Use an adaptive ray epsilon */
		Float rayMinT = ray.mint;
		if (rayMinT == Epsilon)
			rayMinT *= std::max(std::max(std::max(std::abs(ray.o.x),
				std::abs(ray.o.y)), std::abs(ray.o.z)), Epsilon);

		if (rayMinT > mint) mint = rayMinT;
		if (ray.maxt < maxt) maxt = ray.maxt;

		if (EXPECT_TAKEN(maxt > mint)) {
			if (rayIntersectHavran<false>(ray, mint, maxt, its.t, temp)) {
				fillIntersectionRecord<true>(ray, temp, its);
				return true;
			}
		}
	}
	return false;
}

bool ShapeKDTree::rayIntersect(const Ray &ray, Float &t, ConstShapePtr &shape,
		Normal &n, Point2 &uv) const {
	uint8_t temp[MTS_KD_INTERSECTION_TEMP];
	Float mint, maxt;

	t = std::numeric_limits<Float>::infinity();

	++shadowRaysTraced;
	if (m_aabb.rayIntersect(ray, mint, maxt)) {
		/* Use an adaptive ray epsilon */
		Float rayMinT = ray.mint;
		if (rayMinT == Epsilon)
			rayMinT *= std::max(std::max(std::abs(ray.o.x),
				std::abs(ray.o.y)), std::abs(ray.o.z));

		if (rayMinT > mint) mint = rayMinT;
		if (ray.maxt < maxt) maxt = ray.maxt;

		if (EXPECT_TAKEN(maxt > mint)) {
			if (rayIntersectHavran<false>(ray, mint, maxt, t, temp)) {
				const IntersectionCache *cache = reinterpret_cast<const IntersectionCache *>(temp);
				shape = m_shapes[cache->shapeIndex];

				if (m_triangleFlag[cache->shapeIndex]) {
					const TriMesh *trimesh = static_cast<const TriMesh *>(shape);
					const Triangle &tri = trimesh->getTriangles()[cache->primIndex];
					const Point *vertexPositions = trimesh->getVertexPositions();
					const Point2 *vertexTexcoords = trimesh->getVertexTexcoords();
					const uint32_t idx0 = tri.idx[0], idx1 = tri.idx[1], idx2 = tri.idx[2];
					const Point &p0 = vertexPositions[idx0];
					const Point &p1 = vertexPositions[idx1];
					const Point &p2 = vertexPositions[idx2];
					n = normalize(cross(p1-p0, p2-p0));

					if (EXPECT_TAKEN(vertexTexcoords)) {
						const Vector b(1 - cache->u - cache->v, cache->u, cache->v);
						const Point2 &t0 = vertexTexcoords[idx0];
						const Point2 &t1 = vertexTexcoords[idx1];
						const Point2 &t2 = vertexTexcoords[idx2];
						uv = t0 * b.x + t1 * b.y + t2 * b.z;
					} else {
						uv = Point2(0.0f);
					}
				} else {
					/// Uh oh... -- much unnecessary work is done here
					Intersection its;
					its.t = t;
					shape->fillIntersectionRecord(ray,
						reinterpret_cast<const uint8_t*>(temp) + 2*sizeof(IndexType), its);
					n = its.geoFrame.n;
					uv = its.uv;
					if (its.shape)
						shape = its.shape;
				}

				return true;
			}
		}
	}
	return false;
}


bool ShapeKDTree::rayIntersect(const Ray &ray) const {
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
			if (rayIntersectHavran<true>(ray, mint, maxt, t, NULL))
				return true;
	}
	return false;
}

#if defined(MTS_HAS_COHERENT_RT)

/// Ray traversal stack entry for uncoherent ray tracing
struct CoherentKDStackEntry {
	/* Current ray interval */
	RayInterval4 MM_ALIGN16 interval;
	/* Pointer to the far child */
	const ShapeKDTree::KDNode * __restrict node;
};

static StatsCounter coherentPackets("General", "Coherent ray packets");
static StatsCounter incoherentPackets("General", "Incoherent ray packets");

void ShapeKDTree::rayIntersectPacket(const RayPacket4 &packet,
		const RayInterval4 &rayInterval, Intersection4 &its, void *temp) const {
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

	SSEVector itsFound( _mm_cmpgt_ps(interval.mint.ps, interval.maxt.ps));
	SSEVector masked(itsFound);
	if (_mm_movemask_ps(itsFound.ps) == 0xF)
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
				startsAfterSplit = _mm_or_ps(masked.ps,
					_mm_cmplt_ps(t, interval.mint.ps)),
				endsBeforeSplit = _mm_or_ps(masked.ps,
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
			masked.ps = _mm_or_ps(masked.ps,
					_mm_cmpgt_ps(interval.mint.ps, interval.maxt.ps));
			stackIndex++;
		}

		/* Arrived at a leaf node - intersect against primitives */
		const IndexType primStart = currNode->getPrimStart();
		const IndexType primEnd = currNode->getPrimEnd();

		if (EXPECT_NOT_TAKEN(primStart != primEnd)) {
			SSEVector
				searchStart(_mm_max_ps(rayInterval.mint.ps,
					_mm_mul_ps(interval.mint.ps, SSEConstants::om_eps.ps))),
				searchEnd(_mm_min_ps(rayInterval.maxt.ps,
					_mm_mul_ps(interval.maxt.ps, SSEConstants::op_eps.ps)));

			for (IndexType entry=primStart; entry != primEnd; entry++) {
				const TriAccel &kdTri = m_triAccel[m_indices[entry]];
				if (EXPECT_TAKEN(kdTri.k != KNoTriangleFlag)) {
					itsFound.ps = _mm_or_ps(itsFound.ps,
						mitsuba::rayIntersectPacket(kdTri, packet, searchStart.ps, searchEnd.ps, masked.ps, its));
				} else {
					const Shape *shape = m_shapes[kdTri.shapeIndex];

					for (int i=0; i<4; ++i) {
						if (masked.i[i])
							continue;
						Ray ray;
						for (int axis=0; axis<3; axis++) {
							ray.o[axis] = packet.o[axis].f[i];
							ray.d[axis] = packet.d[axis].f[i];
							ray.dRcp[axis] = packet.dRcp[axis].f[i];
						}
						Float t;

						if (shape->rayIntersect(ray, searchStart.f[i], searchEnd.f[i], t,
								reinterpret_cast<uint8_t *>(temp)
								+ i * MTS_KD_INTERSECTION_TEMP + 2*sizeof(IndexType))) {
							its.t.f[i] = t;
							its.shapeIndex.i[i] = kdTri.shapeIndex;
							its.primIndex.i[i] = KNoTriangleFlag;
							itsFound.i[i] = 0xFFFFFFFF;
						}
					}
				}
				searchEnd.ps = _mm_min_ps(searchEnd.ps, its.t.ps);
			}
		}

		/* Abort if the tree has been traversed or if
		   intersections have been found for all four rays */
		if (_mm_movemask_ps(itsFound.ps) == 0xF || --stackIndex < 0)
			break;

		/* Pop from the stack */
		currNode = stack[stackIndex].node;
		interval = stack[stackIndex].interval;
		masked.ps = _mm_or_ps(itsFound.ps,
			_mm_cmpgt_ps(interval.mint.ps, interval.maxt.ps));
	}
}

void ShapeKDTree::rayIntersectPacketIncoherent(const RayPacket4 &packet,
		const RayInterval4 &rayInterval, Intersection4 &its4, void *temp) const {

	++incoherentPackets;
	for (int i=0; i<4; i++) {
		Ray ray;
		Float t;
		for (int axis=0; axis<3; axis++) {
			ray.o[axis] = packet.o[axis].f[i];
			ray.d[axis] = packet.d[axis].f[i];
			ray.dRcp[axis] = packet.dRcp[axis].f[i];
		}
		ray.mint = rayInterval.mint.f[i];
		ray.maxt = rayInterval.maxt.f[i];
		uint8_t *rayTemp = reinterpret_cast<uint8_t *>(temp) + i * MTS_KD_INTERSECTION_TEMP;
		if (ray.mint < ray.maxt && rayIntersectHavran<false>(ray, ray.mint, ray.maxt, t, rayTemp)) {
			const IntersectionCache *cache = reinterpret_cast<const IntersectionCache *>(rayTemp);
			its4.t.f[i] = t;
			its4.shapeIndex.i[i] = cache->shapeIndex;
			its4.primIndex.i[i] = cache->primIndex;
			its4.u.f[i] = cache->u;
			its4.v.f[i] = cache->v;
		}
	}
}

#endif

MTS_IMPLEMENT_CLASS(ShapeKDTree, false, KDTreeBase)
MTS_NAMESPACE_END

#include <mitsuba/render/kdtree.h>

MTS_NAMESPACE_BEGIN

#if defined(MTS_HAS_COHERENT_RT)
static StatsCounter coherentPackets("General", "Coherent ray packets");
static StatsCounter incoherentPackets("General", "Incoherent ray packets");

void KDTree::rayIntersectPacket(const RayPacket4 &packet, 
		const RayInterval4 &rayInterval, Intersection4 &its) const {
	CoherentKDStackEntry MM_ALIGN16 stack[MTS_KD_MAXDEPTH];
	RayInterval4 MM_ALIGN16 interval;

	const KDNode * __restrict currNode = &m_nodes[1];
	int stackIndex = 0;

	++coherentPackets;

	/* First, intersect with the kd-tree AABB to determine
	   the intersection search intervals */
	if (!m_rootBounds.rayIntersectPacket(packet, interval))
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
				currNode = currNode->getOtherChild();
				continue;
			}

			if (EXPECT_TAKEN(_mm_movemask_ps(endsBeforeSplit) == 15)) {
				continue;
			}

			stack[stackIndex].node = currNode->getOtherChild();
			stack[stackIndex].interval.maxt =    interval.maxt;
			stack[stackIndex].interval.mint.ps = _mm_max_ps(t, interval.mint.ps);
			interval.maxt.ps =                   _mm_min_ps(t, interval.maxt.ps);
			masked = _mm_or_ps(masked, 	
					_mm_cmpgt_ps(interval.mint.ps, interval.maxt.ps));
			stackIndex++;
		}

		/* Arrived at a leaf node - intersect against primitives */
		unsigned int primStart = currNode->getPrimStart();
		unsigned int primEnd = currNode->getPrimEnd();

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

			for (unsigned int entry=primStart; entry != primEnd; entry++) {
				const KDTriangle &kdTri = m_triangles[m_indices[entry]];
				if (EXPECT_TAKEN(kdTri.index != KNoTriangleFlag)) {
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

					Normal faceNormal(normalize(cross(v1.v-v0.v, v2.v-v0.v)));
					it.geoFrame.n = faceNormal;
					coordinateSystem(it.geoFrame.n, it.geoFrame.s, it.geoFrame.t);

					/* Intersection refinement step */
					Vector rel = it.p - v0.v;
					Float correction = -dot(rel, faceNormal)/dot(rayD, faceNormal);
					it.t += correction;
					it.p += rayD * correction;

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
	Intersection its;
	unsigned int shapeIndex, primIndex;

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
		if (ray.mint < ray.maxt && rayIntersect(ray, its, ray.mint, ray.maxt, false, shapeIndex, primIndex)) {
			its4.t.f[i] = its.t;
			its4.u.f[i] = its.uv.x;
			its4.v.f[i] = its.uv.y;
			its4.shapeIndex.i[i] = shapeIndex;
			its4.primIndex.i[i] = primIndex;
		}
	}
}

#else // !MTS_HAS_COHERENT_RT

void KDTree::rayIntersectPacket(const Ray *rays, Intersection *its) const {
	rayIntersectPacketIncoherent(rays, its);
}
#endif

MTS_NAMESPACE_END

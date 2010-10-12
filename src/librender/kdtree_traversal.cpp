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

#include <mitsuba/core/triangle.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/kdtree.h>

MTS_NAMESPACE_BEGIN

static StatsCounter raysTraced("General", "Normal rays traced");
static StatsCounter shadowRaysTraced("General", "Shadow rays traced");

bool KDTree::rayIntersect(const Ray &ray, Intersection &its, Float mint, Float maxt, 
		bool shadowRay, unsigned int &shapeIndex, unsigned int &primIndex) const {
	KDStackEntry stack[MTS_KD_MAXDEPTH];
	const KDNode * __restrict farChild,
		         * __restrict currNode = &m_nodes[1];

	static const int prevAxisTable[] = { 2, 0, 1 };
	static const int nextAxisTable[] = { 1, 2, 0 };
	Float t;

#if !defined(MTS_USE_TRIACCEL4)
	Float u, v;
#else
	/* Force these to be aligned on 16-byte boundaries */
	Float _tt[8], _u[8], _v[8];
	Float *tt = STACK_ALIGN16(_tt),
		  *u = STACK_ALIGN16(_u),
		  *v = STACK_ALIGN16(_v);

	/* Move the ray into SSE registers */
	const __m128 
		ray_o = _mm_loadu_ps(&ray.o.x),
		ray_d = _mm_loadu_ps(&ray.d.x);
#endif

	/* Set up the entry point */
	int enPt = 0;
	stack[enPt].t = mint;
	if (mint >= 0.0f)
		stack[enPt].pb = ray(mint);
	else
		stack[enPt].pb = ray.o;

	/* Set up the exit point */
	int exPt = 1;
	stack[exPt].t = maxt;
	stack[exPt].pb = ray(maxt);
	stack[exPt].node = NULL;

	while (currNode != NULL) {
		while (EXPECT_TAKEN(!currNode->isLeaf())) {
			const uint8_t axis = currNode->getAxis();
			const Float splitVal = currNode->getSplit();
			if (stack[enPt].pb[axis] <= splitVal) {
				if (stack[exPt].pb[axis] <= splitVal) {
					/* Cases N1, N2, N3, P5, Z2 and Z3 (see thesis) */
					currNode = currNode->getLeft();
					continue;
				}

				/* stack[exPt].pb[axis] > splitVal. This seems to be
				   a typo in the thesis. */
				if (stack[enPt].pb[axis] == splitVal) {
					/* Case Z1 */
					currNode = currNode->getRight();
					continue;
				}

				/* Case N4 */
				currNode = currNode->getLeft();
				farChild = currNode + 1; // getLeft()
			} else { /* stack[enPt].pb[axis] > splitVal */
				if (splitVal < stack[exPt].pb[axis]) {
					/* Cases P1, P2, P3 and N5 */
					currNode = currNode->getRight();
					continue;
				}
				/* Case P4 */
				farChild = currNode->getLeft();
				currNode = farChild + 1; // getRight()
			}

			/* Cases P4 or N4 - calculate the distance to the split plane */
			t = (splitVal - ray.o[axis]) * ray.dRcp[axis];

			/* Set up a new exit point */
			const int tmp = exPt++;
			if (exPt == enPt) /* Do not overwrite the entry point */
				++exPt;

			const int nextAxis = nextAxisTable[axis];
			const int prevAxis = prevAxisTable[axis];

			stack[exPt].prev = tmp;
			stack[exPt].t = t;
			stack[exPt].node = farChild;
			stack[exPt].pb[axis] = splitVal;
			stack[exPt].pb[nextAxis] = ray.o[nextAxis] + t*ray.d[nextAxis];
			stack[exPt].pb[prevAxis] = ray.o[prevAxis] + t*ray.d[prevAxis];
		}

#if defined(SINGLE_PRECISION)
		const Float eps = 1e-3;
#else
		const Float eps = 1e-5;
#endif
		const Float m_eps = 1-eps, p_eps = 1+eps;

		/* Floating-point arithmetic.. - use both absolute and relative 
		   epsilons when looking for intersections in the subinterval */
		const Float searchStart = std::max(mint, stack[enPt].t * m_eps - eps);
		Float searchEnd   = std::min(maxt, stack[exPt].t * p_eps + eps);

		/* Reached a leaf node */
		bool foundIntersection = false;
		for (unsigned int entry=currNode->getPrimStart(),
				last = currNode->getPrimEnd(); entry != last; entry++) {
#ifdef MTS_USE_TRIACCEL4
			if (m_packedTriangles[entry].rayIntersect(ray_o, ray_d,
				searchStart, searchEnd, tt, u, v, shapeIndex, primIndex)) {
				if (shadowRay)
					return true;
				foundIntersection = true;
				searchEnd = its.t = tt[0];
				its.uv.x = u[0];
				its.uv.y = v[0];
			}
			if (EXPECT_NOT_TAKEN(m_packedTriangles[entry].nonTriFlag)) {
				const TriAccel4 &ta = m_packedTriangles[entry];
				for (int i=0; i<4; ++i) {
					if (ta.index[i] != KNoTriangleFlag)
						continue;
					const Shape *shape = m_shapes[m_packedTriangles[entry].shapeIndex[i]];
					if (shape->rayIntersect(ray, searchStart, searchEnd, t)) {
						if (shadowRay)
							return true;
						foundIntersection = true;
						searchEnd = its.t = t;
						shapeIndex = ta.shapeIndex[i];
						primIndex = ta.index[i];
					}
				}
			}
#elif MTS_USE_TRIACCEL
			const KDTriangle &kdTri = m_triangles[m_indices[entry]];
			if (EXPECT_TAKEN(kdTri.index != KNoTriangleFlag)) {
				if (kdTri.rayIntersect(ray, searchStart, searchEnd, u, v, t)) {
					if (shadowRay)
						return true;
					foundIntersection = true;
					shapeIndex = kdTri.shapeIndex;
					primIndex = kdTri.index;
					searchEnd = its.t = t;
					its.uv.x = u;
					its.uv.y = v;
				}
			} else {
				if (m_shapes[kdTri.shapeIndex]->rayIntersect(ray, searchStart, searchEnd, t)) {
					if (shadowRay)
						return true;
					foundIntersection = true;
					shapeIndex = kdTri.shapeIndex;
					primIndex = kdTri.index;
					searchEnd = its.t = t;
				}
			}
#else
			/* Calculate triangle intersections using the fast algorithm by 
			   Moeller and Trumbore (without having done any pre-computation) */
			const KDTriangle &kdTri = m_triangles[m_indices[entry]];
			if (EXPECT_TAKEN(kdTri.index != KNoTriangleFlag)) {
				const TriMesh *mesh = (const TriMesh *) m_shapes[kdTri.shapeIndex];
				const Triangle &tri = mesh->getTriangles()[kdTri.index];
				if (tri.rayIntersect(mesh->getVertexBuffer(), ray, u, v, t)) {
					if (t >= searchStart && t < searchEnd) {
						if (shadowRay)
							return true;
						searchEnd = its.t = t;
						its.uv.x = u;
						its.uv.y = v;
						const KDTriangle &tri = m_triangles[m_indices[entry]];
						shapeIndex = tri.shapeIndex;
						primIndex = tri.index;
						foundIntersection = true;
					}
				}
			} else {
				if (m_shapes[kdTri.shapeIndex]->rayIntersect(ray, searchStart, searchEnd, t)) {
					if (shadowRay)
						return true;
					foundIntersection = true;
					shapeIndex = kdTri.shapeIndex;
					primIndex = kdTri.index;
					searchEnd = its.t = t;
				}
			}
#endif
		}
		if (foundIntersection) 
			return foundIntersection;

		/* Pop from the stack and advance to the next node on the interval */
		enPt = exPt;
		currNode = stack[exPt].node;
		exPt = stack[enPt].prev;
	}

	return false;
}

bool KDTree::rayIntersect(const Ray &ray, Intersection &its) const {
	Float mint, maxt;

	++raysTraced;
	if (m_rootBounds.rayIntersect(ray, mint, maxt)) {
		its.t = std::numeric_limits<Float>::infinity(); 
		
		/* Adaptive ray epsilon */
		Float rayMinT = ray.mint;
		if (rayMinT == Epsilon)
			rayMinT *= std::max(std::max(std::abs(ray.o.x), 
				std::abs(ray.o.y)), std::abs(ray.o.z));
		
		if (rayMinT > mint)
			mint = rayMinT;
		if (ray.maxt < maxt)
			maxt = ray.maxt;
		if (maxt <= mint)
			return false;

		unsigned int shapeIndex, primIndex;
		if (rayIntersect(ray, its, mint, maxt, false, shapeIndex, primIndex)) {
			const Shape *shape = m_shapes[shapeIndex];
			if (EXPECT_NOT_TAKEN(its.t < mint || its.t > maxt)) // double-check without epsilons
				return false;

			if (EXPECT_TAKEN(primIndex != KNoTriangleFlag)) {
				const TriMesh *triMesh = static_cast<const TriMesh *>(shape);
				const Triangle &t = triMesh->getTriangles()[primIndex];
				const Vertex &v0 = triMesh->getVertexBuffer()[t.idx[0]];
				const Vertex &v1 = triMesh->getVertexBuffer()[t.idx[1]];
				const Vertex &v2 = triMesh->getVertexBuffer()[t.idx[2]];
				its.p = ray(its.t);
				Normal faceNormal(normalize(cross(v1.p-v0.p, v2.p-v0.p)));

#if defined(SINGLE_PRECISION)
				/* Intersection refinement step */
				Float correction = -dot(its.p - v0.p, faceNormal)/dot(ray.d, faceNormal);
				its.t += correction;
				its.p += ray.d * correction;

				if (its.t < mint || its.t > maxt) {
					its.t = std::numeric_limits<Float>::infinity(); 
					return false;
				}
#endif

				const Vector b(1 - its.uv.x - its.uv.y,
					its.uv.x, its.uv.y);

				its.uv = v0.uv * b.x + v1.uv * b.y + v2.uv * b.z;
				its.dpdu = v0.dpdu * b.x + v1.dpdu * b.y + v2.dpdu * b.z;
				its.dpdv = v0.dpdv * b.x + v1.dpdv * b.y + v2.dpdv * b.z;

				its.geoFrame.n = faceNormal;
				its.shape = shape;

				coordinateSystem(its.geoFrame.n, its.geoFrame.s, its.geoFrame.t);
				its.shFrame.n = normalize(v0.n * b.x + v1.n * b.y + v2.n * b.z);
				its.shFrame.s = normalize(its.dpdu - its.shFrame.n 
					* dot(its.shFrame.n, its.dpdu));
				its.shFrame.t = cross(its.shFrame.n, its.shFrame.s);
				its.wi = its.toLocal(-ray.d);
				its.hasUVPartials = false;

#if defined(MTS_HAS_VERTEX_COLORS)
				its.color.fromLinearRGB(
					v0.color[0] * b.x + v1.color[0] * b.y + v2.color[0] * b.z,
					v0.color[1] * b.x + v1.color[1] * b.y + v2.color[1] * b.z,
					v0.color[2] * b.x + v1.color[2] * b.y + v2.color[2] * b.z
				);
#endif
			} else {
				/* Non-triangle shape: intersect again to fill in details */
				if (!shape->rayIntersect(ray, its))
					return false;
			}

			return true;
		}
	}
	its.t = std::numeric_limits<Float>::infinity(); 
	return false;
}

bool KDTree::rayIntersect(const Ray &ray) const {
	Float mint, maxt;

	++shadowRaysTraced;
	Intersection its;

	if (m_rootBounds.rayIntersect(ray, mint, maxt)) {
		its.t = std::numeric_limits<Float>::infinity(); 
		if (ray.mint > mint)
			mint = ray.mint;
		if (ray.maxt < maxt)
			maxt = ray.maxt;
		if (EXPECT_TAKEN(maxt > mint)) {
			unsigned int shapeIndex, primIndex;
			if (rayIntersect(ray, its, mint, maxt, true, shapeIndex, primIndex))
				return true;
		}
	}
	its.t = std::numeric_limits<Float>::infinity(); 
	return false;
}

MTS_NAMESPACE_END

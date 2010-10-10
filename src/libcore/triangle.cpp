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

MTS_NAMESPACE_BEGIN

bool Triangle::rayIntersect(const Vertex *buffer, const Ray &ray,
	Float &u, Float &v, Float &t) const {
	const Point &v0 = buffer[idx[0]].p;
	const Point &v1 = buffer[idx[1]].p;
	const Point &v2 = buffer[idx[2]].p;

	/* find vectors for two edges sharing v[0] */
	Vector edge1 = v1 - v0, edge2 = v2 - v0;

	/* begin calculating determinant - also used to calculate U parameter */
	Vector pvec = cross(ray.d, edge2);

	/* if determinant is near zero, ray lies in plane of triangle */
	Float det = dot(edge1, pvec);

	if (det > -Epsilon && det < Epsilon)
		return false;
	Float inv_det = 1.0f / det;

	/* calculate distance from v[0] to ray origin */
	Vector tvec = ray.o - v0;

	/* calculate U parameter and test bounds */
	u = dot(tvec, pvec) * inv_det;
	if (u < 0.0 || u > 1.0)
		return false;

	/* prepare to test V parameter */
	Vector qvec = cross(tvec, edge1);

	/* calculate V parameter and test bounds */
	v = dot(ray.d, qvec) * inv_det;
	if (v < 0.0 || u + v > 1.0)
		return false;

	/* calculate t, ray intersects triangle */
	t = dot(edge2, qvec) * inv_det;

	return true;
}

Point Triangle::sample(const Vertex *buffer, Normal &normal, 
	const Point2 &sample) const {
	const Point &v0 = buffer[idx[0]].p;
	const Point &v1 = buffer[idx[1]].p;
	const Point &v2 = buffer[idx[2]].p;
	const Normal &n0 = buffer[idx[0]].n;
	const Normal &n1 = buffer[idx[1]].n;
	const Normal &n2 = buffer[idx[2]].n;

	Point2 bary = squareToTriangle(sample);
	Vector sideA = v1 - v0, sideB = v2 - v0;
	Point p = v0 + (sideA * bary.x) + (sideB * bary.y);	
	normal = Normal(normalize(
		n0 * (1.0f - bary.x - bary.y) +
		n1 * bary.x + n2 * bary.y
	));

	return p;
}

Float Triangle::surfaceArea(const Vertex *buffer) const {
	const Point &v0 = buffer[idx[0]].p;
	const Point &v1 = buffer[idx[1]].p;
	const Point &v2 = buffer[idx[2]].p;
	Vector sideA = v1 - v0, sideB = v2 - v0;
	return 0.5f * cross(sideA, sideB).length();
}

#define MAX_VERTS 10

static int sutherlandHodgman(Point *input, int inCount, Point *output, int axis, 
		Float splitPos, bool isMinimum) {
	if (inCount < 3)
		return 0;

	Point cur         = input[0];
	Float sign        = isMinimum ? 1.0f : -1.0f;
	Float distance    = sign * (cur[axis] - splitPos);
	bool  curIsInside = (distance >= 0);
	int   outCount    = 0;

	for (int i=0; i<inCount; ++i) {
		Point next = input[(i+1)%inCount];
		distance = sign * (next[axis] - splitPos);
		bool nextIsInside = (distance >= 0);

		if (curIsInside && nextIsInside) {
			/* Both this and the next vertex are inside, add to the list */
			SAssertEx(outCount + 1 < MAX_VERTS, "Overflow in sutherlandHodgman()!");
			output[outCount++] = next;
		} else if (curIsInside && !nextIsInside) {
			/* Going outside -- add the intersection */
			Float t = (splitPos - cur[axis]) / (next[axis] - cur[axis]);
			SAssertEx(outCount + 1 < MAX_VERTS, "Overflow in sutherlandHodgman()!");
			Point p = cur + (next - cur) * t;
			p[axis] = splitPos; // Avoid roundoff errors
			output[outCount++] = p;
		} else if (!curIsInside && nextIsInside) {
			/* Coming back inside -- add the intersection + next vertex */
			Float t = (splitPos - cur[axis]) / (next[axis] - cur[axis]);
			SAssertEx(outCount + 2 < MAX_VERTS, "Overflow in sutherlandHodgman()!");
			Point p = cur + (next - cur) * t;
			p[axis] = splitPos; // Avoid roundoff errors
			output[outCount++] = p;
			output[outCount++] = next;
		} else {
			/* Entirely outside - do not add anything */
		}
		cur = next;
		curIsInside = nextIsInside;
	}
	return outCount;
}

AABB Triangle::getClippedAABB(const Vertex *buffer, const AABB &aabb) const {
	/* Reserve room for some additional vertices */
	Point vertices1[MAX_VERTS], vertices2[MAX_VERTS];
	int nVertices = 3;

#if defined(MTS_DEBUG_KD)
	AABB origAABB;
#endif


	for (int i=0; i<3; ++i) {
		vertices1[i] = buffer[idx[i]].p;
#if defined(MTS_DEBUG_KD)
		origAABB.expandBy(vertices1[i]);
#endif
	}

	for (int axis=0; axis<3; ++axis) {
		nVertices = sutherlandHodgman(vertices1, nVertices, vertices2, axis, aabb.min[axis], true);
		nVertices = sutherlandHodgman(vertices2, nVertices, vertices1, axis, aabb.max[axis], false);
	}

	AABB result;
	for (int i=0; i<nVertices; ++i) 
		result.expandBy(vertices1[i]);
	result.clip(aabb);

	return result;
}

MTS_NAMESPACE_END

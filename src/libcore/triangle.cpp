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

#include <mitsuba/core/triangle.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

Point Triangle::sample(const Point *positions, const Normal *normals,
		const Point2 *texCoords, Normal &normal, Point2 &uv, const Point2 &sample) const {
	const Point &p0 = positions[idx[0]];
	const Point &p1 = positions[idx[1]];
	const Point &p2 = positions[idx[2]];

	Point2 bary = warp::squareToUniformTriangle(sample);
	Vector sideA = p1 - p0, sideB = p2 - p0;
	Point p = p0 + (sideA * bary.x) + (sideB * bary.y);

	if (normals) {
		const Normal &n0 = normals[idx[0]];
		const Normal &n1 = normals[idx[1]];
		const Normal &n2 = normals[idx[2]];

		normal = Normal(normalize(
			n0 * (1.0f - bary.x - bary.y) +
			n1 * bary.x + n2 * bary.y
		));
	} else {
		normal = Normal(normalize(cross(sideA, sideB)));
	}

	if (texCoords) {
		const Point2 &uv0 = texCoords[idx[0]];
		const Point2 &uv1 = texCoords[idx[1]];
		const Point2 &uv2 = texCoords[idx[2]];

		uv = uv0 * (1.0f - bary.x - bary.y) +
			uv1 * bary.x + uv2 * bary.y;
	} else {
		uv = bary;
	}

	return p;
}

Float Triangle::surfaceArea(const Point *positions) const {
	const Point &p0 = positions[idx[0]];
	const Point &p1 = positions[idx[1]];
	const Point &p2 = positions[idx[2]];
	Vector sideA = p1 - p0, sideB = p2 - p0;
	return 0.5f * cross(sideA, sideB).length();
}

#define MAX_VERTS 10

static int sutherlandHodgman(Point3d *input, int inCount, Point3d *output, int axis,
		double splitPos, bool isMinimum) {
	if (inCount < 3)
		return 0;

	Point3d cur       = input[0];
	double sign       = isMinimum ? 1.0f : -1.0f;
	double distance   = sign * (cur[axis] - splitPos);
	bool  curIsInside = (distance >= 0);
	int   outCount    = 0;

	for (int i=0; i<inCount; ++i) {
		int nextIdx = i+1;
		if (nextIdx == inCount)
			nextIdx = 0;
		Point3d next = input[nextIdx];
		distance = sign * (next[axis] - splitPos);
		bool nextIsInside = (distance >= 0);

		if (curIsInside && nextIsInside) {
			/* Both this and the next vertex are inside, add to the list */
			SAssertEx(outCount + 1 < MAX_VERTS, "Overflow in sutherlandHodgman()!");
			output[outCount++] = next;
		} else if (curIsInside && !nextIsInside) {
			/* Going outside -- add the intersection */
			double t = (splitPos - cur[axis]) / (next[axis] - cur[axis]);
			SAssertEx(outCount + 1 < MAX_VERTS, "Overflow in sutherlandHodgman()!");
			Point3d p = cur + (next - cur) * t;
			p[axis] = splitPos; // Avoid roundoff errors
			output[outCount++] = p;
		} else if (!curIsInside && nextIsInside) {
			/* Coming back inside -- add the intersection + next vertex */
			double t = (splitPos - cur[axis]) / (next[axis] - cur[axis]);
			SAssertEx(outCount + 2 < MAX_VERTS, "Overflow in sutherlandHodgman()!");
			Point3d p = cur + (next - cur) * t;
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

AABB Triangle::getClippedAABB(const Point *positions, const AABB &aabb) const {
	/* Reserve room for some additional vertices */
	Point3d vertices1[MAX_VERTS], vertices2[MAX_VERTS];
	int nVertices = 3;

	/* The kd-tree code will frequently call this function with
	   almost-collapsed AABBs. It's extremely important not to introduce
	   errors in such cases, otherwise the resulting tree will incorrectly
	   remove triangles from the associated nodes. Hence, do the
	   following computation in double precision! */
	for (int i=0; i<3; ++i)
		vertices1[i] = Point3d(positions[idx[i]]);

	for (int axis=0; axis<3; ++axis) {
		nVertices = sutherlandHodgman(vertices1, nVertices, vertices2, axis, aabb.min[axis], true);
		nVertices = sutherlandHodgman(vertices2, nVertices, vertices1, axis, aabb.max[axis], false);
	}

	AABB result;
	for (int i=0; i<nVertices; ++i) {
		for (int j=0; j<3; ++j) {
			double pos = vertices1[i][j];
			result.min[j] = std::min(result.min[j], (Float) math::castflt_down(pos));
			result.max[j] = std::max(result.max[j], (Float) math::castflt_up(pos));
		}
	}
	result.clip(aabb);

	return result;
}

MTS_NAMESPACE_END

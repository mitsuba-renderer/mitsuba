/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#if !defined(__TRIANGLE_H)
#define __TRIANGLE_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN

/** 
 * \brief Simple triangle class including a collection of routines 
 * for analysis and transformation. 
 *
 * Triangles are stored as indices into a vertex array
 * \ingroup libcore
 */
struct MTS_EXPORT_CORE Triangle {
	/// Indices into a vertex buffer
	uint32_t idx[3];

	/// Construct an axis-aligned box, which contains the triangle
	inline AABB getAABB(const Point *positions) const {
		AABB result(positions[idx[0]]);
		result.expandBy(positions[idx[1]]);
		result.expandBy(positions[idx[2]]);
		return result;
	}

	/**
	 * \brief Returns the axis-aligned bounding box of a triangle after it has 
	 * clipped to the extends of another, given AABB. 
	 *
	 * This function uses the Sutherland-Hodgman algorithm to calculate the 
	 * convex polygon, which is created when applying all 6 AABB splitting 
	 * planes to the triangle. Afterwards, the AABB of the newly created 
	 * convex polygon is returned. This function is an important component 
	 * for efficiently creating 'Perfect Split' KD-trees. For more detail, 
	 * see "On building fast kd-Trees for Ray Tracing, and on doing 
	 * that in O(N log N)" by Ingo Wald and Vlastimil Havran
	 */
	AABB getClippedAABB(const Point *positions, const AABB &aabb) const;

	/// Uniformly sample a point on the triangle and return its normal
	Point sample(const Point *positions, const Normal *normals,
			Normal &n, const Point2 &seed) const;

	/// Calculate the surface area of this triangle
	Float surfaceArea(const Point *positions) const;

	/** \brief Ray-triangle intersection test
	 * 
	 * Uses the algorithm presented by Moeller and Trumbore at
	 * http://www.acm.org/jgt/papers/MollerTrumbore97/code.html
	 * Returns true if an intersection has been detected
	 * On success, \a t contains the distance from the ray origin to the
	 * intersection point, and \a u and \a v contain the intersection point in
	 * the local triangle coordinate system
	 */
	FINLINE static bool rayIntersect(const Point &p0, const Point &p1, const Point &p2, 
		const Ray &ray, Float &u, Float &v, Float &t) {
		/* find vectors for two edges sharing v[0] */
		Vector edge1 = p1 - p0, edge2 = p2 - p0;

		/* begin calculating determinant - also used to calculate U parameter */
		Vector pvec = cross(ray.d, edge2);

		/* if determinant is near zero, ray lies in plane of triangle */
		Float det = dot(edge1, pvec);

		if (det > -1e-8f && det < 1e-8f)
			return false;
		Float inv_det = 1.0f / det;

		/* calculate distance from v[0] to ray origin */
		Vector tvec = ray.o - p0;

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

		/* ray intersects triangle -> compute t */
		t = dot(edge2, qvec) * inv_det;

		return true;
	}

	/** \brief Ray-triangle intersection test
	 * 
	 * Uses the algorithm presented by Moeller and Trumbore at
	 * http://www.acm.org/jgt/papers/MollerTrumbore97/code.html
	 * Returns true if an intersection has been detected
	 * On success, \a t contains the distance from the ray origin to the
	 * intersection point, and \a u and \a v contain the intersection point in
	 * the local triangle coordinate system
	 */
	FINLINE bool rayIntersect(const Point *positions, const Ray &ray, Float &u, 
		Float &v, Float &t) const {
		return rayIntersect(
			positions[idx[0]], positions[idx[1]],
			positions[idx[2]], ray, u, v, t);
	}
};

MTS_NAMESPACE_END

#endif /* __TRIANGLE_H */

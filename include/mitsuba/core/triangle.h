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

#if !defined(__TRIANGLE_H)
#define __TRIANGLE_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Mesh vertex data structure. Stores 3D coordinates,
 * vertex normals, UV coordinates and a tangent frame.
 */
struct Vertex {
	Point p;     ///< %Position
	Normal n;    ///< %Normal
	Point2 uv;   ///< %Texture coordinates
	Vector dpdu; ///< Partial derivative of the position with respect to \a u.
	Vector dpdv; ///< Partial derivative of the position with respect to \a v.
#if defined(MTS_HAS_VERTEX_COLORS)
	Float color[3];
#endif

	inline bool operator==(const Vertex &vert) const {
#if defined(MTS_HAS_VERTEX_COLORS)
		return (p == vert.p && n == vert.n && uv == vert.uv 
			     && dpdu == vert.dpdu && dpdv == vert.dpdv
				 && color[0] == vert.color[0]
				 && color[1] == vert.color[1]
				 && color[2] == vert.color[2]);
#else
		return (p == vert.p && n == vert.n && uv == vert.uv 
			     && dpdu == vert.dpdu && dpdv == vert.dpdv);
#endif
	}
	
	inline bool operator!=(const Vertex &vert) const {
		return !operator==(vert);
	}
};

/** 
 * \brief Simple triangle class including a collection of routines 
 * for analysis and transformation. 
 *
 * Triangles are stored as indices into a vertex array
 */
struct MTS_EXPORT_CORE Triangle {
	/// Indices into a vertex buffer
	uint32_t idx[3];

	/// Construct an axis-aligned box, which contains the triangle
	inline AABB getAABB(const Vertex *buffer) const {
		AABB result(buffer[idx[0]].p);
		result.expandBy(buffer[idx[1]].p);
		result.expandBy(buffer[idx[2]].p);
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
	AABB getClippedAABB(const Vertex *buffer, const AABB &aabb) const;

	/** \brief Ray-triangle intersection test
	 * 
	 * Uses the algorithm presented by Moeller and Trumbore at
	 * http://www.acm.org/jgt/papers/MollerTrumbore97/code.html
	 * Returns true if an intersection has been detected
	 * On success, \a t contains the distance from the ray origin to the
	 * intersection point, and \a u and \a v contain the intersection point in
	 * the local triangle coordinate system
	 */
	bool rayIntersect(const Vertex *buffer, const Ray &ray, Float &u, 
		Float &v, Float &t) const;

	/// Uniformly sample a point on the triangle and return its normal
	Point sample(const Vertex *buffer, Normal &n, const Point2 &seed) const;

	/// Calculate the surface area of this triangle
	Float surfaceArea(const Vertex *buffer) const;
};

MTS_NAMESPACE_END

#endif /* __TRIANGLE_H */

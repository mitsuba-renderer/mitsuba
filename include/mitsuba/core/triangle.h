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

#pragma once
#if !defined(__MITSUBA_CORE_TRIANGLE_H_)
#define __MITSUBA_CORE_TRIANGLE_H_

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
     * clipped to the extends of another given AABB.
     *
     * This function uses the Sutherland-Hodgman algorithm to calculate the
     * convex polygon that is created when applying all 6 AABB splitting
     * planes to the triangle. Afterwards, the AABB of the newly created
     * convex polygon is returned. This function is an important component
     * for efficiently creating 'Perfect Split' KD-trees. For more detail,
     * see "On building fast kd-Trees for Ray Tracing, and on doing
     * that in O(N log N)" by Ingo Wald and Vlastimil Havran
     */
    AABB getClippedAABB(const Point *positions, const AABB &aabb) const;

    // Returns the bounding sphere of the triangle
    inline BSphere getBSphere(const Point *positions) const {
        Vector a = (positions[idx[1]] - positions[idx[0]]);
        Vector b = (positions[idx[2]] - positions[idx[0]]);
        Float a2 = dot(a, a);
        Float b2 = dot(b, b);
        Float da = std::sqrt(a2);
        Float db = std::sqrt(b2);
        Vector axb = cross(a, b);
        Float axb2 = dot(axb, axb);
        Float daxb = std::sqrt(axb2);
        return BSphere(positions[idx[0]] + cross(a2 * b - b2 * a, axb) / (2 * axb2),
                       da * db * (a - b).length() / (2 * daxb));
    }

    /// Uniformly sample a point on the triangle and return its normal and UV coordinates
    Point sample(const Point *positions, const Normal *normals,
            const Point2 *texCoords, Normal &n, Point2 &uv,
            const Point2 &seed) const;

    /// Calculate the surface area of this triangle
    Float surfaceArea(const Point *positions) const;

    /** \brief Ray-triangle intersection test
     *
     * Uses the algorithm by Moeller and Trumbore discussed at
     * <tt>http://www.acm.org/jgt/papers/MollerTrumbore97/code.html</tt>.
     *
     * \param p0
     *    Position of the first vertex
     * \param p1
     *    Position of the second vertex
     * \param p2
     *    Position of the third vertex
     * \param ray
     *    The ray segment to be used for the intersection query
     * \param t
     *    Upon success, \a t contains the distance from the ray origin to the
     *    intersection point,
     * \param u
     *   Upon success, \c u will contain the 'U' component of the intersection
     *   in barycentric coordinates
     * \param v
     *   Upon success, \c v will contain the 'V' component of the intersection
     *   in barycentric coordinates
     * \return
     *   \c true if an intersection has been detected
     */
    FINLINE static bool rayIntersect(const Point &p0, const Point &p1, const Point &p2,
        const Ray &ray, Float &u, Float &v, Float &t) {
        /* Find vectors for two edges sharing */
        Vector edge1 = p1 - p0, edge2 = p2 - p0;

        /* Begin calculating determinant - also used to calculate U parameter */
        Vector pvec = cross(ray.d, edge2);

        Float det = dot(edge1, pvec);
        if (det == 0)
            return false;
        Float inv_det = 1.0f / det;

        /* Calculate distance from v[0] to ray origin */
        Vector tvec = ray.o - p0;

        /* Calculate U parameter and test bounds */
        u = dot(tvec, pvec) * inv_det;
        if (u < 0.0 || u > 1.0)
            return false;

        /* Prepare to test V parameter */
        Vector qvec = cross(tvec, edge1);

        /* Calculate V parameter and test bounds */
        v = dot(ray.d, qvec) * inv_det;

        /* Inverted comparison (to catch NaNs) */
        if (v >= 0.0 && u + v <= 1.0) {
            /* ray intersects triangle -> compute t */
            t = dot(edge2, qvec) * inv_det;

            return true;
        }

        return false;
    }

    /** \brief Ray-triangle intersection test
     *
     * Uses the algorithm by Moeller and Trumbore discussed at
     * <tt>http://www.acm.org/jgt/papers/MollerTrumbore97/code.html</tt>.
     *
     * \param positions
     *    Pointer to the vertex positions of the underlying triangle mesh
     * \param index
     *    Index of the triangle that should be intersected
     * \param ray
     *    The ray segment to be used for the intersection query
     * \param t
     *    Upon success, \a t contains the distance from the ray origin to the
     *    intersection point,
     * \param u
     *   Upon success, \c u will contain the 'U' component of the intersection
     *   in barycentric coordinates
     * \param v
     *   Upon success, \c v will contain the 'V' component of the intersection
     *   in barycentric coordinates
     * \return
     *   \c true if an intersection has been detected
     */
    FINLINE bool rayIntersect(const Point *positions, const Ray &ray, Float &u,
        Float &v, Float &t) const {
        return rayIntersect(
            positions[idx[0]], positions[idx[1]],
            positions[idx[2]], ray, u, v, t);
    }
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_TRIANGLE_H_ */

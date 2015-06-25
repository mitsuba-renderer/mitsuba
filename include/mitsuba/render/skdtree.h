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
#if !defined(__MITSUBA_RENDER_SKDTREE_H_)
#define __MITSUBA_RENDER_SKDTREE_H_

#include <mitsuba/render/shape.h>
#include <mitsuba/render/sahkdtree3.h>
#include <mitsuba/render/triaccel.h>

#if defined(MTS_KD_CONSERVE_MEMORY)
#if defined(MTS_HAS_COHERENT_RT)
#error MTS_KD_CONSERVE_MEMORY & MTS_HAS_COHERENT_RT are incompatible
#endif
#endif

#if defined(SINGLE_PRECISION)
/// 64 byte temporary storage for intersection computations
#define MTS_KD_INTERSECTION_TEMP 64
#else
#define MTS_KD_INTERSECTION_TEMP 128
#endif

MTS_NAMESPACE_BEGIN

typedef const Shape * ConstShapePtr;

/**
 * \brief SAH KD-tree acceleration data structure for fast ray-triangle
 * intersections.
 *
 * Implements the construction algorithm for 'perfect split' trees as outlined
 * in the paper "On Bulding fast kd-Trees for Ray Tracing, and on doing that in
 * O(N log N)" by Ingo Wald and Vlastimil Havran. Non-triangle shapes are
 * supported, but most optimizations here target large triangle meshes.
 * For more details regarding the construction algorithm, please refer to
 * the class \ref GenericKDTree.
 *
 * This class offers a choice of two different triangle intersection algorithms:
 * By default, intersections are computed using the "TriAccel" projection with
 * pre-computation method from Ingo Wald's PhD thesis "Realtime Ray Tracing
 * and Interactive Global Illumination". This adds an overhead of 48 bytes per
 * triangle.
 *
 * When compiled with \c MTS_KD_CONSERVE_MEMORY, the Moeller-Trumbore intersection
 * test is used instead, which doesn't need any extra storage. However, it also
 * tends to be quite a bit slower.
 *
 * \sa GenericKDTree
 * \ingroup librender
 */

class MTS_EXPORT_RENDER ShapeKDTree : public SAHKDTree3D<ShapeKDTree> {
	friend class GenericKDTree<AABB, SurfaceAreaHeuristic3, ShapeKDTree>;
	friend class SAHKDTree3D<ShapeKDTree>;
	friend class Instance;
	friend class AnimatedInstance;
	friend class SingleScatter;

public:
	// =============================================================
	//! @{ \name Initialization and tree construction
	// =============================================================
	/// Create an empty kd-tree
	ShapeKDTree();

	/// Add a shape to the kd-tree
	void addShape(const Shape *shape);

	/// Return the list of stored shapes
	inline const std::vector<const Shape *> &getShapes() const { return m_shapes; }

	/**
	 * \brief Return the total number of low-level primitives (triangles
	 * and other low-level primitives)
	 */
	inline SizeType getPrimitiveCount() const {
		return m_shapeMap[m_shapeMap.size()-1];
	}

	/// Return an axis-aligned bounding box containing all primitives
	inline const AABB &getAABB() const { return m_aabb; }

	/// Build the kd-tree (needs to be called before tracing any rays)
	void build();

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Ray tracing routines
	// =============================================================

	/**
	 * \brief Intersect a ray against all primitives stored in the kd-tree
	 * and return detailed intersection information
	 *
	 * \param ray
	 *    A 3-dimensional ray data structure with minimum/maximum
	 *    extent information, as well as a time (which applies when
	 *    the shapes are animated)
	 *
	 * \param its
	 *    A detailed intersection record, which will be filled by the
	 *    intersection query
	 *
	 * \return \c true if an intersection was found
	 */
	bool rayIntersect(const Ray &ray, Intersection &its) const;

	/**
	 * \brief Intersect a ray against all primitives stored in the kd-tree
	 * and return the traveled distance and intersected shape
	 *
	 * This function represents a performance compromise when the
	 * intersected shape must be known, but there is no need for
	 * a detailed intersection record.
	 *
	 * \param ray
	 *    A 3-dimensional ray data structure with minimum/maximum
	 *    extent information, as well as a time (which applies when
	 *    the shapes are animated)
	 *
	 * \param t
	 *    The traveled ray distance will be stored in this parameter

	 * \param shape
	 *    A pointer to the intersected shape will be stored in this
	 *    parameter
	 *
	 * \param n
	 *    The geometric surface normal will be stored in this parameter
	 *
	 * \param uv
	 *    The UV coordinates associated with the intersection will
	 *    be stored here.
	 *
	 * \return \c true if an intersection was found
	 */
	bool rayIntersect(const Ray &ray, Float &t, ConstShapePtr &shape,
		Normal &n, Point2 &uv) const;

	/**
	 * \brief Test a ray for occlusion with respect to all primitives
	 *    stored in the kd-tree.
	 *
	 * This function does not compute a detailed intersection record,
	 * and it never determines the closest intersection, which makes
	 * it quite a bit faster than the other two \c rayIntersect() methods.
	 * However, for this reason, it can only be used to check whether
	 * there is \a any occlusion along a ray or ray segment.
	 *
	 * \param ray
	 *    A 3-dimensional ray data structure with minimum/maximum
	 *    extent information, as well as a time (which applies when
	 *    the shapes are animated)
	 *
	 * \return \c true if there is occlusion
	 */
	bool rayIntersect(const Ray &ray) const;

#if defined(MTS_HAS_COHERENT_RT)
	/**
	 * \brief Intersect four rays with the stored triangle meshes while making
	 * use of ray coherence to do this very efficiently. Requires SSE.
	 */
	void rayIntersectPacket(const RayPacket4 &packet,
		const RayInterval4 &interval, Intersection4 &its, void *temp) const;

	/**
	 * \brief Fallback for incoherent rays
	 * \sa rayIntesectPacket
	 */
	void rayIntersectPacketIncoherent(const RayPacket4 &packet,
		const RayInterval4 &interval, Intersection4 &its, void *temp) const;
#endif
	//! @}
	// =============================================================

	MTS_DECLARE_CLASS()
protected:
	/**
	 * \brief Return the shape index corresponding to a primitive index
	 * seen by the generic kd-tree implementation.
	 *
	 * When this is a triangle mesh, the \a idx parameter is updated to the
	 * triangle index within the mesh.
	 */
	FINLINE IndexType findShape(IndexType &idx) const {
		std::vector<IndexType>::const_iterator it = std::lower_bound(
				m_shapeMap.begin(), m_shapeMap.end(), idx + 1) - 1;
		idx -= *it;
		return (IndexType) (it - m_shapeMap.begin());
	}

 	/// Return the axis-aligned bounding box of a certain primitive
	FINLINE AABB getAABB(IndexType idx) const {
		IndexType shapeIdx = findShape(idx);
		const Shape *shape = m_shapes[shapeIdx];
		if (m_triangleFlag[shapeIdx]) {
			const TriMesh *mesh = static_cast<const TriMesh *>(shape);
			return mesh->getTriangles()[idx].getAABB(mesh->getVertexPositions());
		} else {
			return shape->getAABB();
		}
	}

 	/// Return the AABB of a primitive when clipped to another AABB
	FINLINE AABB getClippedAABB(IndexType idx, const AABB &aabb) const {
		IndexType shapeIdx = findShape(idx);
		const Shape *shape = m_shapes[shapeIdx];
		if (m_triangleFlag[shapeIdx]) {
			const TriMesh *mesh = static_cast<const TriMesh *>(shape);
			return mesh->getTriangles()[idx].getClippedAABB(mesh->getVertexPositions(), aabb);
		} else {
			return shape->getClippedAABB(aabb);
		}
	}

	/// Temporarily holds some intersection information
	struct IntersectionCache {
		SizeType shapeIndex;
		SizeType primIndex;
		Float u, v;
	};

	/**
	 * Check whether a primitive is intersected by the given ray. Some
	 * temporary space is supplied to store data that can later
	 * be used to create a detailed intersection record.
	 */
	FINLINE bool intersect(const Ray &ray, IndexType idx, Float mint,
		Float maxt, Float &t, void *temp) const {
		IntersectionCache *cache =
			static_cast<IntersectionCache *>(temp);

#if defined(MTS_KD_CONSERVE_MEMORY)
		IndexType shapeIdx = findShape(idx);
		if (EXPECT_TAKEN(m_triangleFlag[shapeIdx])) {
			const TriMesh *mesh =
				static_cast<const TriMesh *>(m_shapes[shapeIdx]);
			const Triangle &tri = mesh->getTriangles()[idx];
			Float tempU, tempV, tempT;
			if (tri.rayIntersect(mesh->getVertexPositions(), ray,
						tempU, tempV, tempT)) {
				if (tempT < mint || tempT > maxt)
					return false;
				t = tempT;
				cache->shapeIndex = shapeIdx;
				cache->primIndex = idx;
				cache->u = tempU;
				cache->v = tempV;
				return true;
			}
		} else {
			const Shape *shape = m_shapes[shapeIdx];
			if (shape->rayIntersect(ray, mint, maxt, t,
					reinterpret_cast<uint8_t*>(temp) + 2*sizeof(IndexType))) {
				cache->shapeIndex = shapeIdx;
				cache->primIndex = KNoTriangleFlag;
				return true;
			}
		}
#else
		const TriAccel &ta = m_triAccel[idx];
		if (EXPECT_TAKEN(m_triAccel[idx].k != KNoTriangleFlag)) {
			Float tempU, tempV, tempT;
			if (ta.rayIntersect(ray, mint, maxt, tempU, tempV, tempT)) {
				t = tempT;
				cache->shapeIndex = ta.shapeIndex;
				cache->primIndex = ta.primIndex;
				cache->u = tempU;
				cache->v = tempV;
				return true;
			}
		} else {
			uint32_t shapeIndex = ta.shapeIndex;
			const Shape *shape = m_shapes[shapeIndex];
			if (shape->rayIntersect(ray, mint, maxt, t,
					reinterpret_cast<uint8_t*>(temp) + 2*sizeof(IndexType))) {
				cache->shapeIndex = shapeIndex;
				cache->primIndex = KNoTriangleFlag;
				return true;
			}
		}
#endif
		return false;
	}

	/**
	 * Check whether a primitive is intersected by the given ray. This
	 * version is used for shadow rays, hence no temporary space is supplied.
	 */
	FINLINE bool intersect(const Ray &ray, IndexType idx,
			Float mint, Float maxt) const {
#if defined(MTS_KD_CONSERVE_MEMORY)
		IndexType shapeIdx = findShape(idx);
		if (EXPECT_TAKEN(m_triangleFlag[shapeIdx])) {
			const TriMesh *mesh =
				static_cast<const TriMesh *>(m_shapes[shapeIdx]);
			const Triangle &tri = mesh->getTriangles()[idx];
			Float tempU, tempV, tempT;
			if (tri.rayIntersect(mesh->getVertexPositions(), ray, tempU, tempV, tempT))
				return tempT >= mint && tempT <= maxt;
			return false;
		} else {
			const Shape *shape = m_shapes[shapeIdx];
			return shape->rayIntersect(ray, mint, maxt);
		}
#else
		const TriAccel &ta = m_triAccel[idx];
		uint32_t shapeIndex = ta.shapeIndex;
		const Shape *shape = m_shapes[shapeIndex];
		if (EXPECT_TAKEN(m_triAccel[idx].k != KNoTriangleFlag)) {
			Float tempU, tempV, tempT;
			return ta.rayIntersect(ray, mint, maxt, tempU, tempV, tempT);
		} else {
			return shape->rayIntersect(ray, mint, maxt);
		}
#endif
	}

	/**
	 * \brief After having found a unique intersection, fill a proper record
	 * using the temporary information collected in \ref intersect()
	 */
	template<bool BarycentricPos> FINLINE void fillIntersectionRecord(const Ray &ray,
			const void *temp, Intersection &its) const {
		const IntersectionCache *cache = reinterpret_cast<const IntersectionCache *>(temp);
		const Shape *shape = m_shapes[cache->shapeIndex];
		if (m_triangleFlag[cache->shapeIndex]) {
			const TriMesh *trimesh = static_cast<const TriMesh *>(shape);
			const Triangle &tri = trimesh->getTriangles()[cache->primIndex];
			const Point *vertexPositions = trimesh->getVertexPositions();
			const Normal *vertexNormals = trimesh->getVertexNormals();
			const Point2 *vertexTexcoords = trimesh->getVertexTexcoords();
			const Color3 *vertexColors = trimesh->getVertexColors();
			const TangentSpace *vertexTangents = trimesh->getUVTangents();
			const Vector b(1 - cache->u - cache->v, cache->u, cache->v);

			const uint32_t idx0 = tri.idx[0], idx1 = tri.idx[1], idx2 = tri.idx[2];
			const Point &p0 = vertexPositions[idx0];
			const Point &p1 = vertexPositions[idx1];
			const Point &p2 = vertexPositions[idx2];

			if (BarycentricPos)
				its.p = p0 * b.x + p1 * b.y + p2 * b.z;
			else
				its.p = ray(its.t);

			Vector side1(p1-p0), side2(p2-p0);
			Normal faceNormal(cross(side1, side2));
			Float length = faceNormal.length();
			if (!faceNormal.isZero())
				faceNormal /= length;

			if (EXPECT_NOT_TAKEN(vertexTangents)) {
				const TangentSpace &ts = vertexTangents[cache->primIndex];
				its.dpdu = ts.dpdu;
				its.dpdv = ts.dpdv;
			} else {
				its.dpdu = side1;
				its.dpdv = side2;
			}

			if (EXPECT_TAKEN(vertexNormals)) {
				const Normal
					&n0 = vertexNormals[idx0],
					&n1 = vertexNormals[idx1],
					&n2 = vertexNormals[idx2];

				its.shFrame.n = normalize(n0 * b.x + n1 * b.y + n2 * b.z);

				/* Ensure that the geometric & shading normals face the same direction */
				if (dot(faceNormal, its.shFrame.n) < 0)
					faceNormal = -faceNormal;
			} else {
				its.shFrame.n = faceNormal;
			}
			its.geoFrame = Frame(faceNormal);

			if (EXPECT_TAKEN(vertexTexcoords)) {
				const Point2 &t0 = vertexTexcoords[idx0];
				const Point2 &t1 = vertexTexcoords[idx1];
				const Point2 &t2 = vertexTexcoords[idx2];
				its.uv = t0 * b.x + t1 * b.y + t2 * b.z;
			} else {
				its.uv = Point2(b.y, b.z);
			}

			if (EXPECT_NOT_TAKEN(vertexColors)) {
				const Color3 &c0 = vertexColors[idx0],
							 &c1 = vertexColors[idx1],
							 &c2 = vertexColors[idx2];
				Color3 result(c0 * b.x + c1 * b.y + c2 * b.z);
				its.color.fromLinearRGB(result[0], result[1],
					result[2], Spectrum::EReflectance);
			}

			its.shape = trimesh;
			its.hasUVPartials = false;
			its.primIndex = cache->primIndex;
			its.instance = NULL;
			its.time = ray.time;
		} else {
			shape->fillIntersectionRecord(ray,
				reinterpret_cast<const uint8_t*>(temp) + 2*sizeof(IndexType), its);
		}

		computeShadingFrame(its.shFrame.n, its.dpdu, its.shFrame);
		its.wi = its.toLocal(-ray.d);
	}

	/// Plain shadow ray query (used by the 'instance' plugin)
	inline bool rayIntersect(const Ray &ray, Float _mint, Float _maxt) const {
		Float mint, maxt, tempT = std::numeric_limits<Float>::infinity();
		if (m_aabb.rayIntersect(ray, mint, maxt)) {
			if (_mint > mint) mint = _mint;
			if (_maxt < maxt) maxt = _maxt;

			if (EXPECT_TAKEN(maxt > mint))
				return rayIntersectHavran<true>(ray, mint, maxt, tempT, NULL);
		}
		return false;
	}

	/// Plain intersection query (used by the 'instance' plugin)
	inline bool rayIntersect(const Ray &ray, Float _mint, Float _maxt, Float &t, void *temp) const {
		Float mint, maxt, tempT = std::numeric_limits<Float>::infinity();
		if (m_aabb.rayIntersect(ray, mint, maxt)) {
			if (_mint > mint) mint = _mint;
			if (_maxt < maxt) maxt = _maxt;

			if (EXPECT_TAKEN(maxt > mint)) {
				if (rayIntersectHavran<false>(ray, mint, maxt, tempT, temp)) {
					t = tempT;
					return true;
				}
			}
		}
		return false;
	}

	/// Virtual destructor
	virtual ~ShapeKDTree();
private:
	std::vector<const Shape *> m_shapes;
	std::vector<bool> m_triangleFlag;
	std::vector<IndexType> m_shapeMap;
#if !defined(MTS_KD_CONSERVE_MEMORY)
	TriAccel *m_triAccel;
#endif
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_SKDTREE_H_ */

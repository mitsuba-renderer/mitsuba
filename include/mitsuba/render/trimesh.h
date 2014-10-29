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
#if !defined(__MITSUBA_RENDER_TRIMESH_H_)
#define __MITSUBA_RENDER_TRIMESH_H_

#include <mitsuba/core/triangle.h>
#include <mitsuba/core/pmf.h>
#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Simple tangent space storage for surfaces
 *
 * Note that the \ref dpdu and \ref dpdv vectors are not
 * necessarily orthogonal.
 *
 * \ingroup librender
 * \ingroup libpython
 */
struct TangentSpace {
	/// Position partial with respect to the U parameter of the local chart
	Vector dpdu;

	/// Position partial with respect to the V parameter of the local chart
	Vector dpdv;

	inline TangentSpace() { }
	inline TangentSpace(const Vector &dpdu, const Vector &dpdv)
		: dpdu(dpdu), dpdv(dpdv) { }
	inline TangentSpace(Stream *stream) :
		dpdu(stream), dpdv(stream) {
	}

	inline void serialize(Stream *stream) {
		dpdu.serialize(stream);
		dpdv.serialize(stream);
	}

	inline std::string toString() const {
		std::ostringstream oss;
		oss << "TangentSpace[dpdu=" << dpdu.toString() << ", dpdv=" << dpdv.toString() << "]";
		return oss.str();
	}
};

/** \brief Abstract triangle mesh base class
 * \ingroup librender
 * \ingroup libpython
 */
class MTS_EXPORT_RENDER TriMesh : public Shape {
public:
	/// Create a new, empty triangle mesh with the specified state
	TriMesh(const std::string &name,
			size_t triangleCount, size_t vertexCount,
			bool hasNormals = false,
			bool hasTexcoords = false,
			bool hasVertexColors = false,
			bool flipNormals = false,
			bool faceNormals = false);

	/// Unserialize a triangle mesh
	TriMesh(Stream *stream, InstanceManager *manager);

	/**
	 * \brief Unserialize a triangle mesh
	 *
	 * This is an alternative routine, which only loads triangle data
	 * (no BSDF, Sub-surface integrator, etc.) in a format that will
	 * remain stable as Mitsuba evolves. The files can optionally contain
	 * multiple meshes -- in that case, the specified index determines
	 * which one to load.
	 */
	TriMesh(Stream *stream, int idx = 0);

	// =============================================================
	//! @{ \name General query functions
	// =============================================================

	/// Return the total surface area
	Float getSurfaceArea() const;

	/// Return a bounding box containing the mesh
	AABB getAABB() const;

	/// Return a bounding box containing the mesh
	inline AABB &getAABB() { return m_aabb; }

	/**
	 * \brief Create a triangle mesh approximation of this shape
	 *
	 * Since instances are already triangle meshes, the implementation
	 * just returns a pointer to \a this.
	 */
	ref<TriMesh> createTriMesh();

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Access to the stored triangle mesh
	// =============================================================

	/// Return the number of triangles
	inline size_t getTriangleCount() const { return m_triangleCount; }
	/// Return the number of vertices
	inline size_t getVertexCount() const { return m_vertexCount; }

	/// Return the triangle list (const version)
	inline const Triangle *getTriangles() const { return m_triangles; };
	/// Return the triangle list
	inline Triangle *getTriangles() { return m_triangles; };

	/// Return the vertex positions (const version)
	inline const Point *getVertexPositions() const { return m_positions; };
	/// Return the vertex positions
	inline Point *getVertexPositions() { return m_positions; };

	/// Return the vertex normals (const version)
	inline const Normal *getVertexNormals() const { return m_normals; };
	/// Return the vertex normals
	inline Normal *getVertexNormals() { return m_normals; };
	/// Does the mesh have vertex normals?
	inline bool hasVertexNormals() const { return m_normals != NULL; };

	/// Return the vertex colors (const version)
	inline const Color3 *getVertexColors() const { return m_colors; };
	/// Return the vertex colors
	inline Color3 *getVertexColors() { return m_colors; };
	/// Does the mesh have vertex colors?
	inline bool hasVertexColors() const { return m_colors != NULL; };

	/// Return the vertex texture coordinates (const version)
	inline const Point2 *getVertexTexcoords() const { return m_texcoords; };
	/// Return the vertex texture coordinates
	inline Point2 *getVertexTexcoords() { return m_texcoords; };
	/// Does the mesh have vertex texture coordinates?
	inline bool hasVertexTexcoords() const { return m_texcoords != NULL; };

	/// Return the per-triangle UV tangents (const version)
	inline const TangentSpace *getUVTangents() const { return m_tangents; };
	/// Return the per-triangle UV tangents
	inline TangentSpace *getUVTangents() { return m_tangents; };
	/// Does the mesh have UV tangent information?
	inline bool hasUVTangents() const { return m_tangents != NULL; };

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Sampling routines
	// =============================================================

	/**
	 * \brief Sample a point on the surface of this shape instance
	 * (with respect to the area measure)
	 *
	 * The returned sample density will be uniform over the surface.
	 *
	 * \param pRec
	 *     A position record, which will be used to return the sampled
	 *     position, as well as auxilary information about the sample.
	 *
	 * \param sample
	 *     A uniformly distributed 2D vector
	 */
	void samplePosition(PositionSamplingRecord &pRec,
			const Point2 &sample) const;

	/**
	 * \brief Query the probability density of \ref samplePosition() for
	 * a particular point on the surface.
	 *
	 * This method will generally return the inverse of the surface area.
	 *
	 * \param pRec
	 *     A position record, which will be used to return the sampled
	 *     position, as well as auxilary information about the sample.
	 */

	Float pdfPosition(const PositionSamplingRecord &pRec) const;

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Miscellaneous
	// =============================================================

	/**
	 * \brief Generate per-triangle space basis vectors from
	 * a user-specified set of UV coordinates
	 *
	 * Will throw an exception when no UV coordinates are
	 * associated with the mesh.
	 */
	void computeUVTangents();

	/**
	 * \brief Generate smooth vertex normals?
	 *
	 * \param force
	 *   When this parameter is set to true, the function
	 *   generates normals <em>even</em> when there are
	 *   already existing ones.
	 */
	void computeNormals(bool force = false);

	/**
	 * \brief Rebuild the mesh so that adjacent faces
	 * with a dihedral angle greater than \c maxAngle degrees
	 * are topologically disconnected.
	 *
	 * On the other hand, if the angle is less than \a maxAngle, the code
	 * ensures that the faces  reference the same vertices.
	 * This step is very useful as a pre-process when generating
	 * high-quality smooth shading normals on meshes with creases.
	 * Note: this function is fairly memory intensive and will require
	 * approximately three 3x the storate used by the input mesh.
	 * It will never try to merge vertices with equal positions but
	 * different UV coordinates or vertex colors.
	 */
	void rebuildTopology(Float maxAngle);

	/// Serialize to a file/network stream
	void serialize(Stream *stream, InstanceManager *manager) const;

	/**
	 * \brief Serialize to a file/network stream
	 *
	 * This is an alternative routine, which \a only loads triangle
	 * data (no BSDF, Sub-surface integrator, etc.) in a format that
	 * will remain stable as Mitsuba evolves.
	 */
	void serialize(Stream *stream) const;

	/**
	 * \brief Build a discrete probability distribution
	 * for sampling.
	 *
	 * Called once while loading the scene
	 */
	virtual void configure();

	/**
	 * \brief Return the derivative of the normal vector with
	 * respect to the UV parameterization
	 *
	 * This can be used to compute Gaussian and principal curvatures,
	 * amongst other things.
	 *
	 * \param its
	 *     Intersection record associated with the query
	 * \param dndu
	 *     Parameter used to store the partial derivative of the
	 *     normal vector with respect to \c u
	 * \param dndv
	 *     Parameter used to store the partial derivative of the
	 *     normal vector with respect to \c v
	 * \param shadingFrame
	 *     Specifies whether to compute the derivative of the
	 *     geometric or shading normal of the surface
	 */
	void getNormalDerivative(const Intersection &its,
		Vector &dndu, Vector &dndv, bool shadingFrame) const;

	/**
	 * \brief Return the number of primitives (triangles, hairs, ..)
	 * contributed to the scene by this shape
	 *
	 * Does not include instanced geometry
	 */
	size_t getPrimitiveCount() const;

	/**
	 * \brief Return the number of primitives (triangles, hairs, ..)
	 * contributed to the scene by this shape
	 *
	 * Includes instanced geometry
	 */
	size_t getEffectivePrimitiveCount() const;

	/// Import a shape from the Blender in-memory representation
	static ref<TriMesh> fromBlender(const std::string &name, size_t faceCount, void *facePtr,
		size_t vertexCount, void *vertexPtr, void *uvPtr, void *colPtr, short matNr);

	/// Export a Wavefront OBJ version of this file
	void writeOBJ(const fs::path &path) const;

	/// Export a Stanford PLY version of this file
	void writePLY(const fs::path &path) const;

	/// Return a string representation
	std::string toString() const;

	//! @}
	// =============================================================

	MTS_DECLARE_CLASS()
protected:
	/// Create a new triangle mesh
	TriMesh(const Properties &props);

	/// Virtual destructor
	virtual ~TriMesh();

	/// Load a Mitsuba compressed triangle mesh substream
	void loadCompressed(Stream *stream, int idx = 0);

	/**
	 * \brief Reads the header information of a compressed file, returning
	 * the version ID.
	 *
	 * This function assumes the stream is at the beginning of the compressed
	 * file and leaves the stream located right after the header.
	 */
	static short readHeader(Stream *stream);

	/**
	 * \brief Read the idx-th entry from the offset diccionary at the end of
	 * the stream, which has to be open already, given the file version tag.
	 * This function modifies the position of the stream.
	 */
	static size_t readOffset(Stream *stream, short version, int idx);

	/**
	 * \brief Read the entirety of the end-of-file offset dictionary from the
	 * already open stream, replacing the contents of the input vector.
	 * If the file is not large enough the function returns -1
	 * and does not modify the vector.
	 * This function modifies the position of the stream.
	 */
	 static int readOffsetDictionary(Stream *stream, short version,
		 std::vector<size_t>& outOffsets);

	/// Prepare internal tables for sampling uniformly wrt. area
	void prepareSamplingTable();
protected:
	AABB m_aabb;
	Triangle *m_triangles;
	Point *m_positions;
	Normal *m_normals;
	Point2 *m_texcoords;
	TangentSpace *m_tangents;
	Color3 *m_colors;
	size_t m_triangleCount;
	size_t m_vertexCount;
	bool m_flipNormals;
	bool m_faceNormals;

	/* Surface and distribution -- generated on demand */
	DiscreteDistribution m_areaDistr;
	Float m_surfaceArea;
	Float m_invSurfaceArea;
	ref<Mutex> m_mutex;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_TRIMESH_H_ */

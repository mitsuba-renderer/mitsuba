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

#if !defined(__TRIMESH_H)
#define __TRIMESH_H

#include <mitsuba/core/triangle.h>
#include <mitsuba/core/pdf.h>
#include <mitsuba/render/shape.h>
#include <boost/filesystem/fstream.hpp>

namespace fs = boost::filesystem;

MTS_NAMESPACE_BEGIN

/**
 * \brief Simple tangent space storage for surfaces
 * \ingroup librender
 */
struct TangentSpace {
	/// Position partials wrt. the UV parameterization
	Vector dpdu;
	Vector dpdv;

	inline TangentSpace() { }
	inline TangentSpace(Stream *stream) :
		dpdu(stream), dpdv(stream) {
	}

	inline void serialize(Stream *stream) {
		dpdu.serialize(stream);
		dpdv.serialize(stream);
	}
};

/** \brief Abstract triangle mesh base class
 * \ingroup librender
 */
class MTS_EXPORT_RENDER TriMesh : public Shape {
public:
	/// Create a new, empty triangle mesh with the specified state
	TriMesh(const std::string &name, 
			size_t triangleCount, size_t vertexCount,
			bool hasNormals, bool hasTexcoords, 
			bool hasVertexColors, bool flipNormals = false,
			bool faceNormals = false);

	/// Unserialize a triangle mesh
	TriMesh(Stream *stream, InstanceManager *manager);

	/**
	 * Unserialize a triangle mesh - this is an alternative
	 * routine, which only loads triangle data (no BSDF,
	 * Sub-surface integrator, etc.) in a format that
	 * will remain stable as Mitsuba evolves. The files
	 * can optionally contain multiple meshes -- in that case,
	 * the specified index determines which one to load.
	 */
	TriMesh(Stream *stream, int idx = 0);
	
	// =============================================================
	//! @{ \name General query functions
	// =============================================================

	/// Return the name of this mesh
	virtual std::string getName() const;

	/// Return the total surface area
	virtual Float getSurfaceArea() const;

	/// Return a bounding box containing the mesh
	virtual AABB getAABB() const;

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
	inline const Spectrum *getVertexColors() const { return m_colors; };
	/// Return the vertex colors
	inline Spectrum *getVertexColors() { return m_colors; };
	/// Does the mesh have vertex colors?
	inline bool hasVertexColors() const { return m_colors != NULL; };

	/// Return the vertex texture coordinates (const version)
	inline const Point2 *getVertexTexcoords() const { return m_texcoords; };
	/// Return the vertex texture coordinates
	inline Point2 *getVertexTexcoords() { return m_texcoords; };
	/// Does the mesh have vertex texture coordinates?
	inline bool hasVertexTexcoords() const { return m_texcoords != NULL; };

	/// Return the vertex tangents (const version)
	inline const TangentSpace *getVertexTangents() const { return m_tangents; };
	/// Return the vertex tangents
	inline TangentSpace *getVertexTangents() { return m_tangents; };
	/// Does the mesh have vertex tangents?
	inline bool hasVertexTangents() const { return m_tangents != NULL; };

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Sampling routines
	// =============================================================

	/// Sample a point on the mesh
	Float sampleArea(ShapeSamplingRecord &sRec, const Point2 &sample) const;

	/**
	 * \brief Return the probability density of sampling the 
	 * given point using \ref sampleArea()
	 */
	Float pdfArea(const ShapeSamplingRecord &sRec) const;
	
	//! @}
	// =============================================================
	
	// =============================================================
	//! @{ \name Miscellaneous
	// =============================================================

	/**
	 * \brief Generate tangent space basis vectors. 
	 *
	 * Returns true upon success. A failure is caused by
	 * missing texture coordinates
	 */
	bool computeTangentSpaceBasis();

	/// Generate surface normals
	void computeNormals();

	/**
	 * \brief Regenerate the triangulation so that adjacent faces
	 * with an angle greater than \a maxAngle degrees become disconnected.
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

	/// Export an Wavefront OBJ version of this file
	void writeOBJ(const fs::path &path) const;
	
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
protected:
	std::string m_name;
	AABB m_aabb;
	DiscretePDF m_areaPDF;
	Triangle *m_triangles;
	Point *m_positions;
	Normal *m_normals;
	Point2 *m_texcoords;
	TangentSpace *m_tangents;
	Spectrum *m_colors;
	size_t m_triangleCount;
	size_t m_vertexCount;
	bool m_flipNormals;
	bool m_faceNormals;
	Float m_surfaceArea;
	Float m_invSurfaceArea;
};

MTS_NAMESPACE_END

#endif /* __TRIMESH_H */

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

#if !defined(__TRIMESH_H)
#define __TRIMESH_H

#include <mitsuba/core/triangle.h>
#include <mitsuba/core/pdf.h>
#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

/** \brief Abstract triangle mesh base class
 */
class MTS_EXPORT_RENDER TriMesh : public Shape {
public:
	/// Create a new, empty triangle mesh with the given triangle and vertex count
	TriMesh(size_t triangleCount, size_t vertexCount);

	/// Create a new, empty triangle mesh with the specified data
	TriMesh(const std::string &name, Transform worldToObject,
		Triangle *triangles, size_t triangleCount, 
		Vertex *vertexBuffer, size_t vertexCount);

	/// Unserialize a triangle mesh
	TriMesh(Stream *stream, InstanceManager *manager);

	/**
	 * Unserialize a triangle mesh - this is an alternative
	 * routine, which only loads triangle data (no BSDF,
	 * Sub-surface integrator, etc.) in a format that
	 * will remain stable as mitsuba evolves.
	 */
	TriMesh(Stream *stream);

	/// Return the triangle list
	inline const Triangle *getTriangles() const { return m_triangles; };
	
	/// Return the triangle list
	inline Triangle *getTriangles() { return m_triangles; };

	/// Return the number of triangles
	inline size_t getTriangleCount() const { return m_triangleCount; }

	/// Return the vertex buffer
	inline const Vertex *getVertexBuffer() const { return m_vertexBuffer; };
	
	/// Return the vertex buffer
	inline Vertex *getVertexBuffer() { return m_vertexBuffer; };

	/// Return the number of vertices
	inline size_t getVertexCount() const { return m_vertexCount; }

	/// Sample a point on the mesh
	Float sampleArea(ShapeSamplingRecord &sRec, const Point2 &sample) const;

	/// Generate tangent space basis vectors
	void calculateTangentSpaceBasis(bool hasNormals, bool hasTexCoords, bool complain = true);

	/// Serialize to a file/network stream
	void serialize(Stream *stream, InstanceManager *manager) const;

	/**
	 * Serialize to a file/network stream - this is an alternative
	 * routine, which only loads triangle data (no BSDF,
	 * Sub-surface integrator, etc.) in a format that
	 * will remain stable as mitsuba evolves.
	 */
	void serialize(Stream *stream) const;

	/// Build a discrete probability distribution for sampling. Called once after parsing
	virtual void configure();

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Create a new triangle mesh
	TriMesh(const Properties &props);

	/// Virtual destructor
	virtual ~TriMesh();
protected:
	DiscretePDF m_areaPDF;
	Triangle *m_triangles;
	size_t m_triangleCount;
	Vertex *m_vertexBuffer;
	size_t m_vertexCount;
	bool m_flipNormals;
};

MTS_NAMESPACE_END

#endif /* __TRIMESH_H */

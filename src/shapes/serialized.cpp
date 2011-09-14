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

#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{serialized}{Serialized mesh loader}
 * \order{11}
 * \parameters{
 *     \parameter{filename}{\String}{
 *	     Filename of the gemoetry file that should be loaded
 *	   }
 *     \parameter{faceNormals}{\Boolean}{
 *       When set to \code{true}, Mitsuba will use face normals when rendering
 *       the object, which will give it a faceted apperance. \default{\code{false}}
 *	   }
 *     \parameter{flipNormals}{\Boolean}{
 *       Optional flag to flip all normals. \default{\code{false}, i.e.
 *       the normals are left unchanged}.
 *	   }
 *     \parameter{toWorld}{\Transform}{
 *	      Specifies an optional linear object-to-world transformation.
 *        Note that non-uniform scales are not permitted!
 *        \default{none (i.e. object space $=$ world space)}
 *     }
 * }
 * This plugin represents the most space and time-efficient way
 * of getting geometry into Mitsuba. It uses a highly efficient
 * lossless compressed format for geometry storage. The format will
 * be explained on this page in a subsequent revision of the
 * documentation.
 */
class SerializedMesh : public TriMesh {
public:
	SerializedMesh(const Properties &props) : TriMesh(props) {
		fs::path filePath = Thread::getThread()->getFileResolver()->resolve(
			props.getString("filename"));

		/* Object-space -> World-space transformation */
		Transform objectToWorld = props.getTransform("toWorld", Transform());

		/// When the file contains multiple meshes, this index specifies which one to load 
		int shapeIndex = props.getInteger("shapeIndex", 0);

		m_name = (props.getID() != "unnamed") ? props.getID() 
			: formatString("%s@%i", filePath.stem().c_str(), shapeIndex); 

		/* Load the geometry */
		Log(EInfo, "Loading shape %i from \"%s\" ..", shapeIndex, filePath.leaf().c_str());
		ref<FileStream> stream = new FileStream(filePath, FileStream::EReadOnly);
		stream->setByteOrder(Stream::ELittleEndian);
		ref<TriMesh> mesh = new TriMesh(stream, shapeIndex);
		m_triangleCount = mesh->getTriangleCount();
		m_vertexCount = mesh->getVertexCount();

		m_positions = new Point[m_vertexCount];
		memcpy(m_positions, mesh->getVertexPositions(), sizeof(Point) * m_vertexCount);

		if (mesh->hasVertexNormals()) {
			m_normals = new Normal[m_vertexCount];
			memcpy(m_normals, mesh->getVertexNormals(), sizeof(Normal) * m_vertexCount);
		}

		if (mesh->hasVertexTexcoords()) {
			m_texcoords = new Point2[m_vertexCount];
			memcpy(m_texcoords, mesh->getVertexTexcoords(), sizeof(Point2) * m_vertexCount);
		}

		if (mesh->hasVertexColors()) {
			m_colors = new Spectrum[m_vertexCount];
			memcpy(m_colors, mesh->getVertexColors(), sizeof(Spectrum) * m_vertexCount);
		}

		m_triangles = new Triangle[m_triangleCount];
		memcpy(m_triangles, mesh->getTriangles(), sizeof(Triangle) * m_triangleCount);

		if (!objectToWorld.isIdentity()) {
			for (size_t i=0; i<m_vertexCount; ++i)
				m_positions[i] = objectToWorld(m_positions[i]);
			if (m_normals) {
				for (size_t i=0; i<m_vertexCount; ++i) 
					m_normals[i] = objectToWorld(m_normals[i]);
			}
		}
		if (objectToWorld.det3x3() < 0) {
			for (size_t i=0; i<m_triangleCount; ++i) {
				Triangle &t = m_triangles[i];
				std::swap(t.idx[0], t.idx[1]);
			}
		}
	}

	SerializedMesh(Stream *stream, InstanceManager *manager) : TriMesh(stream, manager) {
	}

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(SerializedMesh, false, TriMesh)
MTS_EXPORT_PLUGIN(SerializedMesh, "Serialized mesh loader");
MTS_NAMESPACE_END

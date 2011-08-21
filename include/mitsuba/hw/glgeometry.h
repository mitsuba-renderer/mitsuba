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

#if !defined(__GLGEOMETRY_H)
#define __GLGEOMETRY_H

#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/hw/gpugeometry.h>

MTS_NAMESPACE_BEGIN

/** \brief OpenGL-based GPUGeometry implementation
 * \ingroup libhw
 */
class MTS_EXPORT_HW GLGeometry : public GPUGeometry {
	friend class GLRenderer;
public:
	/// Create a new GLGeometry instance for the given mesh
	GLGeometry(const TriMesh *mesh);

	/// Upload the geometry object
	void init();

	/// Refresh (re-upload) the geometry object
	void refresh();

	/// Bind the geometry object
	void bind();

	/// Unbind the geometry object
	void unbind();

	/// Free the geometry object from GPU memory
	void cleanup();
	
	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GLGeometry();
protected:
	GLuint m_vertexID;
	GLuint m_indexID;
	GLuint64 m_vertexAddr;
	GLuint64 m_indexAddr;
	GLuint m_vertexSize;
	GLuint m_indexSize;
	int m_stride;
};

MTS_NAMESPACE_END

#endif /* __GLGEOMETRY_H */

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
#if !defined(__MITSUBA_HW_GLGEOMETRY_H_)
#define __MITSUBA_HW_GLGEOMETRY_H_

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
    enum EIdentifier {
        EVertexID = 0,
        EIndexID = 1
    };

    GLuint m_id[2];
    GLuint64 m_addr[2];
    GLuint m_size[2];
    int m_stride;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_GLGEOMETRY_H_ */

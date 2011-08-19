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

#if !defined(__GLRENDERER_H)
#define __GLRENDERER_H

#include <mitsuba/hw/renderer.h>
#if defined(__OSX__)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <OpenGL/glext.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#endif

struct GLEWContextStruct;

/// Get the per-thread context for GLEW-MX
extern MTS_EXPORT_HW GLEWContextStruct *glewGetContext();

MTS_NAMESPACE_BEGIN

/**
 * \brief Specifies the maximum number of triangles that can
 * be sent to the GPU in one batch. This is useful to
 * avoid freezing the OS when dealing with huge inputs.
 */
#define MTS_GL_MAX_QUEUED_TRIS 500000

/**
 * \brief OpenGL implementation of the \ref Renderer interface
 * \ingroup libhw
 */
class MTS_EXPORT_HW GLRenderer : public Renderer {
public:
	/// Construct a new OpenGL rendering interface
	GLRenderer(Session *session);

	/// Initialize the renderer
	virtual void init(Device *device, Renderer *other = NULL);

	/**
	 * \brief Reconfigure the renderer for a certain device
	 * (e.g. after a resize event)
	 */
	void reconfigure(const Device *device);

	/// Shut the renderer down
	virtual void shutdown();

	/// Create a new GPU texture object
	GPUTexture *createGPUTexture(const std::string &name,
		Bitmap *bitmap = NULL);

	/// Create a new GPU geometry object
	GPUGeometry *createGPUGeometry(const TriMesh *mesh);

	/// Create a new GPU program object
	GPUProgram *createGPUProgram(const std::string &name);

	/// Create a new synchronization object
	GPUSync *createGPUSync();

	/// Clear the viewport
	void clear();

	/// Configure the camera
	void setCamera(const ProjectiveCamera *pCamera);

	/// Configure the camera (supports a pixel offset)
	void setCamera(const ProjectiveCamera *pCamera,
		const Point2 &jitter);

	/// Configure the camera (manual)
	void setCamera(const Matrix4x4 &proj, const Matrix4x4 &view);

	/// Set up the renderer for drawing triangle geometry
	void beginDrawingMeshes(bool transmitOnlyPositions = false); 

	/// Send a triangle mesh to the renderer
	void drawTriMesh(const TriMesh *mesh); 

	/// Draw a quad using the given texture
	void blitTexture(const GPUTexture *texture,
		bool flipVertically = false, 
		bool centerHoriz = true, bool centerVert = true,
		const Vector2i &offset = Vector2i(0, 0));

	/// Clean up the renderer after drawing triangle geometry
	void endDrawingMeshes(); 

	/**
	 * \brief Quickly draw all geometry that has been registered 
	 * with the renderer.
	 *
	 * Only transmits positions, hence this is mainly useful for
	 * shadow mapping.
	 */
	void drawAll(const std::vector<std::pair<const GPUGeometry *, Transform> > &geo);

	/// Blit a screen-sized quad
	void blitQuad(bool flipVertically);

	/**
	 * Draw a line of text on the screen. The coordinates are specified
	 * in pixel coordinates, where the upper left corner is the origin
	 */
	void drawText(const Point2i &pos, 
			const Font *font, const std::string &text);

	/// Set the size of point primitives
	void setPointSize(Float size);

	/// Draw a point
	void drawPoint(const Point &p);

	/// Draw a line between two specified points
	void drawLine(const Point &a, const Point &b);

	/// Draw an ellipse with the specified center and axes
	void drawEllipse(const Point &center, 
			const Vector &axis1, const Vector &axis2);

	/// Draw a wire-frame axis-aligned box
	void drawAABB(const AABB &aabb);

	/// Set a depth offset for shadow mapping (0 to disable)
	void setDepthOffset(Float value);

	/// Set the currently active blending mode
	void setBlendMode(EBlendMode mode);

	/// Set the currently active culling mode
	void setCullMode(ECullMode mode);

	/// Activate or deactivate the writing of depth information
	void setDepthMask(bool value);

	/// Activate or deactivate depth testing
	void setDepthTest(bool value);

	/// Activate or deactivate the writing of color information
	void setColorMask(bool value);

	/// Set the current fixed-function pipeline color
	void setColor(const Spectrum &spec, Float alpha = 1.0f);

	/// Push a view transformation onto the matrix stack
	void pushTransform(const Transform &trafo);
	
	/// Pop the last view transformation from the matrix stack
	void popTransform();

	/// Flush outstanding rendering commands
	void flush();

	/// Completely finish outstanding rendering commands
	void finish();

	/// Check for any error indications
	void checkError(bool onlyWarn = true);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GLRenderer();
protected:
	bool m_transmitOnlyPositions;
	bool m_normalsEnabled;
	bool m_texcoordsEnabled;
	bool m_tangentsEnabled;
	bool m_colorsEnabled;
	size_t m_queuedTriangles;
	int m_stride;
};

MTS_NAMESPACE_END

#endif /* __GLRENDERER_H */

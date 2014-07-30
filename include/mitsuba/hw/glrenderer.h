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
#if !defined(__MITSUBA_HW_GLRENDERER_H_)
#define __MITSUBA_HW_GLRENDERER_H_

#include <mitsuba/hw/renderer.h>

struct GLEWContextStruct;

/// Get the per-thread context for GLEW-MX
extern MTS_EXPORT_HW GLEWContextStruct *glewGetContext();

MTS_NAMESPACE_BEGIN

/**
 * \brief Specifies the maximum number of triangles that can
 * be sent to the GPU in one batch. This is useful to
 * avoid freezing the OS when dealing with huge inputs.
 */
#define MTS_GL_MAX_QUEUED_TRIS 250000

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
	GPUGeometry *createGPUGeometry(const Shape *shape);

	/// Create a new GPU program object
	GPUProgram *createGPUProgram(const std::string &name);

	/// Create a new synchronization object
	GPUSync *createGPUSync();

	/// Clear the viewport
	void clear();

	/// Configure the camera
	void setCamera(const ProjectiveCamera *pCamera,
		const Point2 &apertureSample = Point2(0.5f),
		const Point2 &aaSample = Point2(0.5f),
		Float timeSample = 0.5f);

	/// Configure the camera (manual)
	void setCamera(const Matrix4x4 &proj, const Matrix4x4 &view);

	/// Directly set the modelview or projection matrix
	void setMatrix(EMatrixType type, const Matrix4x4 &value);

	/// Fetch the currently set modelview or projection matrix
	Matrix4x4 getMatrix(EMatrixType type) const;

	/// Set up the renderer for drawing triangle geometry
	void beginDrawingMeshes(bool transmitOnlyPositions = false);

	/// Send a triangle mesh to the renderer
	void drawMesh(const TriMesh *geo);

	/// Send a triangle mesh to the renderer
	void drawMesh(const GPUGeometry *geo);

	/// Clean up the renderer after drawing triangle geometry
	void endDrawingMeshes();

	/**
	 * \brief Quickly draw all geometry that has been registered
	 * with the renderer.
	 *
	 * Only transmits positions, hence this is mainly useful for
	 * shadow mapping.
	 */
	void drawAll(const std::vector<TransformedGPUGeometry> &allGeometry);

	/// Draw a quad using the given texture
	void blitTexture(const GPUTexture *texture,
		bool flipVertically = false,
		bool centerHoriz = true, bool centerVert = true,
		const Vector2i &offset = Vector2i(0, 0));

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

	/// Draw a point (2D)
	void drawPoint(const Point2 &p);

	/// Draw a point (2D, integer coordinates)
	void drawPoint(const Point2i &p);

	/// Draw a line between two specified points (2D)
	void drawLine(const Point2 &a, const Point2 &b);

	/// Draw a line between two specified points (2D, integer coordinates)
	void drawLine(const Point2i &a, const Point2i &b);

	/// Draw a rectangle between two specified points (2D)
	void drawRectangle(const Point2 &a, const Point2 &b);

	/// Draw a rectangle between two specified points (2D, integer coordinates)
	void drawRectangle(const Point2i &a, const Point2i &b);

	/// Draw a filled rectangle between two specified points (2D)
	void drawFilledRectangle(const Point2 &a, const Point2 &b);

	/// Draw a filled rectangle between two specified points (2D, integer coordinates)
	void drawFilledRectangle(const Point2i &a, const Point2i &b);

	/// Draw an ellipse with the specified center and axes
	void drawEllipse(const Point &center,
			const Vector &axis1, const Vector &axis2);

	/// Draw a wire-frame axis-aligned box
	void drawAABB(const AABB &aabb);

	/// Set the currently active blending mode
	void setBlendMode(EBlendMode mode);

	/// Set the currently active culling mode
	void setCullMode(ECullMode mode);

	/// Activate or deactivate the writing of depth information
	void setDepthMask(bool value);

	/// Activate or deactivate depth testing
	void setDepthTest(bool value);

	/// Set the current fixed-function pipeline color
	void setColor(const Color3 &color, Float alpha = 1.0f);

	/// Set the current fixed-function pipeline color
	void setColor(const Spectrum &spec, Float alpha = 1.0f);

	/// Set the depth value that is written by \ref clear()
	void setClearDepth(Float depth);

	/// Set the color value that is written by \ref clear()
	void setClearColor(const Color3 &color);

	/// Clear the view and projection transformations
	void clearTransforms();

	/// Flush outstanding rendering commands
	void flush();

	/// Completely finish outstanding rendering commands
	void finish();

	/// Check for any error indications
	void checkError(bool onlyWarn = true);

	/**
	 * \brief Send a debug string to the rendering backend
	 *
	 * This is mainly useful when an OpenGL trace is captured
	 * by a tool such as 'apitrace'.
	 */
	void debugString(const std::string &text);

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

#endif /* __MITSUBA_HW_GLRENDERER_H_ */

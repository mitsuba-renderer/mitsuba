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

#if !defined(__RENDERER_H)
#define __RENDERER_H

#include <mitsuba/hw/device.h>
#include <mitsuba/render/camera.h>
#include <mitsuba/render/trimesh.h>

MTS_NAMESPACE_BEGIN

class GPUTexture;
class GPUGeometry;
class GPUProgram;
class GPUSync;
class Bitmap;
class Font;

/**
 * \brief Helper class, which documents the capabilities of a renderer implementation
 * \ingroup libhw
 */
class MTS_EXPORT_HW RendererCapabilities : public Object {
public:
	enum ECapability {
		EShadingLanguage = 0,
		ERenderToTexture,
		EBufferBlit,
		EFloatingPointBuffer,
		EFloatingPointTextures,
		EMultisampleRenderToTexture,
		EVertexBufferObjects,
		EGeometryShaders,
		ESyncObjects,
		EBindless,
		ECapabilityCount
	};

	inline RendererCapabilities() {
		memset(m_capabilities, 0, sizeof(m_capabilities));
	}

	inline void setSupported(ECapability cap, bool supported) {
		m_capabilities[cap] = supported;
	}

	inline bool isSupported(ECapability cap) const { 
		return m_capabilities[cap];
	}

	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	bool m_capabilities[ECapabilityCount];
};

/** \brief Abstract renderer implementation
 * \ingroup libhw
 */
class MTS_EXPORT_HW Renderer : public Object {
public:
	/* Possible blending modes */
	enum EBlendMode {
		EBlendNone,          // Blending turned off
		EBlendAlpha,         // Normal alpha blending
		EBlendAdditive       // Additive blending
	};

	/* Possible culling modes */
	enum ECullMode {
		ECullNone,     // No culling
		ECullFront,    // Front-face culling
		ECullBack      // Back-face culling
	};

	// Instantiate a new renderer using the appropriate implementation.
	static Renderer *create(Session *session);

	/// Return the renderer's capabilities
	inline const RendererCapabilities *getCapabilities() const {
		return m_capabilities.get();
	}

	/**
	 * Initialize the renderer. Optionally, an existing renderer instance 
	 * can be provided as a second argument -- this establishes a link 
	 * between them to permit sharing of textures, programs, etc.
	 */
	virtual void init(Device *device, Renderer *other = NULL);

	/**
	 * \brief Reconfigure the renderer for a certain device
	 * (e.g. after a resize event)
	 */
	virtual void reconfigure(const Device *device) = 0;

	/// Shut the renderer down
	virtual void shutdown();

	/// Create a new GPU texture object
	virtual GPUTexture *createGPUTexture(const std::string &name,
		Bitmap *bitmap = NULL) = 0;

	/// Create a new GPU geometry object
	virtual GPUGeometry *createGPUGeometry(const TriMesh *mesh) = 0;

	/// Create a new GPU program object
	virtual GPUProgram *createGPUProgram(const std::string &name) = 0;

	/// Create a new synchronization object
	virtual GPUSync *createGPUSync() = 0;

	/// Clear the viewport
	virtual void clear() = 0;

	/// Configure the camera
	virtual void setCamera(const ProjectiveCamera *camera) = 0;

	/// Configure the camera (supports a pixel offset)
	virtual void setCamera(const ProjectiveCamera *camera,
		const Point2 &jitter) = 0;

	/// Configure the camera (manual)
	virtual void setCamera(const Matrix4x4 &proj, const Matrix4x4 &view) = 0;

	/// Set up the renderer for drawing triangle geometry
	virtual void beginDrawingMeshes(bool transmitOnlyPositions = false) = 0; 

	/// Send a triangle mesh to the renderer
	virtual void drawTriMesh(const TriMesh *mesh) = 0; 

	/// Clean up the renderer after drawing triangle geometry
	virtual void endDrawingMeshes() = 0; 

	/**
	 * \brief Quickly draw all geometry that has been registered 
	 * with the renderer.
	 *
	 * Only transmits positions, hence this is mainly useful for
	 * shadow mapping.
	 */
	virtual void drawAll(const std::vector<std::pair<const GPUGeometry *, Transform> > &geo) = 0;

	/// Draw a quad using the given texture
	virtual void blitTexture(const GPUTexture *texture,
		bool flipVertically = false, 
		bool centerHoriz = true, bool centerVert = true,
		const Vector2i &offset = Vector2i(0, 0)) = 0;

	/// Blit a screen-sized quad
	virtual void blitQuad(bool flipVertically) = 0;

	/**
	 * Draw a line of text on the screen. The coordinates are specified
	 * in pixel coordinates, where the upper left corner is the origin
	 */
	virtual void drawText(const Point2i &pos, 
			const Font *font, const std::string &text) = 0;

	/// Set the size of point primitives
	virtual void setPointSize(Float size) = 0;

	/// Draw a point
	virtual void drawPoint(const Point &p) = 0;

	/// Draw a line between two specified points
	virtual void drawLine(const Point &a, const Point &b) = 0;

	/// Draw an ellipse with the specified center and axes
	virtual void drawEllipse(const Point &center, 
			const Vector &axis1, const Vector &axis2) = 0;

	/// Draw a wire-frame axis-aligned box
	virtual void drawAABB(const AABB &aabb) = 0;

	/// Set a depth offset for shadow mapping (0 to disable)
	virtual void setDepthOffset(Float value) = 0;

	/// Set the currently active blending mode
	virtual void setBlendMode(EBlendMode mode) = 0;

	/// Set the currently active culling mode
	virtual void setCullMode(ECullMode mode) = 0;

	/// Activate or deactivate depth testing
	virtual void setDepthTest(bool value) = 0;

	/// Activate or deactivate the writing of depth information
	virtual void setDepthMask(bool value) = 0;

	/// Activate or deactivate the writing of color information
	virtual void setColorMask(bool value) = 0;

	/// Set the current fixed-function pipeline color
	virtual void setColor(const Spectrum &spec, Float alpha = 1.0f) = 0;

	/// Push a view transformation onto the matrix stack
	virtual void pushTransform(const Transform &trafo) = 0;
	
	/// Pop the last view transformation from the matrix stack
	virtual void popTransform() = 0;

	/// Flush outstanding rendering commands
	virtual void flush() = 0;

	/// Completely finish outstanding rendering commands
	virtual void finish() = 0;

	/// Check for any error indications
	virtual void checkError(bool onlyWarn = true) = 0;

	/**
	 * Register a shader with this renderer. Increases the
	 * reference count if it already exists
	 */
	Shader *registerShaderForResource(const HWResource *res);

	/// Look up a shader by the associated HWResource object
	Shader *getShaderForResource(const HWResource *res);

	/**
	 * Decrease the reference count of a shader. Deletes
	 * it when zero is reached
	 */
	void unregisterShaderForResource(const HWResource *res);

	/**
	 * Register a triangle mesh with the renderer. This
	 * will transfer the associated geometry to the GPU,
	 * which accelerates later calls to drawTriMesh()
	 */
	GPUGeometry *registerGeometry(const TriMesh *mesh);

	/// Unregister a triangle mesh from the renderer.
	void unregisterGeometry(const TriMesh *mesh);

	/// Set the log level
	inline void setLogLevel(ELogLevel logLevel) { m_logLevel = logLevel; }

	/// Set the log level for warnings
	inline void setWarnLogLevel(ELogLevel logLevel) { m_warnLogLevel = logLevel; }

	MTS_DECLARE_CLASS()
protected:
	/// Construct a new OpenI rendering interface
	Renderer(Session *session);

	/// Virtual destructor
	virtual ~Renderer();
protected:
	struct ShaderRecord {
		int refCount;
		Shader *shader;
	};

	ref<Session> m_session;
	ref<Device> m_device;
	ref<RendererCapabilities> m_capabilities;
	std::map<const HWResource*, ShaderRecord> m_shaders;
	std::map<const TriMesh*, GPUGeometry *> m_geometry;
	bool m_initialized, m_borrowed;
	std::string m_driverVendor;
	std::string m_driverRenderer;
	std::string m_driverVersion;
	ELogLevel m_logLevel, m_warnLogLevel;
};

MTS_NAMESPACE_END

#endif /* __RENDERER_H */

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

MTS_NAMESPACE_BEGIN

class MTS_EXPORT_HW GLRenderer : public Renderer {
public:
	/// Construct a new OpenGL rendering interface
	GLRenderer(Session *session);

	/// Initialize the renderer
	virtual void init(Device *device, Renderer *other = NULL);

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
	void setCamera(const Matrix4x4 *proj, const Matrix4x4 *view);

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

	/// Fast: draw all geometry that has been registered with the renderer (for shadow mapping)
	void drawAll();

	/// Blit a screen-sized quad
	void blitQuad(bool flipVertically);

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
	void setColor(const Spectrum &spec);

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
};

MTS_NAMESPACE_END

#endif /* __GLRENDERER_H */

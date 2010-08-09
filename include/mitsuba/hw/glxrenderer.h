#if !defined(__GLXRENDERER_H)
#define __GLXRENDERER_H

#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/hw/glxdevice.h>

MTS_NAMESPACE_BEGIN

/** \brief GLX (XFree86) renderer
 */
class MTS_EXPORT_HW GLXRenderer : public GLRenderer {
public:
	/// Create a new renderer
	GLXRenderer(X11Session *session);

	/// Return the rendering context
	inline GLXContext getGLXContext() { return m_context; }
	
	/// Initialize the renderer
	void init(Device *device, Renderer *other = NULL);

	/// Shut the renderer down
	void shutdown();

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GLXRenderer();
private:
    GLXContext m_context;
};

MTS_NAMESPACE_END

#endif /* __GLXRENDERER_H */

#if !defined(__WGLRENDERER_H)
#define __WGLRENDERER_H

#include <mitsuba/hw/wgldevice.h>
#include <mitsuba/hw/glrenderer.h>

MTS_NAMESPACE_BEGIN

/** \brief Windows (WGL) renderer implementation
 */
class MTS_EXPORT_HW WGLRenderer : public GLRenderer {
public:
	/// Create a new renderer
	WGLRenderer(WGLSession *session);

	/// Return the rendering context
	inline HGLRC getWGLContext() { return m_context; }

	/// Initialize the renderer
	void init(Device *device, Renderer *other = NULL);

	/// Shut the renderer down
	void shutdown();

	/// Lookup an OpenGL extension
	virtual void *lookupExtension(const std::string &name) const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~WGLRenderer();
private:
	HGLRC m_context;
};

MTS_NAMESPACE_END

#endif /* __WGLRENDERER_H */

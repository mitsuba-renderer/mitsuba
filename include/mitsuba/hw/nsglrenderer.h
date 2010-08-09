#if !defined(__NSGLRENDERER_H)
#define __NSGLRENDERER_H

#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/hw/nsgldevice.h>

#ifdef __OBJC__
#include <Cocoa/Cocoa.h>
#endif

MTS_NAMESPACE_BEGIN

/** \brief A MacOS X (NSGL) OpenGL Renderer
 */
class MTS_EXPORT_HW NSGLRenderer : public GLRenderer {
public:
	/// Create a new renderer
	NSGLRenderer(NSGLSession *session);

	/// Initialize the renderer
	void init(Device *device, Renderer *other = NULL);

	/// Shut the renderer down
	void shutdown();

	/// Return the rendering context
	inline void *getContext() { return m_context; }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~NSGLRenderer();

	/// Lookup an OpenGL extension
	void *lookupExtension(const std::string &name) const;
private:
#ifdef __OBJC__
	NSOpenGLContext *m_context;
#else
	void *m_context;
#endif
};

MTS_NAMESPACE_END

#endif /* __NSGLRENDERER_H */

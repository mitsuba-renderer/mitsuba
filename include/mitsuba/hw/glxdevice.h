#if !defined(__GLXDEVICE_H)
#define __GLXDEVICE_H

#include <mitsuba/hw/x11device.h>

MTS_NAMESPACE_BEGIN

/** \brief X Windows OpenGL-capable (GLX) device
 */
class MTS_EXPORT_HW GLXDevice : public X11Device {
public:
	/// Create a new device
	GLXDevice(X11Session *session);

	/// Flip the buffers (when using double-buffering)
	void flip();

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GLXDevice();

	/// Create a visual
	virtual XVisualInfo* createVisual();
};

MTS_NAMESPACE_END

#endif /* __GLXDEVICE_H */

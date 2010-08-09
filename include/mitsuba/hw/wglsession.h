#if !defined(__WGLSESSION_H)
#define __WGLSESSION_H

#include <mitsuba/hw/session.h>

MTS_NAMESPACE_BEGIN

/** \brief Windows (WGL) windowing environment session
 */
class MTS_EXPORT_HW WGLSession : public Session {
	friend class WGLDevice;
	friend class WGLRenderer;
public:
	/// Create a new session
	WGLSession();

	/// Initialize the session
	void init();

	/// Shut the session down
	void shutdown();

	/// Process all events and call event callbacks
	void processEvents();

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~WGLSession();
protected:
	HINSTANCE m_hinstance;
	std::string m_wndClassName;
};

MTS_NAMESPACE_END

#endif /* __WGLSESSION_H */

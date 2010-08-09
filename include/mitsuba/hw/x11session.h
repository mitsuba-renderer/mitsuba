#if !defined(__X11SESSION_H)
#define __X11SESSION_H

#include <mitsuba/hw/session.h>
#include <GL/glx.h>
#include <X11/extensions/xf86vmode.h>

MTS_NAMESPACE_BEGIN

/** \brief X Window System (X11R6) session
 */
class MTS_EXPORT_HW X11Session : public Session {
	friend class X11Device;
	friend class GLXDevice;
	friend class GLXRenderer;
public:
	/// Create a new session
	X11Session();

	/// Set the display name (eg. "localhost:0.0")
	void setDisplayName(const std::string &displayname);
	
	/// Initialize the session
	void init();

	/// Shut the session down
	void shutdown();

	/// Process all events and call event callbacks
	void processEvents();

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~X11Session();
protected:
	std::string m_displayName;
	Display *m_display;
	Window m_root;
	int m_screen;
	bool m_hasVidMode;
	bool m_hasGLX;
};

MTS_NAMESPACE_END

#endif /* __X11SESSION_H */

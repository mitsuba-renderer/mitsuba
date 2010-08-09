#if !defined(__X11DEVICE_H)
#define __X11DEVICE_H

#include <mitsuba/hw/device.h>
#include <mitsuba/hw/x11session.h>

MTS_NAMESPACE_BEGIN

/** \brief X Window System (X11R6) device / software surface
 */
class MTS_EXPORT_HW X11Device : public Device {
	friend class X11Session;
	friend class GLXRenderer;
public:
	/// Create a new device
	X11Device(X11Session *session);

	/// Initialize the device
	void init(Device *other = NULL);

	/// Shut the device down
	void shutdown();

	/// Flip the buffers
	virtual void flip();

	/// Only applies in windowed mode
	void setVisible(bool enabled);

	/// Only applies in windowed mode
	void setPosition(const Point2i &position);

	/// Set the window title
	void setTitle(const std::string &title);
	
	/// Display the X11 cursor?
	void showCursor(bool enabled);

	/// Set the cursor grab state
	void setGrab(bool grab);

	/// Move the mouse to another position
	void warpMouse(const Point2i &position);

	/// Handle an X11 event. Called by the session
	void processEvent(const XEvent &event);

	/// Associate a renderer with this device
	void makeCurrent(Renderer *renderer);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~X11Device();

	/// Create a visual
	virtual XVisualInfo* createVisual();
private:
	/// Translate X11 mouse events
	void translateMouse(const XEvent &xEvent, DeviceEvent &event);

	/// Translate X11 keyboard events
	bool translateKeyboard(const XEvent &xEvent, DeviceEvent &event);

	/// Return the X11 window
	inline Window getWindow() { return m_window; }

	/// Return the X11 visual
	inline XVisualInfo *getVisual() { return m_visinfo; }
protected:
	Window m_window;
	XVisualInfo *m_visinfo;
	XF86VidModeModeInfo m_previousMode;
	Atom m_deleteWindow;
	bool m_visible;
	Cursor m_cursor;
    int m_keymap[256];
	Point2i m_mouse;
	int m_modifierState;
	int m_buttonState;
	bool m_grab;
};

MTS_NAMESPACE_END

#endif /* __X11DEVICE_H */

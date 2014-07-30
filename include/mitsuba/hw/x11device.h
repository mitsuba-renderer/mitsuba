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

#if !defined(__MITSUBA_HW_X11DEVICE_H_)
#define __MITSUBA_HW_X11DEVICE_H_

#include <mitsuba/hw/device.h>
#include <mitsuba/hw/x11session.h>

MTS_NAMESPACE_BEGIN

/** \brief X Window System (X11R6) device / software surface
 * \ingroup libhw
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

#endif /* __MITSUBA_HW_X11DEVICE_H_ */

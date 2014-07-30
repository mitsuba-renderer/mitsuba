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

#if !defined(__MITSUBA_HW_NSGLDEVICE_H_)
#define __MITSUBA_HW_NSGLDEVICE_H_

#include <mitsuba/hw/device.h>
#include <mitsuba/hw/nsglsession.h>

MTS_NAMESPACE_BEGIN
class NSGLDevice;
MTS_NAMESPACE_END

#ifdef __OBJC__
#include <Cocoa/Cocoa.h>

@interface CustomView : NSView < NSWindowDelegate > {
	mitsuba::NSGLDevice *m_device;
	unsigned int m_keymap[0xFF];
	unsigned int m_modifiers;
	unsigned int m_buttonMask;
	bool m_mouseInWindow;
	bool m_firstMouseMotion;
	bool m_focused;
	bool m_ignoreNextMouseEvent;
}

- (BOOL) focused;
- (unsigned int) extractModifiers: (unsigned int) modifiers;
- (void) setDevice: (mitsuba::NSGLDevice *) device;
- (void) ignoreNextMouseEvent;
- (void) ignoreFirstMouseMotion;
@end

#endif

MTS_NAMESPACE_BEGIN

/** \brief A MacOS X (NSGL) device
 */
class MTS_EXPORT_HW NSGLDevice : public Device {
public:
	/// Create a new device
	NSGLDevice(NSGLSession *session);

	/// Initialize the device
	void init(Device *other = NULL);

	/// Shut the device down
	void shutdown();

	/// Flip the buffers
	void flip();

	/// Only applies in windowed mode
	void setVisible(bool enabled);

	/// Only applies in windowed mode
	void setPosition(const Point2i &position);

	/// Set the window title
	void setTitle(const std::string &title);

	/// Display the NSGL cursor?
	void showCursor(bool enabled);

	/// Move the mouse to another position
	void warpMouse(const Point2i &position);

	/// Associate a renderer with this device
	void makeCurrent(Renderer *renderer);

	/// Set the cursor grab state
	void setGrab(bool grab);

	// *************************************
	// ************* INTERNAL **************
	// *************************************

	/// Is the NSGL cursor shown?
	inline bool getCursor() const { return m_cursor; }

	/// Is the mouse inside the window?
	bool isMouseInWindow();

	/// Push a device event onto the stack
	void pushEvent(const DeviceEvent &event);

	/**
	 * Deliver all events which have been
	 * received asynchronously
	 */
	void processEvents();

	/// Return the NSWindow
	inline void *getWindow() { return m_window; }

	/// Return the pixel format
	inline void *getPixelFormat() { return m_fmt; }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~NSGLDevice();
private:
#ifdef __OBJC__
	NSWindow *m_window;
	CustomView *m_view;
	NSOpenGLPixelFormat *m_fmt;
	NSOpenGLContext *m_currentContext;
#else
	void *m_window;
	void *m_view;
	void *m_fmt;
	void *m_currentContext;
#endif
	bool m_visible;
	bool m_cursor;
	std::vector<DeviceEvent> m_deviceEvents;
	ref<Mutex> m_mutex;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_NSGLDEVICE_H_ */

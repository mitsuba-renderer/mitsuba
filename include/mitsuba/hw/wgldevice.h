#if !defined(__WGLDEVICE_H)
#define __WGLDEVICE_H

#include <mitsuba/hw/device.h>
#include <mitsuba/hw/wglsession.h>

MTS_NAMESPACE_BEGIN

/** \brief Windows (WGL) device implementation
 */
class MTS_EXPORT_HW WGLDevice : public Device {
public:
	/// Create a new device
	WGLDevice(WGLSession *session);

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
	
	/// Display the WGL cursor?
	void showCursor(bool enabled);

	/// Move the mouse to another position
	void warpMouse(const Point2i &position);

	/// Set the cursor grab state
	void setGrab(bool grab);

	/// Associate a renderer with this device
	void makeCurrent(Renderer *renderer);

	/// WNDPROC for event handling
	static LONG WINAPI WndProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam);

	/// Return the device context
	inline HDC getDeviceContext() const { return m_hdc; }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~WGLDevice();

	/// Initialize the pixel format
	void initPixelFormat(HWND hWnd);

	/// Translate a WIN32 key event
	bool translateKey(WPARAM vkey, LPARAM lParam, DeviceEvent &event);

	/// Translate a WIN32 mouse event
	bool translateMouse(UINT uMsg, WPARAM wParam, DeviceEvent &event);	
protected:
	ref<WGLDevice> m_parent;
	HWND m_hwnd;
	HDC m_hdc;
	PIXELFORMATDESCRIPTOR m_pfd;
	bool m_visible;
	Point2i m_mouse;
	int m_modifierState;
	int m_buttonState;
	int m_special[256];
	char m_std[256];
	bool m_leftShift;
	bool m_rightShift;
	bool m_grab;
	bool m_cursor;
	bool m_mouseInWindow;
	int m_pf;
};

MTS_NAMESPACE_END

#endif /* __WGLDEVICE_H */

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

#include <mitsuba/hw/wglrenderer.h>

#if defined(__OSX__)
#include <OpenGL/gl.h>
#else
#if defined(_WIN32)
# include <windows.h>
#endif
#include <GL/gl.h>
#endif

#define REPEATED_KEYMASK (1<<30)
#define DEFINE_SPECIAL(wsym, sym) m_special[wsym] = sym

#if !defined(WGL_SAMPLE_BUFFERS_ARB)
#define WGL_SAMPLE_BUFFERS_ARB  0x2041
#define WGL_SAMPLES_ARB     0x2042
#endif

#ifndef WGL_ARB_pixel_format
#define WGL_ARB_pixel_format 1
#define WGL_DRAW_TO_WINDOW_ARB                                  0x2001
#define WGL_ACCELERATION_ARB                                    0x2003
#define WGL_SUPPORT_OPENGL_ARB                                  0x2010
#define WGL_COLOR_BITS_ARB                                      0x2014
#define WGL_RED_BITS_ARB                                        0x2015
#define WGL_GREEN_BITS_ARB                                      0x2017
#define WGL_BLUE_BITS_ARB                                       0x2019
#define WGL_ALPHA_BITS_ARB                                      0x201B
#define WGL_DEPTH_BITS_ARB                                      0x2022
#define WGL_STENCIL_BITS_ARB                                    0x2023
#define WGL_DOUBLE_BUFFER_ARB                                   0x2011
#define WGL_FULL_ACCELERATION_ARB                               0x2027

typedef BOOL (APIENTRY * wglChoosePixelFormatARBProc) (HDC hdc, const int *piAttribIList,
        const FLOAT *pfAttribFList, UINT nMaxFormats, int *piFormats, UINT *nNumFormats);
#endif

MTS_NAMESPACE_BEGIN

/**
 * Needed because the WndProc will be called (during initialization)
 * without the pointer to the WGLDevice
 */
static WGLDevice *__global_workaround = NULL;

LONG WINAPI WGLDevice::WndProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
    LRESULT retval = 1;

    /* See comment above */
    WGLDevice *device = reinterpret_cast<WGLDevice *>(GetWindowLongPtr(hWnd, 0));
    if (device == NULL)
        device = __global_workaround;
    PAINTSTRUCT ps;

    DeviceEvent deviceEvent;
    switch (uMsg) {
        case WM_CREATE:
            device->initPixelFormat(hWnd);
            break;
        case WM_SIZE: {
                Vector2i size = Vector2i(LOWORD(lParam), HIWORD(lParam));
                if (device->m_size != size) {
                    device->m_size = size;
                    deviceEvent.setType(EResizeEvent);
                    device->fireDeviceEvent(deviceEvent);
                }
            }
            break;
        case WM_PAINT:
            BeginPaint(hWnd, &ps);
            if (device->getDoubleBuffer())
                SwapBuffers(device->m_hdc);
            EndPaint(hWnd, &ps);
            break;
        case WM_CLOSE:
            deviceEvent.setType(EQuitEvent);
            device->fireDeviceEvent(deviceEvent);
            break;
        case WM_MOUSELEAVE:
            device->m_buttonState = 0;
            device->m_mouseInWindow = false;
            break;
        case WM_KILLFOCUS:
            deviceEvent.setType(ELoseFocusEvent);
            device->fireDeviceEvent(deviceEvent);
            device->m_modifierState = 0;
            device->m_buttonState = 0;
            break;
        case WM_SETFOCUS:
            deviceEvent.setType(EGainFocusEvent);
            device->fireDeviceEvent(deviceEvent);
            device->m_modifierState = 0;
            device->m_buttonState = 0;
            break;
        case WM_SYSKEYUP:
        case WM_KEYUP:
            deviceEvent.setType(EKeyUpEvent);
            if (device->translateKey(wParam, lParam, deviceEvent)) {
                int special = deviceEvent.getKeyboardSpecial();

                if (special == EKeyLShift || special == EKeyRShift) {
                    device->m_modifierState = device->m_modifierState & (~EShiftModifier);
                } else if (special == EKeyLAlt || special == EKeyRAlt) {
                    device->m_modifierState = device->m_modifierState & (~EAltModifier);
                } else if (special == EKeyLControl || special == EKeyRControl) {
                    device->m_modifierState = device->m_modifierState & (~EControlModifier);
                }

                deviceEvent.setKeyboardModifiers(device->m_modifierState);
                device->fireDeviceEvent(deviceEvent);
            }
            break;
        case WM_SYSKEYDOWN:
        case WM_KEYDOWN:
            if (lParam & REPEATED_KEYMASK)
                break;
            deviceEvent.setType(EKeyDownEvent);
            if (device->translateKey(wParam, lParam, deviceEvent)) {
                int special = deviceEvent.getKeyboardSpecial();

                if (special == EKeyLShift || special == EKeyRShift) {
                    device->m_modifierState |= EShiftModifier;
                } else if (special == EKeyLAlt || special == EKeyRAlt) {
                    device->m_modifierState |= EAltModifier;
                } else if (special == EKeyLControl || special == EKeyRControl) {
                    device->m_modifierState |= EControlModifier;
                }

                deviceEvent.setKeyboardModifiers(device->m_modifierState);
                device->fireDeviceEvent(deviceEvent);
            }
            break;
        case WM_MOUSEWHEEL:
        case WM_MOUSEMOVE:
        case WM_LBUTTONDOWN:
        case WM_MBUTTONDOWN:
        case WM_RBUTTONDOWN:
        case WM_LBUTTONUP:
        case WM_MBUTTONUP:
        case WM_RBUTTONUP:
            if (!device->m_mouseInWindow) {
                device->m_mouseInWindow = true;
                TRACKMOUSEEVENT tme;
                tme.cbSize = sizeof(TRACKMOUSEEVENT);
                tme.dwFlags = TME_LEAVE;
                tme.hwndTrack = device->m_hwnd;
                if (!TrackMouseEvent(&tme)) {
                    Log(EError, "Cannot track mouse leave events!");
                }
            }
            deviceEvent.setMousePosition(Point2i(LOWORD(lParam), HIWORD(lParam)));
            if (device->translateMouse(uMsg, wParam, deviceEvent)) {
                if (deviceEvent.getType() == EMouseButtonDownEvent) {
                    if (deviceEvent.getMouseButton() != Device::EWheelUpButton
                        && deviceEvent.getMouseButton() != Device::EWheelDownButton)
                        device->m_buttonState |= deviceEvent.getMouseButton();
                } else if (deviceEvent.getType() == EMouseButtonUpEvent) {
                    device->m_buttonState &= ~deviceEvent.getMouseButton();
                } else if (deviceEvent.getType() == EMouseMotionEvent
                        && device->m_buttonState != 0) {
                    deviceEvent.setMouseButton(device->m_buttonState);
                    deviceEvent.setType(EMouseDragEvent);
                }

                device->fireDeviceEvent(deviceEvent);
            }
            break;
        default:
            retval = DefWindowProc(hWnd, uMsg, wParam, lParam);
            break;
    };

    return (LONG) retval;
}


WGLDevice::WGLDevice(WGLSession *session)
 : Device(session),
 m_hwnd(NULL),
 m_hdc(NULL),
 m_leftShift(false),
 m_rightShift(false),
 m_mouseInWindow(false) {
    m_title = "Mitsuba [wgl]";

    char c;
    for (int i=0; i<256; i++) {
        m_special[i] = ENoSpecial;
        m_std[i] = '\0';
    }

    DEFINE_SPECIAL(VK_BACK, EKeyBackspace);
    DEFINE_SPECIAL(VK_TAB, EKeyTab);
    DEFINE_SPECIAL(VK_CLEAR, EKeyClear);
    DEFINE_SPECIAL(VK_RETURN, EKeyReturn);
    DEFINE_SPECIAL(VK_PAUSE, EKeyPause);
    DEFINE_SPECIAL(VK_ESCAPE, EKeyEscape);
    DEFINE_SPECIAL(VK_DELETE, EKeyDelete);
    DEFINE_SPECIAL(VK_NUMPAD0, EKeyKeyPad0);
    DEFINE_SPECIAL(VK_NUMPAD1, EKeyKeyPad1);
    DEFINE_SPECIAL(VK_NUMPAD2, EKeyKeyPad2);
    DEFINE_SPECIAL(VK_NUMPAD3, EKeyKeyPad3);
    DEFINE_SPECIAL(VK_NUMPAD4, EKeyKeyPad4);
    DEFINE_SPECIAL(VK_NUMPAD5, EKeyKeyPad5);
    DEFINE_SPECIAL(VK_NUMPAD6, EKeyKeyPad6);
    DEFINE_SPECIAL(VK_NUMPAD7, EKeyKeyPad7);
    DEFINE_SPECIAL(VK_NUMPAD8, EKeyKeyPad8);
    DEFINE_SPECIAL(VK_NUMPAD9, EKeyKeyPad9);
    DEFINE_SPECIAL(VK_DECIMAL, EKeyKeyPadPeriod);
    DEFINE_SPECIAL(VK_DIVIDE, EKeyKeyPadDivide);
    DEFINE_SPECIAL(VK_MULTIPLY, EKeyKeyPadMultiply);
    DEFINE_SPECIAL(VK_SUBTRACT, EKeyKeyPadMinus);
    DEFINE_SPECIAL(VK_ADD, EKeyKeyPadPlus);
    DEFINE_SPECIAL(VK_LEFT, EKeyLeft);
    DEFINE_SPECIAL(VK_RIGHT, EKeyRight);
    DEFINE_SPECIAL(VK_UP, EKeyUp);
    DEFINE_SPECIAL(VK_DOWN, EKeyDown);
    DEFINE_SPECIAL(VK_INSERT, EKeyInsert);
    DEFINE_SPECIAL(VK_HOME, EKeyHome);
    DEFINE_SPECIAL(VK_END, EKeyEnd);
    DEFINE_SPECIAL(VK_PRIOR, EKeyPageUp);
    DEFINE_SPECIAL(VK_NEXT, EKeyPageDown);
    DEFINE_SPECIAL(VK_F1, EKeyF1);
    DEFINE_SPECIAL(VK_F2, EKeyF2);
    DEFINE_SPECIAL(VK_F3, EKeyF3);
    DEFINE_SPECIAL(VK_F4, EKeyF4);
    DEFINE_SPECIAL(VK_F5, EKeyF5);
    DEFINE_SPECIAL(VK_F6, EKeyF6);
    DEFINE_SPECIAL(VK_F7, EKeyF7);
    DEFINE_SPECIAL(VK_F8, EKeyF8);
    DEFINE_SPECIAL(VK_F9, EKeyF9);
    DEFINE_SPECIAL(VK_F10, EKeyF10);
    DEFINE_SPECIAL(VK_F11, EKeyF11);
    DEFINE_SPECIAL(VK_F12, EKeyF12);
    DEFINE_SPECIAL(VK_F13, EKeyF13);
    DEFINE_SPECIAL(VK_F14, EKeyF14);
    DEFINE_SPECIAL(VK_F15, EKeyF15);
    DEFINE_SPECIAL(VK_NUMLOCK, EKeyNumLock);
    DEFINE_SPECIAL(VK_CAPITAL, EKeyCapsLock);
    DEFINE_SPECIAL(VK_SCROLL, EKeyScrollLock);
    DEFINE_SPECIAL(VK_RMENU, EKeyLAlt);
    DEFINE_SPECIAL(VK_LMENU, EKeyRAlt);
    DEFINE_SPECIAL(VK_LCONTROL, EKeyLControl);
    DEFINE_SPECIAL(VK_RCONTROL, EKeyRControl);

    for (c='0'; c<='9'; c++)
        m_std[(unsigned char) c] = c;

    for (c='A'; c<='Z'; c++)
        m_std[(unsigned char) c] = c;
}

WGLDevice::~WGLDevice() {
    if (m_initialized)
        shutdown();
}


void WGLDevice::init(Device *other) {
    Device::init(other);

    Log(EDebug, "Initializing WGL device");

    WGLSession *session = static_cast<WGLSession *>(getSession());
    m_parent = static_cast<WGLDevice *>(other);

    __global_workaround = this;
    m_pf = -1;
    for (int i=0; i<2; i++) {
        /* Do this twice, the first pass creates a dummy device &
           rendering context in order to get hold of a proper pixel format
           (windows is soo bugged!), the second pass creates the actual
           window with the matching pixel format
        */
        int extra = 0;
        if (m_resizeAllowed)
            extra = WS_SIZEBOX | WS_MAXIMIZEBOX;

        m_hwnd = CreateWindow(
            session->m_wndClassName.c_str(),
            m_title.c_str(),
            m_fullscreen ? WS_POPUP : (WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX | extra),
            CW_USEDEFAULT,
            CW_USEDEFAULT,
            getSize().x, getSize().y,
            NULL, NULL,
            session->m_hinstance,
            NULL);
    }
    __global_workaround = NULL;

    if (!m_hwnd)
        Log(EError, "Unable to create the window");
    SetWindowLongPtr(m_hwnd, 0, (LONG_PTR) this);
    ShowWindow(m_hwnd, SW_HIDE);

    /* Switch to fullscreen */
    if (m_fullscreen) {
        DEVMODE dm;

        ZeroMemory(&dm, sizeof(dm));
        dm.dmSize = sizeof(dm);
        dm.dmBitsPerPel = GetDeviceCaps(m_hdc, BITSPIXEL);
        dm.dmPelsWidth = m_size.x;
        dm.dmPelsHeight = m_size.y;

        dm.dmFields = DM_PELSWIDTH | DM_PELSHEIGHT | DM_BITSPERPEL;
        if (!ChangeDisplaySettings(&dm, CDS_FULLSCREEN) == DISP_CHANGE_SUCCESSFUL)
            Log(EError, "Could not switch to fullscreen!");
    }

    if (m_fullscreen || m_center) {
        m_position.x = (GetSystemMetrics(SM_CXSCREEN) - m_size.x) / 2;
        m_position.y = (GetSystemMetrics(SM_CYSCREEN) - m_size.y) / 2;
    }

    SetWindowPos(m_hwnd, HWND_TOP, m_position.x, m_position.y, m_size.x, m_size.y, SWP_NOCOPYBITS);
    UpdateWindow(m_hwnd);

    m_mouse = Point2i(-1, -1);
    m_visible = false;
    m_initialized = true;
    m_grab = false;
    m_modifierState = 0;
    m_buttonState = 0;
    m_cursor = true;

    setTitle(m_title);
}

void WGLDevice::setTitle(const std::string &title) {
    Device::setTitle(title);

    if (m_initialized) {
        std::string finalTitle;

        if (m_showFPS && m_fps != 0) {
            finalTitle = formatString("%s - %i FPS", title.c_str(), m_fps);
        } else {
            finalTitle = title;
        }

        SetWindowText(m_hwnd, finalTitle.c_str());
    }
}

bool WGLDevice::translateKey(WPARAM wParam, LPARAM lParam, DeviceEvent &event) {
    BYTE allKeys[256];
    GetKeyboardState(allKeys);
    LRESULT ret = ToAscii((UINT) wParam, (UINT) ((lParam >> 16) & 0x00FF), allKeys, (WORD *) event.getKeyboardInterpreted(), 0);
    event.getKeyboardInterpreted()[ret] = '\0';
    event.setKeyboardSpecial(ENoSpecial);
    event.setKeyboardKey('\0');

    if (wParam < 256) {
        event.setKeyboardSpecial(m_special[wParam]);
        event.setKeyboardKey(m_std[wParam]);
    }

    switch (wParam) {
        case VK_RETURN: {
                /* Determine which return key it was */
                if (lParam & 0x01000000)
                    event.setKeyboardSpecial(EKeyKeyPadEnter);
                else
                    event.setKeyboardSpecial(EKeyReturn);
            }
            break;
        case VK_SHIFT: {
                /* Determine which shift key it was */
                bool lshift = (GetKeyState(VK_LSHIFT) & 0x8000) != 0;
                bool rshift = (GetKeyState(VK_RSHIFT) & 0x8000) != 0;

                if (m_leftShift != lshift) {
                    m_leftShift = lshift;
                    event.setKeyboardSpecial(EKeyLShift);
                } else if (m_rightShift != rshift) {
                    m_rightShift = rshift;
                    event.setKeyboardSpecial(EKeyRShift);
                }
            }
            break;
        case VK_CONTROL: {
                if (lParam & 0x01000000)
                    event.setKeyboardSpecial(EKeyRControl);
                else
                    event.setKeyboardSpecial(EKeyLControl);

                /* Filter out Alt-Gr key presses */
                DWORD msg_time = GetMessageTime();
                MSG next_msg;
                if (PeekMessage(&next_msg, NULL, 0, 0, PM_NOREMOVE)) {
                    if (next_msg.message == WM_KEYDOWN || next_msg.message == WM_SYSKEYDOWN
                        || next_msg.message == WM_KEYUP || next_msg.message == WM_SYSKEYUP) {
                        if (next_msg.wParam == VK_MENU && (next_msg.lParam & 0x01000000)
                            && next_msg.time == msg_time) {
                            event.setKeyboardSpecial(ENoSpecial);
                        }
                    }
                }
            }
            break;
        case VK_MENU:
                if (lParam & 0x01000000)
                    event.setKeyboardSpecial(EKeyRAlt);
                else
                    event.setKeyboardSpecial(EKeyLAlt);
            break;
    }

    if (event.getKeyboardSpecial() != ENoSpecial || event.getKeyboardKey() != '\0') {
        return true;
    } else if (ret == 1) {
        event.setKeyboardKey(event.getKeyboardInterpreted()[0]);
        return true;
    }
    return false;
}


bool WGLDevice::translateMouse(UINT uMsg, WPARAM wParam, DeviceEvent &event) {
    if (uMsg == WM_MOUSEWHEEL) {
        event.setMousePosition(m_mouse);
    }

    if (m_mouse != Point2i(-1, -1)) {
        event.setMouseRelative(event.getMousePosition() - m_mouse);
    } else {
        event.setMouseRelative(Vector2i(0, 0));
    }

    m_mouse = event.getMousePosition();

    switch (uMsg) {
        case WM_MOUSEMOVE:
            event.setType(EMouseMotionEvent);
            event.setMouseButton(ENoButton);
            if (m_grab)
                warpMouse(Point2i(getSize().x / 2,
                    getSize().y/2));
            break;
        case WM_LBUTTONDOWN:
            event.setType(EMouseButtonDownEvent);
            event.setMouseButton(ELeftButton);
            break;
        case WM_MBUTTONDOWN:
            event.setType(EMouseButtonDownEvent);
            event.setMouseButton(EMiddleButton);
            break;
        case WM_RBUTTONDOWN:
            event.setType(EMouseButtonDownEvent);
            event.setMouseButton(ERightButton);
            break;
        case WM_LBUTTONUP:
            event.setType(EMouseButtonUpEvent);
            event.setMouseButton(ELeftButton);
            break;
        case WM_MBUTTONUP:
            event.setType(EMouseButtonUpEvent);
            event.setMouseButton(EMiddleButton);
            break;
        case WM_RBUTTONUP:
            event.setType(EMouseButtonUpEvent);
            event.setMouseButton(ERightButton);
            break;
        case WM_MOUSEWHEEL: {
                int move = (short) HIWORD(wParam);
                event.setType(EMouseButtonDownEvent);

                if (move > 0)
                    event.setMouseButton(EWheelUpButton);
                else if (move < 0)
                    event.setMouseButton(EWheelDownButton);
                else
                    return false;
            }
            break;
        default:
            return false;
    }

    return true;
}

void WGLDevice::initPixelFormat(HWND hWnd) {
    m_hdc = GetDC(hWnd);

    if (m_hdc == 0)
        Log(EError, "Could not acquire a device context!");

    if (m_parent && m_pf == -1) {
        m_pf = GetPixelFormat(m_parent->m_hdc);
        if (m_pf == 0)
            Log(EError, "GetPixelFormat() failed: %s", lastErrorText().c_str());
        if (!DescribePixelFormat(m_parent->m_hdc, m_pf, sizeof(m_pfd), &m_pfd))
            Log(EError, "DescribePixelFormat() failed: %s", lastErrorText().c_str());
    }

    if (m_pf != -1) {
        if (!SetPixelFormat(m_hdc, m_pf, &m_pfd))
            Log(EError, "Could not set the pixel format");
        return;
    }

    /* First create a dummy pixel format and rendering context */
    ZeroMemory(&m_pfd, sizeof(m_pfd));
    m_pfd.nSize = sizeof(m_pfd);
    m_pfd.nVersion = 1;
    m_pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL;
    if (m_doubleBuffer)
        m_pfd.dwFlags |= PFD_DOUBLEBUFFER;
    m_pfd.iPixelType = PFD_TYPE_RGBA;
    m_pfd.cColorBits = GetDeviceCaps(m_hdc, BITSPIXEL);
    m_pfd.cBlueBits = m_blueBits;
    m_pfd.cRedBits = m_redBits;
    m_pfd.cGreenBits = m_greenBits;
    m_pfd.cDepthBits = m_depthBits;
    m_pfd.cAlphaBits = m_alphaBits;
    m_pfd.cStencilBits = m_stencilBits;

    m_pf = ChoosePixelFormat(m_hdc, &m_pfd);
    if (!m_pf)
        Log(EError, "Could not find a matching pixel format");

    if (!SetPixelFormat(m_hdc, m_pf, &m_pfd))
        Log(EError, "Could not set the pixel format");

    HGLRC context = wglCreateContext(m_hdc);
    if (context == 0)
        Log(EError, "Could not create a temporary WGL context");

    if (!wglMakeCurrent(m_hdc, context))
        Log(EError, "wglMakeCurrent failed!");

    wglChoosePixelFormatARBProc wglChoosePixelFormatARB =
        (wglChoosePixelFormatARBProc) wglGetProcAddress("wglChoosePixelFormatARB");

    if (wglChoosePixelFormatARB == NULL)
        Log(EError, "Could not get a function pointer to wglChoosePixelFormatARB!");

    /* Delete the rendering context */
    wglDeleteContext(context);

    int count=0, attribs[32];
    UINT numFormats;

    /* Lookup the actual pixel format */
    attribs[count++] = WGL_DRAW_TO_WINDOW_ARB; attribs[count++] = GL_TRUE;
    attribs[count++] = WGL_SUPPORT_OPENGL_ARB; attribs[count++] = GL_TRUE;
    attribs[count++] = WGL_ACCELERATION_ARB; attribs[count++] = WGL_FULL_ACCELERATION_ARB;
    attribs[count++] = WGL_COLOR_BITS_ARB; attribs[count++] = GetDeviceCaps(m_hdc, BITSPIXEL);
    attribs[count++] = WGL_RED_BITS_ARB; attribs[count++] = m_redBits;
    attribs[count++] = WGL_BLUE_BITS_ARB; attribs[count++] = m_blueBits;
    attribs[count++] = WGL_GREEN_BITS_ARB; attribs[count++] = m_greenBits;
    attribs[count++] = WGL_ALPHA_BITS_ARB; attribs[count++] = m_alphaBits;
    attribs[count++] = WGL_DEPTH_BITS_ARB; attribs[count++] = m_depthBits;
    attribs[count++] = WGL_STENCIL_BITS_ARB; attribs[count++] = m_stencilBits;
    attribs[count++] = WGL_DOUBLE_BUFFER_ARB; attribs[count++] = m_doubleBuffer ? GL_TRUE : GL_FALSE;

    if (m_fsaa > 1) {
        attribs[count++] = WGL_SAMPLE_BUFFERS_ARB; attribs[count++] = GL_TRUE;
        attribs[count++] = WGL_SAMPLES_ARB; attribs[count++] = m_fsaa;
    }
    attribs[count++] = 0; attribs[count++] = 0;

    float fAttributes[] = {0,0};
    bool valid = wglChoosePixelFormatARB(m_hdc,attribs,fAttributes,1,&m_pf,&numFormats);

    if (!(valid && numFormats >= 1))
        Log(EError, "Could not find a matching pixel format!");

    /* Self destruction .. */
    ReleaseDC(hWnd, m_hdc);
    DestroyWindow(hWnd);
}

void WGLDevice::setVisible(bool visible) {
    Assert(m_initialized);

    if (visible)
        ShowWindow(m_hwnd, SW_SHOWNORMAL);
    else
        ShowWindow(m_hwnd, SW_HIDE);
}

void WGLDevice::flip() {
    Assert(m_initialized);

    Device::flip();

    glFinish();

    if (m_doubleBuffer)
        SwapBuffers(m_hdc);
}

void WGLDevice::setPosition(const Point2i &position) {
    Assert(m_initialized);

    SetWindowPos(m_hwnd, HWND_TOP,
        getPosition().x, getPosition().y,
        0, 0, SWP_NOSIZE | SWP_NOZORDER);
}

void WGLDevice::makeCurrent(Renderer *pRenderer) {
    Assert(m_initialized);
    WGLRenderer *renderer = static_cast<WGLRenderer *>(pRenderer);
    if (!wglMakeCurrent(m_hdc, renderer->getWGLContext()))
        Log(EError, "wglMakeCurrent failed!");
}

void WGLDevice::showCursor(bool pShowCursor) {
    Assert(m_initialized);
    if (m_cursor != pShowCursor)
        ::ShowCursor(pShowCursor);
    m_cursor = pShowCursor;
}

void WGLDevice::warpMouse(const Point2i &position) {
    Assert(m_initialized);

    POINT pt;
    MSG  msg;
    pt.x = position.x;
    pt.y = position.y;
    m_mouse = position;
    ClientToScreen(m_hwnd, &pt);
    SetCursorPos(pt.x, pt.y);
    // Not so nice but does the job
    PeekMessage(&msg, m_hwnd, WM_MOUSEFIRST, WM_MOUSELAST, PM_REMOVE);
}

void WGLDevice::setGrab(bool grab) {
    Assert(m_initialized);
    m_grab = grab;
    showCursor(!grab);
}

void WGLDevice::shutdown() {
    Device::shutdown();
    Log(EDebug, "Shutting down WGL device");
    if (m_fullscreen) {
        ChangeDisplaySettings(NULL, 0);
        ShowWindow(m_hwnd, SW_HIDE);
    }
    DestroyWindow(m_hwnd);
    ReleaseDC(m_hwnd, m_hdc);
    m_initialized = false;
}

MTS_IMPLEMENT_CLASS(WGLDevice, false, Device)
MTS_NAMESPACE_END

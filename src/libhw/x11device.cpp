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

#include <mitsuba/hw/glxrenderer.h>
#include <mitsuba/hw/x11device.h>
#include <X11/XKBlib.h>

#define DEFINE_KEY(xsym, sym) m_keymap[xsym & 0xFF] = sym

MTS_NAMESPACE_BEGIN

X11Device::X11Device(X11Session *session)
 : Device(session), m_visinfo(NULL), m_visible(false) {
    m_title = "Mitsuba [x11]";
    for (int i=0; i<256; i++)
        m_keymap[i] = ENoSpecial;

    /* X11 <-> Mitsuba Keyboard mappings */
    DEFINE_KEY(XK_BackSpace, EKeyBackspace);
    DEFINE_KEY(XK_Tab, EKeyTab);
    DEFINE_KEY(XK_Clear, EKeyClear);
    DEFINE_KEY(XK_Return, EKeyReturn);
    DEFINE_KEY(XK_Linefeed, EKeyReturn);
    DEFINE_KEY(XK_Pause, EKeyPause);
    DEFINE_KEY(XK_Escape, EKeyEscape);
    DEFINE_KEY(XK_Delete, EKeyDelete);
    DEFINE_KEY(XK_KP_0, EKeyKeyPad0);
    DEFINE_KEY(XK_KP_1, EKeyKeyPad1);
    DEFINE_KEY(XK_KP_2, EKeyKeyPad2);
    DEFINE_KEY(XK_KP_3, EKeyKeyPad3);
    DEFINE_KEY(XK_KP_4, EKeyKeyPad4);
    DEFINE_KEY(XK_KP_5, EKeyKeyPad5);
    DEFINE_KEY(XK_KP_6, EKeyKeyPad6);
    DEFINE_KEY(XK_KP_7, EKeyKeyPad7);
    DEFINE_KEY(XK_KP_8, EKeyKeyPad8);
    DEFINE_KEY(XK_KP_9, EKeyKeyPad9);
    DEFINE_KEY(XK_KP_Insert, EKeyKeyPad0);
    DEFINE_KEY(XK_KP_End, EKeyKeyPad1);
    DEFINE_KEY(XK_KP_Down, EKeyKeyPad2);
    DEFINE_KEY(XK_KP_Page_Down, EKeyKeyPad3);
    DEFINE_KEY(XK_KP_Left, EKeyKeyPad4);
    DEFINE_KEY(XK_KP_Begin, EKeyKeyPad5);
    DEFINE_KEY(XK_KP_Right, EKeyKeyPad6);
    DEFINE_KEY(XK_KP_Home, EKeyKeyPad7);
    DEFINE_KEY(XK_KP_Up, EKeyKeyPad8);
    DEFINE_KEY(XK_KP_Page_Up, EKeyKeyPad9);
    DEFINE_KEY(XK_KP_Delete, EKeyKeyPadPeriod);
    DEFINE_KEY(XK_KP_Decimal, EKeyKeyPadPeriod);
    DEFINE_KEY(XK_KP_Divide, EKeyKeyPadDivide);
    DEFINE_KEY(XK_KP_Multiply, EKeyKeyPadMultiply);
    DEFINE_KEY(XK_KP_Subtract, EKeyKeyPadMinus);
    DEFINE_KEY(XK_KP_Add, EKeyKeyPadPlus);
    DEFINE_KEY(XK_KP_Enter, EKeyKeyPadEnter);
    DEFINE_KEY(XK_KP_Equal, EKeyKeyPadEquals);
    DEFINE_KEY(XK_Up, EKeyUp);
    DEFINE_KEY(XK_Down, EKeyDown);
    DEFINE_KEY(XK_Right, EKeyRight);
    DEFINE_KEY(XK_Left, EKeyLeft);
    DEFINE_KEY(XK_Insert, EKeyInsert);
    DEFINE_KEY(XK_Home, EKeyHome);
    DEFINE_KEY(XK_End, EKeyEnd);
    DEFINE_KEY(XK_Page_Up, EKeyPageUp);
    DEFINE_KEY(XK_Page_Down, EKeyPageDown);
    DEFINE_KEY(XK_F1, EKeyF1);
    DEFINE_KEY(XK_F2, EKeyF2);
    DEFINE_KEY(XK_F3, EKeyF3);
    DEFINE_KEY(XK_F4, EKeyF4);
    DEFINE_KEY(XK_F5, EKeyF5);
    DEFINE_KEY(XK_F6, EKeyF6);
    DEFINE_KEY(XK_F7, EKeyF7);
    DEFINE_KEY(XK_F8, EKeyF8);
    DEFINE_KEY(XK_F9, EKeyF9);
    DEFINE_KEY(XK_F10, EKeyF10);
    DEFINE_KEY(XK_F11, EKeyF11);
    DEFINE_KEY(XK_F12, EKeyF12);
    DEFINE_KEY(XK_F13, EKeyF13);
    DEFINE_KEY(XK_F14, EKeyF14);
    DEFINE_KEY(XK_F15, EKeyF15);
    DEFINE_KEY(XK_Num_Lock, EKeyNumLock);
    DEFINE_KEY(XK_Caps_Lock, EKeyCapsLock);
    DEFINE_KEY(XK_Scroll_Lock, EKeyScrollLock);
    DEFINE_KEY(XK_Shift_L, EKeyLShift);
    DEFINE_KEY(XK_Shift_R, EKeyRShift);
    DEFINE_KEY(XK_Meta_L, EKeyLMeta);
    DEFINE_KEY(XK_Meta_R, EKeyRMeta);
    DEFINE_KEY(XK_Alt_L, EKeyLAlt);
    DEFINE_KEY(XK_Alt_R, EKeyRAlt);
    DEFINE_KEY(XK_Control_L, EKeyLControl);
    DEFINE_KEY(XK_Control_R, EKeyRControl);
}

X11Device::~X11Device() {
    Log(EDebug, "Destroying X11 device");
    if (m_initialized)
        shutdown();
}

namespace algo {
    struct vsync_sort {
        /**
         * Approximate the vertical retrace speed of the mode info structure
         * (taken from bzflag)
         */
        inline int getVSync(XF86VidModeModeInfo *mode) {
            return ((int) (0.5f + (1000.0f * mode->dotclock) / (mode->htotal * mode->vtotal)));
        }

        bool operator() (XF86VidModeModeInfo *mode1, XF86VidModeModeInfo *mode2) {
            return getVSync(mode1) > getVSync(mode2);
        }
    };
}

XVisualInfo *X11Device::createVisual() {
    XVisualInfo visTempl;
    int visCount;

    X11Session *session = static_cast<X11Session *>(getSession());
    /* Search for visuals on the current screen matching the requrested color depth */
    visTempl.screen = session->m_screen;
    visTempl.depth = m_redBits + m_blueBits + m_greenBits + m_alphaBits;
    XVisualInfo *visinfo = XGetVisualInfo(session->m_display, VisualScreenMask | VisualDepthMask, &visTempl, &visCount);

    if (visinfo == NULL)
        Log(EError, "Could not find a matching visual!");

    return visinfo;
}

void X11Device::init(Device *other) {
    Device::init(other);

    Log(EDebug, "Initializing X11 device");
    X11Session *session = static_cast<X11Session *>(getSession());

    /* Find matching visuals */
    m_visinfo = createVisual();

    /* Fullscreen mode selection */
    std::vector<XF86VidModeModeInfo *> modeList;

    if (m_fullscreen) {
        if (!session->m_hasVidMode)
            Log(EError, "VidMode extension is required for fullscreen display");

        /* Retrieve a list of all video modes */
        int modeCount;
        XF86VidModeModeInfo** modes;
        XF86VidModeGetAllModeLines(session->m_display, session->m_screen, &modeCount, &modes);

        /* Save the current mode */
        m_previousMode = *modes[0];

        /* Find the best matching screen resolution */
        for (int i=0; i<modeCount; ++i) {
            if (modes[i]->hdisplay == m_size.x &&
                modes[i]->vdisplay == m_size.y)
                modeList.push_back(modes[i]);
        }
        /* Release the memory held by the resolution list */
        XFree(modes);

        std::sort(modeList.begin(), modeList.end(), algo::vsync_sort());

        if (modeList.size() == 0) {
            Log(EWarn, "No matching fullscreen resolution found, using windowed mode!");
            m_fullscreen = false;
        }
    }

    /* Create the window's attributes */
    XSetWindowAttributes x11attr;
    x11attr.background_pixel = x11attr.border_pixel =
        BlackPixel(session->m_display, session->m_screen);

    /* Create a colormap for the window */
    Colormap colormap = XCreateColormap(session->m_display, session->m_root, m_visinfo->visual, AllocNone);
    x11attr.colormap = colormap;

    if (m_fullscreen) {
        /* Switch to the best matching resolution */
        XF86VidModeSwitchToMode(session->m_display, session->m_screen, modeList[0]);
        XF86VidModeSetViewPort(session->m_display, session->m_screen, 0, 0);

        x11attr.override_redirect = True;

        m_window = XCreateWindow(session->m_display, session->m_root,
            0, 0, /* x,y position*/
            m_size.x, m_size.y,
            0, m_visinfo->depth, /* border width, color depth */
            InputOutput, m_visinfo->visual,
            CWBackPixel | CWBorderPixel | CWColormap |
            CWOverrideRedirect, &x11attr);

        XWarpPointer(session->m_display, None, m_window, 0, 0, 0, 0, 0, 0);
        XMapRaised(session->m_display, m_window);

        /* Grab the pointer & keyboard */
        XGrabKeyboard(session->m_display, m_window, True,
            GrabModeAsync, GrabModeAsync, CurrentTime);
        XGrabPointer(session->m_display, m_window, True,
            ButtonPressMask, GrabModeAsync, GrabModeAsync,
            m_window, None, CurrentTime);

        m_visible = true;
    } else {
        /* Center the window if needed */
        if (m_center) {
            m_position.x = (DisplayWidth(session->m_display, session->m_screen) - m_size.x) / 2;
            m_position.y = (DisplayHeight(session->m_display, session->m_screen) - m_size.y) / 2;
        }

        /* Create the X window */
        unsigned long mask = CWBackPixel | CWBorderPixel | CWColormap;
        m_window = XCreateWindow(session->m_display, session->m_root,
                m_position.x, m_position.y, m_size.x, m_size.y, 0, m_visinfo->depth,
                InputOutput, m_visinfo->visual, mask, &x11attr);

        if (!m_window)
            Log(EError, "Could not create the window");

        /* Make the window non-resizable */
        XSizeHints *hints = XAllocSizeHints();
        hints->width = m_size.x;
        hints->height = m_size.y;

        if (m_resizeAllowed) {
            hints->min_width = hints->min_height = 10;
            hints->max_width = hints->max_height = INT_MAX;
        } else {
            hints->min_width = hints->max_width = m_size.x;
            hints->min_height = hints->max_height = m_size.y;
        }

        hints->x = m_position.x; hints->y = m_position.y;
        hints->flags = PMaxSize | PMinSize | USSize | USPosition;
        XSetNormalHints(session->m_display, m_window, hints);
        XFree(hints);

        /* Set input hints */
        XWMHints *wmHints = XAllocWMHints();
        wmHints->input = True;
        wmHints->flags = InputHint;
        XSetWMHints(session->m_display, m_window, wmHints);
        XFree(wmHints);

        /* Make the window closeable */
        m_deleteWindow = XInternAtom(session->m_display, "WM_DELETE_WINDOW", False);
        XSetWMProtocols(session->m_display, m_window, &m_deleteWindow, 1);
    }

    /* Mark events types in which we are interested */
    XSelectInput(session->m_display, m_window, FocusChangeMask | KeyPressMask | KeyReleaseMask
            | ButtonPressMask | ButtonReleaseMask | PointerMotionMask | StructureNotifyMask);

    /* Initialize member variables */
    m_cursor = None;
    m_mouse = Point2i(-1, -1);
    m_modifierState = 0;
    m_buttonState = 0;
    m_grab = false;

    m_initialized = true;
    setTitle(m_title);
}

void X11Device::setTitle(const std::string &title) {
    X11Session *session = static_cast<X11Session *>(getSession());

    Device::setTitle(title);
    if (m_initialized) {
        std::string finalTitle;

        if (m_showFPS && m_fps != 0) {
            finalTitle = formatString("%s - %i FPS", title.c_str(), m_fps);
        } else {
            finalTitle = title;
        }

        XStoreName(session->m_display, m_window, finalTitle.c_str());
        XFlush(session->m_display);
    }
}

void X11Device::setPosition(const Point2i &position) {
    Assert(m_initialized);

    X11Session *session = static_cast<X11Session *>(getSession());

    Device::setPosition(position);
    if (m_initialized && !m_fullscreen) {
        XMoveWindow(session->m_display, m_window, m_position.x, m_position.y);
        XFlush(session->m_display);
    }
}

void X11Device::setVisible(bool visible) {
    Assert(m_initialized);

    if (visible && !m_visible) {
        X11Session *session = static_cast<X11Session *>(getSession());
        XMapRaised(session->m_display, m_window);
        XSync(session->m_display, 0);
        m_visible = true;
    } else if (!visible && m_visible) {
        X11Session *session = static_cast<X11Session *>(getSession());
        XUnmapWindow(session->m_display, m_window);
        XSync(session->m_display, 0);
        m_visible = false;
    }
}

void X11Device::warpMouse(const Point2i &position) {
    Assert(m_initialized);

    X11Session *session = static_cast<X11Session *>(getSession());
    XEvent event;
    XWarpPointer(session->m_display, None, m_window, 0, 0, 0, 0, position.x, position.y);
    XSync(session->m_display, False);
    /* Remove the caused event from the queue */
    XCheckTypedWindowEvent(session->m_display, m_window, MotionNotify, &event);
    m_mouse = Point2i(position.x, position.y);
}

void X11Device::showCursor(bool enabled) {
    X11Session *session = static_cast<X11Session *>(getSession());
    if (enabled) {
        if (m_cursor != None) {
            XFreeCursor(session->m_display, m_cursor);
            m_cursor = None;
        }
        XUndefineCursor(session->m_display, m_window);
        XSync(session->m_display, False);
    } else {
        if (m_cursor == None) {
            /* Create a transparent cursor */
            char bm[] = { 0, 0, 0, 0, 0, 0, 0, 0 };
            Pixmap pix = XCreateBitmapFromData(session->m_display, m_window, bm, 8, 8);
            XColor black;

            memset(&black, 0, sizeof(XColor));
            black.flags = DoRed | DoGreen | DoBlue;
            m_cursor = XCreatePixmapCursor(session->m_display, pix, pix, &black, &black, 0, 0);
            XFreePixmap(session->m_display, pix);
        }
        XDefineCursor(session->m_display, m_window, m_cursor);
        XSync(session->m_display, False);
    }
}

void X11Device::shutdown() {
    X11Session *session = static_cast<X11Session *>(getSession());
    Log(EDebug, "Shutting down X11 device");
    Device::shutdown();
    setVisible(false);
    XDestroyWindow(session->m_display, m_window);
    XFree(m_visinfo);

    if (m_fullscreen) {
        /* Switch back to the previous screen resolution */
        XF86VidModeSwitchToMode(session->m_display, session->m_screen, &m_previousMode);
        XF86VidModeSetViewPort(session->m_display, session->m_screen, 0, 0);
    }

    /* In case auto_repeat was left on */
    XKeyboardState xkbs;
    XAutoRepeatOn(session->m_display);
    XGetKeyboardControl(session->m_display, &xkbs);

    if (!xkbs.global_auto_repeat)
        Log(EWarn, "Unable to restore the keyboard auto-repeat flag");

    m_initialized = false;
}

void X11Device::flip() {
    Assert(m_initialized);
}

void X11Device::processEvent(const XEvent &event) {
    DeviceEvent deviceEvent(ENoEvent);
    X11Session *session = static_cast<X11Session *>(getSession());

    if (m_callbacks.size() == 0)
        return;

    switch (event.type) {
    case ClientMessage:
        /* The window close button pressed */
        if ((event.xclient.format == 32) && ((unsigned) event.xclient.data.l[0] == m_deleteWindow)) {
            deviceEvent.setType(EQuitEvent);
        }
        break;
    case FocusIn:
        /* Deactivate auto-repeat */
        XAutoRepeatOff(session->m_display);
        deviceEvent.setType(EGainFocusEvent);
        m_modifierState = 0;
        m_buttonState = 0;
        break;
    case FocusOut:
        /* Reactivate auto-repeat */
        XAutoRepeatOn(session->m_display);
        deviceEvent.setType(ELoseFocusEvent);
        m_modifierState = 0;
        m_buttonState = 0;
        break;
    case ButtonPress:
        deviceEvent.setType(EMouseButtonDownEvent);
        translateMouse(event, deviceEvent);
        m_buttonState |= deviceEvent.getMouseButton();
        break;
    case ButtonRelease:
        deviceEvent.setType(EMouseButtonUpEvent);
        translateMouse(event, deviceEvent);
        m_buttonState &= ~deviceEvent.getMouseButton();
        break;
    case MotionNotify: {
        deviceEvent.setType(m_buttonState == 0 ? EMouseMotionEvent : EMouseDragEvent);
        translateMouse(event, deviceEvent);
        deviceEvent.setMouseButton(m_buttonState);
        if (m_grab)
            warpMouse(Point2i(getSize().x / 2, getSize().y/2));
        int xpos = deviceEvent.getMousePosition().x;
        int ypos = deviceEvent.getMousePosition().y;
        if (xpos > m_size.x || xpos < 0 || ypos > m_size.y || ypos < 0)
            return;
        }
        break;
    case KeyPress:
        if (translateKeyboard(event, deviceEvent)) {
            deviceEvent.setType(EKeyDownEvent);

            int special = deviceEvent.getKeyboardSpecial();

            /* Update the current modifier state */
            if (special == EKeyLShift || special == EKeyRShift) {
                m_modifierState |= EShiftModifier;
            } else if (special == EKeyLAlt || special == EKeyRAlt) {
                m_modifierState |= EAltModifier;
            } else if (special == EKeyLControl || special == EKeyRControl) {
                m_modifierState |= EControlModifier;
            }

            deviceEvent.setKeyboardModifiers(m_modifierState);
        }
        break;
    case KeyRelease:
        if (translateKeyboard(event, deviceEvent)) {
            deviceEvent.setType(EKeyUpEvent);

            int special = deviceEvent.getKeyboardSpecial();

            /* Update the current modifier state */
            if (special == EKeyLShift || special == EKeyRShift) {
                m_modifierState = m_modifierState & (~EShiftModifier);
            } else if (special == EKeyLAlt || special == EKeyRAlt) {
                m_modifierState = m_modifierState & (~EAltModifier);
            } else if (special == EKeyLControl || special == EKeyRControl) {
                m_modifierState = m_modifierState & (~EControlModifier);
            }

            deviceEvent.setKeyboardModifiers(m_modifierState);
        }
        break;
    case MapNotify:
    case UnmapNotify:
        m_modifierState = 0;
        break;
    case ConfigureNotify: {
            Vector2i size(event.xconfigure.width, event.xconfigure.height);
            if (m_size != size) {
                m_size = size;
                deviceEvent.setType(EResizeEvent);
            }
        }
        break;
    case ReparentNotify:
    case Expose:
        break;
    default:
        Log(EWarn, "Unknown event %i received", event.type);
    }
    if (deviceEvent.getType() != ENoEvent)
        fireDeviceEvent(deviceEvent);
}


void X11Device::translateMouse(const XEvent &xEvent, DeviceEvent &event) {
    event.setMousePosition(Point2i(xEvent.xbutton.x, xEvent.xbutton.y));

    /* Calculate relative coordinates */
    if (m_mouse.x != -1 && m_mouse.y != -1) {
        event.setMouseRelative(Vector2i(event.getMousePosition() - m_mouse));
    } else {
        event.setMouseRelative(Vector2i(0,0));
    }

    m_mouse = event.getMousePosition();

    if (xEvent.xbutton.button == Button1)
        event.setMouseButton(ELeftButton);
    else if (xEvent.xbutton.button == Button2)
        event.setMouseButton(EMiddleButton);
    else if (xEvent.xbutton.button == Button3)
        event.setMouseButton(ERightButton);
    else if (xEvent.xbutton.button == Button4)
        event.setMouseButton(EWheelUpButton);
    else if (xEvent.xbutton.button == Button5)
        event.setMouseButton(EWheelDownButton);
    else
        event.setMouseButton(ENoButton);
}

bool X11Device::translateKeyboard(const XEvent &xEvent, DeviceEvent &event) {
    X11Session *session = static_cast<X11Session *>(getSession());

    KeySym sym = XkbKeycodeToKeysym(session->m_display,
        xEvent.xkey.keycode, 0, 0);

    /* Default: 0 */
    event.setKeyboardKey(0);
    event.setKeyboardSpecial(0);

    /* Retrieve an 'interpreted' version of the key event */
    int noc = XLookupString(const_cast<XKeyEvent *>(&xEvent.xkey), event.getKeyboardInterpreted(),
        15, NULL, NULL);
    event.getKeyboardInterpreted()[noc] = 0;

    if (sym != 0) {
        int type = sym >> 8;
        if (type >= 0x00 && type <= 0x0D) {
            /* Alphanumeric key */

            char key = sym & 0xFF;
            if (key >= 'A' && key <= 'Z') {
                /* Convert to lower-case */
                key += ('a' - 'A');
            }

            event.setKeyboardKey(key);
        } else if (type == 0xFE) {
            /* Ignore */
            return false;
        } else if (type == 0xFF) {
            /* Lookup special key */
            event.setKeyboardSpecial(m_keymap[sym & 0xFF]);
            if (event.getKeyboardSpecial()  == ENoSpecial)
                return false;
        } else {
            Log(EWarn, "Unknown X11 keysym: 0x%x", (int) sym);
            return false;
        }
        return true;
    }
    return false;
}

void X11Device::setGrab(bool grab) {
    Assert(m_initialized);

    m_grab = grab;
    showCursor(!grab);
}

void X11Device::makeCurrent(Renderer *renderer) {
    Assert(m_initialized);

    X11Session *session = static_cast<X11Session *>(getSession());
    GLXRenderer *glxRenderer = static_cast<GLXRenderer *>(renderer);
    int retval;
    if (glxRenderer == NULL)
        retval = glXMakeCurrent(session->m_display, None, NULL);
    else
        retval = glXMakeCurrent(session->m_display, m_window, glxRenderer->getGLXContext());

    if (retval != True)
        Log(EError, "Error in glXMakeCurrent - unable to activate the rendering context");
}

MTS_IMPLEMENT_CLASS(X11Device, false, Device)
MTS_NAMESPACE_END

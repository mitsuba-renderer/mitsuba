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

#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/hw/nsgldevice.h>
#include <mitsuba/hw/nsglkeys.h>
#include <mitsuba/hw/nsglrenderer.h>

using namespace mitsuba;

#define DEFINE_KEY(xsym, sym) m_keymap[xsym] = Device::sym

@implementation CustomView
- (id) initWithFrame: (NSRect) frame {
    self = [super initWithFrame: frame];
    if (self) {
        m_modifiers = 0;
        m_firstMouseMotion = true;
        m_focused = false;
        m_mouseInWindow = false;
        m_ignoreNextMouseEvent = false;

        /* Cocoa <-> Mitsuba keyboard mappings */
        for (uint i=0; i<0xFF; i++)
            m_keymap[i] = Device::ENoSpecial;
        DEFINE_KEY(QZ_BACKSPACE, EKeyBackspace);
        DEFINE_KEY(QZ_TAB, EKeyTab);
        DEFINE_KEY(QZ_RETURN, EKeyReturn);
        DEFINE_KEY(QZ_PAUSE, EKeyPause);
        DEFINE_KEY(QZ_ESCAPE, EKeyEscape);
        DEFINE_KEY(QZ_DELETE, EKeyDelete);
        DEFINE_KEY(QZ_KP0, EKeyKeyPad0);
        DEFINE_KEY(QZ_KP1, EKeyKeyPad1);
        DEFINE_KEY(QZ_KP2, EKeyKeyPad2);
        DEFINE_KEY(QZ_KP3, EKeyKeyPad3);
        DEFINE_KEY(QZ_KP4, EKeyKeyPad4);
        DEFINE_KEY(QZ_KP5, EKeyKeyPad5);
        DEFINE_KEY(QZ_KP6, EKeyKeyPad6);
        DEFINE_KEY(QZ_KP7, EKeyKeyPad7);
        DEFINE_KEY(QZ_KP8, EKeyKeyPad8);
        DEFINE_KEY(QZ_KP9, EKeyKeyPad9);
        DEFINE_KEY(QZ_KP_PERIOD, EKeyKeyPadPeriod);
        DEFINE_KEY(QZ_KP_DIVIDE, EKeyKeyPadDivide);
        DEFINE_KEY(QZ_KP_MULTIPLY, EKeyKeyPadMultiply);
        DEFINE_KEY(QZ_KP_MINUS, EKeyKeyPadMinus);
        DEFINE_KEY(QZ_KP_PLUS, EKeyKeyPadPlus);
        DEFINE_KEY(QZ_KP_ENTER, EKeyKeyPadEnter);
        DEFINE_KEY(QZ_KP_EQUALS, EKeyKeyPadEquals);
        DEFINE_KEY(QZ_UP, EKeyUp);
        DEFINE_KEY(QZ_DOWN, EKeyDown);
        DEFINE_KEY(QZ_RIGHT, EKeyRight);
        DEFINE_KEY(QZ_LEFT, EKeyLeft);
        DEFINE_KEY(QZ_INSERT, EKeyInsert);
        DEFINE_KEY(QZ_HOME, EKeyHome);
        DEFINE_KEY(QZ_END, EKeyEnd);
        DEFINE_KEY(QZ_PAGEUP, EKeyPageUp);
        DEFINE_KEY(QZ_PAGEDOWN, EKeyPageDown);
        DEFINE_KEY(QZ_F1, EKeyF1);
        DEFINE_KEY(QZ_F2, EKeyF2);
        DEFINE_KEY(QZ_F3, EKeyF3);
        DEFINE_KEY(QZ_F4, EKeyF4);
        DEFINE_KEY(QZ_F5, EKeyF5);
        DEFINE_KEY(QZ_F6, EKeyF6);
        DEFINE_KEY(QZ_F7, EKeyF7);
        DEFINE_KEY(QZ_F8, EKeyF8);
        DEFINE_KEY(QZ_F9, EKeyF9);
        DEFINE_KEY(QZ_F10, EKeyF10);
        DEFINE_KEY(QZ_F11, EKeyF11);
        DEFINE_KEY(QZ_F12, EKeyF12);
        DEFINE_KEY(QZ_NUMLOCK, EKeyNumLock);
        DEFINE_KEY(QZ_CAPSLOCK, EKeyCapsLock);
        DEFINE_KEY(QZ_SCROLLOCK, EKeyScrollLock);
        DEFINE_KEY(QZ_LSHIFT, EKeyLShift);
        DEFINE_KEY(QZ_RSHIFT, EKeyRShift);
        DEFINE_KEY(QZ_LMETA, EKeyLMeta);
        DEFINE_KEY(QZ_RMETA, EKeyRMeta);
        DEFINE_KEY(QZ_LALT, EKeyLAlt);
        DEFINE_KEY(QZ_RALT, EKeyRAlt);
        DEFINE_KEY(QZ_LCTRL, EKeyLControl);
        DEFINE_KEY(QZ_RCTRL, EKeyRControl);
    }
    return self;
}

- (void) setDevice: (NSGLDevice *) device {
    m_device = device;
}

- (void) ignoreNextMouseEvent {
    m_ignoreNextMouseEvent = true;
}

- (void) ignoreFirstMouseMotion {
    m_firstMouseMotion = false;
}

- (BOOL) focused {
    return m_focused;
}

- (BOOL) acceptsFirstResponder {
    return YES;
}

- (uint) extractModifiers: (uint) modifiers {
    uint result = 0;

    if (modifiers & NSAlphaShiftKeyMask)
        result |= 0x10;
    if (modifiers & NSShiftKeyMask)
        result |= Device::EShiftModifier;
    if (modifiers & NSControlKeyMask)
        result |= Device::EControlModifier;
    if (modifiers & NSAlternateKeyMask)
        result |= Device::EAltModifier;
    if (modifiers & NSCommandKeyMask)
        result |= Device::EMetaModifier;
    return result;
}

- (void) windowDidBecomeKey: (NSNotification *) notification {
    bool cursorInWindow = m_device->isMouseInWindow();
    m_focused = true;
    m_device->pushEvent(DeviceEvent(Device::EGainFocusEvent));
    m_buttonMask = 0;
    if (!m_device->getCursor() && cursorInWindow) {
        m_mouseInWindow = true;
        [NSCursor hide];
    }
}

- (void) windowDidResignKey: (NSNotification *) notification {
    m_focused = false;
    m_device->pushEvent(DeviceEvent(Device::ELoseFocusEvent));
    m_buttonMask = 0;
    if (!m_device->getCursor()) {
        [NSCursor unhide];
    }
}

- (void) handleEvent: (NSEvent *) event {
    DeviceEvent deviceEvent(Device::ENoEvent);
    uint type = [event type];

    switch (type) {
        case NSFlagsChanged: {
                const uint list1[] = { Device::EShiftModifier, Device::EControlModifier, Device::EAltModifier, Device::EMetaModifier, 0x10};
                const uint list2[] = { Device::EKeyLShift, Device::EKeyLControl, Device::EKeyLAlt, Device::EKeyLMeta, Device::EKeyCapsLock};
                uint newModifiers = [self extractModifiers: [event modifierFlags]];
                uint difference = m_modifiers ^ newModifiers;
                m_modifiers = newModifiers;

                for (uint i=0; i<5; i++) {
                    if ((difference & list1[i]) != 0) {
                        deviceEvent.setType((m_modifiers & list1[i]) == 0 ? Device::EKeyUpEvent : Device::EKeyDownEvent);
                        deviceEvent.getKeyboardInterpreted()[0] = '\0';
                        deviceEvent.setKeyboardSpecial(list2[i]);
                        deviceEvent.setKeyboardKey('\0');
                        deviceEvent.setKeyboardModifiers(m_modifiers);
                        m_device->pushEvent(deviceEvent);
                    }
                }
            }
            return;
        case NSLeftMouseUp:
        case NSRightMouseUp:
        case NSOtherMouseUp:
            deviceEvent.setType(Device::EMouseButtonUpEvent);
        case NSLeftMouseDown:
        case NSRightMouseDown:
        case NSOtherMouseDown: {
                NSPoint location = [((NSWindow *) m_device->getWindow()) mouseLocationOutsideOfEventStream];
                if (deviceEvent.getType() == Device::ENoEvent)
                    deviceEvent.setType(Device::EMouseButtonDownEvent);
                deviceEvent.setMousePosition(Point2i((int) location.x, m_device->getSize().y - (int) location.y));
                deviceEvent.setMouseRelative(Vector2i(0, 0));
                uint buttonNumber = [event buttonNumber];
                uint buttonMask = 0;
                if (buttonNumber == 0)
                    buttonMask = Device::ELeftButton;
                else if (buttonNumber == 1)
                    buttonMask = Device::ERightButton;
                else if (buttonNumber == 2)
                    buttonMask = Device::EMiddleButton;
                else
                    return;
                if (deviceEvent.getType() == Device::EMouseButtonDownEvent)
                    m_buttonMask |= buttonMask;
                else
                    m_buttonMask &= ~buttonMask;
                deviceEvent.setMouseButton(buttonMask);
            }
            break;
        case NSScrollWheel: {
                NSPoint location = [((NSWindow *) m_device->getWindow()) mouseLocationOutsideOfEventStream];
                float deltaX = [event deltaX], deltaY = [event deltaY];
                if (deltaX > 0 || deltaY > 0)
                    deviceEvent.setMouseButton(Device::EWheelUpButton);
                else
                    deviceEvent.setMouseButton(Device::EWheelDownButton);
                deviceEvent.setMousePosition(Point2i((int) location.x, m_device->getSize().y - (int) location.y));
                deviceEvent.setMouseRelative(Vector2i(0, 0));
                deviceEvent.setType(Device::EMouseButtonDownEvent);
                m_device->pushEvent(deviceEvent);
                deviceEvent.setType(Device::EMouseButtonUpEvent);
                m_device->pushEvent(deviceEvent);
                return;
            }
            break;
        case NSLeftMouseDragged:
        case NSRightMouseDragged:
        case NSOtherMouseDragged:
        case NSMouseMoved: {
                if (m_ignoreNextMouseEvent) {
                    m_ignoreNextMouseEvent = false;
                    return;
                }
                NSPoint location = [((NSWindow *) m_device->getWindow()) mouseLocationOutsideOfEventStream];
                Point2i absolute((int) location.x, m_device->getSize().y - (int) location.y);

                bool cursorInWindow = m_device->isMouseInWindow();
                if (m_buttonMask == 0)
                    deviceEvent.setType(Device::EMouseMotionEvent);
                else
                    deviceEvent.setType(Device::EMouseDragEvent);
                deviceEvent.setMousePosition(absolute);
                deviceEvent.setMouseButton(m_buttonMask);
                if (m_firstMouseMotion) {
                    deviceEvent.setMouseRelative(Vector2i(absolute.x, absolute.y));
                    m_firstMouseMotion = false;
                } else {
                    deviceEvent.setMouseRelative(Vector2i((int) [event deltaX], (int) [event deltaY]));
                }

                if (!cursorInWindow && m_focused && !m_device->getCursor()) {
                    if (m_mouseInWindow) {
                        m_mouseInWindow = false;
                        [NSCursor unhide];
                    }
                } else if (cursorInWindow && m_focused && !m_device->getCursor()) {
                    if (!m_mouseInWindow) {
                        m_mouseInWindow = true;
                        [NSCursor hide];
                    }
                }

                if (cursorInWindow)
                    m_mouseInWindow = true;

                if (absolute.x > m_device->getSize().x || absolute.x < 0
                    || absolute.y > m_device->getSize().y || absolute.y< 0)
                    return;
            }
            break;
        case NSKeyUp:
            deviceEvent.setType(Device::EKeyUpEvent);
        case NSKeyDown: {
            if (deviceEvent.getType() == Device::ENoEvent)
                deviceEvent.setType(Device::EKeyDownEvent);
            if ([event isARepeat])
                return;

            NSString *characters = [event characters];
            uint count = [characters length];

            if (count == 0 || count == 1) {
                unsigned char scanCode = [event keyCode];
                uint special = m_keymap[scanCode];

                strncpy(deviceEvent.getKeyboardInterpreted(), [characters UTF8String], 15);

                if (special != 0) {
                    deviceEvent.setKeyboardKey('\0');
                    deviceEvent.setKeyboardSpecial(m_keymap[scanCode]);
                    deviceEvent.setKeyboardModifiers([self extractModifiers: [event modifierFlags]]);
                } else {
                    if (count > 0)
                        deviceEvent.setKeyboardKey([[event charactersIgnoringModifiers] UTF8String][0]);
                    else
                        deviceEvent.setKeyboardKey('\0');
                    deviceEvent.setKeyboardSpecial(Device::ENoSpecial);
                    deviceEvent.setKeyboardModifiers(0);
                }
            }
        }
    }

    if (deviceEvent.getType() != Device::ENoEvent)
        m_device->pushEvent(deviceEvent);
}

- (void) mouseDown: (NSEvent *) event {
    [self handleEvent: event];
}

- (void) mouseUp: (NSEvent *) event {
    [self handleEvent: event];
}

- (void) mouseMoved: (NSEvent *) event {
    [self handleEvent: event];
}

- (void) mouseDragged: (NSEvent *) event {
    [self handleEvent: event];
}

- (void) rightMouseDown: (NSEvent *) event {
    [self handleEvent: event];
}

- (void) rightMouseUp: (NSEvent *) event {
    [self handleEvent: event];
}

- (void) rightMouseMoved: (NSEvent *) event {
    [self handleEvent: event];
}

- (void) rightMouseDragged: (NSEvent *) event {
    [self handleEvent: event];
}

- (void) otherMouseDown: (NSEvent *) event {
    [self handleEvent: event];
}

- (void) otherMouseUp: (NSEvent *) event {
    [self handleEvent: event];
}

- (void) otherMouseMoved: (NSEvent *) event {
    [self handleEvent: event];
}

- (void) otherMouseDragged: (NSEvent *) event {
    [self handleEvent: event];
}

- (void) scrollWheel: (NSEvent *) event {
    [self handleEvent: event];
}

- (void) keyDown: (NSEvent *) event {
    [self handleEvent: event];
}

- (void) keyUp: (NSEvent *) event {
    [self handleEvent: event];
}

- (void) flagsChanged: (NSEvent *) event {
    [self handleEvent: event];
}

- (BOOL) windowShouldClose: (id) sender {
    m_device->pushEvent(DeviceEvent(Device::EQuitEvent));
    return NO;
}

@end

MTS_NAMESPACE_BEGIN

NSGLDevice::NSGLDevice(NSGLSession *session)
 : Device(session), m_visible(false), m_cursor(true) {
    m_title = "Mitsuba [nsgl]";
}

NSGLDevice::~NSGLDevice() {
    Log(EDebug, "Destroying NSGL device");
    if (m_initialized)
        shutdown();
}

void NSGLDevice::init(Device *other) {
    Device::init(other);

    __mts_init_cocoa();

    NSOpenGLPixelFormatAttribute attribs[32];
    uint i=0;

    Device::init();

    m_currentContext = NULL;

    Log(EDebug, "Initializing NSGL device");

    NSRect contentRect = NSMakeRect(m_position.x, m_position.y,
        m_size.x, m_size.y);

    /* Protect the event queue */
    m_mutex = new Mutex();

    /* Create the device window */
    m_window = [[NSWindow alloc] initWithContentRect: contentRect
        styleMask: NSTitledWindowMask | NSClosableWindowMask | NSMiniaturizableWindowMask
        backing: NSBackingStoreBuffered defer: NO];
    if (m_window == nil)
        Log(EError, "Could not create window");

    if (m_center)
        [m_window center];

    /* Create a sub-view as drawing destination and in order to catch events */
    m_view = [[CustomView alloc] initWithFrame: contentRect];
    if (m_view == nil)
        Log(EError, "Could not create view");

    [m_view setDevice: this];
    [[m_window contentView] addSubview: m_view];
    [m_window setDelegate: m_view];
    [m_window setAcceptsMouseMovedEvents: YES];

    /* Pixel format setup */
    AssertEx(m_redBits == m_blueBits || m_redBits == m_greenBits, "NSGL does not support individual color depths");

    attribs[i++] = NSOpenGLPFAColorSize;
    attribs[i++] = (NSOpenGLPixelFormatAttribute) m_redBits;

    attribs[i++] = NSOpenGLPFAAlphaSize;
    attribs[i++] = (NSOpenGLPixelFormatAttribute) m_alphaBits;

    attribs[i++] = NSOpenGLPFADepthSize;
    attribs[i++] = (NSOpenGLPixelFormatAttribute) m_depthBits;

    attribs[i++] = NSOpenGLPFAStencilSize;
    attribs[i++] = (NSOpenGLPixelFormatAttribute) m_stencilBits;

    if (m_doubleBuffer) {
        attribs[i++] = NSOpenGLPFADoubleBuffer;
    }

    if (m_fsaa > 1) {
        attribs[i++] = NSOpenGLPFASampleBuffers; attribs[i++] = (NSOpenGLPixelFormatAttribute) 1;
        attribs[i++] = NSOpenGLPFASamples; attribs[i++] = (NSOpenGLPixelFormatAttribute) m_fsaa;
    }

    attribs[i++] = NSOpenGLPFANoRecovery; // Never switch renderers
    attribs[i++] = NSOpenGLPFAWindow;
    attribs[i++] = (NSOpenGLPixelFormatAttribute) 0;
    attribs[i++] = (NSOpenGLPixelFormatAttribute) 0;

    m_fmt = [[NSOpenGLPixelFormat alloc] initWithAttributes: attribs];

    if (m_fmt == nil)
        Log(EError, "Could not create OpenGL pixel format!");

    m_initialized = true;
    setTitle(m_title);
}

void NSGLDevice::shutdown() {
    Device::shutdown();
    Log(EDebug, "Shutting down NSGL device");
    setVisible(false);
    [m_fmt release];
    [m_view release];
    [m_window release];
    m_initialized = false;
}


void NSGLDevice::setTitle(const std::string &title) {
    std::string finalTitle;

    Assert(m_initialized);
    if (m_showFPS && m_fps != 0) {
        finalTitle = formatString("%s - %i FPS", title.c_str(), m_fps);
    } else {
        finalTitle = title;
    }

    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];

    NSString *string = [NSString stringWithUTF8String: finalTitle.c_str()];
    [m_window setTitle: string];
    Device::setTitle(title);

    [pool release];
}

void NSGLDevice::setPosition(const Point2i &position) {
    Assert(m_initialized);

    NSPoint point = NSMakePoint(position.x,
        CGDisplayPixelsHigh(kCGDirectMainDisplay) - position.y);

    [m_window setFrameTopLeftPoint: point];
    Device::setPosition(position);
}

void NSGLDevice::setVisible(bool visible) {
    Assert(m_initialized);

    if (visible && !m_visible) {
        [m_window makeKeyAndOrderFront: nil];
        m_visible = true;
    } else if (!visible && m_visible) {
        [m_window orderOut: nil];
        m_visible = false;
    }
}

void NSGLDevice::warpMouse(const Point2i &position) {
    Assert(m_initialized);

    NSRect frame = [m_window frame];

    CGPoint point;
    point.x = position.x + frame.origin.x;
    point.y = position.y + frame.origin.y;

    /* Avoids cursor freezing */
    [m_view ignoreNextMouseEvent];
    CGEventSourceRef evsrc =
        CGEventSourceCreate(kCGEventSourceStateCombinedSessionState);
    CGEventSourceSetLocalEventsSuppressionInterval(evsrc, 0.0);
    CGWarpMouseCursorPosition(point);
    CFRelease(evsrc);
}

void NSGLDevice::setGrab(bool grab) {
    Assert(m_initialized);

    if (grab) {
        warpMouse(Point2i(getSize().x / 2, getSize().y/2));
        CGAssociateMouseAndMouseCursorPosition(false);
    } else {
        CGAssociateMouseAndMouseCursorPosition(true);
    }
    [m_view ignoreFirstMouseMotion];
    showCursor(!grab);
}

void NSGLDevice::showCursor(bool enabled) {
    Assert(m_initialized);

    if (!m_cursor && enabled) {
        [NSCursor unhide];
        m_cursor = true;
    } else if (m_cursor && !enabled) {
        m_cursor = false;
        if (isMouseInWindow() && [m_view focused]) {
            [NSCursor hide];
        }
    }
}

void NSGLDevice::flip() {
    Assert(m_initialized);

    Device::flip();

    if (m_doubleBuffer) {
        Assert(m_currentContext != NULL);
        [m_currentContext flushBuffer];
    }
}

void NSGLDevice::pushEvent(const DeviceEvent &event) {
    m_mutex->lock();
    m_deviceEvents.push_back(event);
    m_mutex->unlock();
}

void NSGLDevice::processEvents() {
    Assert(m_initialized);

    m_mutex->lock();
    for (std::vector<DeviceEvent>::iterator it = m_deviceEvents.begin(); it!=m_deviceEvents.end(); ++it)
        fireDeviceEvent(*it);
    m_deviceEvents.clear();
    m_mutex->unlock();
}

bool NSGLDevice::isMouseInWindow() {
    return m_fullscreen || NSPointInRect([m_window mouseLocationOutsideOfEventStream], [m_view frame]);
}

void NSGLDevice::makeCurrent(Renderer *renderer) {
    Assert(m_initialized);

    m_currentContext = static_cast<NSOpenGLContext *>(static_cast<NSGLRenderer *>(renderer)->getContext());
    [m_currentContext setView: m_view];
    [m_currentContext makeCurrentContext];
}

MTS_IMPLEMENT_CLASS(NSGLDevice, false, Device)
MTS_NAMESPACE_END

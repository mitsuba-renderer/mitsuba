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

#include <mitsuba/hw/device.h>
#include <mitsuba/core/timer.h>
#if defined(WIN32)
#include <mitsuba/hw/wgldevice.h>
#elif defined(__OSX__)
#include <mitsuba/hw/nsgldevice.h>
#else
#include <mitsuba/hw/glxdevice.h>
#endif

MTS_NAMESPACE_BEGIN

Device::Device(Session *name) {
    m_initialized = false;
    m_redBits = m_greenBits = m_blueBits = 8;
    m_alphaBits = 0;
    m_depthBits = 16;
    m_stencilBits = 8;
    m_doubleBuffer = true;
    m_fullscreen = false;
    m_fsaa = 1;
    m_size = Vector2i(640, 480);
    m_position = Point2i(0, 0);
    m_center = true;
    m_session = name;
    m_showFPS = true;
    m_fps = 0;
    m_fpsCounter = 0;
    m_lastTime = 0;
    m_timer = new Timer();
    m_resizeAllowed = true;
}

Device::~Device() {
}

Device *Device::create(Session *name) {
#if defined(WIN32)
    return new WGLDevice(static_cast<WGLSession *>(name));
#elif defined(__OSX__)
    return new NSGLDevice(static_cast<NSGLSession *>(name));
#else
    return new GLXDevice(static_cast<X11Session *>(name));
#endif
}

void Device::setPosition(const Point2i &position) {
    m_position = position;
}

void Device::setTitle(const std::string &title) {
    m_title = title;
}

void Device::setSize(const Vector2i &dimension) {
    Assert(!m_initialized);
    m_size = dimension;
}

void Device::setRedBits(int redBits) {
    Assert(!m_initialized);
    m_redBits = redBits;
}

void Device::setGreenBits(int greenBits) {
    Assert(!m_initialized);
    m_greenBits = greenBits;
}

void Device::setBlueBits(int blueBits) {
    Assert(!m_initialized);
    m_blueBits = blueBits;
}

void Device::setFSAA(int fsaa) {
    Assert(!m_initialized);
    m_fsaa = fsaa;
}

void Device::setColorBits(int colorBits) {
    Assert(!m_initialized);
    m_redBits = m_greenBits = m_blueBits = colorBits;
}

void Device::setAlphaBits(int alphaBits) {
    Assert(!m_initialized);
    m_alphaBits = alphaBits;
}

void Device::setDepthBits(int depthBits) {
    Assert(!m_initialized);
    m_depthBits = depthBits;
}

void Device::setStencilBits(int stencilBits) {
    Assert(!m_initialized);
    m_stencilBits = stencilBits;
}

void Device::setDoubleBuffer(bool doubleBuffer) {
    Assert(!m_initialized);
    m_doubleBuffer = doubleBuffer;
}

void Device::setFullscreen(bool fullscreen) {
    Assert(!m_initialized);
    m_fullscreen = fullscreen;
}

void Device::setResizeAllowed(bool resizeAllowed) {
    Assert(!m_initialized);
    m_resizeAllowed = resizeAllowed;
}

void Device::setCenter(bool center) {
    Assert(!m_initialized);
    m_center = center;
}

void Device::init(Device *other) {
    Assert(!m_initialized);
    m_session->m_devices.push_back(this);
}

void Device::shutdown() {
    Assert(m_initialized);
    m_session->m_devices.erase(
        std::remove(m_session->m_devices.begin(), m_session->m_devices.end(), this),
        m_session->m_devices.end()
    );
}

void Device::addCallback(DeviceEventListener *callback) {
    m_callbacks.push_back(callback);
}

void Device::removeCallback(DeviceEventListener *callback) {
    m_callbacks.remove(callback);
}

void Device::fireDeviceEvent(const DeviceEvent &event) {
    for (std::list<DeviceEventListener *>::iterator it = m_callbacks.begin();
            it != m_callbacks.end(); ++it) {
        (*it)->deviceEventOccurred(event);
    }
}

void Device::flip() {
    m_fpsCounter++;
    int currentTime = m_timer->getMilliseconds();
    if (currentTime - m_lastTime > 1000) {
        m_lastTime = currentTime;
        m_fps = m_fpsCounter;
        m_fpsCounter = 0;
        if (m_showFPS)
            setTitle(m_title);
    }
}

std::string DeviceEvent::toString() const {
    std::ostringstream oss;
    oss << "DeviceEvent[type=";
    switch (m_type) {
        case Device::EQuitEvent:
            oss << "quit";
            break;
        case Device::EGainFocusEvent:
            oss << "gainFocus";
            break;
        case Device::ELoseFocusEvent:
            oss << "loseFocus";
            break;
        case Device::EKeyDownEvent:
            oss << "keyDown, key='" << (m_keyboard.key == 0 ? ' ' : m_keyboard.key) << "', special="
                << (int) m_keyboard.special << ", modifiers=" << (int) m_keyboard.modifiers
                << ", interpreted='" << m_keyboard.interpreted << "'";
            break;
        case Device::EKeyUpEvent:
            oss << "keyUp, key='" << (m_keyboard.key == 0 ? ' ' : m_keyboard.key) << "', special="
                << (int) m_keyboard.special << ", modifiers=" << (int) m_keyboard.modifiers
                << ", interpreted='" << m_keyboard.interpreted << "'";
            break;
        case Device::EMouseMotionEvent:
            oss << "mouseMotion, x=" << m_mouse.x << ", y=" << m_mouse.y
                << ", xrel=" << m_mouse.xrel << ", yrel=" << m_mouse.yrel
                << ", button=" << (int) m_mouse.button;
            break;
        case Device::EMouseEnterEvent:
            oss << "mouseEnter, x=" << m_mouse.x << ", y=" << m_mouse.y
                << ", xrel=" << m_mouse.xrel << ", yrel=" << m_mouse.yrel
                << ", button=" << (int) m_mouse.button;
            break;
        case Device::EMouseLeaveEvent:
            oss << "mouseLeave, x=" << m_mouse.x << ", y=" << m_mouse.y
                << ", xrel=" << m_mouse.xrel << ", yrel=" << m_mouse.yrel
                << ", button=" << (int) m_mouse.button;
            break;
        case Device::EMouseDragEvent:
            oss << "mouseDrag, x=" << m_mouse.x << ", y=" << m_mouse.y
                << ", xrel=" << m_mouse.xrel << ", yrel=" << m_mouse.yrel
                << ", button=" << (int) m_mouse.button;
            break;
        case Device::EMouseBeginDragEvent:
            oss << "mouseBeginDrag, x=" << m_mouse.x << ", y=" << m_mouse.y
                << ", xrel=" << m_mouse.xrel << ", yrel=" << m_mouse.yrel
                << ", button=" << (int) m_mouse.button;
            break;
        case Device::EMouseEndDragEvent:
            oss << "mouseEndDrag, x=" << m_mouse.x << ", y=" << m_mouse.y
                << ", xrel=" << m_mouse.xrel << ", yrel=" << m_mouse.yrel
                << ", button=" << (int) m_mouse.button;
            break;
        case Device::EMouseButtonDownEvent:
            oss << "mouseButtonDown, x=" << m_mouse.x << ", y=" << m_mouse.y
                << ", xrel=" << m_mouse.xrel << ", yrel=" << m_mouse.yrel
                << ", button=" << (int) m_mouse.button;
            break;
        case Device::EMouseButtonUpEvent:
            oss << "mouseButtonUp, x=" << m_mouse.x << ", y=" << m_mouse.y
                << ", xrel=" << m_mouse.xrel << ", yrel=" << m_mouse.yrel
                << ", button=" << (int) m_mouse.button;
            break;
        case Device::EMouseDoubleClickEvent:
            oss << "mouseDoubleClick, x=" << m_mouse.x << ", y=" << m_mouse.y
                << ", xrel=" << m_mouse.xrel << ", yrel=" << m_mouse.yrel
                << ", button=" << (int) m_mouse.button;
            break;
        case Device::ENoEvent:
            oss << "none";
            break;
        default:
            oss << "unknown";
            break;
    }
    oss << "]";
    return oss.str();
}

MTS_IMPLEMENT_CLASS(Device, true, Object)
MTS_NAMESPACE_END

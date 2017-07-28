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

#pragma once
#if !defined(__MITSUBA_HW_DEVICE_H_)
#define __MITSUBA_HW_DEVICE_H_

#include <mitsuba/hw/session.h>
#include <list>

MTS_NAMESPACE_BEGIN

class Renderer;

/** \brief The device event structure encapsulates event
 * information such as mouse movement or key presses
 * \ingroup libhw
 */
struct MTS_EXPORT_HW DeviceEvent {
public:
    /// Default constructor
    inline DeviceEvent() { }

    /// Type constructor
    inline DeviceEvent(unsigned short type) { m_type = type; }

    /// Return the event type
    inline unsigned short getType() const {  return m_type; }

    /// Set the event type
    inline void setType(int type) { m_type = type; }

    /// Get the mouse position vector
    inline void setMousePosition(const Point2i &position) { m_mouse.x = position.x; m_mouse.y = position.y; }

    /// Set the mouse position vector
    inline Point2i getMousePosition() const { return Point2i(m_mouse.x, m_mouse.y); }

    /// Get the relative mouse movement vector
    inline void setMouseRelative(const Vector2i &relative) { m_mouse.xrel = relative.x; m_mouse.yrel = relative.y; }

    /// Set the relative mouse movement vector
    inline Vector2i getMouseRelative() const { return Vector2i(m_mouse.xrel, m_mouse.yrel); }

    /// Set the enum of the pressed mouse button
    inline void setMouseButton(unsigned short button) { m_mouse.button = button; }

    /// Get the enum of the pressed mouse button
    inline unsigned short getMouseButton() const { return m_mouse.button; }

    /// Set the pressed keyboard key (latin1)
    inline void setKeyboardKey(char key) { m_keyboard.key = key; }

    /// Get the pressed keyboard key (latin1)
    inline char getKeyboardKey() const { return m_keyboard.key; }

    /// Set the pressed keyboard key special identifier enum (see Device::ESpecialKeys)
    inline void setKeyboardSpecial(unsigned short special) { m_keyboard.special = special; }

    /// Get the pressed keyboard key special identifier enum (see Device::ESpecialKeys)
    inline unsigned short getKeyboardSpecial() const { return m_keyboard.special; }

    /// Set the keyboard modifiers (see Device::EKeyboardModifiers)
    inline void setKeyboardModifiers(unsigned short modifiers) { m_keyboard.modifiers = modifiers; }

    /// Get the keyboard modifiers (see Device::EKeyboardModifiers)
    inline unsigned short getKeyboardModifiers() const { return m_keyboard.modifiers; }

    /// Get the interpreted keypress data
    inline const char* getKeyboardInterpreted() const { return m_keyboard.interpreted; }

    /// Get the interpreted keypress data
    inline char* getKeyboardInterpreted() { return m_keyboard.interpreted; }

    /// Set the event source
    inline void setActionSource(Object *source) { m_action.source = source; }

    /// Get the event source
    inline Object *getActionSource() { return m_action.source; }

    /// Return a string representation
    std::string toString() const;
private:
    unsigned short m_type;

    union {
        struct {
            unsigned short special, modifiers;
            char key, interpreted[16];
        } m_keyboard;

        struct {
            int x, y, xrel, yrel;
            unsigned short button;
        } m_mouse;

        struct {
            Object *source;
        } m_action;
    };
};

/** \brief Abstract device event callback
 * \ingroup libhw
 */
class MTS_EXPORT_HW DeviceEventListener {
public:
    /** \brief Called when a device event occurs
     * \param event The event data structure
     * \return True if the result has been handled, false otherwise
     */
    virtual bool deviceEventOccurred(const DeviceEvent &event) = 0;
protected:
    /// Virtual destructor
    virtual ~DeviceEventListener() { }
};

/** \brief An abstract drawing device
 * \ingroup libhw
 */
class MTS_EXPORT_HW Device : public Object {
public:
    /// Device event types
    enum EEventType {
        ENoEvent = 0x0000,
        EQuitEvent = 0x0001,
        EKeyDownEvent = 0x0002,
        EKeyUpEvent = 0x0004,
        EMouseMotionEvent = 0x0008,
        EMouseDragEvent = 0x0010,
        EMouseButtonDownEvent = 0x0020,
        EMouseButtonUpEvent = 0x0040,
        EMouseEnterEvent = 0x0080,
        EMouseLeaveEvent = 0x0100,
        EMouseBeginDragEvent = 0x0200,
        EMouseEndDragEvent = 0x0400,
        EMouseDoubleClickEvent = 0x0800,
        EGainFocusEvent = 0x1000,
        ELoseFocusEvent = 0x2000,
        EResizeEvent = 0x4000
    };

    /// Device keyboard event modifiers
    enum EKeyboardModifiers {
        EShiftModifier = 0x01,
        EControlModifier = 0x02,
        EAltModifier = 0x04,
        EMetaModifier = 0x08
    };

    /// Device keyboard event modifiers
    enum EMouseButton {
        ENoButton = 0x0,
        ELeftButton = 0x01,
        EMiddleButton = 0x02,
        ERightButton = 0x04,
        EWheelUpButton = 0x08,
        EWheelDownButton = 0x10
    };

    /// Device special keys
    enum ESpecialKeys {
        ENoSpecial = 0,
        EKeyEscape,
        EKeyF1,
        EKeyF2,
        EKeyF3,
        EKeyF4,
        EKeyF5,
        EKeyF6,
        EKeyF7,
        EKeyF8,
        EKeyF9,
        EKeyF10,
        EKeyF11,
        EKeyF12,
        EKeyF13,
        EKeyF14,
        EKeyF15,
        EKeyBackspace,
        EKeyTab,
        EKeyClear,
        EKeyReturn,
        EKeyPause,
        EKeyInsert,
        EKeyDelete,
        EKeyUp,
        EKeyDown,
        EKeyLeft,
        EKeyRight,
        EKeyHome,
        EKeyEnd,
        EKeyPageUp,
        EKeyPageDown,
        EKeyNumLock,
        EKeyCapsLock,
        EKeyScrollLock,
        EKeyLShift,
        EKeyRShift,
        EKeyLAlt,
        EKeyRAlt,
        EKeyLMeta,
        EKeyRMeta,
        EKeyLControl,
        EKeyRControl,
        EKeyKeyPad0,
        EKeyKeyPad1,
        EKeyKeyPad2,
        EKeyKeyPad3,
        EKeyKeyPad4,
        EKeyKeyPad5,
        EKeyKeyPad6,
        EKeyKeyPad7,
        EKeyKeyPad8,
        EKeyKeyPad9,
        EKeyKeyPadPeriod,
        EKeyKeyPadDivide,
        EKeyKeyPadMultiply,
        EKeyKeyPadMinus,
        EKeyKeyPadPlus,
        EKeyKeyPadEnter,
        EKeyKeyPadEquals,
        EKeyLastSpecialKey
    };

    /// Construct a new device using the appropriate implementation
    static Device *create(Session *session);

    /// Return the dimension of the device
    inline Vector2i getSize() const { return m_size; }

    /// Set the dimension of the device
    void setSize(const Vector2i &dimension);

    /// Return the aspect ratio of the device
    inline Float getAspect() const { return (Float) m_size.x / (Float) m_size.y; }

    /// Return the position of the device
    inline Point2i getPosition() const { return m_position; }

    /// Set the position of the device
    virtual void setPosition(const Point2i &position);

    /// Set the FSAA sample count, do this before Init()
    void setFSAA(int fsaa);

    /// Return the FSAA sample count
    inline int getFSAA() const { return m_fsaa; }

    /// Only applies to devices, which are UI windows
    virtual void setVisible(bool visible) = 0;

    /** \brief A convenience method.
     *
     * Sets the amount of bits for the red, green and
     * blue components
     */
    void setColorBits(int colorBits);

    /// Set the amount of bits for the red component
    void setRedBits(int redBits);

    /// Return the amount of bits for the red component
    inline int getRedBits() const { return m_redBits; }

    /// Set the amount of bits for the green component
    void setGreenBits(int greenBits);

    /// Return the amount of bits for the green component
    inline int getGreenBits() const { return m_greenBits; }

    /// Set the amount of bits for the blue component
    void setBlueBits(int blueBits);

    /// Return the amount of bits for the blue component
    inline int getBlueBits() const { return m_blueBits; }

    /// Set the amount of bits for the alpha component
    void setAlphaBits(int alphaBits);

    /// Return the amount of bits for the alpha component
    inline int getAlphaBits() const { return m_alphaBits; }

    /// Set the amount of bits for the depth component
    void setDepthBits(int depthBits);

    /// Return the amount of bits for the depth component
    inline int getDepthBits() const { return m_depthBits; }

    /// Set the amount of bits for the stencil component
    void setStencilBits(int stencilBits);

    /// Return the amount of bits for the stencil component
    inline int getStencilBits() const { return m_stencilBits; }

    /// Define whether to enable double buffering
    void setDoubleBuffer(bool doubleBuffer);

    // Return whether double buffering is enabled
    inline bool getDoubleBuffer() const { return m_doubleBuffer; }

    // Define whether to enable full screen drawing
    void setFullscreen(bool fullscreen);

    /// Return whether full screen drawing is enabled
    inline bool getFullscreen() const { return m_fullscreen; }

    // Specify whether resizing the window is allowed
    void setResizeAllowed(bool resizeAllowed);

    /// Return whether it is possible to resize the window
    inline bool isResizeAllowed() const { return m_resizeAllowed; }

    /// Define whether to enable window centering
    void setCenter(bool center);

    /// Return whether window centering is enabled
    inline bool getCenter() const { return m_center; }

    /// Define whether to show the frames per second
    inline void setShowFPS(bool showFPS) { m_showFPS = showFPS; }

    /// Return whether to show the frames per second
    inline bool getShowFPS() const { return m_showFPS; }

    /// Return the frames per second (0 if no data is available)
    inline int getFPS() const { return m_fps; }

    /// Set the x window position
    void setXPos(int xpos);

    /// Set the window title
    virtual void setTitle(const std::string &title);

    /// Return the window title
    inline const std::string &getTitle() const { return m_title; }

    /// Get the session
    inline const Session *getSession() const { return m_session.get(); }

    /// Get the session
    inline Session *getSession() { return m_session; }

    /**
     * Initialize the renderer. Optionally, an existing device instance
     * can be provided as a second argument -- this is primarily meant
     * to create a device that will be able to support a shared context
     * with another device.
     */
    virtual void init(Device *other = NULL);

    /// Shut the device down
    virtual void shutdown();

    /// Add an event callback to the device
    void addCallback(DeviceEventListener *callback);

    /// Remove an event callback from the device
    void removeCallback(DeviceEventListener *callback);

    /// Associate the renderer with this device
    virtual void makeCurrent(Renderer *renderer) = 0;

    /// Flip the buffers (when using double buffering)
    virtual void flip();

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~Device();

    /// Create a new device
    Device(Session *session);

    /** \brief Send a device event using
     * the registered callbacks
     */
    void fireDeviceEvent(const DeviceEvent &event);
protected:
    ref<Session> m_session;
    ref<Timer> m_timer;
    Vector2i m_size;
    Point2i m_position;
    int m_fsaa;
    int m_redBits, m_greenBits, m_blueBits;
    int m_alphaBits, m_depthBits, m_stencilBits;
    bool m_doubleBuffer, m_initialized, m_fullscreen;
    bool m_center, m_showFPS, m_resizeAllowed;
    int m_fpsCounter, m_fps, m_lastTime;
    std::string m_title;
    std::list<DeviceEventListener *> m_callbacks;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_DEVICE_H_ */

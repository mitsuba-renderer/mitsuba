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

#if !defined(__MITSUBA_HW_X11SESSION_H_)
#define __MITSUBA_HW_X11SESSION_H_

#include <mitsuba/hw/session.h>
#include <GL/glx.h>
#include <X11/extensions/xf86vmode.h>

MTS_NAMESPACE_BEGIN

/** \brief X Window System (X11R6) session
 * \ingroup libhw
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

    /**
     * \brief Process all events and call event callbacks.
     *
     * This function will run until the \c stop parameter is set
     * to \c true from within an event callback.
     */
    void processEventsBlocking(bool &stop);

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

#endif /* __MITSUBA_HW_X11SESSION_H_ */

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
#if !defined(__MITSUBA_HW_SESSION_H_)
#define __MITSUBA_HW_SESSION_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

class Device;

/** \brief Abstract windowing environment session
 * \ingroup libhw
 */
class MTS_EXPORT_HW Session : public Object {
    friend class Device;
public:
    /// Create a new session using the appropriate implementation
    static Session *create();

    /// Initialize the session
    virtual void init();

    /// Shut the session down
    virtual void shutdown();

    /// Process all events and call event callbacks
    virtual void processEvents() = 0;

    /**
     * \brief Process all events and call event callbacks.
     *
     * This function will run until the \c stop parameter is set
     * to \c true from within an event callback.
     */
    virtual void processEventsBlocking(bool &stop) = 0;

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~Session() { }

    /// Create a new session
    Session();
protected:
    bool m_initialized;
    std::vector<Device *> m_devices;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_SESSION_H_ */

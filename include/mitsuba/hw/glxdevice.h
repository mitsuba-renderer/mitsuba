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

#if !defined(__MITSUBA_HW_GLXDEVICE_H_)
#define __MITSUBA_HW_GLXDEVICE_H_

#include <mitsuba/hw/x11device.h>

MTS_NAMESPACE_BEGIN

/** \brief X Windows OpenGL-capable (GLX) device
 * \ingroup libhw
 */
class MTS_EXPORT_HW GLXDevice : public X11Device {
public:
    /// Create a new device
    GLXDevice(X11Session *session);

    /// Flip the buffers (when using double-buffering)
    void flip();

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~GLXDevice();

    /// Create a visual
    virtual XVisualInfo* createVisual();
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_GLXDEVICE_H_ */

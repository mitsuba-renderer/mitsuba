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

#if !defined(__MITSUBA_HW_NSGLRENDERER_H_)
#define __MITSUBA_HW_NSGLRENDERER_H_

#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/hw/nsgldevice.h>

#ifdef __OBJC__
#include <Cocoa/Cocoa.h>
#endif

MTS_NAMESPACE_BEGIN

/** \brief A MacOS X (NSGL) OpenGL Renderer
 * \ingroup libhw
 */
class MTS_EXPORT_HW NSGLRenderer : public GLRenderer {
public:
    /// Create a new renderer
    NSGLRenderer(NSGLSession *session);

    /// Initialize the renderer
    void init(Device *device, Renderer *other = NULL);

    /// Shut the renderer down
    void shutdown();

    /// Return the rendering context
    inline void *getContext() { return m_context; }

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~NSGLRenderer();

    /// Lookup an OpenGL extension
    void *lookupExtension(const std::string &name) const;
private:
#ifdef __OBJC__
    NSOpenGLContext *m_context;
#else
    void *m_context;
#endif
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_NSGLRENDERER_H_ */

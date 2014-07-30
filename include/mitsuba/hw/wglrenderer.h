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
#if !defined(__MITSUBA_HW_WGLRENDERER_H_)
#define __MITSUBA_HW_WGLRENDERER_H_

#include <mitsuba/hw/wgldevice.h>
#include <mitsuba/hw/glrenderer.h>

MTS_NAMESPACE_BEGIN

/** \brief Windows (WGL) renderer implementation
 * \ingroup libhw
 */
class MTS_EXPORT_HW WGLRenderer : public GLRenderer {
public:
	/// Create a new renderer
	WGLRenderer(WGLSession *session);

	/// Return the rendering context
	inline HGLRC getWGLContext() { return m_context; }

	/// Initialize the renderer
	void init(Device *device, Renderer *other = NULL);

	/// Shut the renderer down
	void shutdown();

	/// Lookup an OpenGL extension
	virtual void *lookupExtension(const std::string &name) const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~WGLRenderer();
private:
	HGLRC m_context;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_WGLRENDERER_H_ */

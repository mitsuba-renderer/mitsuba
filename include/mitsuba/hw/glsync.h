/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#if !defined(__GLSYNC_H)
#define __GLSYNC_H

#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/hw/gpusync.h>

MTS_NAMESPACE_BEGIN

/** \brief OpenGL-based GPUSync implementation
 * \ingroup libhw
 */
class MTS_EXPORT_HW GLSync : public GPUSync {
public:
	/// Allocate memory for a new synchronization object
	GLSync();

	/// Create the synchronization object on the GL
	void init();

	/// Wait on the fence (blocking)
	void wait();

	/// Enqueue a wait command, but do not block
	void enqueueWait();

	/// Remove the synchronization object
	void cleanup();

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	~GLSync();
protected:
	GLsync m_sync;
};

MTS_NAMESPACE_END

#endif /* __GLSYNC_H */

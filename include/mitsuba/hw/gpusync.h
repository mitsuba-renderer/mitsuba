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
#if !defined(__MITSUBA_HW_GPUSYNC_H_)
#define __MITSUBA_HW_GPUSYNC_H_

#include <mitsuba/hw/renderer.h>

MTS_NAMESPACE_BEGIN

/** \brief Abstract GPU synchronization object implementing
 * a memory fence operation.
 * \ingroup libhw
 */
class MTS_EXPORT_HW GPUSync : public Object {
public:
	/// Allocate memory for a new synchronization object
	GPUSync();

	/// Create the synchronization object on the GPU
	virtual void init() = 0;

	/// Wait on the fence (blocking)
	virtual void wait() = 0;

	/// Enqueue a wait command, but do not block
	virtual void enqueueWait() = 0;

	/// Remove the synchronization object
	virtual void cleanup() = 0;

	/// Return a string representation of this class
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GPUSync();
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_GPUSYNC_H_ */

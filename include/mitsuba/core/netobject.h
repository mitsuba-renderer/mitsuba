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
#if !defined(__MITSUBA_CORE_NETOBJECT_H_)
#define __MITSUBA_CORE_NETOBJECT_H_

#include <mitsuba/core/cobject.h>

MTS_NAMESPACE_BEGIN

/** \brief Abstract interface for objects that reference shared network
 * resources.
 *
 * When a networked object is serialized as part of a parallel process
 * executed on multiple machines, the object is first given the
 * opportunity to bind named resources to the process (by a call to
 * \ref bindUsedResources()). These will then be distributed to all
 * participating compute servers. Once unserialized on the remote side,
 * \ref wakeup() is called to let the object re-associate with the
 * shared resources.
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE NetworkedObject : public ConfigurableObject {
public:
	/// Bind any used resources to the process \a proc
	virtual void bindUsedResources(ParallelProcess *proc) const;

	/// Retrieve any required resources
	virtual void wakeup(ConfigurableObject *parent,
		std::map<std::string, SerializableObject *> &params);

	/// Serialize this object to a stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~NetworkedObject() { }

	/// Constructor
	inline NetworkedObject(const Properties &props) : ConfigurableObject(props) { }

	/// Unserialize a configurable object
	inline NetworkedObject(Stream *stream, InstanceManager *manager)
	 : ConfigurableObject(stream, manager) {
	}
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_NETOBJECT_H_ */

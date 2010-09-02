/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__NETOBJECT_H)
#define __NETOBJECT_H

#include <mitsuba/core/sched.h>

MTS_NAMESPACE_BEGIN

/** \brief Interface of an abstract object referencing
 * globally shared resources. When it is serialized for use in a
 * parallel process executed on several machines, the object
 * is first given the opportunity to bind named resources to 
 * the process, which will then be distributed to all participating 
 * compute servers. Once unserialized on the remote side, 
 * <tt>wakeup</tt> is called to let the object re-associate
 * with the shared resources.
 */
class MTS_EXPORT_CORE NetworkedObject : public ConfigurableObject {
public:
	/// Constructor
	inline NetworkedObject(const Properties &props) : ConfigurableObject(props) { }

	/// Unserialize a configurable object
	inline NetworkedObject(Stream *stream, InstanceManager *manager) 
	 : ConfigurableObject(stream, manager) {
	}

	virtual void bindUsedResources(ParallelProcess *proc) const;
	virtual void wakeup(std::map<std::string, SerializableObject *> &params);

	/// Serialize this object to a stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~NetworkedObject() { }
};

MTS_NAMESPACE_END

#endif /* __NETOBJECT_H */

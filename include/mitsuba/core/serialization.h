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
#if !defined(__MITSUBA_CORE_SERIALIZATION_H_)
#define __MITSUBA_CORE_SERIALIZATION_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/** \brief Base class of all reference-counted objects with serialization support
 *
 * To support unserialization from a stream, the implementation should use one of the
 * RTTI macros \ref MTS_IMPLEMENT_CLASS_S or \ref MTS_IMPLEMENT_CLASS_IS.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE SerializableObject : public Object {
public:
	/// Unserialize a serializable object
	SerializableObject(Stream *stream, InstanceManager *manager);

	/// Serialize this object to a stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const = 0;

	MTS_DECLARE_CLASS()
protected:
	/// Construct a serializable object
	inline SerializableObject() { }

	/// Virtual deconstructor
	virtual ~SerializableObject() { }
};

/** \brief Coordinates the serialization and unserialization of object graphs
 *
 * When serializaing a complicated object graph to a binary data stream,
 * the instance manager annotates the data stream to avoid serializing
 * objects twice or becoming stuck in a cyclic dependency. This allows
 * arbitrary connected graphs to be serialized.
 *
 * Similarly when unserializing a stream, it ensures that the resulting
 * object graph has the same structure.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE InstanceManager : public Object {
	friend class SerializableObject;
public:
	/// \brief Construct a new instance manager
	InstanceManager();

	/// Retrieve an instance from the given stream
	SerializableObject *getInstance(Stream *stream);

	/// Store an instance to the given stream
	void serialize(Stream *stream, const SerializableObject *inst);

	MTS_DECLARE_CLASS()
private:
	/// Virtual destructor
	virtual ~InstanceManager();

	/// Called from the unserialization constructor of SerializableObject
	void registerInstance(SerializableObject *object);
private:
	unsigned int m_counter, m_lastID;
	std::vector<SerializableObject *> m_fullyAllocated;
	std::map<unsigned int, SerializableObject *> m_idToObj;
	std::map<const SerializableObject *, unsigned int> m_objToId;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_SERIALIZATION_H_ */

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

#include <mitsuba/core/serialization.h>

MTS_NAMESPACE_BEGIN

// #define DEBUG_SERIALIZATION 1

InstanceManager::InstanceManager() : m_counter(0) {
#ifdef DEBUG_SERIALIZATION
    Log(EDebug, "Creating an instance manger");
#endif
}

InstanceManager::~InstanceManager() {
#ifdef DEBUG_SERIALIZATION
    Log(EDebug, "Destroying an instance manager");
#endif
    for (std::vector<SerializableObject *>::iterator it = m_fullyAllocated.begin();
        it!= m_fullyAllocated.end(); ++it) {
        (*it)->decRef();
    }
}

SerializableObject *InstanceManager::getInstance(Stream *stream) {
    m_lastID = stream->readUInt();
    if (m_lastID == 0) {
        return NULL;
    } else if (m_idToObj.find(m_lastID) != m_idToObj.end()) {
        return m_idToObj[m_lastID];
    } else {
        SerializableObject *object = NULL;
        std::string className = stream->readString();
#ifdef DEBUG_SERIALIZATION
        Log(EDebug, "Unserializing a class of type '%s'", className.c_str());
#endif
        const Class *theClass = Class::forName(className);
        if (theClass == NULL)
            Log(EError, "Class with name '%s' not found!", className.c_str());
        try {
            object = static_cast<SerializableObject *>
                (theClass->unserialize(stream, this));
        } catch (std::exception &e) {
            Log(EError, "Encountered an exception while unserializing an "
                "instance of \"%s\": \"%s\"!", className.c_str(), e.what());
        }
        m_fullyAllocated.push_back(object);
        object->incRef();
        return object;
    }
}

void InstanceManager::registerInstance(SerializableObject *object) {
    m_idToObj[m_lastID] = object;
}

void InstanceManager::serialize(Stream *stream, const SerializableObject *inst) {
    if (inst == NULL) {
        stream->writeUInt(0);
    } else if (m_objToId.find(inst) != m_objToId.end()) {
        stream->writeUInt(m_objToId[inst]);
    } else {
#ifdef DEBUG_SERIALIZATION
        Log(EDebug, "Serializing a class of type '%s'", inst->getClass()->getName().c_str());
#endif
        stream->writeUInt(++m_counter);
        stream->writeString(inst->getClass()->getName());
        m_objToId[inst]=m_counter;
        inst->serialize(stream, this);
    }
}

SerializableObject::SerializableObject(Stream *stream, InstanceManager *manager) {
    manager->registerInstance(this);
}

MTS_IMPLEMENT_CLASS(SerializableObject, true, Object)
MTS_IMPLEMENT_CLASS(InstanceManager, false, Object)
MTS_NAMESPACE_END

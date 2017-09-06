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

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

typedef std::map<std::string, Class *, std::less<std::string> > ClassMap;
static ClassMap *__classes;
bool Class::m_isInitialized = false;
typedef Object *(*inst_function)();
typedef Object *(*unSer_function)(Stream *, InstanceManager *);

Class::Class(const std::string &name, bool abstract, const std::string &superClassName, void *instPtr, void *unSerPtr)
 : m_name(name), m_abstract(abstract), m_superClass(NULL), m_superClassName(superClassName), m_instPtr(instPtr), m_unSerPtr(unSerPtr) {
    if (__classes == NULL)
        __classes = new ClassMap();

    (*__classes)[name] = this;
}

const Class *Class::forName(const char *name) {
    if (name != NULL) {
        if (__classes->find(name) != __classes->end())
            return (*__classes)[name];
    }

    return NULL;
}

const Class *Class::forName(const std::string &name) {
    if (__classes->find(name) != __classes->end())
        return (*__classes)[name];

    return NULL;
}

bool Class::derivesFrom(const Class *theClass) const {
    const Class *mClass = this;

    while (mClass != NULL) {
        if (theClass == mClass)
            return true;
        mClass = mClass->getSuperClass();
    }

    return false;
}

void Class::initializeOnce(Class *theClass) {
    const std::string &base = theClass->m_superClassName;

    if (base != "") {
        if (__classes->find(base) != __classes->end()) {
            theClass->m_superClass = (*__classes)[base];
        } else {
            std::cerr << "Critical error during the static RTTI initialization: " << endl
                << "Could not locate the base class '" << base << "' while initializing '"
                << theClass->getName() << "'!" << endl;
            exit(-1);
        }
    }
}

void Class::staticInitialization() {
    for (ClassMap::iterator it = __classes->begin();
            it != __classes->end(); ++it)
        initializeOnce(it->second);
    m_isInitialized = true;
}

Object *Class::instantiate() const {
    if (m_instPtr == NULL) {
        SLog(EError, "RTTI error: An attempt to instantiate a "
            "class lacking the instantiation feature occurred (%s)!",
            getName().c_str());
    }
    return ((inst_function) m_instPtr)();
}

Object *Class::unserialize(Stream *stream, InstanceManager *manager) const {
    if (m_unSerPtr == NULL) {
        SLog(EError, "RTTI error: An attempt to instantiate a "
            "class lacking the unserialization feature occurred (%s)!",
            getName().c_str());
    }
    return ((unSer_function) m_unSerPtr)(stream, manager);
}

void Class::staticShutdown() {
    for (ClassMap::iterator it = __classes->begin(); it!= __classes->end(); ++it)
        delete (*it).second;
    delete __classes;
    __classes = NULL;
    m_isInitialized = false;
}

MTS_NAMESPACE_END

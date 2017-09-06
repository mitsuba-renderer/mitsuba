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

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/lock.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/cobject.h>
#include <mitsuba/core/version.h>

#if !defined(__WINDOWS__)
# include <dlfcn.h>
#else
# include <windows.h>
#endif

MTS_NAMESPACE_BEGIN

// -----------------------------------------------------------------------
//  Abstract plugin module implementation
// -----------------------------------------------------------------------

namespace {
    typedef void *(*CreateInstanceFunc)(const Properties &props);
    typedef void *(*CreateUtilityFunc)();
    typedef char *(*GetDescriptionFunc)();
}

struct Plugin::PluginPrivate {
#if defined(__WINDOWS__)
    HMODULE handle;
#else
    void *handle;
#endif
    const std::string shortName;
    const fs::path path;
    bool isUtility;
    GetDescriptionFunc getDescription;
    CreateInstanceFunc createInstance;
    CreateUtilityFunc createUtility;

    PluginPrivate(const std::string &sn, const fs::path &p)
    : shortName(sn), path(p) {}
};

Plugin::Plugin(const std::string &shortName, const fs::path &path)
 : d(new PluginPrivate(shortName, path)) {
#if defined(__WINDOWS__)
    d->handle = LoadLibraryW(path.c_str());
    if (!d->handle) {
        SLog(EError, "Error while loading plugin \"%s\": %s",
                d->path.string().c_str(), lastErrorText().c_str());
    }
#else
    d->handle = dlopen(path.string().c_str(), RTLD_LAZY | RTLD_LOCAL);
    if (!d->handle) {
        SLog(EError, "Error while loading plugin \"%s\": %s",
            d->path.string().c_str(), dlerror());
    }
#endif
    try {
        d->getDescription = (GetDescriptionFunc) getSymbol("GetDescription");
    } catch (...) {
#if defined(__WINDOWS__)
        FreeLibrary(d->handle);
#else
        dlclose(d->handle);
#endif
        throw;
    }

    d->createInstance = NULL;
    d->createUtility = NULL;
    d->isUtility = false;

    if (hasSymbol("CreateUtility")) {
        d->createUtility = (CreateUtilityFunc) getSymbol("CreateUtility");
        d->isUtility = true;
    } else {
        d->createInstance = (CreateInstanceFunc) getSymbol("CreateInstance");
    }
    Statistics::getInstance()->logPlugin(shortName, getDescription());

    /* New classes must be registered within the class hierarchy */
    Class::staticInitialization();
}

bool Plugin::hasSymbol(const std::string &sym) const {
#if defined(__WINDOWS__)
    void *ptr = GetProcAddress(d->handle, sym.c_str());
#else
    void *ptr = dlsym(d->handle, sym.c_str());
#endif
    return ptr != NULL;
}

void *Plugin::getSymbol(const std::string &sym) {
#if defined(__WINDOWS__)
    void *data = GetProcAddress(d->handle, sym.c_str());
    if (!data) {
        SLog(EError, "Could not resolve symbol \"%s\" in \"%s\": %s",
            sym.c_str(), d->path.string().c_str(), lastErrorText().c_str());
    }
#else
    void *data = dlsym(d->handle, sym.c_str());
    if (!data) {
        SLog(EError, "Could not resolve symbol \"%s\" in \"%s\": %s",
            sym.c_str(), d->path.string().c_str(), dlerror());
    }
#endif
    return data;
}

ConfigurableObject *Plugin::createInstance(const Properties &props) const {
    return (ConfigurableObject *) d->createInstance(props);
}

Utility *Plugin::createUtility() const {
    return (Utility *) d->createUtility();
}

bool Plugin::isUtility() const {
    return d->isUtility;
}

std::string Plugin::getDescription() const {
    return d->getDescription();
}

const fs::path& Plugin::getPath() const {
    return d->path;
}

const std::string& Plugin::getShortName() const {
    return d->shortName;
}

Plugin::~Plugin() {
#if defined(__WINDOWS__)
    FreeLibrary(d->handle);
#else
    dlclose(d->handle);
#endif
}

// -----------------------------------------------------------------------
//  Plugin manager
// -----------------------------------------------------------------------

ref<PluginManager> PluginManager::m_instance = NULL;

PluginManager::PluginManager() {
    m_mutex = new Mutex();
}

PluginManager::~PluginManager() {
    /* Release the memory used by plugin modules */
    for (std::map<std::string, Plugin *>::iterator it = m_plugins.begin();
        it != m_plugins.end(); ++it) {
        delete (*it).second;
    }
}

ConfigurableObject *PluginManager::createObject(const Class *classType,
    const Properties &props) {
    ConfigurableObject *object;

    {
        LockGuard lock(m_mutex);
        ensurePluginLoaded(props.getPluginName());
        object = m_plugins[props.getPluginName()]->createInstance(props);
    }
    if (!object->getClass()->derivesFrom(classType))
        Log(EError, "Type mismatch when loading plugin \"%s\": Expected "
        "an instance of \"%s\"", props.getPluginName().c_str(), classType->getName().c_str());
    if (object->getClass()->isAbstract())
        Log(EError, "Error when loading plugin \"%s\": Identifies itself as an abstract class",
        props.getPluginName().c_str());
    return object;
}

ConfigurableObject *PluginManager::createObject(const Properties &props) {
    ConfigurableObject *object;

    {
        LockGuard lock(m_mutex);
        ensurePluginLoaded(props.getPluginName());
        object = m_plugins[props.getPluginName()]->createInstance(props);
    }
    if (object->getClass()->isAbstract())
        Log(EError, "Error when loading plugin \"%s\": Identifies itself as an abstract class",
        props.getPluginName().c_str());
    return object;
}

std::vector<std::string> PluginManager::getLoadedPlugins() const {
    std::vector<std::string> list;
    LockGuard lock(m_mutex);
    for (std::map<std::string, Plugin *>::const_iterator it = m_plugins.begin();
        it != m_plugins.end(); ++it) {
        list.push_back((*it).first);
    }
    return list;
}

void PluginManager::ensurePluginLoaded(const std::string &name) {
    /* Plugin already loaded? */
    if (m_plugins[name] != NULL)
        return;

    /* Build the full plugin file name */
    fs::path shortName = fs::path("plugins") / name;
#if defined(__WINDOWS__)
    shortName.replace_extension(".dll");
#elif defined(__OSX__)
    shortName.replace_extension(".dylib");
#else
    shortName.replace_extension(".so");
#endif

    const FileResolver *resolver = Thread::getThread()->getFileResolver();
    fs::path path = resolver->resolve(shortName);

    if (fs::exists(path)) {
        Log(EInfo, "Loading plugin \"%s\" ..", shortName.string().c_str());
        m_plugins[name] = new Plugin(shortName.string(), path);
        return;
    }

    /* Plugin not found! */
    Log(EError, "Plugin \"%s\" not found!", name.c_str());
}

void PluginManager::staticInitialization() {
    m_instance = new PluginManager();
}

void PluginManager::staticShutdown() {
    m_instance = NULL;
}

Version::Version(const std::string &versionString) {
    std::vector<std::string> tokens = tokenize(trim(versionString), ".");
    if (tokens.size() != 3)
        SLog(EError, "Unable to parse version string \"%s\"!", versionString.c_str());
    char *end_ptr = NULL;
    m_major = strtol(tokens[0].c_str(), &end_ptr, 10);
    if (*end_ptr != '\0')
        SLog(EError, "Unable to parse the major program version \"%i\"!", tokens[0].c_str());
    m_minor = strtol(tokens[1].c_str(), &end_ptr, 10);
    if (*end_ptr != '\0')
        SLog(EError, "Unable to parse the minor program version \"%i\"!", tokens[1].c_str());
    m_release = strtol(tokens[2].c_str(), &end_ptr, 10);
    if (*end_ptr != '\0')
        SLog(EError, "Unable to parse the release program version \"%i\"!", tokens[2].c_str());
}

std::string Version::toString() const {
    return formatString("%i.%i.%i", m_major, m_minor, m_release);
}

std::string Version::toStringComplete() const {
    std::ostringstream oss;

    oss << m_major << "." << m_minor << "." << m_release << " (";
#if defined(__WINDOWS__)
    oss << "Windows, ";
#elif defined(__LINUX__)
    oss << "Linux, ";
#elif defined(__OSX__)
    oss << "Mac OS, ";
#else
    oss << "Unknown, ";
#endif

#if defined(__64BIT__)
    oss << "64 bit)";
#else
    oss << "32 bit)";
#endif
    return oss.str();
}

MTS_IMPLEMENT_CLASS(PluginManager, false, Object)
MTS_NAMESPACE_END

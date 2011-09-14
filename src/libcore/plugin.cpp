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

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/lock.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/cobject.h>
#include <mitsuba/core/version.h>

#if !defined(WIN32)
#include <dlfcn.h>
#endif

MTS_NAMESPACE_BEGIN

// -----------------------------------------------------------------------
//  Abstract plugin module implementation
// -----------------------------------------------------------------------

Plugin::Plugin(const std::string &shortName, const fs::path &path) 
 : m_shortName(shortName), m_path(path) {
#if defined(WIN32)
	m_handle = LoadLibrary(path.file_string().c_str());
	if (!m_handle) {
		SLog(EError, "Error while loading plugin \"%s\": %s", 
				m_path.file_string().c_str(), lastErrorText().c_str());
	}
#else
	m_handle = dlopen(path.file_string().c_str(), RTLD_LAZY | RTLD_LOCAL);
	if (!m_handle) {
		SLog(EError, "Error while loading plugin \"%s\": %s",
			m_path.file_string().c_str(), dlerror());
	}
#endif
	try {
		m_getDescription = (GetDescriptionFunc) getSymbol("GetDescription");
	} catch (...) {
#if defined(WIN32)
		FreeLibrary(m_handle);
#else
		dlclose(m_handle);
#endif
		throw;
	}

	m_createInstance = NULL;
	m_createUtility = NULL;
	m_isUtility = false;

	if (hasSymbol("CreateUtility")) {
		m_createUtility = (CreateUtilityFunc) getSymbol("CreateUtility");
		m_isUtility = true;
	} else {
		m_createInstance = (CreateInstanceFunc) getSymbol("CreateInstance");
	}
	Statistics::getInstance()->logPlugin(shortName, getDescription());

	/* New classes must be registered within the class hierarchy */
	Class::staticInitialization();
}

bool Plugin::hasSymbol(const std::string &sym) const {
#if defined(WIN32)
	void *ptr = GetProcAddress(m_handle, sym.c_str());
#else
	void *ptr = dlsym(m_handle, sym.c_str());
#endif
	return ptr != NULL;
}

void *Plugin::getSymbol(const std::string &sym) {
#if defined(WIN32)
	void *data = GetProcAddress(m_handle, sym.c_str());
	if (!data) {
		SLog(EError, "Could not resolve symbol \"%s\" in \"%s\": %s",
			sym.c_str(), m_path.file_string().c_str(), lastErrorText().c_str());
	}
#else
	void *data = dlsym(m_handle, sym.c_str());
	if (!data) {
		SLog(EError, "Could not resolve symbol \"%s\" in \"%s\": %s",
			sym.c_str(), m_path.file_string().c_str(), dlerror());
	}
#endif
	return data;
}

ConfigurableObject *Plugin::createInstance(const Properties &props) const {
	return (ConfigurableObject *) m_createInstance(props);
}

Utility *Plugin::createUtility() const {
	return (Utility *) m_createUtility();
}

std::string Plugin::getDescription() const {
	return m_getDescription();
}

Plugin::~Plugin() {
#if defined(WIN32)
	FreeLibrary(m_handle);
#else
	dlclose(m_handle);
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

	m_mutex->lock();
	try {
		ensurePluginLoaded(props.getPluginName());
		object = m_plugins[props.getPluginName()]->createInstance(props);
	} catch (std::runtime_error &e) {
		m_mutex->unlock();
		throw e;
	} catch (std::exception &e) {
		m_mutex->unlock();
		throw e;
	}
	m_mutex->unlock();
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

	m_mutex->lock();
	try {
		ensurePluginLoaded(props.getPluginName());
		object = m_plugins[props.getPluginName()]->createInstance(props);
	} catch (std::runtime_error &e) {
		m_mutex->unlock();
		throw e;
	} catch (std::exception &e) {
		m_mutex->unlock();
		throw e;
	}
	m_mutex->unlock();
	if (object->getClass()->isAbstract())
		Log(EError, "Error when loading plugin \"%s\": Identifies itself as an abstract class",
		props.getPluginName().c_str());
	return object;
}

std::vector<std::string> PluginManager::getLoadedPlugins() const {
	std::vector<std::string> list;
	m_mutex->lock();
	for (std::map<std::string, Plugin *>::const_iterator it = m_plugins.begin();
		it != m_plugins.end(); ++it) {
		list.push_back((*it).first);
	}
	m_mutex->unlock();
	return list;
}

void PluginManager::ensurePluginLoaded(const std::string &name) {
	/* Plugin already loaded? */
	if (m_plugins[name] != NULL)
		return;

	/* Build the full plugin file name */
#if defined(WIN32)
	std::string shortName = std::string("plugins/") + name + std::string(".dll");
#elif defined(__OSX__)
	std::string shortName = std::string("plugins/") + name + std::string(".dylib");
#else
	std::string shortName = std::string("plugins/") + name + std::string(".so");
#endif
	const FileResolver *resolver = Thread::getThread()->getFileResolver();
	fs::path path = resolver->resolve(shortName);

	if (fs::exists(path)) {
		Log(EInfo, "Loading plugin \"%s\" ..", shortName.c_str());
		m_plugins[name] = new Plugin(shortName, path.file_string());
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

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fresolver.h>
#include <stdexcept>

#if !defined(WIN32)
#include <dlfcn.h>
#endif

MTS_NAMESPACE_BEGIN

// -----------------------------------------------------------------------
//  Abstract plugin module implementation
// -----------------------------------------------------------------------

Plugin::Plugin(const std::string &shortName, const std::string &path) 
 : m_shortName(shortName), m_path(path) {
#if defined(WIN32)
	m_handle = LoadLibrary(path.c_str());
	if (!m_handle) {
		SLog(EError, "Error while loading plugin \"%s\": %s", m_path.c_str(),
				lastErrorText().c_str());
	}
#else
	m_handle = dlopen(path.c_str(), RTLD_LAZY | RTLD_LOCAL);
	if (!m_handle) {
		SLog(EError, "Error while loading plugin \"%s\": %s", m_path.c_str(),
				dlerror());
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

	try {
		m_createInstance = (CreateInstanceFunc) getSymbol("CreateInstance");
	} catch (const std::exception &ex) {
		m_createUtility = (CreateUtilityFunc) getSymbol("CreateUtility");
		m_isUtility = true;
	}
	Statistics::getInstance()->logPlugin(shortName, getDescription());

	/* New classes must be registered within the class hierarchy */
	Class::staticInitialization();
}

void *Plugin::getSymbol(const std::string &sym) {
#if defined(WIN32)
	void *data = GetProcAddress(m_handle, sym.c_str());
	if (!data) {
		SLog(EError, "Could not resolve symbol \"%s\" in \"%s\": %s",
			sym.c_str(), m_path.c_str(), lastErrorText().c_str());
	}
#else
	void *data = dlsym(m_handle, sym.c_str());
	if (!data) {
		SLog(EError, "Could not resolve symbol \"%s\" in \"%s\": %s",
			sym.c_str(), m_path.c_str(), dlerror());
	}
#endif
	return data;
}

ConfigurableObject *Plugin::createInstance(const Properties &props) const {
	return (ConfigurableObject *) m_createInstance(props);
}

Utility *Plugin::createUtility(UtilityServices *us) const {
	return (Utility *) m_createUtility(us);
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
	FileResolver *resolver = FileResolver::getInstance();
	std::string fullName = resolver->resolve(shortName);

	if (FileStream::exists(fullName)) {
		Log(EInfo, "Loading plugin \"%s\" ..", shortName.c_str());
		m_plugins[name] = new Plugin(shortName, fullName);
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

MTS_IMPLEMENT_CLASS(PluginManager, false, Object)
MTS_NAMESPACE_END

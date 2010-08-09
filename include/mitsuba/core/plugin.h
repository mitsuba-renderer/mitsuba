#if !defined(__PLUGIN_H)
#define __PLUGIN_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

class MTS_EXPORT_CORE Plugin {
	typedef void *(*CreateInstanceFunc)(const Properties &props);
	typedef char *(*GetDescriptionFunc)();
public:
	/// Load a plugin from the supplied path
	Plugin(const std::string &shortName, const std::string &path);

	/// Virtual destructor
	virtual ~Plugin();

	/// Return an instance of the class implemented by this plugin
	ConfigurableObject *createInstance(const Properties &props) const;

	/// Return a description of this plugin
	std::string getDescription() const;
	
	/// Return the path of this plugin
	inline const std::string &getPath() const { return m_path; }
	
	/// Return a short name of this plugin
	inline const std::string &getShortName() const { return m_shortName; }
protected:
	/// Resolve the given symbol and return a pointer
	void *getSymbol(const std::string &sym);
private:
#if defined(WIN32)
	HMODULE m_handle;
#else
	void *m_handle;
#endif
	std::string m_shortName;
	std::string m_path;
	CreateInstanceFunc m_createInstance;
	GetDescriptionFunc m_getDescription;
};

class FileResolver;

class MTS_EXPORT_CORE PluginManager : public Object {
public:
	/// Return the global plugin manager
	inline static PluginManager *getInstance() {
		return m_instance;
	}

	/// Ensure that a plugin is loaded and ready
	void ensurePluginLoaded(const std::string &name);

	/// Return the list of loaded plugins
	std::vector<std::string> getLoadedPlugins() const;

	/// Instantiate an object
	ConfigurableObject *createObject(
		const Class *classType,
		const Properties &props
	);

	MTS_DECLARE_CLASS()
protected:
	PluginManager();
	
	/// Destruct and unload all plugins
	~PluginManager();
private:
	std::map<std::string, Plugin *> m_plugins;
	mutable ref<Mutex> m_mutex;
	static ref<PluginManager> m_instance;
};

MTS_NAMESPACE_END

#endif /* __PLUGIN_H */

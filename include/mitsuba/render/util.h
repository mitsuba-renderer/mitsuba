#if !defined(__UTILITY_H)
#define __UTILITY_H

#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/**
 * Contains functionality provided to utility plugins,
 * such as loading a scene from an XML file.
 */
class MTS_EXPORT_RENDER UtilityServices {
public:
	virtual Scene *loadScene(const std::string &filename) = 0;
};

/** \brief Abstract utility class -- can be used to implement
 * loadable utility plugins that perform various actions. They
 * can be started using the 'mtsutil' launcher.
 */
class MTS_EXPORT_RENDER Utility : public Object {
public:
	inline Utility(UtilityServices *services) 
		: m_utilityServices(services) { }

	/**
	 * Run the utility. The supplied <tt>argc</tt>
	 * and <tt>argv</tt> parameters contain any 
	 * extra arguments passed to mtsutil. The value
	 * returned here will be used as the return value of the
	 * 'mtsutil' process.
	 */
	virtual int run(int argc, char **argv) = 0;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Utility() { }

	/// Load a scene
	inline Scene *loadScene(const std::string &fname) {
		return m_utilityServices->loadScene(fname);
	}
private:
	UtilityServices *m_utilityServices;
};

#define MTS_EXPORT_UTILITY(name, descr) \
	extern "C" { \
		void MTS_EXPORT *CreateUtility(UtilityServices *us) { \
			return new name(us); \
		} \
		const char MTS_EXPORT *GetDescription() { \
			return descr; \
		} \
	}

MTS_NAMESPACE_END

#endif /* __UTILITY_H */

#if !defined(__SHANDLER_H)
#define __SHANDLER_H

#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/sax/AttributeList.hpp>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/core/plugin.h>
#include <stack>
#include <map>

XERCES_CPP_NAMESPACE_USE
using namespace mitsuba;

/**
 * XML parser for mitsuba scene files. Uses Xerces-C and SAX
 */
class SceneHandler : public HandlerBase {
public:
	SceneHandler();
	SceneHandler(const std::map<std::string, std::string> &params,
		bool isIncludedFile = false);
	virtual ~SceneHandler();

	// -----------------------------------------------------------------------
	//  Implementation of the SAX DocumentHandler interface
	// -----------------------------------------------------------------------
	virtual void startDocument();
	virtual void endDocument();
	virtual void startElement(
		const XMLCh* const name,
		AttributeList& attributes
	);
	virtual void endElement(const XMLCh* const name);
	virtual void characters(const XMLCh* const chars, const unsigned int length);

	inline const Scene *getScene() const { return m_scene.get(); }
	inline Scene *getScene() { return m_scene; }

	inline std::string transcode(const XMLCh * const xmlName) const {
		char *value = XMLString::transcode(xmlName);
		std::string result(value);
		XMLString::release(&value);
		return result;
	}

	// -----------------------------------------------------------------------
	//  Implementation of the SAX ErrorHandler interface
	// -----------------------------------------------------------------------
	void warning(const SAXParseException& exc);
	void error(const SAXParseException& exc);
	void fatalError(const SAXParseException& exc);
private:
	struct ParseContext {
		inline ParseContext(ParseContext *_parent)
		 : parent(_parent) {
		}

		ParseContext *parent;
		Properties properties;
		std::map<std::string, std::string> attributes;
		std::vector<std::pair<std::string, ConfigurableObject *> > children;
	};

	ref<Scene> m_scene;
	std::map<std::string, ConfigurableObject *> m_objects;
	std::map<std::string, std::string> m_params;
	PluginManager *m_pluginManager;
	std::stack<ParseContext> m_context;
	Transform m_transform;
	bool m_isIncludedFile;
};

#endif /* __SHANDLER_H */

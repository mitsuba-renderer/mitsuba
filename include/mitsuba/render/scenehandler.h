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

#if !defined(__SHANDLER_H)
#define __SHANDLER_H

#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/sax/AttributeList.hpp>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/version.h>
#include <stack>
#include <map>

XERCES_CPP_NAMESPACE_BEGIN
class SAXParser;
XERCES_CPP_NAMESPACE_END

XERCES_CPP_NAMESPACE_USE
MTS_NAMESPACE_BEGIN

/**
 * \brief This exception is thrown when attempting to load an outdated file
 * \ingroup librender
 */
class VersionException : public std::runtime_error {
public:
	inline VersionException(const std::string &str, const Version &version)
		: std::runtime_error(str), m_version(version) { }

	inline const Version &getVersion() const { return m_version; }
private:
	Version m_version;
};

/**
 * \brief XML parser for Mitsuba scene files. To be used with the
 * SAX interface of Xerces-C++.
 *
 * \remark In the Python bindings, only the static function
 *         \ref loadScene() is exposed.
 * \ingroup librender
 * \ingroup libpython
 */
class MTS_EXPORT_RENDER SceneHandler : public HandlerBase {
public:
	typedef std::map<std::string, ConfigurableObject *> NamedObjectMap;
	typedef std::map<std::string, std::string, SimpleStringOrdering> ParameterMap;

	SceneHandler(const SAXParser *parser, const ParameterMap &params,
			NamedObjectMap *objects = NULL, bool isIncludedFile = false);
	virtual ~SceneHandler();

	/// Convenience method -- load a scene from a given filename 
	static ref<Scene> loadScene(const fs::path &filename,
		const ParameterMap &params= ParameterMap());

	/// Initialize Xerces-C++ (needs to be called once at program startup)
	static void staticInitialization();

	/// Free the memory taken up by staticInitialization()
	static void staticShutdown();

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

	// -----------------------------------------------------------------------
	//  Implementation of the SAX ErrorHandler interface
	// -----------------------------------------------------------------------
	void warning(const SAXParseException& exc);
	void error(const SAXParseException& exc);
	void fatalError(const SAXParseException& exc);
protected:
	inline std::string transcode(const XMLCh * const xmlName) const {
		char *value = XMLString::transcode(xmlName);
		std::string result(value);
		XMLString::release(&value);
		return result;
	}
	Float parseFloat(const std::string &name, const std::string &str,
			Float defVal = -1) const;

	void clear();

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

	const SAXParser *m_parser;
	ref<Scene> m_scene;
	ParameterMap m_params;
	NamedObjectMap *m_namedObjects;
	PluginManager *m_pluginManager;
	std::stack<ParseContext> m_context;
	Transform m_transform;
	bool m_isIncludedFile;
};

MTS_NAMESPACE_END

#endif /* __SHANDLER_H */

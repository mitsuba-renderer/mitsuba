/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

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
MTS_NAMESPACE_BEGIN

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

MTS_NAMESPACE_END

#endif /* __SHANDLER_H */

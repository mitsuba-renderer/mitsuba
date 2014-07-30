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

#include <mitsuba/core/platform.h>

// Mitsuba's "Assert" macro conflicts with Xerces' XSerializeEngine::Assert(...).
// This becomes a problem when using a PCH which contains mitsuba/core/logger.h
#if defined(Assert)
# undef Assert
#endif
#include <xercesc/parsers/SAXParser.hpp>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scenehandler.h>

MTS_NAMESPACE_BEGIN

ref<Scene> Utility::loadScene(const fs::path &filename,
		const ParameterMap &params) {
	return SceneHandler::loadScene(filename, params);
}

ref<Scene> Utility::loadSceneFromString(const std::string &content,
		const ParameterMap &params) {
	return SceneHandler::loadSceneFromString(content, params);
}

MTS_IMPLEMENT_CLASS(Utility, true, Object)
MTS_NAMESPACE_END

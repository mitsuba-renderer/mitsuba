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

#pragma once
#if !defined(__MITSUBA_RENDER_UTIL_H_)
#define __MITSUBA_RENDER_UTIL_H_

#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/** \brief Abstract utility class -- can be used to implement
 * loadable utility plugins that perform various actions. They
 * can be started using the 'mtsutil' launcher.
 * \ingroup librender
 */
class MTS_EXPORT_RENDER Utility : public Object {
public:
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
	typedef std::map<std::string, std::string, SimpleStringOrdering> ParameterMap;

	/// Virtual destructor
	virtual ~Utility() { }

	/// Load a scene from an external file
	ref<Scene> loadScene(const fs::path &fname,
		const ParameterMap &params= ParameterMap());

	/// Load a scene from a string
	ref<Scene> loadSceneFromString(const std::string &content,
		const ParameterMap &params= ParameterMap());
};

#define MTS_DECLARE_UTILITY() \
	MTS_DECLARE_CLASS()

#define MTS_EXPORT_UTILITY(name, descr) \
	MTS_IMPLEMENT_CLASS(name, false, Utility) \
	extern "C" { \
		void MTS_EXPORT *CreateUtility() { \
			return new name(); \
		} \
		const char MTS_EXPORT *GetDescription() { \
			return descr; \
		} \
	}

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_UTIL_H_ */

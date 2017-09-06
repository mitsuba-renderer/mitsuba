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
#if !defined(__MITSUBA_CORE_PLUGIN_H_)
#define __MITSUBA_CORE_PLUGIN_H_

#include <mitsuba/mitsuba.h>
#include <boost/filesystem.hpp>
#include <boost/scoped_ptr.hpp>

MTS_NAMESPACE_BEGIN

/**
 * \brief Abstract plugin class -- represents loadable configurable objects
 * and utilities.
 *
 * Please see the \ref ConfigurableObject and \ref Utility classes for
 * details.
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE Plugin {
public:
    /// Load a plugin from the supplied path
    Plugin(const std::string &shortName, const fs::path &path);

    /// Virtual destructor
    virtual ~Plugin();

    /// Is this a configurable object plugin or an utility plugin?
    bool isUtility() const;

    /// Return an instance of the class implemented by this plugin
    ConfigurableObject *createInstance(const Properties &props) const;

    /// Return an utility instance (if this is an utility plugin)
    Utility *createUtility() const;

    /// Return a description of this plugin
    std::string getDescription() const;

    /// Return the path of this plugin
    const fs::path &getPath() const;

    /// Return a short name of this plugin
    const std::string &getShortName() const;
protected:
    /// Resolve the given symbol and return a pointer
    void *getSymbol(const std::string &sym);
    /// Check whether a certain symbol is provided by the plugin
    bool hasSymbol(const std::string &sym) const;
private:
    struct PluginPrivate;
    boost::scoped_ptr<PluginPrivate> d;
};

/**
 * \brief The plugin manager is responsible for resolving and
 * loading external plugins.
 *
 * Ordinarily, this class will be used by making repeated calls to
 * the \ref createObject() methods. The generated instances are then
 * assembled into a final object graph, such as a scene. One such
 * examples is the \ref SceneHandler class, which parses an XML
 * scene file by esentially translating the XML elements into calls
 * to \ref createObject().
 *
 * Since this kind of construction method can be tiresome when
 * dynamically building scenes from Python, this class has an
 * additional Python-only method \c create(), which works as follows:
 *
 * \code
 * from mitsuba.core import *
 *
 * pmgr = PluginManager.getInstance()
 * camera = pmgr.create({
 *     "type" : "perspective",
 *     "toWorld" : Transform.lookAt(
 *         Point(0, 0, -10),
 *         Point(0, 0, 0),
 *         Vector(0, 1, 0)
 *     ),
 *     "film" : {
 *         "type" : "ldrfilm",
 *         "width" : 1920,
 *         "height" : 1080
 *     }
 * })
 * \endcode
 *
 * The above snippet constructs a \ref Camera instance from a
 * plugin named \c perspective.so/dll/dylib and adds a child object
 * named \c film, which is a \ref Film instance loaded from the
 * plugin \c ldrfilm.so/dll/dylib. By the time the function
 * returns, the object hierarchy has already been assembled, and the
 * \ref ConfigurableObject::configure() methods of every object
 * has been called.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
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

    /**
     * \brief Instantiate a plugin, verify its type,
     * and return the newly created instance.
     *
     * \param classType Expected type of the plugin. An
     *    exception will be thrown if it turns out not
     *    to derive from this class.
     * \param props A \ref Properties instance containing
     *    all information required to find and construct
     *    the plugin.
     */
    ConfigurableObject *createObject(
        const Class *classType,
        const Properties &props
    );

    /**
     * \brief Instantiate a plugin and return the new
     * instance (without verifying its type).
     *
     * \param props A \ref Properties instance containing
     *    all information required to find and construct
     *    the plugin.
     */
    ConfigurableObject *createObject(
        const Properties &props
    );

    /// Initializes the global plugin manager instance
    static void staticInitialization();

    /// Free the memory taken by staticInitialization()
    static void staticShutdown();

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

#endif /* __MITSUBA_CORE_PLUGIN_H_ */

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
#if !defined(__MITSUBA_CORE_COBJECT_H_)
#define __MITSUBA_CORE_COBJECT_H_

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/serialization.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

/** \brief Generic serializable object, which supports construction
* from a Properties instance.
 *
 * All plugins in Mitsuba derive from ConfigurableObject. This mechanism
 * lets them accept parameters specified in an external XML file. Additionally,
 * they can have child objects, which correspond to nested instantiation
 * requests in the XML file.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE ConfigurableObject : public SerializableObject {
public:
    /**
     * \brief Notify the \ref ConfigurableObject instance about
     * its parent object
     *
     * The default implementation does nothing.
     */
    virtual void setParent(ConfigurableObject *parent);

    /// Add a child (default implementation throws an error)
    virtual void addChild(const std::string &name, ConfigurableObject *child);

    /// Add an unnamed child
    inline void addChild(ConfigurableObject *child) { addChild("", child); }

    /** \brief Configure the object (called \a once after construction
       and addition of all child \ref ConfigurableObject instances)) */
    virtual void configure();

    /// Serialize this object to a binary data stream
    virtual void serialize(Stream *stream, InstanceManager *manager) const;

    /// Return the identifier associated with this instance (or "unnamed")
    inline const std::string &getID() const { return m_properties.getID(); }

    /// Set the identifier associated with this instance
    inline void setID(const std::string &name) { m_properties.setID(name); }

    /**
     * \brief Return the properties object that was originally used to
     * create this instance
     *
     * This feature mainly of use for editors and other graphical
     * user interfaces, which present the properties of an object
     * in some form.
     */
    inline const Properties &getProperties() const { return m_properties; }

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~ConfigurableObject() { }

    /// Construct a configurable object
    inline ConfigurableObject(const Properties &props)
        : SerializableObject(), m_properties(props) { }

    /// Unserialize a configurable object
    ConfigurableObject(Stream *stream, InstanceManager *manager);
protected:
    Properties m_properties;
};

/** \brief This macro creates the binary interface, which Mitsuba
 * requires to load a plugin.
 *
 * \ingroup libcore
 */
#define MTS_EXPORT_PLUGIN(name, descr) \
    extern "C" { \
        void MTS_EXPORT *CreateInstance(const Properties &props) { \
            return new name(props); \
        } \
        const char MTS_EXPORT *GetDescription() { \
            return descr; \
        } \
    }

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_COBJECT_H_ */

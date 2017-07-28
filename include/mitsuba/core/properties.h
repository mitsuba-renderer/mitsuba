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
#if !defined(__MITSUBA_CORE_PROPERTIES_H_)
#define __MITSUBA_CORE_PROPERTIES_H_

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/transform.h>

MTS_NAMESPACE_BEGIN

struct PropertyElement;

/** \brief Associative parameter map for constructing
 * subclasses of \ref ConfigurableObject.
 *
 * Note that the Python bindings for this class do not implement
 * the various type-dependent getters and setters. Instead, they
 * are accessed just like a normal Python map, e.g:
 *
 * \code
 * myProps = mitsuba.core.Properties("pluginName")
 * myProps["stringProperty"] = "hello"
 * myProps["spectrumProperty"] = mitsuba.core.Spectrum(1.0)
 * \endcode
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE Properties {
public:
    /// Supported types of properties
    enum EPropertyType {
        /// Boolean value (true/false)
        EBoolean = 0,
        /// 64-bit signed integer
        EInteger,
        /// Floating point value
        EFloat,
        /// 3D point
        EPoint,
        /// 3D vector
        EVector,
        /// 4x4 transform for homogeneous coordinates
        ETransform,
        /// An animated 4x4 transformation
        EAnimatedTransform,
        /// Discretized color spectrum
        ESpectrum,
        /// Arbitrary-length string
        EString,
        /// Arbitrary data (pointer+size)
        EData
    };

    /// Simple pointer-size pair for passing arbitrary data (e.g. between plugins)
    struct Data {
        uint8_t *ptr;
        size_t size;

        inline bool operator==(const Data &d) const {
            return ptr == d.ptr && size == d.size;
        }

        inline bool operator!=(const Data &d) const {
            return !operator==(d);
        }
    };

    /// Construct an empty property container
    Properties();

    /// Construct an empty property container and set the plugin name
    Properties(const std::string &pluginName);

    /// Copy constructor
    Properties(const Properties &props);

    /// Release all memory
    ~Properties();

    /// Set the associated plugin name
    inline void setPluginName(const std::string &name) { m_pluginName = name; }
    /// Get the associated plugin name
    inline const std::string &getPluginName() const { return m_pluginName; }

    /// Returns the associated identifier (or the string "unnamed")
    inline const std::string &getID() const { return m_id; }
    /// Set the associated identifier
    inline void setID(const std::string &id) { m_id = id; }

    /// Set a boolean value
    void setBoolean(const std::string &name, const bool &value, bool warnDuplicates = true);
    /// Get an boolean value
    bool getBoolean(const std::string &name) const;
    /// Get an boolean value (with default);
    bool getBoolean(const std::string &name, const bool &defVal) const;

    /// Set an integer value
    void setInteger(const std::string &name, const int &value, bool warnDuplicates = true);
    /// Get an integer value
    int getInteger(const std::string &name) const;
    /// Get an integer value (with default);
    int getInteger(const std::string &name, const int &defVal) const;

    /// Set an integer value
    void setLong(const std::string &name, const int64_t &value, bool warnDuplicates = true);
    /// Get an integer value
    int64_t getLong(const std::string &name) const;
    /// Get an integer value (with default);
    int64_t getLong(const std::string &name, const int64_t &defVal) const;

    /// Set a size value
    void setSize(const std::string &name, const size_t &value, bool warnDuplicates = true);
    /// Get a size value
    size_t getSize(const std::string &name) const;
    /// Get an size value (with default);
    size_t getSize(const std::string &name, const size_t &defVal) const;

    /// Set a single precision floating point value
    void setFloat(const std::string &name, const Float &value, bool warnDuplicates = true);
    /// Get a single precision floating point value
    Float getFloat(const std::string &name) const;
    /// Get a single precision floating point value (with default)
    Float getFloat(const std::string &name, const Float &defVal) const;

    /// Set an arbitrary data value
    void setData(const std::string &name, const Data &value, bool warnDuplicates = true);
    /// Get an arbitrary data value
    Data getData(const std::string &name) const;
    /// Get an arbitrary data value (with default)
    Data getData(const std::string &name, const Data &defVal) const;

    /// Set a linear transformation
    void setTransform(const std::string &name, const Transform &value, bool warnDuplicates = true);
    /// Get a linear transformation
    Transform getTransform(const std::string &name) const;
    /// Get a linear transformation (with default)
    Transform getTransform(const std::string &name, const Transform &defVal) const;

    /// Set an animated linear transformation
    void setAnimatedTransform(const std::string &name, const AnimatedTransform *value, bool warnDuplicates = true);
    /// Get an animated linear transformation
    ref<const AnimatedTransform> getAnimatedTransform(const std::string &name) const;
    /// Get an animated linear transformation (with default)
    ref<const AnimatedTransform> getAnimatedTransform(const std::string &name, const AnimatedTransform *defVal) const;
    /// Get an animated linear transformation (with default)
    ref<const AnimatedTransform> getAnimatedTransform(const std::string &name, const Transform &defVal) const;

    /// Set a spectral power distribution
    void setSpectrum(const std::string &name, const Spectrum &value, bool warnDuplicates = true);
    /// Get a spectral power distribution
    Spectrum getSpectrum(const std::string &name) const;
    /// Get a spectral power distribution (with default)
    Spectrum getSpectrum(const std::string &name, const Spectrum &defVal) const;

    /// Set a 3d point
    void setPoint(const std::string &name, const Point &value, bool warnDuplicates = true);
    /// Get a 3d point
    Point getPoint(const std::string &name) const;
    /// Get a 3d point (with default)
    Point getPoint(const std::string &name, const Point &defVal) const;

    /// Set a 3d vector
    void setVector(const std::string &name, const Vector &value, bool warnDuplicates = true);
    /// Get a 3d vector
    Vector getVector(const std::string &name) const;
    /// Get a 3d vector (with default)
    Vector getVector(const std::string &name, const Vector &defVal) const;

    /// Set a string
    void setString(const std::string &name, const std::string &value, bool warnDuplicates = true);
    /// Get a string
    std::string getString(const std::string &name) const;
    /// Get a string (with default)
    std::string getString(const std::string &name, const std::string &defVal) const;

    /// Return one of the parameters (converting it to a string if necessary)
    std::string getAsString(const std::string &name) const;
    /// Return one of the parameters (converting it to a string if necessary, with default value)
    std::string getAsString(const std::string &name, const std::string &defVal) const;

    /// Copy an attribute from another Properties object and potentially rename it
    void copyAttribute(const Properties &properties,
        const std::string &sourceName, const std::string &targetName);

    /// Store an array containing the names of all stored properties
    void putPropertyNames(std::vector<std::string> &results) const;

    /// Return an array containing the names of all stored properties
    inline std::vector<std::string> getPropertyNames() const {
        std::vector<std::string> results;
        putPropertyNames(results);
        return results;
    }

    /// Manually mark a certain property as queried
    void markQueried(const std::string &name) const;

    /// Check if a certain property was queried
    bool wasQueried(const std::string &name) const;

    /// Verify if a value with the specified name exists
    bool hasProperty(const std::string &name) const;

    /**
     * \brief Remove a property with the specified name
     * \return \c true upon success
     */
    bool removeProperty(const std::string &name);

    /// Return the property of a type
    EPropertyType getType(const std::string &name) const;

    /// Return the list of un-queried attributed
    std::vector<std::string> getUnqueried() const;

    /// Assignment operator
    void operator=(const Properties &props);

    /// Equality comparison operator
    bool operator==(const Properties &props) const;

    /// Inequality comparision operator
    inline bool operator!=(const Properties &props) const {
        return !operator==(props);
    }

    /// Merge a properties record into the current one
    void merge(const Properties &props);

    /// Return a string representation
    std::string toString() const;
private:
    std::map<std::string, PropertyElement> *m_elements;
    std::string m_pluginName, m_id;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_PROPERTIES_H_ */

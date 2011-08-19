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

#if !defined(__PROPERTIES_H)
#define __PROPERTIES_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/transform.h>
#include <boost/variant.hpp>

MTS_NAMESPACE_BEGIN

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
	};

	/// Construct an empty property container
	Properties() : m_id("unnamed") { }

	/// Construct an empty property container and set the plugin name
	Properties(const std::string &pluginName) : m_pluginName(pluginName), m_id("unnamed") { }

	/// Set the associated plugin name
	inline void setPluginName(const std::string &name) { m_pluginName = name; }
	/// Get the associated plugin name
	inline const std::string &getPluginName() const { return m_pluginName; }
	
	/// Returns the associated ID (or the string "unnamed")
	inline const std::string &getID() const { return m_id; }
	/// Set the associated ID
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

	/// Store an array containing the names of all stored properties
	inline void putNames(std::vector<std::string> &results) const {
		for (std::map<std::string, Element>::const_iterator it = m_elements.begin();
			it != m_elements.end(); ++it) 
			results.push_back((*it).first);
	}
	
	/// Return an array containing the names of all stored properties
	inline std::vector<std::string> getNames() const {
		std::vector<std::string> results;
		putNames(results);
		return results;
	}

	/// Manually mark a certain property as queried
	void markQueried(const std::string &name) const;

	/// Check if a certain property was queried
	bool wasQueried(const std::string &name) const;

	/// Verify if a value with the specified name exists
	bool hasProperty(const std::string &name) const;

	/// Return the property of a type
	EPropertyType getType(const std::string &name) const;

	/// Return the list of un-queried attributed
	std::vector<std::string> getUnqueried() const;

	/// Return a string representation
	std::string toString() const;
private:
	/// \cond
	typedef boost::variant<
		bool, int64_t, Float, Point, Vector, Transform,
		Spectrum, std::string, Data> ElementData;

	struct Element {
		ElementData data;
		mutable bool queried;
	};
	/// \endcond

	std::map<std::string, Element> m_elements;
	std::string m_pluginName, m_id;
};

MTS_NAMESPACE_END

#endif /* __PROPERTIES_H */

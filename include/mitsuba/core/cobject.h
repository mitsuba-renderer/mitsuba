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

#if !defined(__COBJECT_H)
#define __COBJECT_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/serialization.h>

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
	/// Set the parent object
	virtual void setParent(ConfigurableObject *parent);

	/// Return the parent object
	inline ConfigurableObject *getParent() { return m_parent; }

	/// Return the parent object (const version)
	inline const ConfigurableObject *getParent() const { return m_parent; }

	/// Add a child (default implementation throws an error)
	virtual void addChild(const std::string &name, ConfigurableObject *child);

	/// Add an unnamed child
	inline void addChild(ConfigurableObject *child) { addChild("", child); }

	/** \brief Configure the object (called \a once after construction
	   and addition of all child ConfigurableObjects) */
	virtual void configure();

	/// Serialize this object to a binary data stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~ConfigurableObject() { }
	
	/// Construct a configurable object
	inline ConfigurableObject(const Properties &props) 
		: SerializableObject(), m_parent(NULL) { }
	
	/// Unserialize a configurable object
	ConfigurableObject(Stream *stream, InstanceManager *manager);
protected:
	ConfigurableObject *m_parent;
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

#endif /* __COBJECT_H */

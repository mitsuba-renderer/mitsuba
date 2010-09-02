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

#if !defined(__COBJECT_H)
#define __COBJECT_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/serialization.h>

MTS_NAMESPACE_BEGIN

/** \brief Interface of a generic serializable object that can be 
 * configured by an instance of the Properties class.
 */
class MTS_EXPORT_CORE ConfigurableObject : public SerializableObject {
public:
	/// Constructor
	inline ConfigurableObject(const Properties &props) : SerializableObject(), m_configured(false), m_parent(NULL) { }

	/// Unserialize a configurable object
	ConfigurableObject(Stream *stream, InstanceManager *manager);

	/// Set the parent object
	virtual void setParent(ConfigurableObject *parent);

	/// Get the parent object
	inline ConfigurableObject *getParent() { return m_parent; }
	
	/// Get the parent object
	inline const ConfigurableObject *getParent() const { return m_parent; }

	/// Add a child (default implementation throws an error)
	virtual void addChild(const std::string &name, ConfigurableObject *child);

	/** \brief Configure the object (called _once_ after construction
	   and addition of all child ConfigurableObjects. */
	virtual void configure();

	/// Serialize this object to a stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~ConfigurableObject() { }
protected:
	bool m_configured;
	ConfigurableObject *m_parent;
};

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

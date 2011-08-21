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

#include <mitsuba/core/properties.h>
#include <mitsuba/core/netobject.h>

MTS_NAMESPACE_BEGIN

#define DEFINE_PROPERTY_ACCESSOR(Type, BaseType, TypeName, ReadableName) \
	void Properties::set##TypeName(const std::string &name, const Type &value, bool warnDuplicates) { \
		if (hasProperty(name) && warnDuplicates) \
			SLog(EWarn, "Property \"%s\" has already been specified!", name.c_str()); \
		m_elements[name].data = (BaseType) value; \
		m_elements[name].queried = false; \
	} \
	\
	Type Properties::get##TypeName(const std::string &name) const { \
		std::map<std::string, Element>::const_iterator it = m_elements.find(name); \
		if (it == m_elements.end()) \
			SLog(EError, "Property \"%s\" missing", name.c_str()); \
		const BaseType *result = boost::get<BaseType>(&it->second.data); \
		if (!result) \
			SLog(EError, "The property \"%s\" has the wrong type (expected <" #ReadableName ">). The " \
					"complete property record is :\n%s", name.c_str(), toString().c_str()); \
		it->second.queried = true; \
		return (Type) *result; \
	} \
	\
	Type Properties::get##TypeName(const std::string &name, const Type &defVal) const { \
		std::map<std::string, Element>::const_iterator it = m_elements.find(name); \
		if (it == m_elements.end()) \
			return defVal; \
		const BaseType *result = boost::get<BaseType>(&it->second.data); \
		if (!result) \
			SLog(EError, "The property \"%s\" has the wrong type (expected <" #ReadableName ">). The " \
					"complete property record is :\n%s", name.c_str(), toString().c_str()); \
		it->second.queried = true; \
		return (Type) *result; \
	}

DEFINE_PROPERTY_ACCESSOR(bool, bool, Boolean, bool)
DEFINE_PROPERTY_ACCESSOR(int64_t, int64_t, Long, integer)
DEFINE_PROPERTY_ACCESSOR(int, int64_t, Integer, integer)
DEFINE_PROPERTY_ACCESSOR(size_t, int64_t, Size, integer)
DEFINE_PROPERTY_ACCESSOR(Float, Float, Float, float)
DEFINE_PROPERTY_ACCESSOR(Point, Point, Point, point)
DEFINE_PROPERTY_ACCESSOR(Vector, Vector, Vector, vector)
DEFINE_PROPERTY_ACCESSOR(Transform, Transform, Transform, transform)
DEFINE_PROPERTY_ACCESSOR(Spectrum, Spectrum, Spectrum, spectrum)
DEFINE_PROPERTY_ACCESSOR(std::string, std::string, String, string)
DEFINE_PROPERTY_ACCESSOR(Properties::Data, Properties::Data, Data, data)

class type_visitor : public boost::static_visitor<Properties::EPropertyType> {
public:
	Properties::EPropertyType operator()(const bool &) const             { return Properties::EBoolean; }
	Properties::EPropertyType operator()(const int64_t &) const          { return Properties::EInteger; }
	Properties::EPropertyType operator()(const Float &) const            { return Properties::EFloat; }
	Properties::EPropertyType operator()(const Point &) const            { return Properties::EPoint; }
	Properties::EPropertyType operator()(const Vector &) const           { return Properties::EVector; }
	Properties::EPropertyType operator()(const Transform &) const        { return Properties::ETransform; }
	Properties::EPropertyType operator()(const Spectrum &) const         { return Properties::ESpectrum; }
	Properties::EPropertyType operator()(const std::string &) const      { return Properties::EString; }
	Properties::EPropertyType operator()(const Properties::Data &) const { return Properties::EData; }
};

bool Properties::hasProperty(const std::string &name) const {
	return m_elements.find(name) != m_elements.end();
}

std::vector<std::string> Properties::getUnqueried() const {
	std::map<std::string, Element>::const_iterator it = m_elements.begin();
	std::vector<std::string> result;

	for (; it != m_elements.end(); ++it) {
		if (!(*it).second.queried)
			result.push_back((*it).first);
	}

	return result;
}

Properties::EPropertyType Properties::getType(const std::string &name) const {
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if (it == m_elements.end())
		SLog(EError, "Property \"%s\" has not been specified!", name.c_str());
	
	type_visitor myVisitor;
	return boost::apply_visitor(myVisitor, it->second.data);
}

std::string Properties::toString() const {
	std::map<std::string, Element>::const_iterator it = m_elements.begin();
	std::ostringstream oss;

	oss << "Properties[" << endl
		<< "  pluginName = \"" << m_pluginName << "\"," << endl
		<< "  elements = {" << endl;
	while (it != m_elements.end()) {
		oss << "    \"" << (*it).first << "\" -> ";
		const ElementData &data = (*it).second.data;
		EPropertyType type = boost::apply_visitor(type_visitor(), data);
		switch (type) {
			case EBoolean:
				oss << (boost::get<bool>(data) ? "true" : "false");
				break;
			case EInteger:
				oss << boost::get<int64_t>(data);
				break;
			case EFloat:
				oss << boost::get<Float>(data);
				break;
			case EPoint:
				oss << boost::get<Point>(data).toString();
				break;
			case ETransform:
				oss << boost::get<Transform>(data).toString();
				break;
			case ESpectrum:
				oss << boost::get<Spectrum>(data).toString();
				break;
			case EString:
				oss << "\"" << boost::get<std::string>(data) << "\"";
				break;
			case EData:
				oss << boost::get<Data>(data).ptr << " (size=" 
					<< boost::get<Data>(data).size << ")";
				break;
			default:
				SLog(EError, "Encountered an unknown property type!");
		}
		if (++it != m_elements.end())
			oss << ",";
		oss << endl;
	}
	oss << "  }" << endl
		<< "]" << endl;
	return oss.str();
}

void Properties::markQueried(const std::string &name) const {
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if (it == m_elements.end())
		return;
	it->second.queried = true;
}

bool Properties::wasQueried(const std::string &name) const {
	std::map<std::string, Element>::const_iterator it = m_elements.find(name);
	if (it == m_elements.end())
		SLog(EError, "Could not find parameter \"%s\"!", name.c_str());
	return it->second.queried;
}

ConfigurableObject::ConfigurableObject(Stream *stream, InstanceManager *manager) 
 : SerializableObject(stream, manager) {
	m_parent = static_cast<ConfigurableObject *>(manager->getInstance(stream));
}

void ConfigurableObject::setParent(ConfigurableObject *parent) {
	m_parent = parent;
}

void ConfigurableObject::configure() {
}

void ConfigurableObject::serialize(Stream *stream, InstanceManager *manager) const {
	if (!getClass()->isSerializable())
		Log(EError, "Error: trying to serialize an instance of type '%s', which does "
			"not have full serialization support!", getClass()->getName().c_str());
	manager->serialize(stream, m_parent);
}

void ConfigurableObject::addChild(const std::string &name, ConfigurableObject *child) {
	SLog(EError, "ConfigurableObject::addChild(\"%s\", %s) not implemented in \"%s\"", 
		name.c_str(), child->toString().c_str(), toString().c_str());
}

void NetworkedObject::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
}
	
void NetworkedObject::bindUsedResources(ParallelProcess *proc) const {
}

void NetworkedObject::wakeup(std::map<std::string, SerializableObject *> &params) {
}

MTS_IMPLEMENT_CLASS(ConfigurableObject, true, SerializableObject)
MTS_IMPLEMENT_CLASS(NetworkedObject, true, ConfigurableObject)
MTS_NAMESPACE_END

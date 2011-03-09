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

#include <mitsuba/core/properties.h>
#include <mitsuba/render/luminaire.h>
#include <mitsuba/render/medium.h>

MTS_NAMESPACE_BEGIN

Luminaire::Luminaire(const Properties &props)
 : ConfigurableObject(props) {
	m_type = 0;

	// Transformation from the luminaire's local coordinates to world coordiantes
	m_luminaireToWorld = props.getTransform("toWorld", Transform());

	// Importance sampling weight (used by the luminaire sampling code in \ref Scene)
	m_samplingWeight = props.getFloat("samplingWeight", 1.0f);

	m_name = props.getID();

	m_worldToLuminaire = m_luminaireToWorld.inverse();
	AssertEx(!m_luminaireToWorld.hasScale(), "toWorld transformation can't have scale factors!");
	m_intersectable = false;
}

Luminaire::Luminaire(Stream *stream, InstanceManager *manager)
 : ConfigurableObject(stream, manager) {
	m_medium = static_cast<Medium *>(manager->getInstance(stream));
	m_samplingWeight = stream->readFloat();
	m_type = (EType) stream->readInt();
	m_intersectable = stream->readBool();
	m_worldToLuminaire = Transform(stream);
	m_name = stream->readString();
	m_luminaireToWorld = m_worldToLuminaire.inverse();
}

Luminaire::~Luminaire() {
}

void Luminaire::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
	manager->serialize(stream, m_medium.get());
	stream->writeFloat(m_samplingWeight);
	stream->writeInt(m_type);
	stream->writeBool(m_intersectable);
	m_worldToLuminaire.serialize(stream);
	stream->writeString(m_name);
}

bool Luminaire::createEmissionRecord(EmissionRecord &eRec, const Ray &ray) const {
	Log(EError, "createEmissionRecord(): Not implemented!");
	return false;
}

void Luminaire::preprocess(const Scene *scene) {
}

bool Luminaire::isBackgroundLuminaire() const {
	return false;
}

Spectrum Luminaire::Le(const Ray &ray) const {
	Log(EError, "Luminaire::Le(const Ray &) is not implemented!");
	return Spectrum(0.0f);
}
	
Spectrum Luminaire::Le(const ShapeSamplingRecord &sRec, const Vector &d) const {
	Log(EError, "Luminaire::Le(const ShapeSamplingRecord sRec&, const Vector &) is not implemented!");
	return Spectrum(0.0f);
}

std::string EmissionRecord::toString() const {
	std::ostringstream oss;
	oss << "EmissionRecord[" << std::endl
		<< "  sRec = " << indent(sRec.toString()) << "," << std::endl
		<< "  d = " << d.toString() << "," << std::endl
		<< "  value = " << value.toString() << "," << std::endl
		<< "  pdfArea = " << pdfArea << "," << std::endl
		<< "  pdfDir = " << pdfDir << "," << std::endl
		<< "  luminaire = " << ((luminaire == NULL ) ? "null" : indent(((Object *) luminaire)->toString()).c_str()) << std::endl
		<< "]";
	return oss.str();
}

std::string LuminaireSamplingRecord::toString() const {
	std::ostringstream oss;
	oss << "LuminaireSamplingRecord[" << std::endl
		<< "  sRec = " << indent(sRec.toString()) << "," << std::endl
		<< "  d = " << d.toString() << "," << std::endl
		<< "  pdf = " << pdf << "," << std::endl
		<< "  value = " << value.toString() << "," << std::endl
		<< "  luminaire = " << ((luminaire == NULL ) ? "null" : indent(((Object *) luminaire)->toString()).c_str()) << std::endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(Luminaire, true, ConfigurableObject)
MTS_NAMESPACE_END

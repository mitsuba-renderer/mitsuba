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

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/luminaire.h>
#include <mitsuba/render/medium.h>

MTS_NAMESPACE_BEGIN

Shape::Shape(const Properties &props) 
 : ConfigurableObject(props), m_occluder(false) { }

Shape::Shape(Stream *stream, InstanceManager *manager) 
 : ConfigurableObject(stream, manager) {
	m_bsdf = static_cast<BSDF *>(manager->getInstance(stream));
	m_subsurface = static_cast<Subsurface *>(manager->getInstance(stream));
	m_luminaire = static_cast<Luminaire *>(manager->getInstance(stream));
	m_interiorMedium = static_cast<Medium *>(manager->getInstance(stream));
	m_exteriorMedium = static_cast<Medium *>(manager->getInstance(stream));
	m_occluder = stream->readBool();
}

Shape::~Shape() { }


void Shape::configure() {
	if ((hasSubsurface() || isLuminaire()) && m_bsdf == NULL) {
		/* Light source & no BSDF -> set an all-absorbing BSDF to turn
		   the shape into an occluder. This is needed for the path
		   tracer implementation to work correctly. */
		Properties props("diffuse");
		props.setSpectrum("reflectance", Spectrum(0.0f));
		addChild(static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), props)));
	}
}
	
bool Shape::isCompound() const {
	return false;
}
	
std::string Shape::getName() const {
	return "Unnamed";
}

Shape *Shape::getElement(int i) {
	return NULL;
}
	
AABB Shape::getClippedAABB(const AABB &box) const {
	AABB result = getAABB();
	result.clip(box);
	return result;
}

Float Shape::sampleSolidAngle(ShapeSamplingRecord &sRec, 
		const Point &from, const Point2 &sample) const {
	/* Turns the area sampling routine into one that samples wrt. solid angles */
	Float pdfArea = sampleArea(sRec, sample);
	Vector lumToPoint = from - sRec.p;
	Float distSquared = lumToPoint.lengthSquared(), dp = dot(lumToPoint, sRec.n);
	if (dp > 0)
		return pdfArea * distSquared * std::sqrt(distSquared) / dp;
	else
		return 0.0f;
}

Float Shape::pdfSolidAngle(const ShapeSamplingRecord &sRec, const Point &from) const {
	/* Turns the area sampling routine into one that samples wrt. solid angles */
	Vector lumToPoint = from - sRec.p;
	Float distSquared = lumToPoint.lengthSquared();
	Float invDP = std::max((Float) 0, std::sqrt(distSquared) / dot(lumToPoint, sRec.n));
	return pdfArea(sRec) * distSquared * invDP;
}

void Shape::addChild(const std::string &name, ConfigurableObject *child) {
	const Class *cClass = child->getClass();
	if (cClass->derivesFrom(MTS_CLASS(BSDF))) {
		m_bsdf = static_cast<BSDF *>(child);
		m_occluder = true;
	} else if (cClass->derivesFrom(MTS_CLASS(Luminaire))) {
		Assert(m_luminaire == NULL);
		m_luminaire = static_cast<Luminaire *>(child);
		if (m_luminaire && m_exteriorMedium)
			m_luminaire->setMedium(m_exteriorMedium);
	} else if (cClass->derivesFrom(MTS_CLASS(Subsurface))) {
		Assert(m_subsurface == NULL);
		if (m_interiorMedium != NULL)
			Log(EError, "Shape \"%s\" has both an interior medium "
				"and a subsurface integrator -- please choose one or the other!", getName().c_str());
		m_subsurface = static_cast<Subsurface *>(child);
	} else if (cClass->derivesFrom(MTS_CLASS(Medium))) {
		if (name == "interior") {
			Assert(m_interiorMedium == NULL);
			if (m_subsurface != NULL)
				Log(EError, "Shape \"%s\" has both an interior medium "
					"and a subsurface integrator -- please choose one or the other!", getName().c_str());
			m_interiorMedium = static_cast<Medium *>(child);
		} else if (name == "exterior") {
			Assert(m_exteriorMedium == NULL);
			m_exteriorMedium = static_cast<Medium *>(child);
			if (m_luminaire)
				m_luminaire->setMedium(m_exteriorMedium);
		} else {
			Log(EError, "Shape: Invalid medium child (must be named "
				"'interior' or 'exterior')!");
		}
	} else {
		ConfigurableObject::addChild(name, child);
	}
}

const KDTreeBase<AABB> *Shape::getKDTree() const {
	return NULL;
}

void Shape::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
	manager->serialize(stream, m_bsdf.get());
	manager->serialize(stream, m_subsurface.get());
	manager->serialize(stream, m_luminaire.get());
	manager->serialize(stream, m_interiorMedium.get());
	manager->serialize(stream, m_exteriorMedium.get());
	stream->writeBool(m_occluder);
}
	
Float Shape::getSurfaceArea() const {
	Log(EError, "%s::getSurfaceArea(): Not implemented!",
			getClass()->getName().c_str());
	return 0.0f;
}

bool Shape::rayIntersect(const Ray &ray, Float mint, 
		Float maxt, Float &t, void *temp) const {
	Log(EError, "%s::rayIntersect(): Not implemented!",
			getClass()->getName().c_str());
	return false;
}

bool Shape::rayIntersect(const Ray &ray, Float mint, 
		Float maxt) const {
	Log(EError, "%s::rayIntersect(): Not implemented!",
			getClass()->getName().c_str());
	return false;
}


void Shape::fillIntersectionRecord(const Ray &ray, 
		const void *temp, Intersection &its) const {
	Log(EError, "%s::fillIntersectionRecord(): Not implemented!",
			getClass()->getName().c_str());
}

Float Shape::sampleArea(ShapeSamplingRecord &sRec, 
		const Point2 &sample) const {
	Log(EError, "%s::sampleArea(): Not implemented!",
			getClass()->getName().c_str());
	return 0.0f;
}

Float Shape::pdfArea(const ShapeSamplingRecord &sRec) const {
	Log(EError, "%s::pdfArea(): Not implemented!",
			getClass()->getName().c_str());
	return 0.0f;
}

ref<TriMesh> Shape::createTriMesh() {
	return NULL;
}

std::string ShapeSamplingRecord::toString() const {
	std::ostringstream oss;
	oss << "ShapeSamplingRecord[" << std::endl
		<< "  p = " << p.toString() << "," << std::endl
		<< "  n = " << n.toString() << std::endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(Shape, true, ConfigurableObject)
MTS_NAMESPACE_END

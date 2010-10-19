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

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/luminaire.h>

MTS_NAMESPACE_BEGIN

Shape::Shape(const Properties &props) 
 : ConfigurableObject(props) { }

Shape::Shape(Stream *stream, InstanceManager *manager) 
 : ConfigurableObject(stream, manager) {
	m_bsdf = static_cast<BSDF *>(manager->getInstance(stream));
	m_subsurface = static_cast<Subsurface *>(manager->getInstance(stream));
	m_luminaire = static_cast<Luminaire *>(manager->getInstance(stream));
}

Shape::~Shape() {
}


void Shape::configure() {
	/* Ensure that there is at least some default BSDF */
	if (m_bsdf == NULL) {
		m_bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(BSDF::m_theClass, Properties("lambertian")));
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
	if (dp > 0) {
		return pdfArea * distSquared * std::sqrt(distSquared) / dp;
	} else {
		return 0.0f;
	}
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
	if (cClass->derivesFrom(BSDF::m_theClass)) {
		m_bsdf = static_cast<BSDF *>(child);
	} else if (cClass->derivesFrom(Luminaire::m_theClass)) {
		Assert(m_luminaire == NULL);
		m_luminaire = static_cast<Luminaire *>(child);
	} else if (cClass->derivesFrom(Subsurface::m_theClass)) {
		Assert(m_subsurface == NULL);
		m_subsurface = static_cast<Subsurface *>(child);
	} else {
		Log(EError, "Shape: Invalid child node!");
	}
}

void Shape::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
	manager->serialize(stream, m_bsdf.get());
	manager->serialize(stream, m_subsurface.get());
	manager->serialize(stream, m_luminaire.get());
}

bool Shape::rayIntersect(const Ray &ray, Float mint, 
		Float maxt, Float &t, void *temp) const {
	Log(EError, "%s::rayIntersect(): Not implemented!",
			getClass()->getName().c_str());
	return false;
}

void Shape::fillIntersectionRecord(const Ray &ray, Float t, 
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

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

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sensor.h>

MTS_NAMESPACE_BEGIN

Shape::Shape(const Properties &props)
 : ConfigurableObject(props) {
	m_name = props.getID();
}

Shape::Shape(Stream *stream, InstanceManager *manager)
 : ConfigurableObject(stream, manager) {
	m_name = stream->readString();
	m_bsdf = static_cast<BSDF *>(manager->getInstance(stream));
	m_subsurface = static_cast<Subsurface *>(manager->getInstance(stream));
	m_emitter = static_cast<Emitter *>(manager->getInstance(stream));
	m_sensor = static_cast<Sensor *>(manager->getInstance(stream));
	m_interiorMedium = static_cast<Medium *>(manager->getInstance(stream));
	m_exteriorMedium = static_cast<Medium *>(manager->getInstance(stream));
}

Shape::~Shape() { }

void Shape::configure() {
	if (m_bsdf == NULL) {
		ref<BSDF> bsdf = NULL;
		if (isEmitter() || isSensor() || hasSubsurface()) {
			/* Light source / sensor and no BSDF! -> set an all-absorbing BSDF */
			Properties props("diffuse");
			props.setSpectrum("reflectance", Spectrum(0.0f));
			bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(BSDF), props));
		} else if (!isMediumTransition()) {
			/* A surface without BSDF, which is not a medium
			   transition/sensor/emitter/subsurface emitter doesn't make
			   much sense. Assign it a 0.5 Lambertian BRDF for convenience */
			Properties props("diffuse");
			props.setSpectrum("reflectance", Spectrum(0.5f));
			bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(BSDF), props));
		} else {
			/* Assign a "null" BSDF */
			bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(BSDF), Properties("null")));
		}
		bsdf->configure();
		addChild(bsdf);
	}

	if ((m_bsdf->getType() & BSDF::ENull) && (isEmitter() || isSensor() || hasSubsurface()))
		Log(EError, "Shape \"%s\" has an index-matched BSDF and an "
			"emitter/sensor/subsurface attachment. This is not allowed!", getName().c_str());
}

void Shape::adjustTime(Intersection &its, Float time) const {
	its.time = time;
	/* Do nothing else by default */
}

bool Shape::isCompound() const {
	return false;
}

std::string Shape::getName() const {
	return m_name;
}

Shape *Shape::getElement(int i) {
	return NULL;
}

AABB Shape::getClippedAABB(const AABB &box) const {
	AABB result = getAABB();
	result.clip(box);
	return result;
}

void Shape::sampleDirect(DirectSamplingRecord &dRec,
			const Point2 &sample) const {
	/* Piggyback on sampleArea() */
	samplePosition(dRec, sample);

	dRec.d = dRec.p - dRec.ref;

	Float distSquared = dRec.d.lengthSquared();
	dRec.dist = std::sqrt(distSquared);
	dRec.d /= dRec.dist;
	Float dp = absDot(dRec.d, dRec.n);
	dRec.pdf *= dp != 0 ? (distSquared / dp) : 0.0f;
	dRec.measure = ESolidAngle;
}

Float Shape::pdfDirect(const DirectSamplingRecord &dRec) const {
	Float pdfPos = pdfPosition(dRec);

	if (dRec.measure == ESolidAngle)
		return pdfPos * (dRec.dist * dRec.dist) / absDot(dRec.d, dRec.n);
	else if (dRec.measure == EArea)
		return pdfPos;
	else
		return 0.0f;
}

void Shape::addChild(const std::string &name, ConfigurableObject *child) {
	const Class *cClass = child->getClass();
	if (cClass->derivesFrom(MTS_CLASS(BSDF))) {
		m_bsdf = static_cast<BSDF *>(child);
	} else if (cClass->derivesFrom(MTS_CLASS(Emitter))) {
		Emitter *emitter = static_cast<Emitter *>(child);
		if (m_emitter != NULL)
			Log(EError, "Tried to attach multiple emitters to a shape!");
		if (emitter) {
			if (!emitter->isOnSurface())
				Log(EError, "Tried to attach an incompatible emitter to a surface!");
			if (m_exteriorMedium)
				emitter->setMedium(m_exteriorMedium);
		}
		m_emitter = emitter;
	} else if (cClass->derivesFrom(MTS_CLASS(Sensor))) {
		Sensor *sensor = static_cast<Sensor *>(child);
		if (m_sensor != NULL)
			Log(EError, "Tried to attach multiple sensors to a shape!");
		if (sensor) {
			if (!sensor->isOnSurface())
				Log(EError, "Tried to attach an incompatible sensor to a surface!");
			if (m_exteriorMedium)
				sensor->setMedium(m_exteriorMedium);
		}
		m_sensor = sensor;
	} else if (cClass->derivesFrom(MTS_CLASS(Subsurface))) {
		Assert(m_subsurface == NULL);
		if (m_interiorMedium != NULL)
			Log(EError, "Shape \"%s\" has both an interior medium "
				"and a subsurface scattering model -- please choose one or the other!", getName().c_str());
		m_subsurface = static_cast<Subsurface *>(child);
	} else if (cClass->derivesFrom(MTS_CLASS(Medium))) {
		if (name == "interior") {
			Assert(m_interiorMedium == NULL || m_interiorMedium == child);
			if (m_subsurface != NULL)
				Log(EError, "Shape \"%s\" has both an interior medium "
					"and a subsurface scattering model -- please choose one or the other!", getName().c_str());
			m_interiorMedium = static_cast<Medium *>(child);
		} else if (name == "exterior") {
			Assert(m_exteriorMedium == NULL || m_exteriorMedium == child);
			m_exteriorMedium = static_cast<Medium *>(child);
			if (m_emitter)
				m_emitter->setMedium(m_exteriorMedium);
			if (m_sensor)
				m_sensor->setMedium(m_exteriorMedium);
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
	stream->writeString(m_name);
	manager->serialize(stream, m_bsdf.get());
	manager->serialize(stream, m_subsurface.get());
	manager->serialize(stream, m_emitter.get());
	manager->serialize(stream, m_sensor.get());
	manager->serialize(stream, m_interiorMedium.get());
	manager->serialize(stream, m_exteriorMedium.get());
}

Float Shape::getSurfaceArea() const { NotImplementedError("getSurfaceArea"); }
bool Shape::rayIntersect(const Ray &ray, Float mint,
		Float maxt, Float &t, void *temp) const { NotImplementedError("rayIntersect"); }
bool Shape::rayIntersect(const Ray &ray, Float mint,
		Float maxt) const { NotImplementedError("rayIntersect"); }

void Shape::fillIntersectionRecord(const Ray &ray,
		const void *temp, Intersection &its) const {
	NotImplementedError("fillIntersectionRecord"); }

void Shape::getCurvature(const Intersection &its, Float &H, Float &K,
		bool shadingFrame) const {
	Vector dndu, dndv;
	getNormalDerivative(its, dndu, dndv, shadingFrame);

	/* Compute the coefficients of the first and second fundamental form */
	Float
		E =  dot(its.dpdu, its.dpdu),
		F =  dot(its.dpdu, its.dpdv),
		G =  dot(its.dpdv, its.dpdv),
		e = -dot(its.dpdu, dndu),
		f = -dot(its.dpdv, dndu),
		g = -dot(its.dpdv, dndv),
		invDenom = 1.0f / (E*G - F*F);

	K = (e*g - f*f) * invDenom;
	H = .5f*(e*G - 2*f*F + g*E) * invDenom;
}

void Shape::getNormalDerivative(const Intersection &its,
		Vector &dndu, Vector &dndv, bool shadingFrame) const {
	NotImplementedError("getNormalDerivative");
}

void Shape::samplePosition(PositionSamplingRecord &pRec, const Point2 &sample) const {
	NotImplementedError("samplePosition");
}

Float Shape::pdfPosition(const PositionSamplingRecord &pRec) const {
	NotImplementedError("pdfPosition");
}

void Shape::copyAttachments(Shape *shape) {
	m_bsdf = shape->getBSDF();
	m_emitter = shape->getEmitter();
	m_sensor = shape->getSensor();
	m_subsurface = shape->getSubsurface();
	m_interiorMedium = shape->getInteriorMedium();
	m_exteriorMedium = shape->getInteriorMedium();
}

ref<TriMesh> Shape::createTriMesh() {
	return NULL;
}

MTS_IMPLEMENT_CLASS(Shape, true, ConfigurableObject)
MTS_NAMESPACE_END

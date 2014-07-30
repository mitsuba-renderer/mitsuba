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

#include <mitsuba/render/sensor.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/core/track.h>
#include <mitsuba/core/plugin.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

Sensor::Sensor(const Properties &props)
 : AbstractEmitter(props) {
	m_shutterOpen = props.getFloat("shutterOpen", 0.0f);
	Float shutterClose = props.getFloat("shutterClose", 0.0f);
	m_shutterOpenTime = shutterClose - m_shutterOpen;

	if (m_shutterOpenTime < 0)
		Log(EError, "Shutter opening time must be less than "
			         "or equal to the shutter closing time!");

	if (m_shutterOpenTime == 0)
		m_type |= EDeltaTime;
}

Sensor::Sensor(Stream *stream, InstanceManager *manager)
 : AbstractEmitter(stream, manager) {
	m_film = static_cast<Film *>(manager->getInstance(stream));
	m_sampler = static_cast<Sampler *>(manager->getInstance(stream));
	m_shutterOpen = stream->readFloat();
	m_shutterOpenTime = stream->readFloat();
}

Sensor::~Sensor() {
}

void Sensor::serialize(Stream *stream, InstanceManager *manager) const {
	AbstractEmitter::serialize(stream, manager);
	manager->serialize(stream, m_film.get());
	manager->serialize(stream, m_sampler.get());
	stream->writeFloat(m_shutterOpen);
	stream->writeFloat(m_shutterOpenTime);
}

void Sensor::setShutterOpenTime(Float time) {
	m_shutterOpenTime = time;
	if (m_shutterOpenTime == 0)
		m_type |= EDeltaTime;
	else
		m_type &= ~EDeltaTime;
}


Spectrum Sensor::eval(const Intersection &its, const Vector &d,
		Point2 &samplePos) const {
	Log(EError, "%s::eval(const Intersection &, const Vector &, "
		"Point2&) is not implemented!", getClass()->getName().c_str());
	return Spectrum(0.0f);
}

bool Sensor::getSamplePosition(const PositionSamplingRecord &pRec,
		const DirectionSamplingRecord &dRec, Point2 &samplePosition) const {
	Log(EError, "%s::getSamplePosition(const PositionSamplingRecord &, "
		"const DirectionSamplingRecord &, Point2&) is not implemented!",
		getClass()->getName().c_str());
	return false;
}

void Sensor::configure() {
	if (m_film == NULL) {
		/* Instantiate an EXR film by default */
		m_film = static_cast<Film*> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Film), Properties("hdrfilm")));
		m_film->configure();
	}

	if (m_sampler == NULL) {
		/* No sampler has been selected - load an independent filter with 4 samples/pixel by default */
		Properties props("independent");
		props.setInteger("sampleCount", 4);
		m_sampler = static_cast<Sampler *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Sampler), props));
		m_sampler->configure();
	}

	m_aspect = m_film->getSize().x /
	   (Float) m_film->getSize().y;

	m_resolution = Vector2(m_film->getCropSize());
	m_invResolution = Vector2(
		(Float) 1 / m_resolution.x,
		(Float) 1 / m_resolution.y);
}

Spectrum Sensor::sampleRayDifferential(RayDifferential &ray,
		const Point2 &samplePosition,
		const Point2 &apertureSample,
		Float timeSample) const {
	Spectrum result = sampleRay(ray, samplePosition,
		apertureSample, timeSample);

	/* Sample a ray for X+1 */
	Ray tempRay;
	sampleRay(tempRay, samplePosition + Vector2(1, 0),
			apertureSample, timeSample);
	ray.rxOrigin = tempRay.o;
	ray.rxDirection = tempRay.d;

	/* Sample a ray for Y+1 */
	sampleRay(tempRay, samplePosition + Vector2(0, 1),
			apertureSample, timeSample);
	ray.ryOrigin = tempRay.o;
	ray.ryDirection = tempRay.d;
	ray.hasDifferentials = true;

	return result;
}

Float Sensor::pdfTime(const Ray &ray, EMeasure measure) const {
	if (ray.time < m_shutterOpen || ray.time > m_shutterOpenTime + m_shutterOpenTime)
		return 0.0f;

	if (m_shutterOpenTime == 0 && measure == EDiscrete)
		return 1.0f;
	else if (m_shutterOpenTime > 0 && measure == ELength)
		return 1.0f / m_shutterOpenTime;
	else
		return 0.0f;
}

void Sensor::addChild(const std::string &name, ConfigurableObject *child) {
	if (child->getClass()->derivesFrom(MTS_CLASS(Sampler))) {
		m_sampler = static_cast<Sampler *>(child);
	} else if (child->getClass()->derivesFrom(MTS_CLASS(Film))) {
		m_film = static_cast<Film *>(child);
	} else {
		AbstractEmitter::addChild(name, child);
	}
}

ProjectiveCamera::ProjectiveCamera(const Properties &props) : Sensor(props) {
	/* Distance to the near clipping plane */
	m_nearClip = props.getFloat("nearClip", 1e-2f);
	/* Distance to the far clipping plane */
	m_farClip = props.getFloat("farClip", 1e4f);
	/* Distance to the focal plane */
	m_focusDistance = props.getFloat("focusDistance", m_farClip);

	if (m_nearClip <= 0)
		Log(EError, "The 'nearClip' parameter must be greater than zero!");
	if (m_nearClip >= m_farClip)
		Log(EError, "The 'nearClip' parameter must be smaller than 'farClip'.");

	m_type |= EProjectiveCamera;
}

ProjectiveCamera::ProjectiveCamera(Stream *stream, InstanceManager *manager)
	: Sensor(stream, manager) {
	m_nearClip = stream->readFloat();
	m_farClip = stream->readFloat();
	m_focusDistance = stream->readFloat();
}

void ProjectiveCamera::serialize(Stream *stream, InstanceManager *manager) const {
	Sensor::serialize(stream, manager);
	stream->writeFloat(m_nearClip);
	stream->writeFloat(m_farClip);
	stream->writeFloat(m_focusDistance);
}

void ProjectiveCamera::setFocusDistance(Float focusDistance) {
	if (m_focusDistance != focusDistance) {
		m_focusDistance = focusDistance;
		m_properties.setFloat("focusDistance", focusDistance, false);
	}
}

void ProjectiveCamera::setNearClip(Float nearClip) {
	if (m_nearClip != nearClip) {
		m_nearClip = nearClip;
		m_properties.setFloat("nearClip", nearClip, false);
	}
}

void ProjectiveCamera::setFarClip(Float farClip) {
	if (m_farClip != farClip) {
		m_farClip = farClip;
		m_properties.setFloat("farClip", farClip, false);
	}
}

ProjectiveCamera::~ProjectiveCamera() {
}

void ProjectiveCamera::setWorldTransform(const Transform &trafo) {
	m_worldTransform = new AnimatedTransform(trafo);
	m_properties.setTransform("toWorld", trafo, false);
}

void ProjectiveCamera::setWorldTransform(AnimatedTransform *trafo) {
	m_worldTransform = trafo;
	m_properties.setAnimatedTransform("toWorld", trafo, false);
}

PerspectiveCamera::PerspectiveCamera(const Properties &props)
	: ProjectiveCamera(props), m_xfov(0.0f) {
	props.markQueried("fov");
	props.markQueried("fovAxis");
	props.markQueried("focalLength");

	if (m_properties.hasProperty("fov") && m_properties.hasProperty("focalLength"))
		Log(EError, "Please specify either a focal length ('focalLength') or a "
			"field of view ('fov')!");
}

PerspectiveCamera::PerspectiveCamera(Stream *stream, InstanceManager *manager)
	: ProjectiveCamera(stream, manager), m_xfov(0.0f) {
	setXFov(stream->readFloat());
}

PerspectiveCamera::~PerspectiveCamera() {
}

void PerspectiveCamera::configure() {
	ProjectiveCamera::configure();
	if (m_xfov != 0)
		return;

	if (m_properties.hasProperty("fov")) {
		Float fov = m_properties.getFloat("fov");

		std::string fovAxis =
			boost::to_lower_copy(m_properties.getString("fovAxis", "x"));

		if (fovAxis == "smaller")
			fovAxis = m_aspect > 1 ? "y" : "x";
		else if (fovAxis == "larger")
			fovAxis = m_aspect > 1 ? "x" : "y";

		if (fovAxis == "x")
			setXFov(fov);
		else if (fovAxis == "y")
			setYFov(fov);
		else if (fovAxis == "diagonal")
			setDiagonalFov(fov);
		else
			Log(EError, "The 'fovAxis' parameter must be set "
				"to one of 'smaller', 'larger', 'diagonal', 'x', or 'y'!");
	} else {
		std::string f = m_properties.getString("focalLength", "50mm");
		if (boost::ends_with(f, "mm"))
			f = f.substr(0, f.length()-2);

		char *end_ptr = NULL;
		Float value = (Float) strtod(f.c_str(), &end_ptr);
		if (*end_ptr != '\0')
			SLog(EError, "Could not parse the focal length (must be of the form "
				"<x>mm, where <x> is a positive integer)!");

		m_properties.removeProperty("focalLength");
		setDiagonalFov(2 * 180/M_PI* std::atan(std::sqrt((Float) (36*36+24*24)) / (2*value)));
	}
}

void PerspectiveCamera::serialize(Stream *stream, InstanceManager *manager) const {
	ProjectiveCamera::serialize(stream, manager);
	stream->writeFloat(m_xfov);
}

void PerspectiveCamera::setXFov(Float xfov) {
	if (xfov <= 0 || xfov >= 180)
		Log(EError, "The horizontal field of view must be "
			"in the interval (0, 180)!");
	if (xfov != m_xfov) {
		m_xfov = xfov;
		m_properties.setFloat("fov", xfov, false);
		m_properties.setString("fovAxis", "x", false);
	}
}

void PerspectiveCamera::setYFov(Float yfov) {
	setXFov(radToDeg(2*std::atan(
		std::tan(0.5f * degToRad(yfov)) * m_aspect)));
}

void PerspectiveCamera::setDiagonalFov(Float dfov) {
	Float diagonal = 2 * std::tan(0.5f * degToRad(dfov));
	Float width = diagonal / std::sqrt(1.0f + 1.0f / (m_aspect*m_aspect));
	setXFov(radToDeg(2*std::atan(width*0.5f)));
}


Float PerspectiveCamera::getYFov() const {
	return radToDeg(2*std::atan(
		std::tan(0.5f * degToRad(m_xfov)) / m_aspect));
}

Float PerspectiveCamera::getDiagonalFov() const {
	Float width = std::tan(0.5f * degToRad(m_xfov));
	Float diagonal = width * std::sqrt(1.0f + 1.0f / (m_aspect*m_aspect));
	return radToDeg(2*std::atan(diagonal));
}

MTS_IMPLEMENT_CLASS(PerspectiveCamera, true, ProjectiveCamera)
MTS_IMPLEMENT_CLASS(ProjectiveCamera, true, Sensor)
MTS_IMPLEMENT_CLASS(Sensor, true, AbstractEmitter)
MTS_NAMESPACE_END

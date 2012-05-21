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

#include <mitsuba/render/camera.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/core/plugin.h>

MTS_NAMESPACE_BEGIN

Camera::Camera(const Properties &props)
 : ConfigurableObject(props), m_properties(props) {
	m_cameraToWorld = props.getTransform("toWorld", Transform());
	m_shutterOpen = props.getFloat("shutterOpen", 0.0f);
	m_shutterClose = props.getFloat("shutterClose", 0.0f);
	if (m_shutterOpen > m_shutterClose)
		Log(EError, "Shutter opening time must be less than "
			"or equal to the shutter closing time!");
	m_worldToCamera = m_cameraToWorld.inverse();
	m_position = m_cameraToWorld(Point(0,0,0));
	m_shutterOpenTime = m_shutterClose - m_shutterOpen;
}

Camera::Camera(Stream *stream, InstanceManager *manager)
 : ConfigurableObject(stream, manager) {
	m_film = static_cast<Film *>(manager->getInstance(stream));
	m_sampler = static_cast<Sampler *>(manager->getInstance(stream));
	m_medium = static_cast<Medium *>(manager->getInstance(stream));
	m_worldToCamera = Transform(stream);
	m_cameraToWorld = Transform(stream);
	m_shutterOpen = stream->readFloat();
	m_shutterClose = stream->readFloat();
	m_position = m_cameraToWorld(Point(0,0,0));
	m_shutterOpenTime = m_shutterClose - m_shutterOpen;
}

Camera::~Camera() {
}

void Camera::setParent(ConfigurableObject *parent) {
	// Do not store the parent - this is useful in case
	// the camera subtree needs to be serialized by itself.
}

void Camera::configure() {
	if (m_film == NULL) {
		/* Instantiate an EXR film by default */
		m_film = static_cast<Film*> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Film), Properties("exrfilm")));
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
}

Point Camera::getPosition(const Point2 &sample) const {
	return m_position; // default impl.
}

void Camera::setFilm(Film *film) {
	m_film = film;
	configure();
}

void Camera::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
	manager->serialize(stream, m_film.get());
	manager->serialize(stream, m_sampler.get());
	manager->serialize(stream, m_medium.get());
	m_worldToCamera.serialize(stream);
	m_cameraToWorld.serialize(stream);
	stream->writeFloat(m_shutterOpen);
	stream->writeFloat(m_shutterClose);
}

void Camera::generateRayDifferential(const Point2 &sample,
	const Point2 &lensSample, Float timeSample, RayDifferential &ray) const {
	Ray tempRay;
	generateRay(sample, lensSample, timeSample, ray);
	Point2 temp = sample; temp.x += 1; 
	generateRay(temp, lensSample, timeSample, tempRay);
	ray.rxOrigin = tempRay.o; ray.rxDirection = tempRay.d;
	temp = sample; temp.y += 1;
	generateRay(temp, lensSample, timeSample, tempRay);
	ray.ryOrigin = tempRay.o; ray.ryDirection = tempRay.d;
	ray.hasDifferentials = true;
}

void Camera::addChild(const std::string &name, ConfigurableObject *child) {
	if (child->getClass()->derivesFrom(MTS_CLASS(Sampler))) {
		Assert(m_sampler == NULL);
		m_sampler = static_cast<Sampler *>(child);
	} else if (child->getClass()->derivesFrom(MTS_CLASS(Film))) {
		Assert(m_film == NULL);
		m_film = static_cast<Film *>(child);
	} else if (child->getClass()->derivesFrom(MTS_CLASS(Medium))) {
		Assert(m_medium == NULL);
		m_medium = static_cast<Medium *>(child);
	} else if (child->getClass()->derivesFrom(MTS_CLASS(Medium))) {
		Assert(m_medium == NULL);
		m_medium = static_cast<Medium *>(child);
	} else {
		Log(EError, "Camera: Invalid child node!");
	}
}

bool Camera::positionToSample(const Point &p, Point2 &sample) const {
	Log(EError, "%s::positionToSample(): not implemented!", 
			getClass()->getName().c_str());
	return false;
}

Float Camera::areaDensity(const Point2 &p) const {
	Log(EError, "%s::areaDensity(): not implemented!", 
			getClass()->getName().c_str());
	return 0.0f;
}

bool Camera::needsLensSample() const {
	return false;
}

ProjectiveCamera::ProjectiveCamera(Stream *stream, InstanceManager *manager)
 : Camera(stream, manager) {
	m_cameraToScreen = Transform(stream);
	m_cameraToScreenGL = Transform(stream);
	m_nearClip = stream->readFloat();
	m_farClip = stream->readFloat();
	m_aspect = stream->readFloat();
}
	
ProjectiveCamera::ProjectiveCamera(const Properties &props) : Camera(props) {
	/* Near clipping plane distance */
	m_nearClip = props.getFloat("nearClip", 1e-2f);
	/* Far clipping plane distance */
	m_farClip = props.getFloat("farClip", 1e4f);

	if (m_nearClip <= 0)
		Log(EError, "The 'nearClip' parameter must be greater than zero!");
	if (m_nearClip >= m_farClip)
		Log(EError, "The 'nearClip' parameter must be smaller than 'farClip'.");
}

void ProjectiveCamera::configure() {
	Camera::configure();
	m_aspect = (Float) m_film->getSize().x / (Float) m_film->getSize().y;
}

void ProjectiveCamera::serialize(Stream *stream, InstanceManager *manager) const {
	Camera::serialize(stream, manager);

	m_cameraToScreen.serialize(stream);
	m_cameraToScreenGL.serialize(stream);
	stream->writeFloat(m_nearClip);
	stream->writeFloat(m_farClip);
	stream->writeFloat(m_aspect);
}

PerspectiveCamera::PerspectiveCamera(const Properties &props)
 : ProjectiveCamera(props) {
	/* Field of view of the camera (in degrees) */
	m_fov = props.getFloat("fov", 45);
	/* Specifies which side of the image plane should cover the
	   field of view specified in the <tt>fov</tt> parameter
	*/
	m_mapSmallerSide = props.getBoolean("mapSmallerSide", true);
	/* Distance to the focal plane */
	m_focusDepth = props.getFloat("focusDepth", m_farClip);
	/* World-space aperture radius */
	m_apertureRadius = props.getFloat("apertureRadius", 0.0f);
}
	
PerspectiveCamera::PerspectiveCamera(Stream *stream, InstanceManager *manager)
 : ProjectiveCamera(stream, manager) {
	m_fov = stream->readFloat();
	m_apertureRadius = stream->readFloat();
	m_focusDepth = stream->readFloat();
	m_mapSmallerSide = stream->readBool();
}

void PerspectiveCamera::serialize(Stream *stream, InstanceManager *manager) const {
	ProjectiveCamera::serialize(stream, manager);

	stream->writeFloat(m_fov);
	stream->writeFloat(m_apertureRadius);
	stream->writeFloat(m_focusDepth);
	stream->writeBool(m_mapSmallerSide);
}

void PerspectiveCamera::configure() {
	ProjectiveCamera::configure();

	if ((m_aspect >= 1.0f) ^ !m_mapSmallerSide) {
		m_yfov = m_fov;
		m_xfov = radToDeg(2 * std::atan(std::tan(degToRad(m_fov)/2) * m_aspect));
	} else {
		m_xfov = m_fov;
		m_yfov = radToDeg(2 * std::atan(std::tan(degToRad(m_fov)/2) / m_aspect));
	}

	m_imagePlaneSize.x = 2 * std::tan(degToRad(m_xfov)/2);
	m_imagePlaneSize.y = 2 * std::tan(degToRad(m_yfov)/2);
	m_imagePlanePixelSize.x = m_imagePlaneSize.x / getFilm()->getSize().x;
	m_imagePlanePixelSize.y = m_imagePlaneSize.y / getFilm()->getSize().y;
	m_imagePlaneInvArea = 1 / (m_imagePlaneSize.x * m_imagePlaneSize.y);
}

Float PerspectiveCamera::importance(const Point2 &p) const {
	Float x = (p.x * m_imagePlanePixelSize.x) - .5f * m_imagePlaneSize.x;
	Float y = (p.y * m_imagePlanePixelSize.y) - .5f * m_imagePlaneSize.y;

	return std::pow(1 + x*x+y*y, (Float) (3.0/2.0)) * m_imagePlaneInvArea;
}

Float PerspectiveCamera::importance(const Vector &v) const {
	Vector localV;
	m_worldToCamera(v, localV);
	if (localV.z <= 0.0f) 
		return 0.0f;
	Float invZ = 1.0f / localV.z;
	localV.x *= invZ; localV.y *= invZ;
	if (2*std::abs(localV.x)>m_imagePlaneSize.x 
	 || 2*std::abs(localV.y)>m_imagePlaneSize.y) 
		return 0.0f;
	return std::pow(1 + localV.x*localV.x+localV.y*localV.y, 
		(Float) (3.0/2.0)) * m_imagePlaneInvArea;
}

bool PerspectiveCamera::needsLensSample() const {
	return m_apertureRadius > 0.0f;
}

Vector2 PerspectiveCamera::getImagePlaneSize(Float dist) const {
	return m_imagePlaneSize * dist;
}

MTS_IMPLEMENT_CLASS(Camera, true, ConfigurableObject)
MTS_IMPLEMENT_CLASS(ProjectiveCamera, true, Camera)
MTS_IMPLEMENT_CLASS(PerspectiveCamera, true, ProjectiveCamera)
MTS_NAMESPACE_END

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

#if !defined(__CAMERA_H)
#define __CAMERA_H

#include <mitsuba/render/film.h>

MTS_NAMESPACE_BEGIN

/** \brief Abstract camera base class. A camera turns a sample on
 * the image plane into a 3D ray. It uses two supporting child 
 * objects: a \re Sampler and a \ref Film instance.
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER Camera : public ConfigurableObject {
public:
	// =============================================================
	//! @{ \name Ray generation
	// =============================================================

	/// Create a ray from the given sample
	virtual void generateRay(const Point2 &sample, 
		const Point2 &lensSample,
		Float timeSample, Ray &ray) const = 0;

	/// Create ray differentials from the given sample
	void generateRayDifferential(const Point2 &sample, 
		const Point2 &lensSample, Float timeSample, 
		RayDifferential &ray) const;
	
	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Other query functions
	// =============================================================

	/**
	 * Turn a world-space position into fractional pixel coordinates.
	 * Returns false if the computed position is not visible through 
	 * the film's crop window
	 */
	virtual bool positionToSample(const Point &p, Point2 &sample) const;

	/// Does generateRay() expect a proper lens sample?
	virtual bool needsLensSample() const;
	
	/// Does generateRay() expect a proper time sample?
	inline bool needsTimeSample() const { return m_shutterOpenTime > 0; }

	/// Return the time value of the shutter opening event
	inline Float getShutterOpen() const { return m_shutterOpen; }

	/// Return the length, for which the shutter remains open
	inline Float getShutterOpenTime() const { return m_shutterOpenTime; }

	/// Return the time value of the shutter closing event
	inline Float getShutterClose() const { return m_shutterClose; }

	//! @}
	// =============================================================
	
	// =============================================================
	//! @{ \name Other query functions
	// =============================================================

	/// Return the camera position (approximate in the case of finite sensor area)
	inline const Point &getPosition() const { return m_position; }

	/// Return the camera position at the provided screen sample
	virtual Point getPosition(const Point2 &sample) const;

	/// Return the view transformation
	inline const Transform &getViewTransform() const { return m_worldToCamera; }

	/// Return the inverse view transformation
	inline const Transform &getInverseViewTransform() const { return m_cameraToWorld; }

	/// Set the view transformation
	inline void setViewTransform(const Transform &pTransform) {
		m_worldToCamera = pTransform;
		m_cameraToWorld = m_worldToCamera.inverse();
		m_position = m_cameraToWorld(Point(0,0,0));
	}
	
	/// Set the inverse view transformation
	inline void setInverseViewTransform(const Transform &pTransform) {
		m_cameraToWorld = pTransform;
		m_worldToCamera = m_cameraToWorld.inverse();
		m_position = m_cameraToWorld(Point(0,0,0));
	}

	/// Return the image plane normal
	inline Normal getImagePlaneNormal() const {
		return Normal(normalize(m_cameraToWorld(Vector(0, 0, 1))));
	}

	/**
	 * Calculate the pixel area density at a position on the image plane.
	 * Returns zero for cameras with an infinitesimal sensor (e.g. pinhole cameras).
	 */
	virtual Float areaDensity(const Point2 &p) const;

	//! @}
	// =============================================================
	
	// =============================================================
	//! @{ \name Miscellaneous
	// =============================================================

	/**
	 * \brief Return the camera's sampler.
	 *
	 * This is the 'root' sampler, which will later be cloned a 
	 * number of times to provide each participating worker thread 
	 * with its own instance (see \ref Scene::getSampler()). 
	 * Therefore, this sampler should never be used for anything 
	 * except creating clones.
	 */
	inline Sampler *getSampler() { return m_sampler; }

	/**
	 * \brief Return the camera's sampler (const version).
	 *
	 * This is the 'root' sampler, which will later be cloned a 
	 * number of times to provide each participating worker thread 
	 * with its own instance (see \ref Scene::getSampler()). 
	 * Therefore, this sampler should never be used for anything 
	 * except creating clones.
	 */
	inline const Sampler *getSampler() const { return m_sampler.get(); }

	/// Set the camera's film
	void setFilm(Film *film);

	/// Return the camera's film
	inline Film *getFilm() { return m_film; }

	/// Return the camera's film
	inline const Film *getFilm() const { return m_film.get(); }

	/// Return a pointer to the medium that surrounds the camera
	inline Medium *getMedium() { return m_medium; }
	
	/// Return a pointer to the medium that surrounds the camera (const version)
	inline const Medium *getMedium() const { return m_medium.get(); }

	/// Set the medium that surrounds the camera
	inline void setMedium(Medium *medium) { m_medium = medium; }

	/// Add a child ConfigurableObject
	virtual void addChild(const std::string &name, ConfigurableObject *child);
	/// Add an unnamed child
	inline void addChild(ConfigurableObject *child) { addChild("", child); }

	/// Serialize this camera to a binary data stream	
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	// Set the parent node
	void setParent(ConfigurableObject *parent);

	/// Return the properties of this camera 
	inline const Properties &getProperties() const { return m_properties; }

	/** \brief Configure the object (called _once_ after construction
	   and addition of all child ConfigurableObjects. */
	virtual void configure();

	//! @}
	// =============================================================

	MTS_DECLARE_CLASS()
protected:
	/// Create a new camera
	Camera(const Properties &props);

	/// Unserialize a new camera
	Camera(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~Camera();
protected:
	ref<Film> m_film;
	ref<Sampler> m_sampler;
	Transform m_worldToCamera, m_cameraToWorld;
	Point m_position;
	Properties m_properties;
	Float m_shutterOpen, m_shutterClose, m_shutterOpenTime;
	ref<Medium> m_medium;
};
 

/**
 * Projective camera base class
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER ProjectiveCamera : public Camera {
public:
	/// Return the projection transformation
	inline const Transform &getProjectionTransform() const { return m_cameraToScreen; }

	/// Return the projection transformation (using GL clip space coordinates)
	inline const Transform &getGLProjectionTransform() const { return m_cameraToScreenGL; }

	/// Return a projection transformation that includes a pixel offset (using GL clip space coordinates)
	virtual Transform getGLProjectionTransform(const Point2 &jitter) const = 0;

	/// Serialize this camera to a binary data stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;
	
	/** \brief Configure the object (called _once_ after construction
	   and addition of all child ConfigurableObjects. */
	virtual void configure();

	/// Return the near clip plane distance
	inline Float getNearClip() const { return m_nearClip; }

	/// Return the far clip plane distance
	inline Float getFarClip() const { return m_farClip; }

	MTS_DECLARE_CLASS()
protected:
	ProjectiveCamera(const Properties &props);
	ProjectiveCamera(Stream *stream, InstanceManager *manager);
	virtual ~ProjectiveCamera() { }
protected:
	Transform m_cameraToScreen, m_cameraToScreenGL;
	Float m_nearClip, m_farClip, m_aspect;
};

/**
 * \brief Base class of all perspective cameras 
 *
 * Provides solid angle computation routines useful 
 * for importance-based integrators.
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER PerspectiveCamera : public ProjectiveCamera {
public:
	/** \brief Configure the object (called _once_ after construction
	   and addition of all child ConfigurableObjects. */
	virtual void configure();

	/// Serialize this camera to a binary data stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Calculate the image plane size for a plane of the given distance
	Vector2 getImagePlaneSize(Float dist) const;

	/**
	 * \brief Calculate the importance of the given image-plane point
	 * expressed in fractional pixel coordinates. 
	 *
	 * Note that this value varies even within the individual pixels 
	 * due to the non-uniformity of rays generated by a strategy that
	 * uniformly samples points on the image plane.
	 */
	Float importance(const Point2 &p) const;

	/// Similar to importanceCamera(), but instead takes a world-space direction
	Float importance(const Vector &v) const;

	/// Return 

	/// Return the field of view along the X axis
	inline Float getXFov() const { return m_xfov; }

	/// Return the field of view along the Y axis
	inline Float getYFov() const { return m_yfov; }

	/// Get the field of view as specified in the scene file
	inline Float getFov() const { return m_fov; }

	/// Set the field of view as specified in the scene file
	void setFov(Float fov) { m_fov = fov; }

	/// Set the aperture radius
	void setApertureRadius(Float apertureRadius) { m_apertureRadius = apertureRadius; }

	/// Return the aperture radius
	inline Float getApertureRadius() const { return m_apertureRadius; }

	/// Set the focal distance
	void setFocusDepth(Float focusDepth) { m_focusDepth = focusDepth; }

	/// Return the focal distance
	inline Float getFocusDepth() const { return m_focusDepth; }

	/// Does generateRay() expect a proper lens sample?
	bool needsLensSample() const;

	MTS_DECLARE_CLASS()
protected:
	PerspectiveCamera(const Properties &props);
	PerspectiveCamera(Stream *stream, InstanceManager *manager);
	virtual ~PerspectiveCamera() { }
protected:
	Float m_fov, m_xfov, m_yfov;
	Vector2 m_imagePlaneSize;
	Vector2 m_imagePlanePixelSize;
	Float m_imagePlaneInvArea;
	Float m_apertureRadius;
	Float m_focusDepth;
	bool m_mapSmallerSide;
};

MTS_NAMESPACE_END

#endif /* __CAMERA_H */

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

#if !defined(__CAMERA_H)
#define __CAMERA_H

#include <mitsuba/render/film.h>

MTS_NAMESPACE_BEGIN

/** \brief Abstract camera base class. A camera turns a sample on
 * the image plane into a 3D ray. For this, it requires two 
 * supporting objects: a <tt>Sampler</tt> and a <tt>Film</tt> instance.
 */
class MTS_EXPORT_RENDER Camera : public ConfigurableObject {
public:
	/// Create a ray from the given sample
	virtual void generateRay(const Point2 &sample, const Point2 &lensSample,
		Ray &ray) const = 0;

	/// Create ray differentials from the given sample
	void generateRayDifferential(const Point2 &sample, 
		const Point2 &lensSample, RayDifferential &ray) const;

	/**
	 * Turn a world-space position into fractional pixel coordinates.
	 * Returns false if the computed position is not visible through 
	 * the film's crop window
	 */
	virtual bool positionToSample(const Point &p, Point2 &sample) const = 0;

	/// Does generateRay() expect a proper lens sample?
	virtual bool needsLensSample() const = 0;

	/// Return the camera position (approximate in the case of finite sensor area)
	inline const Point &getPosition() const { return m_position; }

	/// Return the camera position at the provided screen sample
	virtual Point getPosition(const Point2 &sample) const;

	/**
	 * Calculate the pixel area density at a position on the image plane.
	 * Returns zero for cameras with an infinitesimal sensor (e.g. pinhole cameras).
	 */
	virtual Float areaDensity(const Point2 &p) const = 0;

	/// Return the camera's sampler
	inline Sampler *getSamplerX() { return m_sampler; }

	/**
	 * Return the camera's sampler. This is the 'root' sampler,
	 * which will later be replicated for submission to all
	 * participating workers.
	 */
	inline const Sampler *getSamplerX() const { return m_sampler.get(); }

	/// Return the image plane normal
	inline Normal getImagePlaneNormal() const {
		return Normal(normalize(m_cameraToWorld(Vector(0, 0, 1))));
	}
	
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

	/// Set the camera's film
	void setFilm(Film *film);

	/// Return the camera's film
	inline Film *getFilm() { return m_film; }

	/// Return the camera's film
	inline const Film *getFilm() const { return m_film.get(); }

	/// Add a child ConfigurableObject
	virtual void addChild(const std::string &name, ConfigurableObject *child);

	/// Serialize this camera to disk	
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	// Set the parent node
	void setParent(ConfigurableObject *parent);

	/// Return the properties of this camera 
	inline const Properties &getProperties() const { return m_properties; }

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
};
	
class MTS_EXPORT_RENDER ProjectiveCamera : public Camera {
public:
	/// Return the projection transformation
	inline const Transform &getProjectionTransform() const { return m_cameraToScreen; }

	/// Return the projection transformation (using GL clip space coordinates)
	inline const Transform &getGLProjectionTransform() const { return m_cameraToScreenGL; }

	/// Return a projection transformation that includes a pixel offset (using GL clip space coordinates)
	virtual Transform getGLProjectionTransform(const Point2 &jitter) const = 0;

	/// Serialize this camera to disk
	virtual void serialize(Stream *stream, InstanceManager *manager) const;
	
	/** \brief Configure the object (called _once_ after construction
	   and addition of all child ConfigurableObjects. */
	virtual void configure();

	MTS_DECLARE_CLASS()
protected:
	inline ProjectiveCamera(const Properties &props) : Camera(props) { }
	ProjectiveCamera(Stream *stream, InstanceManager *manager);
	virtual ~ProjectiveCamera() { }
protected:
	Transform m_cameraToScreen, m_cameraToScreenGL;
	Float m_aspect;
};

/**
 * Base class of all pinhole cameras. Provides solid angle computation
 * routines useful for importance-based integrators.
 */
class MTS_EXPORT_RENDER PinholeCamera : public ProjectiveCamera {
public:
	/** \brief Configure the object (called _once_ after construction
	   and addition of all child ConfigurableObjects. */
	virtual void configure();

	/// Serialize this camera to disk
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Calculate the image plane size for a plane of the given distance
	Vector2 getImagePlaneSize(Float dist) const;

	/**
	 * Calculate the solid angle subtended by the area [xs,xe]x[ys,ye] 
	 * on the image plane (in pixel coordinates starting at 0). 
	 */
	Float solidAngle(Float xs, Float xe, Float ys, Float ye) const;

	/**
	 * Calculate the importance of the given image-plane point
	 * expressed in fractional pixel coordinates. Note that this value
	 * varies even within the individual pixels due to the non-uniformity
	 * of rays generated by a strategy that uniformly samples points on 
	 * the image plane.
	 */
	Float importance(const Point2 &p) const;

	/// Similar to importanceCamera(), but instead takes a world-space direction
	Float importance(const Vector &v) const;

	/// Return the field of view along the X axis
	inline Float getXFov() const { return m_xfov; }

	/// Return the field of view along the Y axis
	inline Float getYFov() const { return m_yfov; }

	/// Get the field of view as specified in the scene file
	inline Float getFov() const { return m_fov; }

	/// Set the field of view as specified in the scene file
	void setFov(float fov) { m_fov = fov; }

	MTS_DECLARE_CLASS()
protected:
	PinholeCamera(const Properties &props);
	PinholeCamera(Stream *stream, InstanceManager *manager);
	virtual ~PinholeCamera() { }
protected:
	Float m_fov, m_xfov, m_yfov;
	Vector2 m_imagePlaneSize;
	Vector2 m_imagePlanePixelSize;
	Float m_imagePlaneInvArea;
	bool m_mapSmallerSide;
};

MTS_NAMESPACE_END

#endif /* __CAMERA_H */

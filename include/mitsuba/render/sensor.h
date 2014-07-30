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

#pragma once
#if !defined(__MITSUBA_RENDER_SENSOR_H_)
#define __MITSUBA_RENDER_SENSOR_H_

#include <mitsuba/render/common.h>
#include <mitsuba/render/film.h>
#include <mitsuba/render/emitter.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Abstract sensor interface
 *
 * This class provides an abstract interface to all sensor plugins in Mitsuba.
 * It exposes functions for evaluating and sampling the response function of the
 * sensor, and it allows querying the probability density of the sampling method.
 *
 * Somewhat curiously, the \ref Sensor class derives from \ref AbstractEmitter.
 * The reason for this is that much like radiance, the spectral response of a
 * sensor can be interpreted as emitted quantity named \a importance. The
 * \ref Sensor interface thus inherits almost all of the emitter API and only
 * needs to add a few camera-specific methods on top.
 *
 * The concept of interpreting sensor response as an emitted quantity and
 * the resulting flexibility of being able to dynamically transition between
 * emitter and receiver interpretations of luminaires and sensors is a key
 * insight that enables the construction of powerful bidirectional rendering
 * techniques It is the reason why the API to these components may seem
 * somewhat unorthodox.
 *
 * In Mitsuba, a sensor can be as simple as an irradiance meter that performs a
 * single measurement along a specified ray, but it can also represent sensors
 * that are more commonly used in computer graphics, such as a perspective camera
 * based on the thin lens equation.
 *
 * An important difference between a luminaire and a sensor is that the sensor
 * records spectral measurements to a film, and for that reason it needs a
 * mapping between rays and film pixel coordinates. Apart from that, the
 * interfaces are almost identical.
 *
 * Mitsuba assumes that a sensor always has a form of "shutter", which opens
 * for a certain time, during which the exposure takes place. The sensor
 * itself may also undergo motion while the shutter is open, but a more
 * complicated dependence on time is not allowed.
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER Sensor : public AbstractEmitter {
public:
	/**
	 * \brief This list of flags is used to additionally characterize
	 * and classify the response functions of different types of sensors
	 *
	 * \sa AbstractEmitter::EEmitterType
	 */
	enum ESensorFlags {
		/// Sensor response contains a Dirac delta term with respect to time
		EDeltaTime             = 0x010,

		/// Does the \ref sampleRay() function need an aperture sample?
		ENeedsApertureSample   = 0x020,

		/// Is the sensor a projective camera?
		EProjectiveCamera      = 0x100,

		/// Is the sensor a perspective camera?
		EPerspectiveCamera     = 0x200,

		/// Is the sensor an orthographic camera?
		EOrthographicCamera    = 0x400,

		/// Does the sample given to \ref samplePosition() determine the pixel coordinates
		EPositionSampleMapsToPixels  = 0x1000,

		/// Does the sample given to \ref sampleDirection() determine the pixel coordinates
		EDirectionSampleMapsToPixels = 0x2000
	};

	// =============================================================
	//! @{ \name Additional sensor-related sampling functions
	// =============================================================

 	/**
	 * \brief Importance sample a ray according to the sensor response
	 *
	 * This function combines all three of the steps of sampling a time,
	 * ray position, and direction value. It does not return any auxiliary
	 * sampling information and is mainly meant to be used by unidirectional
	 * rendering techniques.
	 *
	 * Note that this function potentially uses a different sampling
	 * strategy compared to the sequence of running \ref sampleArea()
	 * and \ref sampleDirection(). The reason for this is that it may
	 * be possible to switch to a better technique when sampling both
	 * position and direction at the same time.
	 *
	 * \param ray
	 *    A ray data structure to be populated with a position
	 *    and direction value
	 *
	 * \param samplePosition
	 *    Denotes the desired sample position in fractional pixel
	 *    coordinates relative to the crop window of the underlying
	 *    film.
	 *
	 * \param apertureSample
	 *    A uniformly distributed 2D vector that is used to sample
	 *    a position on the aperture of the sensor if necessary.
	 *    (Any value is valid when \ref needsApertureSample() == \c false)
	 *
	 * \param timeSample
	 *    A uniformly distributed 1D vector that is used to sample
	 *    the temporal component of the emission profile.
	 *    (Or any value when \ref needsTimeSample() == \c false)
	 *
	 * \return
	 *    An importance weight associated with the sampled ray.
	 *    This accounts for the difference between the sensor response
	 *    and the sampling density function.
	 *
	 * \remark
	 *    In the Python API, the signature of this function is
	 *    <tt>spectrum, ray = sensor.sampleRay(samplePosition, apertureSample)</tt>
	 */
	virtual Spectrum sampleRay(Ray &ray,
		const Point2 &samplePosition,
		const Point2 &apertureSample,
		Float timeSample) const = 0;

	/**
	 * \brief Importance sample a ray differential according to the
	 * sensor response
	 *
	 * This function combines all three of the steps of sampling a time,
	 * ray position, and direction value. It does not return any auxiliary
	 * sampling information and is mainly meant to be used by unidirectional
	 * rendering techniques.
	 *
	 * Note that this function potentially uses a different sampling
	 * strategy compared to the sequence of running \ref sampleArea()
	 * and \ref sampleDirection(). The reason for this is that it may
	 * be possible to switch to a better technique when sampling both
	 * position and direction at the same time.
	 *
	 * The default implementation computes differentials using several
	 * internal calls to \ref sampleRay(). Subclasses of the \ref Sensor
	 * interface may optionally provide a more efficient approach.
	 *
	 * \param ray
	 *    A ray data structure to be populated with a position
	 *    and direction value
 	 *
	 * \param samplePosition
	 *    Denotes the desired sample position in fractional pixel
	 *    coordinates relative to the crop window of the underlying
	 *    film.
	 *
	 * \param apertureSample
	 *    A uniformly distributed 2D vector that is used to sample
	 *    a position on the aperture of the sensor if necessary.
	 *    (Any value is valid when \ref needsApertureSample() == \c false)

	 * \param timeSample
	 *    A uniformly distributed 1D vector that is used to sample
	 *    the temporal component of the emission profile.
	 *    (Or any value when \ref needsTimeSample() == \c false)
	 *
	 * \return
	 *    An importance weight associated with the sampled ray.
	 *    This accounts for the difference between the sensor response
	 *    and the sampling density function.
	 *
	 * \remark
	 *    In the Python API, the signature of this function is
	 *    <tt>spectrum, ray = sensor.sampleRayDifferential(samplePosition, apertureSample)</tt>
	 */
	virtual Spectrum sampleRayDifferential(RayDifferential &ray,
		const Point2 &samplePosition,
		const Point2 &apertureSample,
		Float timeSample) const;

	/// Importance sample the temporal part of the sensor response function
	inline Float sampleTime(Float sample) const {
		return m_shutterOpen + m_shutterOpenTime * sample;
	}

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Additional query functions
	// =============================================================

	/**
	 * \brief Return the emitted importance for the given surface intersection
	 *
	 * This is function is used when a sensor has been hit by a
	 * ray in a particle tracing-style integrator, and it subsequently needs to
	 * be queried for the emitted importance along the negative ray direction.
	 *
	 * It efficiently computes the product of \ref evalPosition()
	 * and \ref evalDirection(), though note that it does not include the
	 * cosine foreshortening factor of the latter method.
	 *
	 * This function is provided here as a fast convenience function for
	 * unidirectional rendering techniques that support intersecting the
	 * sensor. The default implementation throws an exception, which
	 * states that the method is not implemented.
	 *
	 * \param its
	 *    An intersect record that specfies the query position
	 *
	 * \param d
	 *    A unit vector, which specifies the query direction
	 *
	 * \param result
	 *    This argument is used to return the 2D sample position
	 *    (i.e. the fractional pixel coordinates) associated
	 *    with the intersection.
	 *
	 * \return
	 *    The emitted importance
	 *
	 * \remark
	 *    In the Python API, the signature of this function is
	 *    <tt>spectrum, samplePos = sensor.eval(its, d)</tt>
	 */
	virtual Spectrum eval(const Intersection &its, const Vector &d,
			Point2 &samplePos) const;

	/**
	 * \brief Return the sample position associated with a given
	 * position and direction sampling record
	 *
	 * \param dRec
	 *    A direction sampling record, which specifies the query direction
	 *
	 * \param pRec
	 *    A position sampling record, which specifies the query position
	 *
	 * \return \c true if the specified ray is visible by the camera
	 *
	 * \remark
	 *    In the Python API, the signature of this function is
	 *    <tt>visible, position = sensor.getSamplePosition(pRec, dRec)</tt>
	 */
	virtual bool getSamplePosition(const PositionSamplingRecord &pRec,
			const DirectionSamplingRecord &dRec, Point2 &position) const;

	/**
	 * \brief Evaluate the temporal component of the sampling density
	 * implemented by the \ref sampleRay() method.
	 */
	Float pdfTime(const Ray &ray, EMeasure measure) const;

	/// Return the time value of the shutter opening event
	inline Float getShutterOpen() const { return m_shutterOpen; }

	/// Set the time value of the shutter opening event
	void setShutterOpen(Float time) { m_shutterOpen = time; }

	/// Return the length, for which the shutter remains open
	inline Float getShutterOpenTime() const { return m_shutterOpenTime; }

	/// Set the length, for which the shutter remains open
	void setShutterOpenTime(Float time);

	/**
	 * \brief Does the method \ref sampleRay() require a uniformly distributed
	 * sample for the time-dependent component?
	 */
	inline bool needsTimeSample() const { return !(m_type & EDeltaTime); }

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Miscellaneous
	// =============================================================

	/**
	 * \brief Does the method \ref sampleRay() require a uniformly
	 * distributed sample for the aperture component?
	 */
	inline bool needsApertureSample() const { return m_type & ENeedsApertureSample; }

	/// Return the \ref Film instance associated with this sensor
	inline Film *getFilm() { return m_film; }

	/// Return the \ref Film instance associated with this sensor (const)
	inline const Film *getFilm() const { return m_film.get(); }

	/// Return the aspect ratio of the sensor and its underlying film
	inline Float getAspect() const { return m_aspect; }

	/**
	 * \brief Return the sensor's sample generator
	 *
	 * This is the \a root sampler, which will later be cloned a
	 * number of times to provide each participating worker thread
	 * with its own instance (see \ref Scene::getSampler()).
	 * Therefore, this sampler should never be used for anything
	 * except creating clones.
	 */
	inline Sampler *getSampler() { return m_sampler; }

	/**
	 * \brief Return the sensor's sampler (const version).
	 *
	 * This is the \a root sampler, which will later be cloned a
	 * number of times to provide each participating worker thread
	 * with its own instance (see \ref Scene::getSampler()).
	 * Therefore, this sampler should never be used for anything
	 * except creating clones.
	 */
	inline const Sampler *getSampler() const { return m_sampler.get(); }

	/// Serialize this sensor to a binary data stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name ConfigurableObject interface
	// =============================================================
	/// Add a child ConfigurableObject
	virtual void addChild(const std::string &name, ConfigurableObject *child);
	/// Add an unnamed child
	inline void addChild(ConfigurableObject *child) { addChild("", child); }

	/** \brief Configure the object (called \a once after construction
	   and addition of all child \ref ConfigurableObject instances). */
	virtual void configure();

	//! @}
	// =============================================================

	MTS_DECLARE_CLASS()
protected:
	/// Construct a new sensor instance
	Sensor(const Properties &props);

	/// Unserialize a sensor instance from a binary data stream
	Sensor(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~Sensor();
protected:
	ref<Film> m_film;
	ref<Sampler> m_sampler;
	Vector2 m_resolution;
	Vector2 m_invResolution;
	Float m_shutterOpen;
	Float m_shutterOpenTime;
	Float m_aspect;
};

/**
 * \brief Projective camera interface
 *
 * This class provides an abstract interface to several types of sensors that
 * are commonly used in computer graphics, such as perspective and orthographic
 * camera models.
 *
 * The interface is meant to be implemented by any kind of sensor, whose
 * world to clip space transformation can be explained using only linear
 * operations on homogeneous coordinates.
 *
 * A useful feature of \ref ProjectiveCamera sensors is that their view can be
 * rendered using the traditional OpenGL pipeline.
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER ProjectiveCamera : public Sensor {
public:
	using Sensor::getWorldTransform;

	/// Return the world-to-view (aka "view") transformation at time \c t
	inline const Transform getViewTransform(Float t) const {
		return getWorldTransform()->eval(t).inverse();
	}

	/// Return the view-to-world transformation at time \c t
	inline const Transform getWorldTransform(Float t) const {
		return getWorldTransform()->eval(t);
	}

	/**
	 * \brief Overwrite the view-to-world transformation
	 * with a static (i.e. non-animated) transformation.
	 */
	void setWorldTransform(const Transform &trafo);

	/**
	 * \brief Overwrite the view-to-world transformation
	 * with an animated transformation
	 */
	void setWorldTransform(AnimatedTransform *trafo);

	/**
	 * \brief Return a projection matrix suitable for rendering the
	 * scene using OpenGL
	 *
	 * For scenes involving a narrow depth of field and antialiasing,
	 * it is necessary to average many separately rendered images using
	 * different pixel offsets and aperture positions.
	 *
	 * \param apertureSample
	 *     Sample for rendering with defocus blur. This should be a
	 *     uniformly distributed random point in [0,1]^2 (or any value
	 *     when \ref needsApertureSample() == \c false)
	 *
	 * \param aaSample
	 *     Sample for antialiasing. This should be a uniformly
	 *     distributed random point in [0,1]^2.
	 */
	virtual Transform getProjectionTransform(const Point2 &apertureSample,
			const Point2 &aaSample) const = 0;

	/// Serialize this camera to a binary data stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Return the near clip plane distance
	inline Float getNearClip() const { return m_nearClip; }

	/// Set the near clip plane distance
	void setNearClip(Float nearClip);

	/// Return the far clip plane distance
	inline Float getFarClip() const { return m_farClip; }

	/// Set the far clip plane distance
	void setFarClip(Float farClip);

	/// Return the distance to the focal plane
	inline Float getFocusDistance() const { return m_focusDistance; }

	/// Set the distance to the focal plane
	void setFocusDistance(Float focusDistance);

	MTS_DECLARE_CLASS()
protected:
	/// Construct a new camera instance
	ProjectiveCamera(const Properties &props);

	/// Unserialize a camera instance from a binary data stream
	ProjectiveCamera(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~ProjectiveCamera();
protected:
	Float m_nearClip;
	Float m_farClip;
	Float m_focusDistance;
};

/**
 * \brief Perspective camera interface
 *
 * This class provides an abstract interface to several types of sensors that
 * are commonly used in computer graphics, such as perspective and orthographic
 * camera models.
 *
 * The interface is meant to be implemented by any kind of sensor, whose
 * world to clip space transformation can be explained using only linear
 * operations on homogeneous coordinates.
 *
 * A useful feature of \ref ProjectiveCamera sensors is that their view can be
 * rendered using the traditional OpenGL pipeline.
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER PerspectiveCamera : public ProjectiveCamera {
public:
	// =============================================================
	//! @{ \name Field of view-related
	// =============================================================

	/// Return the horizontal field of view in degrees
	inline Float getXFov() const { return m_xfov; }

	/// Set the horizontal field of view in degrees
	void setXFov(Float xfov);

	/// Return the vertical field of view in degrees
	Float getYFov() const;

	/// Set the vertical field of view in degrees
	void setYFov(Float yfov);

	/// Return the diagonal field of view in degrees
	Float getDiagonalFov() const;

	/// Set the diagonal field of view in degrees
	void setDiagonalFov(Float dfov);

	//! @}
	// =============================================================

	/** \brief Configure the object (called \a once after construction
	   and addition of all child \ref ConfigurableObject instances). */
	virtual void configure();

	/// Serialize this camera to a binary data stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	MTS_DECLARE_CLASS()
protected:
	/// Construct a new perspective camera instance
	PerspectiveCamera(const Properties &props);

	/// Unserialize a perspective camera instance from a binary data stream
	PerspectiveCamera(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~PerspectiveCamera();
protected:
	Float m_xfov;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_SENSOR_H_ */

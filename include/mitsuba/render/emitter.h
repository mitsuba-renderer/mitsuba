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
#if !defined(__MITSUBA_RENDER_EMITTER_H_)
#define __MITSUBA_RENDER_EMITTER_H_

#include <mitsuba/render/common.h>
#include <mitsuba/core/track.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/cobject.h>
#include <mitsuba/render/shader.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Abstract radiance/importance emitter interface
 *
 * This class implements an abstract interface to all emitters in Mitsuba. The
 * term emitter is interpreted in a loose bidirectional sense here, where
 * luminaires and sensors are both considered to be emitters of radiance and
 * importance, respectively.
 *
 * Subclasses must implement functions for evaluating and sampling the
 * emission profile and furthermore support querying the probability density
 * of the provided sampling technique.
 *
 * Subclasses must also provide a specialized \a direct sampling method
 * (a generalization of direct illumination sampling to both emitters \a and
 * sensors). A direct sampling is given an arbitrary input position in the
 * scene and in turn returns a sampled emitter position and direction, which
 * has a nonzero contribution towards the provided position. The main idea is
 * that direct sampling reduces the underlying space from 4D to 2D, hence it
 * is often possible to use smarter sampling techniques than in the fully
 * general case.
 *
 * Since the emission profile is defined as function over both positions
 * and directions, there are functions to sample and query \a each of the
 * two components separately. Furthermore, there is a convenience function
 * to sample both at the same time, which is mainly used by unidirectional
 * rendering algorithms that do not need this level of flexibility.
 *
 * One underlying assumption of this interface is that position and
 * direction sampling will happen <em>in sequence</em>. This means that
 * the direction sampling step is allowed to statistically condition on
 * properties of the preceding position sampling step.
 *
 * When rendering scenes involving participating media, it is important
 * to know what medium surrounds the sensors and light sources. For
 * this reason, every emitter instance keeps a reference to a medium
 * (or \c NULL when it is surrounded by vacuum).
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER AbstractEmitter : public ConfigurableObject {
public:
    /**
     * \brief Flags used to classify the emission profile of
     * different types of emitters
     */
    enum EEmitterType {
        /// Emission profile contains a Dirac delta term with respect to direction
        EDeltaDirection = 0x01,

        /// Emission profile contains a Dirac delta term with respect to position
        EDeltaPosition  = 0x02,

        /// Is the emitter associated with a surface in the scene?
        EOnSurface      = 0x04
    };

    // =============================================================
    //! @{ \name Sampling interface
    // =============================================================

    /**
     * \brief Importance sample the spatial component of the
     * emission profile.
     *
     * This function takes an uniformly distributed 2D vector
     * and maps it to a position on the surface of the emitter.
     *
     * Some implementations may choose to implement extra functionality
     * based on the value of \c extra: for instance, Sensors
     * (which are a subclass of \ref AbstractEmitter) perform uniform
     * sampling over the entire image plane if <tt>extra == NULL</tt>,
     * but other values, they will restrict sampling to a pixel-sized
     * rectangle with that offset.
     *
     * The default implementation throws an exception.
     *
     * \param pRec
     *    A position record to be populated with the sampled
     *    position and related information
     *
     * \param sample
     *    A uniformly distributed 2D vector (or any value,
     *    when \ref needsPositionSample() == \c false)
     *
     * \param extra
     *    An additional 2D vector provided to the sampling
     *    routine -- its use is implementation-dependent.
     *
     * \return
     *    An importance weight associated with the sampled position.
     *    This accounts for the difference between the spatial part of the
     *    emission profile and the density function.
     */
    virtual Spectrum samplePosition(PositionSamplingRecord &pRec,
        const Point2 &sample, const Point2 *extra = NULL) const;

    /**
     * \brief Conditioned on the spatial component, importance
     * sample the directional part of the emission profile.
     *
     * Some implementations may choose to implement extra functionality
     * based on the value of \c extra: for instance, Sensors
     * (which are a subclass of \ref AbstractEmitter) perform uniform
     * sampling over the entire image plane if <tt>extra == NULL</tt>,
     * but other values, they will restrict sampling to a pixel-sized
     * rectangle with that offset.
     *
     * The default implementation throws an exception.
     *
     * \param dRec
     *    A direction record to be populated with the sampled
     *    direction and related information
     *
     * \param pRec
     *    A position record generated by a preceding call
     *    to \ref samplePosition()
     *
     * \param sample
     *    A uniformly distributed 2D vector (or any value
     *    when \ref needsDirectionSample() == \c false)
     *
     * \return
     *    An importance weight associated with the sampled direction.
     *    This accounts for the difference between the directional part of the
     *    emission profile and the density function.
     */
    virtual Spectrum sampleDirection(
        DirectionSamplingRecord &dRec,
        PositionSamplingRecord &pRec,
        const Point2 &sample,
        const Point2 *extra = NULL) const;

    /**
     * \brief \a Direct sampling: given a reference point in the
     * scene, sample an emitter position that contributes towards it.
     *
     * Given an arbitrary reference point in the scene, this method
     * samples a position on the emitter that has a nonzero contribution
     * towards that point.
     * This can be seen as a generalization of direct illumination sampling
     * so that it works on both luminaires and sensors.
     *
     * Ideally, the implementation should importance sample the product of
     * the emission profile and the geometry term between the reference point
     * and the position on the emitter.
     *
     * The default implementation throws an exception.
     *
     * \param dRec
     *    A direct sampling record that specifies the reference point and
     *    a time value. After the function terminates, it will be
     *    populated with the position sample and related information
     *
     * \param sample
     *    A uniformly distributed 2D vector (or any value
     *    when \ref needsDirectSample() == \c false)
     *
     * \return
     *    An importance weight associated with the sample. Includes
     *    any geometric terms between the emitter and the reference point.
     */
    virtual Spectrum sampleDirect(DirectSamplingRecord &dRec,
            const Point2 &sample) const;

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Query functions for the emission profile and
    //!          sampling density functions
    // =============================================================

    /**
     * \brief Evaluate the spatial component of the emission profile
     *
     * \param pRec
     *    A position sampling record, which specifies the query location
     *
     * \return The component of the emission profile that depends on
     * the position (i.e. emitted power per unit area for luminaires and
     * sensor response, or inverse power per unit area for sensors)
     */
    virtual Spectrum evalPosition(const PositionSamplingRecord &pRec) const;

    /**
     * \brief Evaluate the directional component of the emission profile
     *
     * When querying a smooth (i.e. non-degenerate) component, it already
     * multiplies the result by the cosine foreshortening factor with
     * respect to the outgoing direction.
     *
     * \param dRec
     *    A direction sampling record, which specifies the query direction
     *
     * \param pRec
     *    A position sampling record, which specifies the query position
     *
     * \return The component of the emission profile that depends on
     * the direction (having units of 1/steradian)
     */

    virtual Spectrum evalDirection(const DirectionSamplingRecord &dRec,
            const PositionSamplingRecord &pRec) const;

    /**
     * \brief Evaluate the spatial component of the sampling density
     * implemented by the \ref samplePosition() method
     *
     * \param pRec
     *    A position sampling record, which specifies the query location
     *
     * \return
     *    The area density at the supplied position
     */
    virtual Float pdfPosition(const PositionSamplingRecord &pRec) const;

    /**
     * \brief Evaluate the directional component of the sampling density
     * implemented by the \ref sampleDirection() method
     *
     * \param dRec
     *    A direction sampling record, which specifies the query direction
     *
     * \param pRec
     *    A position sampling record, which specifies the query position
     *
     * \return
     *    The directional density at the supplied position
     */

    virtual Float pdfDirection(const DirectionSamplingRecord &dRec,
        const PositionSamplingRecord &pRec) const;

    /**
     * \brief Evaluate the probability density of the \a direct sampling
     * method implemented by the \ref sampleDirect() method.
     *
     * \param dRec
     *    A direct sampling record, which specifies the query
     *    location. Note that this record need not be completely
     *    filled out. The important fields are \c p, \c n, \c ref,
     *    \c dist, \c d, \c measure, and \c uv.
     *
     * \return
     *    The density expressed with respect to the requested measure
     *    (usually \ref ESolidAngle)
     */
    virtual Float pdfDirect(const DirectSamplingRecord &dRec) const;

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Other query functions
    // =============================================================

    /**
     * \brief Return a listing of classification flags combined
     * using binary OR.
     *
     * \sa EEmitterType
     */
    inline uint32_t getType() const { return m_type; }

    /// Return the local space to world space transformation
    inline const AnimatedTransform *getWorldTransform() const
            { return m_worldTransform.get(); }

    /// Set the local space to world space transformation
    inline void setWorldTransform(AnimatedTransform *trafo)
            { m_worldTransform = trafo; }

    /**
     * \brief Does the method \ref samplePosition() require a uniformly
     * distributed sample for the spatial component?
     */
    inline bool needsPositionSample() const { return !(m_type & EDeltaPosition); }

    /**
     * \brief Does the method \ref sampleDirection() require a uniformly
     * distributed sample for the direction component?
     */
    inline bool needsDirectionSample() const { return !(m_type & EDeltaDirection); }

    /**
     * \brief Does the emitter lie on some kind of surface?
     */
    inline bool isOnSurface() const { return m_type & EOnSurface; }

    /**
     * \brief Does the sensor have a degenerate directional or spatial
     * distribution?
     */
    inline bool isDegenerate() const { return m_type & (EDeltaPosition | EDeltaDirection); }

    /**
     * \brief Does the method \ref sampleDirect() require a uniformly
     * distributed sample?
     *
     * Since sampleDirect() essentially causes a 2D reduction of the
     * sampling domain, this is the case exactly when the original
     * domain was four-dimensionsional.
     */
    inline bool needsDirectSample() const {
        return needsPositionSample() && needsDirectionSample();
    }

    /**
     * \brief Return the measure associated with the \ref sampleDirect()
     * operation
     */
    inline EMeasure getDirectMeasure() const {
        return needsDirectSample() ? ESolidAngle : EDiscrete;
    }

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Miscellaneous
    // =============================================================

    /// Return a pointer to the medium that surrounds the emitter
    inline Medium *getMedium() { return m_medium; }

    /// Return a pointer to the medium that surrounds the emitter (const version)
    inline const Medium *getMedium() const { return m_medium.get(); }

    /// Return the shape, to which the emitter is currently attached
    inline Shape *getShape() { return m_shape; }

    /// Return the shape, to which the emitter is currently attached (const version)
    inline const Shape *getShape() const { return m_shape; }

    /**
     * \brief Create a special shape that represents the emitter
     *
     * Some types of emitters are inherently associated with a surface, yet
     * this surface is not explicitly needed for many kinds of rendering
     * algorithms.
     *
     * An example would be an environment map, where the associated shape
     * is a sphere surrounding the scene. Another example would be a
     * perspective camera with depth of field, where the associated shape
     * is a disk representing the aperture (remember that this class
     * represents emitters in a generalized bidirectional sense, which
     * includes sensors).
     *
     * When this shape is in fact needed by the underlying rendering algorithm,
     * this function can be called to create it. The default implementation
     * simply returns \c NULL.
     *
     * \param scene
     *     A pointer to the associated scene (the created shape is
     *     allowed to depend on it)
     */
    virtual ref<Shape> createShape(const Scene *scene);

    /**
     * \brief Return an axis-aligned box bounding the spatial
     * extents of the emitter
     */
    virtual AABB getAABB() const = 0;

    /// Set the medium that surrounds the emitter
    inline void setMedium(Medium *medium) { m_medium = medium; }

    /// Serialize this emitter to a binary data stream
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

    //! @}
    // =============================================================

    MTS_DECLARE_CLASS()
protected:
    /// Construct a new emitter instance
    AbstractEmitter(const Properties &props);

    /// Unserialize a emitter instance from a binary data stream
    AbstractEmitter(Stream *stream, InstanceManager *manager);

    /// Virtual destructor
    virtual ~AbstractEmitter();
protected:
    ref<const AnimatedTransform> m_worldTransform;
    ref<Medium> m_medium;
    Shape *m_shape;
    uint32_t m_type;
};

/**
 * \brief Abstract radiance emitter interface
 *
 * This class provides an abstract interface to all emitter plugins in Mitsuba.
 * It exposes functions for evaluating and sampling the emission profile, and
 * it allows querying the probability density of the sampling method.
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER Emitter : public AbstractEmitter, public HWResource {
public:
    /**
     * \brief This list of flags is used to additionally characterize
     * and classify the response functions of different types of sensors
     *
     * \sa AbstractEmitter::EEmitterType
     */
    enum EEmitterFlags {
        /// Is this an environment emitter, such as a HDRI map?
        EEnvironmentEmitter = 0x010
    };

    // =================================================================
    //! @{ \name Additional emitter-related sampling and query functions
    // =================================================================

    /**
     * \brief Return the radiant emittance for the given surface intersection
     *
     * This is function is used when an area light source has been hit by a
     * ray in a path tracing-style integrator, and it subsequently needs to
     * be queried for the emitted radiance along the negative ray direction.
     *
     * It efficiently computes the product of \ref evalPosition()
     * and \ref evalDirection() \a divided by the absolute cosine of the
     * angle between \c d and \c its.shFrame.n.
     *
     * This function is provided here as a fast convenience function for
     * unidirectional rendering techniques. The default implementation
     * throws an exception, which states that the method is not implemented.
     *
     * \param its
     *    An intersect record that specfies the query position
     * \param d
     *    A unit vector, which specifies the query direction
     * \return
     *    The radiant emittance
     */
    virtual Spectrum eval(const Intersection &its, const Vector &d) const;

    /**
     * \brief Importance sample a ray according to the emission profile
     *
     * This function combines both steps of choosing a ray origin and
     * direction value. It does not return any auxiliary sampling
     * information and is mainly meant to be used by unidirectional
     * rendering techniques.
     *
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
     * \param spatialSample
     *    Denotes the sample that is used to choose the spatial component
     *
     * \param directionalSample
     *    Denotes the sample that is used to choose the directional component
     *
     * \param time
     *    Scene time value to be associated with the sample
     *
     * \return
     *    An importance weight associated with the sampled ray.
     *    This accounts for the difference between the emission profile
     *    and the sampling density function.
     */
    virtual Spectrum sampleRay(Ray &ray,
        const Point2 &spatialSample,
        const Point2 &directionalSample,
        Float time) const;

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Miscellaneous
    // =============================================================

    /**
     * \brief Return the luminaire's sampling weight
     *
     * This is used by the luminaire importance sampling
     * routines in \ref Scene.
     */
    inline Float getSamplingWeight() const { return m_samplingWeight; }

    /**
     * \brief Return a bitmap representation of the emitter
     *
     * Some types of light sources (projection lights, environment maps)
     * are closely tied to an underlying bitmap data structure. This function
     * can be used to return this information for various purposes.
     *
     * When the class implementing this interface is a bitmap-backed texture,
     * this function directly returns the underlying bitmap. When it is procedural,
     * a bitmap version must first be generated. In this case, the parameter
     * \ref sizeHint is used to control the target size. The default
     * value <tt>-1, -1</tt> allows the implementation to choose a suitable
     * size by itself.
     *
     * \remark The default implementation throws an exception
     */
    virtual ref<Bitmap> getBitmap(const Vector2i &sizeHint = Vector2i(-1, -1)) const;

    /// Serialize this emitter to a binary data stream
    virtual void serialize(Stream *stream, InstanceManager *manager) const;

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Functionality related to environment emitters
    // =============================================================

    /// Is this an environment emitter? (e.g. an HDRI environment map?)
    inline bool isEnvironmentEmitter() const { return m_type & EEnvironmentEmitter; }

    /**
     * \brief Return the radiant emittance from an environment emitter
     *
     * This is function is called by unidirectional rendering techniques
     * (e.g. a path tracer) when no scene object has been intersected, and
     * the scene has been determined to contain an environment emitter.
     *
     * The default implementation throws an exception.
     *
     * \param ray
     *    Specifies the ray direction that should be queried
     */
    virtual Spectrum evalEnvironment(const RayDifferential &ray) const;

    /**
     * \brief Fill out a data record that can be used to query the direct
     * illumination sampling density of an environment emitter.
     *
     * This is function is mainly called by unidirectional rendering
     * techniques (e.g. a path tracer) when no scene object has been
     * intersected, and the (hypothetical) sampling density of the
     * environment emitter needs to be known by a multiple importance
     * sampling technique.
     *
     * The default implementation throws an exception.
     *
     * \param dRec
     *    The direct illumination sampling record to be filled
     *
     * \param ray
     *    Specifies the ray direction that should be queried
     *
     * \return \c true upon success
     */
    virtual bool fillDirectSamplingRecord(
            DirectSamplingRecord &dRec, const Ray &ray) const;

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Miscellaneous
    // =============================================================

    /// Is this a compound emitter consisting of several sub-objects?
    virtual bool isCompound() const;

    /**
     * \brief Return a sub-element of a compound emitter.
     *
     * When expanding emitters, the scene will repeatedly call this
     * function with increasing indices. Returning \a NULL indicates
     * that no more are available.
     */
    virtual Emitter *getElement(size_t index);

    //! @}
    // =============================================================

    MTS_DECLARE_CLASS()
protected:
    /// Construct a new emitter instance
    Emitter(const Properties &props);

    /// Unserialize a emitter instance from a binary data stream
    Emitter(Stream *stream, InstanceManager *manager);

    /// Virtual destructor
    virtual ~Emitter();
protected:
    Float m_samplingWeight;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_EMITTER_H_ */

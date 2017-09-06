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
#if !defined(__MITSUBA_RENDER_SCENE_H_)
#define __MITSUBA_RENDER_SCENE_H_

#include <mitsuba/core/netobject.h>
#include <mitsuba/core/pmf.h>
#include <mitsuba/core/aabb.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/render/skdtree.h>
#include <mitsuba/render/sensor.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/phase.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Principal scene data structure
 *
 * This class holds information on surfaces, emitters and participating media
 * and coordinates rendering jobs. It also provides useful query routines that
 * are mostly used by the \ref Integrator implementations.
 *
 * \ingroup librender
 * \ingroup libpython
 */
class MTS_EXPORT_RENDER Scene : public NetworkedObject {
public:
    // =============================================================
    //! @{ \name Initialization and rendering
    // =============================================================

    /// Construct a new, empty scene (with the default properties)
    Scene();

    /// Construct a new, empty scene
    Scene(const Properties &props);

    /// Create a shallow clone of a scene
    Scene(Scene *scene);

    /// Unserialize a scene from a binary data stream
    Scene(Stream *stream, InstanceManager *manager);

    /**
     * \brief Initialize the scene
     *
     * This function \a must be called before using any
     * of the methods in this class.
     */
    void initialize();

    /**
     *\brief Invalidate the kd-tree
     *
     * This function must be called if, after running \ref initialize(),
     * additional geometry is added to the scene.
     */
    void invalidate();

    /**
     * \brief Initialize the scene for bidirectional rendering algorithms.
     *
     * This ensures that certain "special" shapes (such as the aperture
     * of the sensor) are added to the scene. This function should be called
     * before using any of the methods in this class.
     */
    void initializeBidirectional();

    /**
     * \brief Perform any pre-processing steps before rendering
     *
     * This function should be called after \ref initialize() and
     * before rendering the scene. It might do a variety of things,
     * such as constructing photon maps or executing distributed overture
     * passes.
     *
     * Progress is tracked by sending status messages to a provided
     * render queue (the parameter \c job is required to discern multiple
     * render jobs occurring in parallel).
     *
     * The last three parameters are resource IDs of the associated scene,
     * sensor and sample generator, which have been made available to all
     * local and remote workers.
     *
     * \return \c true upon successful completion.
     */
    bool preprocess(RenderQueue *queue, const RenderJob *job,
        int sceneResID, int sensorResID, int samplerResID);

    /**
     * \brief Render the scene as seen by the scene's main sensor.
     *
     * Progress is tracked by sending status messages to a provided
     * render queue (the parameter \c job is required to discern multiple
     * render jobs occurring in parallel).
     *
     * The last three parameters are resource IDs of the associated scene,
     * sensor and sample generator, which have been made available to all
     * local and remote workers.
     *
     * \return \c true upon successful completion.
     */
    bool render(RenderQueue *queue, const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID);

    /**
     * \brief Perform any post-processing steps after rendering
     *
     * Progress is tracked by sending status messages to a provided
     * render queue (the parameter \c job is required to discern multiple
     * render jobs occurring in parallel).
     *
     * The last three parameters are resource IDs of the associated scene,
     * sensor and sample generator, which have been made available to all
     * local and remote workers.
     */
    void postprocess(RenderQueue *queue, const RenderJob *job,
        int sceneResID, int sensorResID, int samplerResID);

    /// Write out the current (partially rendered) image
    void flush(RenderQueue *queue, const RenderJob *job);

    /**
     * \brief Cancel a running rendering job
     *
     * This function can be called asynchronously, e.g. from a GUI.
     * In this case, \ref render() will quit with a return value of
     * \c false.
     */
    void cancel();

    /// Add a child node to the scene
    void addChild(const std::string &name, ConfigurableObject *child);

    /// Add an unnamed child
    inline void addChild(ConfigurableObject *child) { addChild("", child); }

    /** \brief Configure this object (called \a once after construction
       and addition of all child \ref ConfigurableObject instances).) */
    void configure();

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Ray tracing
    // =============================================================

    /**
     * \brief Intersect a ray against all primitives stored in the scene
     * and return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum
     *    extent information, as well as a time value (which applies
     *    when the shapes are in motion)
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \return \c true if an intersection was found
     */
    inline bool rayIntersect(const Ray &ray, Intersection &its) const {
        return m_kdtree->rayIntersect(ray, its);
    }

    /**
     * \brief Intersect a ray against all primitives stored in the scene
     * and return the traveled distance and intersected shape
     *
     * This function represents a performance improvement when the
     * intersected shape must be known, but there is no need for
     * a detailed intersection record.
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum
     *    extent information, as well as a time value (which applies
     *    when the shapes are in motion)
     *
     * \param t
     *    The traveled ray distance will be stored in this parameter

     * \param shape
     *    A pointer to the intersected shape will be stored in this
     *    parameter
     *
     * \param n
     *    The geometric surface normal will be stored in this parameter
     *
     * \param uv
     *    The UV coordinates associated with the intersection will
     *    be stored here.
     *
     * \return \c true if an intersection was found
     */
    inline bool rayIntersect(const Ray &ray, Float &t,
            ConstShapePtr &shape, Normal &n, Point2 &uv) const {
        return m_kdtree->rayIntersect(ray, t, shape, n, uv);
    }

    /**
     * \brief Intersect a ray against all primitives stored in the scene
     * and \a only determine whether or not there is an intersection.
     *
     * This is by far the fastest ray tracing method. This performance
     * improvement comes with a major limitation though: this function
     * cannot provide any additional information about the detected
     * intersection (not even its position).
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum
     *    extent information, as well as a time value (which applies
     *    when the shapes are in motion)
     *
     * \return \c true if an intersection was found
     */
    inline bool rayIntersect(const Ray &ray) const {
        return m_kdtree->rayIntersect(ray);
    }

    /**
     * \brief Return the transmittance between \c p1 and \c p2 at the
     * specified time.
     *
     * This function is essentially a continuous version of \ref isOccluded(),
     * which additionally accounts for the presence of participating media
     * and surface interactions that attenuate a ray without changing
     * its direction (i.e. geometry with an alpha mask)
     *
     * The implementation correctly handles arbitrary amounts of index-matched
     * medium transitions. The \c interactions parameter can be used to
     * specify a maximum number of possible surface interactions and medium
     * transitions between \c p1 and \c p2. When this number is exceeded,
     * the function returns zero.
     *
     * Note that index-mismatched boundaries (i.e. a transition from air to
     * water) are not supported by this function. The integrator needs to take
     * care of these in some other way.
     *
     * \param p1
     *     Source position
     * \param p2
     *     Target position
     * \param p1OnSurface
     *     Is the source position located on a surface? This information is
     *     necessary to set up the right ray epsilons for the kd-tree traversal
     * \param p2OnSurface
     *     Is the target position located on a surface?
     * \param medium
     *     The medium at \c p1
     * \param interactions
     *    Specifies the maximum permissible number of index-matched medium
     *    transitions or \ref BSDF::ENull scattering events on the way
     *    to the light source. (<tt>interactions<0</tt> means arbitrarily many).
     *    When the function returns a nonzero result, this parameter will
     *    additionally be used to return the actual number of intermediate
     *    interactions.
     * \param time
     *     Associated scene time value for the transmittance computation
     * \param sampler
     *     Optional: A sample generator. This may be used
     *     to compute a random unbiased estimate of the transmission.
     * \return An spectral-valued transmittance value with components
     *     between zero and one.
     */
    Spectrum evalTransmittance(const Point &p1, bool p1OnSurface,
        const Point &p2, bool p2OnSurface, Float time, const Medium *medium,
        int &interactions, Sampler *sampler = NULL) const;

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Ray tracing support for bidirectional algorithms
    // =============================================================

    /**
     * \brief Intersect a ray against all scene primitives \a and
     * "special" primitives, such as the aperture of a sensor.
     *
     * This function does exactly the same thing as \ref rayIntersect,
     * except that it additionally performs intersections against a
     * list of "special" shapes that are intentionally kept outside
     * of the main scene kd-tree (e.g. because they are not static
     * and might change from rendering to rendering). This is needed
     * by some bidirectional techniques that e.g. care about
     * intersections with the sensor aperture.
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum
     *    extent information, as well as a time value (which applies
     *    when the shapes are in motion)
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersectAll(const Ray &ray, Intersection &its) const;

    /**
     * \brief Intersect a ray against all normal and "special" primitives
     * and only return the traveled distance and intersected shape
     *
     * This function represents a performance improvement when the
     * intersected shape must be known, but there is no need for
     * a detailed intersection record.
     *
     * This function does exactly the same thing as \ref rayIntersect,
     * except that it additionally performs intersections against a
     * list of "special" shapes that are intentionally kept outside
     * of the main scene kd-tree (e.g. because they are not static
     * and might change from rendering to rendering). This is needed
     * by some bidirectional techniques that e.g. care about
     * intersections with the sensor aperture.
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum
     *    extent information, as well as a time value (which applies
     *    when the shapes are in motion)
     *
     * \param t
     *    The traveled ray distance will be stored in this parameter

     * \param shape
     *    A pointer to the intersected shape will be stored in this
     *    parameter
     *
     * \param n
     *    The geometric surface normal will be stored in this parameter
     *
     * \param uv
     *    The UV coordinates associated with the intersection will
     *    be stored here.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersectAll(const Ray &ray, Float &t,
            ConstShapePtr &shape, Normal &n, Point2 &uv) const;

    /**
     * \brief Intersect a ray against all normal and "special" primitives
     * and \a only determine whether or not there is an intersection.
     *
     * This is by far the fastest ray tracing method. This performance
     * improvement comes with a major limitation though: this function
     * cannot provide any additional information about the detected
     * intersection (not even its position).
     *
     * This function does exactly the same thing as \ref rayIntersect,
     * except that it additionally performs intersections against a
     * list of "special" shapes that are intentionally kept outside
     * of the main scene kd-tree (e.g. because they are not static
     * and might change from rendering to rendering). This is needed
     * by some bidirectional techniques that e.g. care about
     * intersections with the sensor aperture.
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum
     *    extent information, as well as a time value (which applies
     *    when the shapes are in motion)
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersectAll(const Ray &ray) const;

    /**
     * \brief Return the transmittance between \c p1 and \c p2 at the
     * specified time (and acount for "special" primitives).
     *
     * This function is essentially a continuous version of \ref isOccluded(),
     * which additionally accounts for the presence of participating media
     * and surface interactions that attenuate a ray without changing
     * its direction (i.e. geometry with an alpha mask)
     *
     * The implementation correctly handles arbitrary amounts of index-matched
     * medium transitions. The \c interactions parameter can be used to
     * specify a maximum number of possible surface interactions and medium
     * transitions between \c p1 and \c p2. When this number is exceeded,
     * the function returns zero.
     *
     * Note that index-mismatched boundaries (i.e. a transition from air to
     * water) are not supported by this function. The integrator needs to take
     * care of these in some other way.
     *
     * This function does exactly the same thing as \ref evalTransmittance,
     * except that it additionally performs intersections against a
     * list of "special" shapes that are intentionally kept outside
     * of the main scene kd-tree (e.g. because they are not static
     * and might change from rendering to rendering). This is needed
     * by some bidirectional techniques that care about intersections
     * with the sensor aperture, etc.
     *
     * \param p1
     *     Source position
     * \param p2
     *     Target position
     * \param p1OnSurface
     *     Is the source position located on a surface? This information is
     *     necessary to set up the right ray epsilons for the kd-tree traversal
     * \param p2OnSurface
     *     Is the target position located on a surface?
     * \param medium
     *     The medium at \c p1
     * \param interactions
     *    Specifies the maximum permissible number of index-matched medium
     *    transitions or \ref BSDF::ENull scattering events on the way
     *    to the light source. (<tt>interactions<0</tt> means arbitrarily many).
     *    When the function returns a nonzero result, this parameter will
     *    additionally be used to return the actual number of intermediate
     *    interactions.
     * \param time
     *     Associated scene time value for the transmittance computation
     * \param sampler
     *     Optional: A sample generator. This may be used
     *     to compute a random unbiased estimate of the transmission.
     * \return An spectral-valued transmittance value with components
     *     between zero and one.
     */
    Spectrum evalTransmittanceAll(const Point &p1, bool p1OnSurface,
        const Point &p2, bool p2OnSurface, Float time, const Medium *medium,
        int &interactions, Sampler *sampler = NULL) const;

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Direct sampling techniques
    // =============================================================

    /**
     * \brief Direct illumination sampling routine
     *
     * Given an arbitrary reference point in the scene, this method samples a
     * position on an emitter that has a nonzero contribution towards that point.
     *
     * Ideally, the implementation should importance sample the product of
     * the emission profile and the geometry term between the reference point
     * and the position on the emitter.
     *
     * \param dRec
     *    A direct illumination sampling record that specifies the
     *    reference point and a time value. After the function terminates,
     *    it will be populated with the position sample and related information
     *
     * \param sample
     *    A uniformly distributed 2D vector
     *
     * \param testVisibility
     *    When set to \c true, a shadow ray will be cast to ensure that the
     *    sampled emitter position and the reference point are mutually visible.
     *
     * \return
     *    An importance weight given by the radiance received along
     *    the sampled ray divided by the sample probability.
     */
    Spectrum sampleEmitterDirect(DirectSamplingRecord &dRec,
            const Point2 &sample, bool testVisibility = true) const;

    /**
     * \brief "Direct illumination" sampling routine for the main scene sensor
     *
     * Given an arbitrary reference point in the scene, this method samples a
     * position on an sensor that has a nonzero contribution towards that point.
     * This function can be interpreted as a generalization of a direct
     * illumination sampling strategy to sensors.
     *
     * Ideally, the implementation should importance sample the product of
     * the response profile and the geometry term between the reference point
     * and the position on the emitter.
     *
     * \param dRec
     *    A direct illumination sampling record that specifies the
     *    reference point and a time value. After the function terminates,
     *    it will be populated with the position sample and related information
     *
     * \param sample
     *    A uniformly distributed 2D vector
     *
     * \param testVisibility
     *    When set to \c true, a shadow ray will be cast to ensure that the
     *    sampled sensor position and the reference point are mutually visible.
     *
     * \return
     *    An importance weight given by the importance emitted along
     *    the sampled ray divided by the sample probability.
     */
    Spectrum sampleSensorDirect(DirectSamplingRecord &dRec,
            const Point2 &sample, bool testVisibility = true) const;

    /**
     * \brief Direct illumination sampling with support for participating
     * media (medium variant)
     *
     * Given an arbitrary reference point in the scene, this method samples a
     * position on an emitter that has a nonzero contribution towards that point.
     * In comparison to \ref sampleEmitterDirect, this version also accounts for
     * attenuation by participating media and should be used when \c dRec.p
     * lies \a inside a medium, i.e. \a not on a surface!
     *
     * Ideally, the implementation should importance sample the product of
     * the emission profile and the geometry term between the reference point
     * and the position on the emitter.
     *
     * \param dRec
     *    A direct illumination sampling record that specifies the
     *    reference point and a time value. After the function terminates,
     *    it will be populated with the position sample and related information
     *
     * \param medium
     *    The medium located at the reference point (or \c NULL for vacuum).
     *
     * \param interactions
     *    Specifies the maximum permissible number of index-matched medium
     *    transitions or \ref BSDF::ENull scattering events on the way
     *    to the light source. (<tt>interactions<0</tt> means arbitrarily many).
     *    When the function returns a nonzero result, this parameter will
     *    additionally be used to return the actual number of intermediate
     *    interactions.
     *
     * \param sample
     *    A uniformly distributed 2D vector
     *
     * \param sampler
     *    Optional: a pointer to a sample generator. Some particular
     *    implementations can do a better job at sampling when they have
     *    access to additional random numbers.
     *
     * \return
     *    An importance weight given by the radiance received along
     *    the sampled ray divided by the sample probability.
     */
    Spectrum sampleAttenuatedEmitterDirect(DirectSamplingRecord &dRec,
            const Medium *medium, int &interactions, const Point2 &sample,
            Sampler *sampler = NULL) const;

    /**
     * \brief "Direct illumination" sampling routine for the main scene sensor
     * with support for participating media (medium variant)
     *
     * Given an arbitrary reference point in the scene, this method samples a
     * position on an sensor that has a nonzero response towards that point.
     * In comparison to \ref sampleSensorDirect, this version also accounts for
     * attenuation by participating media and should be used when \c dRec.p
     * lies \a inside a medium, i.e. \a not on a surface!
     * This function can be interpreted as a generalization of a direct
     * illumination sampling strategy to sensors.
     *
     * Ideally, the implementation should importance sample the product of
     * the response profile and the geometry term between the reference point
     * and the position on the sensor.
     *
     * \param dRec
     *    A direct illumination sampling record that specifies the
     *    reference point and a time value. After the function terminates,
     *    it will be populated with the position sample and related information
     *
     * \param medium
     *    The medium located at the reference point (or \c NULL for vacuum).
     *
     * \param interactions
     *    Specifies the maximum permissible number of index-matched medium
     *    transitions or \ref BSDF::ENull scattering events on the way
     *    to the light source. (<tt>interactions<0</tt> means arbitrarily many).
     *    When the function returns a nonzero result, this parameter will
     *    additionally be used to return the actual number of intermediate
     *    interactions.
     *
     * \param sample
     *    A uniformly distributed 2D vector
     *
     * \param sampler
     *    Optional: a pointer to a sample generator. Some particular
     *    implementations can do a better job at sampling when they have
     *    access to additional random numbers.
     *
     * \return
     *    An importance weight given by the radiance received along
     *    the sampled ray divided by the sample probability.
     */
    Spectrum sampleAttenuatedSensorDirect(DirectSamplingRecord &dRec,
            const Medium *medium, int &interactions, const Point2 &sample,
            Sampler *sampler = NULL) const;

    /**
     * \brief Direct illumination sampling with support for participating
     * media (surface variant)
     *
     * Given an arbitrary reference point in the scene, this method samples a
     * position on an emitter that has a nonzero contribution towards that point.
     * In comparison to \ref sampleEmitterDirect, this version also accounts for
     * attenuation by participating media and should be used when the target
     * position lies on a surface.
     *
     * Ideally, the implementation should importance sample the product of
     * the emission profile and the geometry term between the reference point
     * and the position on the emitter.
     *
     * \param dRec
     *    A direct illumination sampling record that specifies the
     *    reference point and a time value. After the function terminates,
     *    it will be populated with the position sample and related information
     *
     * \param its
     *    An intersection record associated with the reference point in
     *    \c dRec. This record is needed to determine the participating
     *    medium between the emitter sample and the reference point
     *    when \c its marks a medium transition.
     *
     * \param medium
     *    The medium located at \c its (or \c NULL for vacuum). When the shape
     *    associated with \c its marks a medium transition, it does not matter
     *    which of the two media is specified.
     *
     * \param interactions
     *    Specifies the maximum permissible number of index-matched medium
     *    transitions or \ref BSDF::ENull scattering events on the way
     *    to the light source. (<tt>interactions<0</tt> means arbitrarily many).
     *    When the function returns a nonzero result, this parameter will
     *    additionally be used to return the actual number of intermediate
     *    interactions.
     *
     * \param sample
     *    A uniformly distributed 2D vector
     *
     * \param sampler
     *    Optional: a pointer to a sample generator. Some particular
     *    implementations can do a better job at sampling when they have
     *    access to additional random numbers.
     *
     * \return
     *    An importance weight given by the radiance received along
     *    the sampled ray divided by the sample probability.
     */
    Spectrum sampleAttenuatedEmitterDirect(DirectSamplingRecord &dRec,
            const Intersection &its, const Medium *medium, int &interactions,
            const Point2 &sample, Sampler *sampler = NULL) const;

    /**
     * \brief "Direct illumination" sampling routine for the main scene sensor
     * with support for participating media (surface variant)
     *
     * Given an arbitrary reference point in the scene, this method samples a
     * position on an sensor that has a nonzero response towards that point.
     * In comparison to \ref sampleSensorDirect, this version also accounts for
     * attenuation by participating media and should be used when the target
     * position lies on a surface.
     *
     * Ideally, the implementation should importance sample the product of
     * the emission profile and the geometry term between the reference point
     * and the position on the sensor.
     *
     * \param dRec
     *    A direct illumination sampling record that specifies the
     *    reference point and a time value. After the function terminates,
     *    it will be populated with the position sample and related information
     *
     * \param its
     *    An intersection record associated with the reference point in
     *    \c dRec. This record is needed to determine the participating
     *    medium between the sensor sample and the reference point
     *    when \c its marks a medium transition.
     *
     * \param medium
     *    The medium located at \c its (or \c NULL for vacuum). When the shape
     *    associated with \c its marks a medium transition, it does not matter
     *    which of the two media is specified.
     *
     * \param interactions
     *    Specifies the maximum permissible number of index-matched medium
     *    transitions or \ref BSDF::ENull scattering events on the way
     *    to the light source. (<tt>interactions<0</tt> means arbitrarily many).
     *    When the function returns a nonzero result, this parameter will
     *    additionally be used to return the actual number of intermediate
     *    interactions.
     *
     * \param sample
     *    A uniformly distributed 2D vector
     *
     * \param sampler
     *    Optional: a pointer to a sample generator. Some particular
     *    implementations can do a better job at sampling when they have
     *    access to additional random numbers.
     *
     * \return
     *    An importance weight given by the radiance received along
     *    the sampled ray divided by the sample probability.
     */
    Spectrum sampleAttenuatedSensorDirect(DirectSamplingRecord &dRec,
            const Intersection &its, const Medium *medium, int &interactions,
            const Point2 &sample, Sampler *sampler = NULL) const;

    /**
     * \brief Evaluate the probability density of the \a direct sampling
     * method implemented by the \ref sampleEmitterDirect() method.
     *
     * \param dRec
     *    A direct sampling record, which specifies the query
     *    location. Note that this record need not be completely
     *    filled out. The important fields are \c p, \c n, \c ref,
     *    \c dist, \c d, \c measure, and \c uv.
     *
     * \param p
     *    The world-space position that would have been passed to \ref
     *    sampleEmitterDirect()
     *
     * \return
     *    The density expressed with respect to the requested measure
     *    (usually \ref ESolidAngle)
     */
    Float pdfEmitterDirect(const DirectSamplingRecord &dRec) const;

    /**
     * \brief Evaluate the probability density of the \a direct sampling
     * method implemented by the \ref sampleSensorDirect() method.
     *
     * \param dRec
     *    A direct sampling record, which specifies the query
     *    location. Note that this record need not be completely
     *    filled out. The important fields are \c p, \c n, \c ref,
     *    \c dist, \c d, \c measure, and \c uv.
     *
     * \param p
     *    The world-space position that would have been passed to \ref
     *    sampleSensorDirect()
     *
     * \return
     *    The density expressed with respect to the requested measure
     *    (usually \ref ESolidAngle)
     */
    Float pdfSensorDirect(const DirectSamplingRecord &dRec) const;

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Emission sampling techniques
    // =============================================================

    /**
     * \brief Sample a position according to the emission profile
     * defined by the emitters in the scene.
     *
     * To sample the directional component, please use the
     * \ref Emitter::sampleDirection() method.
     *
     * \param pRec
     *    A position record to be populated with the sampled
     *    position and related information
     *
     * \param sample
     *    A uniformly distributed 2D vector
     *
     * \return
     *    An importance weight associated with the sampled position.
     *    This accounts for the difference in the spatial part of the
     *    emission profile and the density function.
     */
    Spectrum sampleEmitterPosition(PositionSamplingRecord &pRec,
        const Point2 &sample) const;

    /**
     * \brief Sample a position on the main sensor of the scene.
     *
     * This function is provided here mainly for symmetry
     * with respect to \ref sampleEmitterPosition().
     *
     * To sample the directional component, please use the
     * \ref Sensor::sampleDirection() method.
     *
     * \param pRec
     *    A position record to be populated with the sampled
     *    position and related information
     *
     * \param sample
     *    A uniformly distributed 2D vector
     *
     * \param extra
     *    An additional 2D vector provided to the sampling
     *    routine -- its use is implementation-dependent.
     *
     * \return
     *    An importance weight associated with the sampled position.
     *    This accounts for the difference in the spatial part of the
     *    response profile and the density function.
     */
    inline Spectrum sampleSensorPosition(PositionSamplingRecord &pRec,
        const Point2 &sample, const Point2 *extra = NULL) const {
        pRec.object = m_sensor.get();
        return m_sensor->samplePosition(pRec, sample, extra);
    }

    /**
     * \brief Evaluate the spatial component of the sampling density
     * implemented by the \ref sampleEmitterPosition() method
     *
     * \param pRec
     *    A position sampling record, which specifies the query location
     *
     * \return
     *    The area density at the supplied position
     */
    Float pdfEmitterPosition(const PositionSamplingRecord &pRec) const;

    /**
     * \brief Evaluate the spatial component of the sampling density
     * implemented by the \ref sampleSensorPosition() method
     *
     * \param pRec
     *    A position sampling record, which specifies the query location
     *
     * \return
     *    The area density at the supplied position
     */
    inline Float pdfSensorPosition(const PositionSamplingRecord &pRec) const {
        return m_sensor->pdfPosition(pRec);
    }

    /**
     * \brief Return the discrete probability of choosing a
     * certain emitter in <tt>sampleEmitter*</tt>
     */
    inline Float pdfEmitterDiscrete(const Emitter *emitter) const {
        return emitter->getSamplingWeight() * m_emitterPDF.getNormalization();
    }

    /**
     * \brief Importance sample a ray according to the emission profile
     * defined by the sensors in the scene
     *
     * This function combines both steps of choosing a ray origin and
     * direction value. It does not return any auxiliary sampling
     * information and is mainly meant to be used by unidirectional
     * rendering techniques.
     *
     * Note that this function potentially uses a different sampling
     * strategy compared to the sequence of running \ref sampleEmitterPosition()
     * and \ref Emitter::sampleDirection(). The reason for this is that it may
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
    Spectrum sampleEmitterRay(Ray &ray,
        const Emitter* &emitter,
        const Point2 &spatialSample,
        const Point2 &directionalSample,
        Float time) const;

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Environment emitters
    // =============================================================

    /// Return the scene's environment emitter (if there is one)
    inline const Emitter *getEnvironmentEmitter() const { return m_environmentEmitter.get(); }

    /// Does the scene have a environment emitter?
    inline bool hasEnvironmentEmitter() const { return m_environmentEmitter.get() != NULL; }

    /**
     * \brief Return the environment radiance for a ray that did not intersect
     * any of the scene objects.
     *
     * This is primarily meant for path tracing-style integrators.
     */
    inline Spectrum evalEnvironment(const RayDifferential &ray) const {
        return hasEnvironmentEmitter() ?
            m_environmentEmitter->evalEnvironment(ray) : Spectrum(0.0f);
    }

    /**
     * \brief Return the environment radiance for a ray that did not intersect
     * any of the scene objects. This method additionally considers
     * transmittance by participating media
     *
     * This is primarily meant for path tracing-style integrators.
     */
    inline Spectrum evalAttenuatedEnvironment(const RayDifferential &ray,
            const Medium *medium, Sampler *sampler) const {
        if (!m_environmentEmitter)
            return Spectrum(0.0f);
        Spectrum result = evalEnvironment(ray);
        if (medium)
            result *= medium->evalTransmittance(ray, sampler);
        return result;
    }

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Miscellaneous
    // =============================================================

    /// Return an axis-aligned bounding box containing the whole scene
    inline const AABB &getAABB() const {
        return m_aabb;
    }

    /**
     * \brief Is the main scene sensor degenerate?  (i.e. has it
     * collapsed to a point or line)
     *
     * Note that this function only cares about the spatial component
     * of the sensor -- its value does not depend on whether the directional
     * response function is degenerate.
     */
    inline bool hasDegenerateSensor() const { return m_degenerateSensor; }

    /**
     * \brief Area \a all emitters in this scene degenerate?
     * (i.e. they has collapsed to a point or line)
     *
     * Note that this function only cares about the spatial component
     * of the emitters -- its value does not depend on whether the
     * directional emission profile is degenerate.
     */
    inline bool hasDegenerateEmitters() const { return m_degenerateEmitters; }

    /// Return a bounding sphere containing the whole scene
    inline BSphere getBSphere() const {
        // todo: switch to something smarter at some point
        return m_aabb.getBSphere();
    }

    /// Does the scene contain participating media?
    inline bool hasMedia() const { return !m_media.empty(); }

    /**
     * \brief Set the main scene sensor.
     *
     * Note that the main sensor is not included when this Scene instance
     * is serialized -- the sensor field will be \c NULL after
     * unserialization. This is intentional so that the sensor can
     * be changed without having to re-transmit the whole scene.
     * Hence, it needs to be submitted separately and re-attached
     * on the remote side using \ref setSensor().
     **/
    void setSensor(Sensor *sensor);

    /// \brief Remove a sensor from the scene's sensor list
    void removeSensor(Sensor *sensor);

    /// \brief Add a sensor to the scene's sensor list
    void addSensor(Sensor *sensor);

    /// Return the scene's sensor
    inline Sensor *getSensor() { return m_sensor; }

    /// Return the scene's sensor (const version)
    inline const Sensor *getSensor() const { return m_sensor.get(); }

    /**
     * \brief Return the list of sensors that are specified
     * by the scene.
     *
     * As scene can have multiple sensors -- however, during
     * a rendering, there will always be one "main" sensor that
     * is currently active.
     *
     * \sa getSensor
     */
    inline ref_vector<Sensor> &getSensors() { return m_sensors; }

    /**
     * \brief Return the list of sensors that are specified
     * by the scene (const version)
     *
     * As scene can have multiple sensors -- however, during
     * a rendering, there will always be one "main" sensor that
     * is currently active.
     *
     * \sa getSensor
     */
    inline const ref_vector<Sensor> &getSensors() const { return m_sensors; }

    /**
     * \brief Set the scene's integrator.
     *
     * Note that the integrator is not included when this Scene instance
     * is serialized -- the integrator field will be \c NULL after
     * unserialization. This is intentional so that the integrator can
     * be changed without having to re-transmit the whole scene. Hence,
     * the integrator needs to be submitted separately and re-attached
     * on the remote side using \ref setIntegrator().
     **/
    inline void setIntegrator(Integrator *integrator) { m_integrator = integrator; }

    /// Return the scene's integrator
    inline Integrator *getIntegrator() { return m_integrator; }
    /// Return the scene's integrator (const version)
    inline const Integrator *getIntegrator() const { return m_integrator.get(); }

    /**
     * \brief Set the scene's sampler.
     *
     * Note that the sampler is not included when this Scene instance
     * is serialized -- the sampler field will be \c NULL after
     * unserialization. This is intentional so that the sampler can
     * be changed without having to re-transmit the whole scene.
     * Hence, the sampler needs to be submitted separately
     * and re-attached on the remote side using \ref setSampler().
     **/
    inline void setSampler(Sampler *sampler) { m_sampler = sampler; }

    /**
     * \brief Return the scene's sampler.
     *
     * Note that when rendering using multiple different threads, each
     * thread will be passed a shallow copy of the scene, which has a
     * different sampler instance. This helps to avoid locking/contention
     * issues and ensures that different threads render with different
     * random number sequences. The sampler instance provided here is a
     * clone of the original sampler specified in the sensor.
     */
    inline Sampler *getSampler() { return m_sampler; }
    /// Return the scene's sampler
    inline const Sampler *getSampler() const { return m_sampler.get(); }

    /// Return the scene's film
    inline Film *getFilm() { return m_sensor->getFilm(); }
    /// Return the scene's film
    inline const Film *getFilm() const { return m_sensor->getFilm(); }

    /// Return the scene's kd-tree accelerator
    inline ShapeKDTree *getKDTree() { return m_kdtree; }
    /// Return the scene's kd-tree accelerator
    inline const ShapeKDTree *getKDTree() const { return m_kdtree.get(); }

    /// Return the a list of all subsurface integrators
    inline ref_vector<Subsurface> &getSubsurfaceIntegrators() { return m_ssIntegrators; }
    /// Return the a list of all subsurface integrators
    inline const ref_vector<Subsurface> &getSubsurfaceIntegrators() const { return m_ssIntegrators; }

    /// Return the scene's triangular meshes (a subset of \ref getShapes())
    inline std::vector<TriMesh *> &getMeshes() { return m_meshes; }
    /// Return the scene's triangular meshes (a subset of \ref getShapes())
    inline const std::vector<TriMesh *> &getMeshes() const { return m_meshes; }
    /// Return the scene's normal shapes (including triangular meshes)
    inline ref_vector<Shape> &getShapes() { return m_shapes; }
    /// Return the scene's normal shapes (including triangular meshes)
    inline const ref_vector<Shape> &getShapes() const { return m_shapes; }

    /// Return a set of special shapes related to emitter/sensor geometry in bidirectional renderings
    inline ref_vector<Shape> &getSpecialShapes() { return m_specialShapes; }
    /// Return a set of special shapes related to emitter/sensor geometry in bidirectional renderings
    inline const ref_vector<Shape> &getSpecialShapes() const { return m_specialShapes; }

    /// Return the scene's emitters
    inline ref_vector<Emitter> &getEmitters() { return m_emitters; }
    /// Return the scene's emitters
    inline const ref_vector<Emitter> &getEmitters() const { return m_emitters; }
    /// Return the scene's participating media
    inline ref_vector<Medium> &getMedia() { return m_media; }
    /// Return the scene's participating media
    inline const ref_vector<Medium> &getMedia() const { return m_media; }
    /// Return referenced objects (such as textures, BSDFs)
    inline ref_vector<ConfigurableObject> &getReferencedObjects() { return m_objects; }
    /// Return referenced objects (such as textures, BSDFs)
    inline const ref_vector<ConfigurableObject> &getReferencedObjects() const { return m_objects; }

    /// Return the name of the file containing the original description of this scene
    inline const fs::path &getSourceFile() const { return *m_sourceFile; }
    /// Set the name of the file containing the original description of this scene
    void setSourceFile(const fs::path &name);
    /// Return the render output filename
    inline const fs::path &getDestinationFile() const { return *m_destinationFile; }
    /// Set the render output filename
    void setDestinationFile(const fs::path &name);

    /// Does the destination file already exist?
    inline bool destinationExists() const { return m_sensor->getFilm()->destinationExists(*m_destinationFile); }

    /// Set the block resolution used to split images into parallel workloads
    inline void setBlockSize(uint32_t size) { m_blockSize = size; }
    /// Return the block resolution used to split images into parallel workloads
    inline uint32_t getBlockSize() const { return m_blockSize; }

    /// Serialize the whole scene to a network/file stream
    void serialize(Stream *stream, InstanceManager *manager) const;

    /* NetworkedObject implementation */
    void bindUsedResources(ParallelProcess *proc) const;
    void wakeup(ConfigurableObject *parent,
        std::map<std::string, SerializableObject *> &params);

    /// Return a string representation
    std::string toString() const;

    //! @}
    // =============================================================

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~Scene();

    /// \cond
    /// Add a shape to the scene
    void addShape(Shape *shape);
    /// \endcond
private:
    ref<ShapeKDTree> m_kdtree;
    ref<Sensor> m_sensor;
    ref<Integrator> m_integrator;
    ref<Sampler> m_sampler;
    ref<Emitter> m_environmentEmitter;
    ref_vector<Shape> m_shapes;
    ref_vector<Shape> m_specialShapes;
    ref_vector<Sensor> m_sensors;
    ref_vector<Emitter> m_emitters;
    ref_vector<ConfigurableObject> m_objects;
    ref_vector<NetworkedObject> m_netObjects;
    ref_vector<Subsurface> m_ssIntegrators;
    ref_vector<Medium> m_media;
    std::vector<TriMesh *> m_meshes;
    fs::path *m_sourceFile;
    fs::path *m_destinationFile;
    DiscreteDistribution m_emitterPDF;
    AABB m_aabb;
    uint32_t m_blockSize;
    bool m_degenerateSensor;
    bool m_degenerateEmitters;
};

MTS_NAMESPACE_END

#include <mitsuba/render/records.inl>

#endif /* __MITSUBA_RENDER_SCENE_H_ */

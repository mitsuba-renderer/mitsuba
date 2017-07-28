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
#if !defined(__MITSUBA_RENDER_INTEGRATOR_H_)
#define __MITSUBA_RENDER_INTEGRATOR_H_

#include <mitsuba/core/netobject.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Abstract integrator base-class; does not make any assumptions on
 * how radiance is computed.
 *
 * In Mitsuba, the different rendering techniques are collectively referred to as
 * \a integrators, since they perform integration over a high-dimensional
 * space. Each integrator represents a specific approach for solving
 * the light transport equation---usually favored in certain scenarios, but
 * at the same time affected by its own set of intrinsic limitations.
 * Therefore, it is important to carefully select an integrator based on
 * user-specified accuracy requirements and properties of the scene to be
 * rendered.
 *
 * This is the base class of all integrators; it does not make any assumptions on
 * how radiance is computed, which allows for many different kinds of implementations
 * ranging from software-based path tracing and Markov-Chain based techniques such
 * as Metropolis Light Transport up to hardware-accelerated rasterization.
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER Integrator : public NetworkedObject {
public:
    /**
     * \brief Possibly perform a pre-process task.
     *
     * This function is called automatically before the main rendering process;
     * the default implementation does nothing.
     *
     * The last three parameters are resource IDs of the associated scene,
     * sensor and sample generator, which have been made available to all
     * local and remote workers.i
     */
    virtual bool preprocess(const Scene *scene, RenderQueue *queue,
        const RenderJob *job, int sceneResID, int sensorResID,
        int samplerResID);

    /**
     * \brief Render the scene as seen by the default sensor.
     *
     * Progress is tracked by sending status messages to a provided render queue.
     * The parameter \c job is required to discern multiple render jobs occurring in
     * parallel. The last three parameters are resource IDs of the associated
     * scene, sensor and sample generator, which have been made available to
     * all local and remote workers. Returns true upon successful completion.
     */
    virtual bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
        int sceneResID, int sensorResID, int samplerResID) = 0;

    /**
     * \brief Cancel a running render job
     *
     * This function can be called asynchronously to cancel a running render
     * job. In this case, \ref render() will quit with a return value of
     * \c false.
     */
    virtual void cancel() = 0;

    /**
     * \brief Possibly perform a post-process task.
     *
     * This function is called automatically before the main rendering process;
     * the default implementation does nothing.
     *
     * The last three parameters are resource IDs of the associated scene,
     * sensor and sample generator, which have been made available to all
     * local and remote workers.i
     */
    virtual void postprocess(const Scene *scene, RenderQueue *queue,
        const RenderJob *job, int sceneResID, int sensorResID,
        int samplerResID);

    /**
     * \brief Configure the sample generator for use with this integrator
     *
     * This function is called once after instantiation and can be used to
     * inform the sampler implementation about specific sample requirements
     * of this integrator.
     */
    virtual void configureSampler(const Scene *scene, Sampler *sampler);

    /**
     * \brief Return the nested integrator (if any)
     *
     * When the integrator contains a nested integrator, this function can
     * be used to query for it
     */
    virtual const Integrator *getSubIntegrator(int index) const;

    /// Serialize this integrator to a binary data stream
    void serialize(Stream *stream, InstanceManager *manager) const;

    MTS_DECLARE_CLASS()
protected:
    /// Create a integrator
    Integrator(const Properties &props);

    /// Unserialize an integrator
    Integrator(Stream *stream, InstanceManager *manager);

    /// Virtual destructor
    virtual ~Integrator() { }
};

/**
 * \brief Radiance query record data structure used by \ref SamplingIntegrator
 * \ingroup librender
 */
struct MTS_EXPORT_RENDER RadianceQueryRecord {
public:
    /// List of suported query types. These can be combined by a binary OR.
    enum ERadianceQuery {
        /// Emitted radiance from a luminaire intersected by the ray
        EEmittedRadiance          = 0x0001,

        /// Emitted radiance from a subsurface integrator */
        ESubsurfaceRadiance       = 0x0002,

        /// Direct (surface) radiance */
        EDirectSurfaceRadiance    = 0x0004,

        /*! \brief Indirect (surface) radiance, where the last bounce did not go
            through a Dirac delta BSDF */
        EIndirectSurfaceRadiance  = 0x0008,

        /*! \brief Indirect (surface) radiance, where the last bounce went
           through a Dirac delta BSDF */
        ECausticRadiance          = 0x0010,

        /// In-scattered radiance due to volumetric scattering (direct)
        EDirectMediumRadiance     = 0x0020,

        /// In-scattered radiance due to volumetric scattering (indirect)
        EIndirectMediumRadiance   = 0x0040,

        /// Distance to the next surface intersection
        EDistance                 = 0x0080,

        /*! \brief Store an opacity value, which is equal to 1 when a shape
           was intersected and 0 when the ray passes through empty space.
           When there is a participating medium, it can also take on fractional
           values. */
        EOpacity                  = 0x0100,

        /*! \brief A ray intersection may need to be performed. This can be set to
           zero if the caller has already provided the intersection */
        EIntersection             = 0x0200,

        /* Radiance from volumes */
        EVolumeRadiance           = EDirectMediumRadiance | EIndirectMediumRadiance,

        /// Radiance query without emitted radiance, ray intersection required
        ERadianceNoEmission       = ESubsurfaceRadiance | EDirectSurfaceRadiance
            | EIndirectSurfaceRadiance | ECausticRadiance | EDirectMediumRadiance
            | EIndirectMediumRadiance | EIntersection,

        /// Default radiance query, ray intersection required
        ERadiance                 = ERadianceNoEmission | EEmittedRadiance,

        /// Radiance + opacity
        ESensorRay                = ERadiance | EOpacity
    };

    /// Additional flags that can be specified in the \ref extra field
    enum EExtraFlags {
        /// This is a query by an irradiance cache
        ECacheQuery    = 0x01,
        /// This is a query by an adaptive integrator
        EAdaptiveQuery = 0x02
    };

    /// Construct an invalid radiance query record
    inline RadianceQueryRecord()
     : type(0), scene(NULL), sampler(NULL), medium(NULL),
       depth(0), alpha(0), dist(-1), extra(0) {
    }

    /// Construct a radiance query record for the given scene and sampler
    inline RadianceQueryRecord(const Scene *scene, Sampler *sampler)
     : type(0), scene(scene), sampler(sampler), medium(NULL),
       depth(0), alpha(0), dist(-1), extra(0) {
    }

    /// Copy constructor
    inline RadianceQueryRecord(const RadianceQueryRecord &rRec)
     : type(rRec.type), scene(rRec.scene), sampler(rRec.sampler), medium(rRec.medium),
       depth(rRec.depth), alpha(rRec.alpha), dist(rRec.dist), extra(rRec.extra) {
    }

    /// Begin a new query of the given type
    inline void newQuery(int _type, const Medium *_medium) {
        type = _type;
        medium = _medium;
        depth = 1;
        extra = 0;
        alpha = 1;
    }

    /// Initialize the query record for a recursive query
    inline void recursiveQuery(const RadianceQueryRecord &parent, int _type) {
        type = _type;
        scene = parent.scene;
        sampler = parent.sampler;
        depth = parent.depth+1;
        medium = parent.medium;
        extra = parent.extra;
    }

    /// Initialize the query record for a recursive query
    inline void recursiveQuery(const RadianceQueryRecord &parent) {
        type = parent.type | EIntersection;
        scene = parent.scene;
        sampler = parent.sampler;
        depth = parent.depth+1;
        medium = parent.medium;
        extra = parent.extra;
    }

    /**
     * \brief Search for a ray intersection
     *
     * This function does several things at once: if the
     * intersection has already been provided, it returns.
     *
     * Otherwise, it
     * 1. performs the ray intersection
     * 2. computes the transmittance due to participating media
     *   and stores it in \c transmittance.
     * 3. sets the alpha value (if \c EAlpha is set in \c type)
     * 4. sets the distance value (if \c EDistance is set in \c type)
     * 5. clears the \c EIntersection flag in \c type
     *
     * \return \c true if there is a valid intersection.
     */
    inline bool rayIntersect(const RayDifferential &ray);

    /// Retrieve a 2D sample
    inline Point2 nextSample2D();

    /// Retrieve a 1D sample
    inline Float nextSample1D();

    /// Return a string representation
    std::string toString() const;
public:
    // An asterisk (*) marks entries, which may be overwritten
    // by the callee.

    /// Query type (*)
    int type;

    /// Pointer to the associated scene
    const Scene *scene;

    /// Sample generator
    Sampler *sampler;

    /// Pointer to the current medium (*)
    const Medium *medium;

    /// Current depth value (# of light bounces) (*)
    int depth;

    /// Surface interaction data structure (*)
    Intersection its;

    /// Opacity value of the associated pixel (*)
    Float alpha;

    /**
     * Ray distance to the first surface interaction
     * (if requested by the query type EDistance) (*)
     */
    Float dist;

    /**
     * Internal flag, which can be used to pass additional information
     * amonst recursive calls inside an integrator. The use
     * is dependent on the particular integrator implementation. (*)
     */
    int extra;
};

/** \brief Abstract base class, which describes integrators
 * capable of computing samples of the scene's radiance function.
 * \ingroup librender
 */
class MTS_EXPORT_RENDER SamplingIntegrator : public Integrator {
public:
    /**
     * \brief Sample the incident radiance along a ray. Also requires
     * a radiance query record, which makes this request more precise.
     */
    virtual Spectrum Li(const RayDifferential &ray,
        RadianceQueryRecord &rRec) const = 0;

    /**
     * \brief Estimate the irradiance at a given surface point
     *
     * The default implementation simply samples the hemisphere using
     * cosine-weighted sampling and a configurable number of rays.
     * An integrator such as irradiance caching will provide something
     * smarter.
     *
     * \param scene
     *     Const pointer to the underlying scene
     * \param its
     *     Describes the surface location where the irradiance is to be computed
     * \param medium
     *     Const pointer to the medium that encloses the ray
     *     <tt>(its.p, its.shFrame.n)</tt>. A value of \c NULL corresponds
     *     to vacuum.
     * \param sampler
     *     A pointer to a sample generator
     * \param nSamples
     *     How many samples should be taken
     * \param includeIndirect
     *     Include indirect illumination in the estimate?
     */
    virtual Spectrum E(const Scene *scene, const Intersection &its,
        const Medium *medium, Sampler *sampler, int nSamples,
        bool includeIndirect) const;

    /**
     * \brief Perform the main rendering task
     *
     * The work is automatically parallelized to multiple cores and
     * remote machines. The default implementation uniformly generates
     * samples on the sensor aperture and image plane as specified by
     * the used sampler. The average of the estimated radiance along the
     * associated rays in a pixel region is then taken as an approximation
     * of that pixel's radiance value. For adaptive strategies, have a look at
     * the \c adaptive plugin, which is an extension of this class.
     */
    bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
        int sceneResID, int sensorResID, int samplerResID);

    /**
     * This can be called asynchronously to cancel a running render job.
     * In this case, <tt>render()</tt> will quit with a return value of
     * <tt>false</tt>.
     */
    void cancel();

    /**
     * This method does the main work of <tt>render()</tt> and
     * runs in parallel for a series of image blocks, which are
     * being processed at a time.
     *
     * \param scene
     *    Pointer to the underlying scene
     * \param sensor
     *    Pointer to the sensor used to render the image
     * \param sampler
     *    Pointer to the sampler used to render the image
     * \param block
     *    Pointer to the image block to be filled
     * \param points
     *    Specifies the traversal order, i.e. using a space-filling
     *    curve. To limit the size of the array, it is currently assumed
     *    that the block size is smaller than 256x256
     * \param stop
     *    Reference to a boolean, which will be set to true when
     *    the user has requested that the program be stopped
     */
    virtual void renderBlock(const Scene *scene, const Sensor *sensor,
        Sampler *sampler, ImageBlock *block, const bool &stop,
        const std::vector< TPoint2<uint8_t> > &points) const;

    /**
     * <tt>NetworkedObject</tt> implementation:
     * When a parallel rendering process starts, the integrator is
     * given the opportunity to attach globally shared resources to
     * the process. This is useful for distributing heavy data
     * structures (e.g. photon maps) without having to re-transmit
     * them every time an image is rendered.
     */
    virtual void bindUsedResources(ParallelProcess *proc) const;

    /**
     * <tt>NetworkedObject</tt> implementation:
     * Called once just before this integrator instance is asked
     * to process an image block. In comparison to <tt>preprocess()</tt>
     * this will be executed on _every_ instance of this class, which is
     * useful for connecting to globally shared resources (photon maps,
     * irradiance caches, ..) after having been unserialized on a
     * remote machine. A list of resources bound to the associated
     * parallel process is given as a parameter.
     */
    virtual void wakeup(ConfigurableObject *parent,
        std::map<std::string, SerializableObject *> &params);

    /// Serialize this integrator to a binary data stream
    void serialize(Stream *stream, InstanceManager *manager) const;

    MTS_DECLARE_CLASS()
protected:
    /// Create a integrator
    SamplingIntegrator(const Properties &props);

    /// Unserialize an integrator
    SamplingIntegrator(Stream *stream, InstanceManager *manager);

    /// Virtual destructor
    virtual ~SamplingIntegrator() { }
protected:
    /// Used to temporarily cache a parallel process while it is in operation
    ref<ParallelProcess> m_process;
};

/*
 * \brief Base class of all recursive Monte Carlo integrators, which compute
 * unbiased solutions to the rendering equation (and optionally
 * the radiative transfer equation).
 * \ingroup librender
 */
class MTS_EXPORT_RENDER MonteCarloIntegrator : public SamplingIntegrator {
public:
    /// Serialize this integrator to a binary data stream
    void serialize(Stream *stream, InstanceManager *manager) const;

    MTS_DECLARE_CLASS()
protected:
    /// Create a integrator
    MonteCarloIntegrator(const Properties &props);
    /// Unserialize an integrator
    MonteCarloIntegrator(Stream *stream, InstanceManager *manager);
    /// Virtual destructor
    virtual ~MonteCarloIntegrator() { }
protected:
    int m_maxDepth;
    int m_rrDepth;
    bool m_strictNormals;
    bool m_hideEmitters;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_INTEGRATOR_H_ */

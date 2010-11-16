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

#if !defined(__INTEGRATOR_H)
#define __INTEGRATOR_H

#include <mitsuba/core/netobject.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

/**
 * Abstract integrator base-class. Does not make any assumptions on
 * how radiance is computed. 
 * Amongst other things, this allows the use of hardware-accelerated 
 * rasterization, which directly operates on the camera's film and has 
 * no global knowledge about radiance within the scene. Other possibilities
 * are sampling- or particle tracing-based integrators.
 */
class MTS_EXPORT_RENDER Integrator : public NetworkedObject {
public:
	/**
	 * Possibly perform a pre-process task. The last three parameters are
	 * resource IDs of the associated scene, camera and sample generator,
	 * which have been made available to all local and remote workers.
	 * The default implementation simply returns.
	 */
	virtual bool preprocess(const Scene *scene, RenderQueue *queue, 
		const RenderJob *job, int sceneResID, int cameraResID, 
		int samplerResID);

	/**
	 * Render the scene as seen by the default camera. Progress is tracked
	 * by sending status messages to a provided render queue. The parameter
	 * <tt>job</tt> is required to discern multiple render jobs occurring in 
	 * parallel. The last three parameters are resource IDs of the associated 
	 * scene, camera and sample generator, which have been made available to 
	 * all local and remote workers. Returns true upon successful completion.
	 */
	virtual bool render(Scene *scene, RenderQueue *queue, const RenderJob *job, 
		int sceneResID, int cameraResID, int samplerResID) = 0;

	/**
	 * This can be called asynchronously to cancel a running render job.
	 * In this case, <tt>render()</tt> will quit with a return value of 
	 * <tt>false</tt>.
	 */
	virtual void cancel() = 0;

	/**
	 * Possibly perform a post-process task. The last three parameters are
	 * resource IDs of the associated scene, camera and sample generator,
	 * which have been made available to all local and remote workers.
	 * The default implementation simply returns.
	 */
	virtual void postprocess(const Scene *scene, RenderQueue *queue, 
		const RenderJob *job, int sceneResID, int cameraResID,
		int samplerResID);

	/**
	 * This is called once after instantiation and can be used to
	 * inform the sampler implementation about sample requirements
	 * of this integrator.
	 */
	virtual void configureSampler(Sampler *sampler);

	/**
	 * When the integrator contains a nested integrator, this function can
	 * be used to query for it
	 */
	virtual const Integrator *getSubIntegrator() const;

	/// Serialize this integrator to a binary data stream
	void serialize(Stream *stream, InstanceManager *manager) const;

	/// Return the properties of this integrator
	inline const Properties &getProperties() const { return m_properties; }

	MTS_DECLARE_CLASS()
protected:
	/// Create a integrator
	Integrator(const Properties &props);

	/// Unserialize an integrator
	Integrator(Stream *stream, InstanceManager *manager); 

	/// Virtual destructor
	virtual ~Integrator() { }
protected:
	Properties m_properties;
};

struct MTS_EXPORT_RENDER RadianceQueryRecord {
public:
	/**
	 * List of suported query types. These can be combined by a binary OR.
	 */
	enum ERadianceQuery {
		/* Emitted radiance from a luminaire intersected by the ray */
		EEmittedRadiance     = 0x0001,

		/* Emitted radiance from a subsurface integrator */
		ESubsurfaceRadiance  = 0x0002,

		/* Direct (surface) radiance */
		EDirectRadiance      = 0x0004,

		/* Indirect (surface) radiance, where the last bounce did not go
		   through a Dirac delta BSDF */
		EIndirectRadiance    = 0x0008,

		/* Indirect (surface) radiance, where the last bounce went
		   through a Dirac delta BSDF */
		ECausticRadiance     = 0x0010,

		/* In-scattered radiance due to volumetric scattering (direct) */
		EInscatteredDirectRadiance = 0x0020,

		/* In-scattered radiance due to volumetric scattering (indirect) */
		EInscatteredIndirectRadiance = 0x0040,

		/* Distance to the next surface intersection */
		EDistance            = 0x0080,

		/* Opacity value: 1 when a surface was hit, 0 when the ray leads
		   into empty space. When there is a participating medium,
		   this can also take on fractional values. */
		EOpacity             = 0x0100,

		/* A ray intersection may need to be performed. This can be set to 
		   zero if the caller has already provided the intersection */
		EIntersection        = 0x0200,

		/* Radiance from volumes */
		EVolumeRadiance      = EInscatteredDirectRadiance | EInscatteredIndirectRadiance,

		/* Radiance query without emitted radiance, ray intersection required */
		ERadianceNoEmission  = ESubsurfaceRadiance | EDirectRadiance | EIndirectRadiance
			| ECausticRadiance | EInscatteredDirectRadiance | EInscatteredIndirectRadiance | EIntersection,

		/* Default radiance query, ray intersection required */
		ERadiance = ERadianceNoEmission | EEmittedRadiance,

		/* Radiance + opacity */
		ECameraRay = ERadiance | EOpacity
	};

	/// Construct an invalid radiance query record
	inline RadianceQueryRecord() 
	 : type(0), scene(NULL), sampler(NULL),
	   depth(0), alpha(0), dist(-1), wavelength(-1), extra(0) {
	}

	/// Construct a radiance query record for the given scene and sampler
	inline RadianceQueryRecord(const Scene *scene, Sampler *sampler) 
	 : type(0), scene(scene), sampler(sampler), 
	   depth(0), alpha(0), dist(-1), wavelength(-1), extra(0) {
	}
	
	/// Copy constructor
	inline RadianceQueryRecord(const RadianceQueryRecord &rRec) 
	 : type(rRec.type), scene(rRec.scene), sampler(rRec.sampler), 
	   depth(rRec.depth), alpha(rRec.alpha), dist(rRec.dist),
	   wavelength(rRec.wavelength), extra(rRec.extra) {
	}

	/// Begin a new query of the given type
	inline void newQuery(int _type) {
		type = _type;
		depth = 1;
		wavelength = -1;
		extra = 0;
	}

	/// Initialize the query record for a recursive query
	inline void recursiveQuery(const RadianceQueryRecord &parent, int _type) {
		type = _type;
		scene = parent.scene;
		sampler = parent.sampler;
		depth = parent.depth+1;
		wavelength = parent.wavelength;
		extra = 0;
	}

	/**
	 * Search for a ray intersection. This
	 * does several things at once - if the intersection has 
	 * already been provided, the function returns.
	 * Otherwise, it
	 * - performs the ray intersection
	 * - computes the attenuation due to participating media
	 *   and stores it in <tt>attenuation</tt>.
	 * - sets the alpha value (if <tt>EAlpha</tt> is set in <tt>type</tt>)
	 * - sets the distance value (if <tt>EDistance</tt> is set in <tt>type</tt>)
	 * - clears the <tt>EIntersection</tt> flag in <tt>type</tt>
	 * Returns true if there is a valid intersection.
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

	/// Current depth value (# of light bounces) (*)
	int depth;

	/// Surface interaction data structure (*)
	Intersection its;

	/// Attenuation along the current ray (*)
	Spectrum attenuation;

	/// Opacity value of the associated pixel (*)
	Float alpha;

	/**
	 * Ray distance to the first surface interaction 
	 * (if requested by the query type EDistance) (*)
	 */
	Float dist;

	/**
	 * In some cases, the integrator may be forced to restrict
	 * radiance computations to one wavelength (e.g. when intersecting
	 * a dielectric material with dispersion or while path
	 * tracing through a highly scattering medium with a non-constant
	 * scattering coefficient). This attribute is used to store the
	 * chosen wavelength. (*)
	 */
	int wavelength;

	/**
	 * Internal flag, which can be used to pass additional information 
	 * amonst recursive calls inside an integrator. The use
	 * is dependent on the particular integrator implementation. (*)
	 */
	int extra;
};

/** \brief Abstract base class, which describes integrators
 * capable of computing samples of the scene's radiance function.
 */
class MTS_EXPORT_RENDER SampleIntegrator : public Integrator {
public:
	/**
	 * Sample the incident radiance along a ray. Also requires
	 * a radiance query record, which makes this request more precise.
	 */
	virtual Spectrum Li(const RayDifferential &ray,
		RadianceQueryRecord &rRec) const = 0;

	/**
	 * Estimate the irradiance at a given surface point. The
	 * default implementation simply samples the hemisphere using
	 * cosine-weighted sampling and a configurable number of rays.
	 */
	virtual Spectrum E(const Scene *scene, const Point &p, const
		Normal &n, Float time, Sampler *sampler) const; 

	/**
	 * Perform the main rendering task. The work is automatically
	 * parallelized to multiple cores and remote machines. The default
	 * implementation uniformly generates samples on the camera lens 
	 * and image plane as specified by the used sampler. The average 
	 * of the estimated radiance along the associated rays in a pixel
	 * region is then taken as an approximation of that pixel's 
	 * radiance value. For adaptive strategies, have a look at the 
	 * <tt>errctrl</tt> plugin, which is an extension of this class.
	 */
	bool render(Scene *scene, RenderQueue *queue, const RenderJob *job, 
		int sceneResID, int cameraResID, int samplerResID);

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
	 */
	virtual void renderBlock(const Scene *scene, const Camera *camera, 
		Sampler *sampler, ImageBlock *block, const bool &stop) const;

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
	virtual void wakeup(std::map<std::string, SerializableObject *> &params);

	/// Serialize this integrator to disk
	void serialize(Stream *stream, InstanceManager *manager) const;

	MTS_DECLARE_CLASS()
protected:
	/// Create a integrator
	SampleIntegrator(const Properties &props);

	/// Unserialize an integrator
	SampleIntegrator(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~SampleIntegrator() { }
protected:
	/// Used to temporarily cache a parallel process while it is in operation
	ref<ParallelProcess> m_process;
	unsigned int m_irrSamples;
	bool m_irrIndirect;
};

/*
 * \brief Base class of all recursive Monte Carlo integrators, which compute
 * unbiased solutions to the rendering equation (and optionally
 * the radiative transfer equation).
 */
class MTS_EXPORT_RENDER MonteCarloIntegrator : public SampleIntegrator {
public:
	/// Serialize this integrator to disk
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
};

MTS_NAMESPACE_END

#endif /* __INTEGRATOR_H */

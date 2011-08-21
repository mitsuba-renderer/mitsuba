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

#if !defined(__SCENE_H)
#define __SCENE_H

#include <mitsuba/core/netobject.h>
#include <mitsuba/core/pdf.h>
#include <mitsuba/core/aabb.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/render/skdtree.h>
#include <mitsuba/render/camera.h>
#include <mitsuba/render/luminaire.h>
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
 * This class holds information on surfaces, luminaires and participating media
 * and coordinates rendering jobs. It also provides useful query routines that 
 * are mostly used by the \ref Integrator implementations.
 *
 * \ingroup librender
 * \ingroup libpython
 */
class MTS_EXPORT_RENDER Scene : public NetworkedObject {
public:
	/**
	 * When this scene is used as a test case, the following enumeration
	 * gives more information on the performed type of test.
	 */
	enum ETestType {
		/* No test */
		ENone = 0,

		/* Perform a statistical test using Student's T statistic */
		ETTest,

		/* Compare the relative error against a specified upper bound */
		ERelativeError
	};

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
	 * \brief Initialize the scene. This function \a must be called 
	 * before using any of the methods in this class.
	 */
	void initialize();

	/**
	 * Pre-process step - should be called after initialize() and 
	 * before rendering the scene. This might do a variety of things, 
	 * such as constructing photon maps or executing distributed overture 
	 * passes. Progress is tracked by sending status messages to a provided 
	 * render queue. The parameter \c job is required to discern 
	 * multiple render jobs occurring in parallel. The last three parameters 
	 * are resource IDs of the associated scene, camera and sample generator,
	 * which have been made available to  all local and remote workers.
	 * Returns true upon successful completion.
	 */
	bool preprocess(RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID);

	/**
	 * Render the scene as seen by the scene's main camera. Progress is tracked
	 * by sending status messages to a provided render queue. The parameter
	 * \c job is required to discern multiple render jobs occurring in 
	 * parallel. The last three parameters are resource IDs of the associated 
	 * scene, camera and sample generator, which have been made available to 
	 * all local and remote workers. Returns true upon successful completion.
	 */
	bool render(RenderQueue *queue, const RenderJob *job,
			int sceneResID, int cameraResID, int samplerResID);

	/// Post-process step after rendering. Parameters are explained above
	void postprocess(RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID);

	/// Write out the current (partially rendered) image
	void flush();

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

	/** \brief Configure this object (called _once_ after construction
	   and addition of all child ConfigurableObjects.) */
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
	 *    extent information, as well as a time (which applies when
	 *    the shapes are animated)
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
	 * and return detailed intersection information
	 *
	 * This variant also computes the attenuation until the next
	 * surface intersection (or until infinity if there is no
	 * intersection). Index-matched medium transitions are transparently
	 * handled and do not result in an intersection.
	 *
	 * \param ray
	 *    A 3-dimensional ray data structure with minimum/maximum
	 *    extent information, as well as a time (which applies when
	 *    the shapes are animated)
	 *
	 * \param medium
	 *    Initial medium containing the ray \c ray
	 *
	 * \param its
	 *    A detailed intersection record, which will be filled by the
	 *    intersection query
	 *
	 * \param indexMatchedMediumTransition
	 *    This parameter is used to return whether or not this routine
	 *    passed through an index-matched medium-to-medium transition 
	 *    while computing the attenuation.
	 *
	 * \param transmittance
	 *    Transmittance from the ray origin to the returned intersection
	 *    (or infinity if there is no intersection)
	 *
	 * \return \c true if an intersection was found
	 */
	bool attenuatedRayIntersect(const Ray &ray, const Medium * medium,
		Intersection &its, bool &indexMatchedMediumTransition,
		Spectrum &transmittance, Sampler *sampler = NULL) const;

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
	 *    extent information, as well as a time (which applies when
	 *    the shapes are animated)
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
	 * \return \c true if an intersection was found
	 */
	inline bool rayIntersect(const Ray &ray, Float &t, 
			ConstShapePtr &shape, Normal &n) const {
		return m_kdtree->rayIntersect(ray, t, shape, n);
	}

	/**
	 * \brief Test for occlusion between \c p1 and \c p2 at the
	 * specified time
	 *
	 * \return \c true if an occluder is located on the line segment
	 * between \c p1 and \c p2.
	 */
	inline bool isOccluded(const Point &p1, const Point &p2, Float time) const {
		Ray ray(p1, p2-p1, time);
		ray.mint = ShadowEpsilon;
		ray.maxt = 1-ShadowEpsilon;
		return m_kdtree->rayIntersect(ray);
	}

	/**
	 * \brief Return the transmittance between \c p1 and \c p2 at
	 * the specified time.
	 *
	 * This function is essentially a continuous version of \ref isOccluded,
	 * which accounts for the presence of participating media. 
	 *
	 * The implementation correctly handles arbitrary amounts of index-matched
	 * medium transitions, which means that it won't just stop and return zero
	 * after encountering one. Mismatched boundaries need to be handled
	 * differently though.
	 *
	 * \return An spectral-valued transmittance value with components
	 * between zero and one.
	 */

	Spectrum getTransmittance(const Point &p1, const Point &p2, Float time, 
		const Medium *medium, Sampler *sampler = NULL) const;
	
	/// Return an axis-aligned bounding box containing the whole scene
	inline const AABB &getAABB() const {
		return m_aabb;
	}
	
	/// Return a bounding sphere containing the whole scene
	inline const BSphere &getBSphere() const {
		return m_bsphere;
	}
	
	//! @}
	// =============================================================
	
	// =============================================================
	//! @{ \name Luminaire query & sampling functions
	// =============================================================

	/**
	 * Sample a visible point on a luminaire (ideally uniform wrt. solid angle at \c p)
	 * \param p
	 *     An arbitrary point in the scene
	 * \param time
	 *    Associated time value -- this is needed to check the visibility when
	 *    objects are potentially moving over time
	 * \param lRec
	 *    A luminaire sampling record, which will hold information such as the
	 *    probability density, etc.
	 * \param testVisibility
	 *    If this is true, a shadow-ray will be cast to ensure that no surface
	 *    blocks the path lRec.sRec.p <-> p.
	 * \return
     *    \c true if sampling was successful
	 */
	bool sampleLuminaire(const Point &p, Float time, LuminaireSamplingRecord &lRec,
			const Point2 &sample, bool testVisibility = true) const;

	/**
	 * \brief Convenience method - similar to \ref sampleLuminaire(), but also attenuates
	 * \c lRec.value by the integrated extinction coefficient along the connection path.
	 * A visibility test is always performed. This function is meant to be used when 
	 * \c p does not lie on a surface.
	 *
	 * \param p
	 *     An arbitrary point in the scene
	 * \param time
	 *    Associated time value -- this is needed to check the visibility when
	 *    objects are potentially moving over time
	 * \param medium
	 *    The medium located at the ray origin (or \c NULL for vacuum).
	 * \param lRec
	 *    A luminaire sampling record, which will hold information such 
	 *    as the probability density, etc.
	 * \return
     *    \c true if sampling was successful
	 */
	bool sampleAttenuatedLuminaire(const Point &p, Float time, 
		const Medium *medium, LuminaireSamplingRecord &lRec, 
		const Point2 &sample, Sampler *sampler = NULL) const;

	/**
	 * \brief Convenience method - similar to \ref sampleLuminaire(), but also attenuates
	 * \c lRec.value by the integrated extinction coefficient along the connection path.
	 * A visibility test is always performed.
	 *
	 * \param its
	 *     An arbitrary intersection
	 * \param medium
	 *    The medium located at the ray origin (or \c NULL for vacuum).
	 * \param lRec
	 *    A luminaire sampling record, which will hold information such 
	 *    as the probability density, etc.
	 * \return
     *    \c true if sampling was successful
	 */
	bool sampleAttenuatedLuminaire(const Intersection &its, 
		const Medium *medium, LuminaireSamplingRecord &lRec, 
		const Point2 &sample, Sampler *sampler = NULL) const;

	/**
	 * \brief Return the probability density associated with the sample \c lRec 
	 * when calling \ref sampleLuminaire() with the point \c p.
	 *
	 * \param p
	 *     An arbitrary point in the scene
	 * \param lRec
	 *     Luminaire sampling record, whose associated density is to
	 *     be determined
	 * \param delta
	 *     When this parameter is set to true, only components with a
	 *     Dirac delta density are queried. Otherwise, they are left out.
	 */
	Float pdfLuminaire(const Point &p, const LuminaireSamplingRecord &lRec,
			bool delta = false) const;

	/**
	 * \brief Sample a particle leaving one of the scene's luminaires and 
	 * return a ray describing its path as well as record containing detailed
	 * probability density information. 
	 *
	 * Two uniformly distributed 2D samples are required. This method does
	 * exactly the same as calling \c sampleEmissionArea and 
	 * \c eRec.luminaire->sampleEmissionDirection(..) in sequence, 
	 * modulating \c eRec.Le by the return value of the latter and dividing
	 * by the product of the spatial and directional sampling densities.
	 */

	void sampleEmission(EmissionRecord &lRec, Point2 &sample1, Point2 &sample2) const;

	/**
	 * \brief Sample only the spatial part of the emission sampling strategy
	 * implemented in \c sampleEmission.
	 *
	 * An examplary use of this method is bidirectional path tracing or MLT, 
	 * where the area and direction sampling steps take place in different 
	 * vertices.
	 *
	 * After the function call terminates, the area density as well as the
	 * spatially dependent emittance component will be stored in \c eRec.
	 */
	void sampleEmissionArea(EmissionRecord &lRec, Point2 &sample) const;

	/**
	 * \brief Given an emitted particle, populate the emission record with the 
	 * relevant probability densities.
	 *
	 * When \c delta is set to true, only components with a Dirac delta density
	 * are considered in the query. Otherwise, they are left out.
	 */
	void pdfEmission(EmissionRecord &eRec, bool delta) const;

	/**
	 * \brief Return the background radiance for a ray that did not intersect
	 * any of the scene objects.
	 *
	 * This is primarily meant for path tracing-style integrators.
	 */
	inline Spectrum LeBackground(const Ray &ray) const {
		return hasBackgroundLuminaire() ? m_backgroundLuminaire->Le(ray) : Spectrum(0.0f);
	}

	/**
	 * \brief Return the background radiance for a ray that did not intersect
	 * any of the scene objects. This method additionally considers
	 * transmittance by participating media
	 *
	 * This is primarily meant for path tracing-style integrators.
	 */
	inline Spectrum LeAttenuatedBackground(const Ray &ray, const Medium *medium, Sampler *sampler) const {
		if (!m_backgroundLuminaire)
			return Spectrum(0.0f);
		Spectrum result = LeBackground(ray);
		if (medium)
			result *= medium->getTransmittance(ray, sampler);
		return result;
	}
	
	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Miscellaneous
	// =============================================================

	/// Does the scene contain participating media?
	inline bool hasMedia() const { return !m_media.empty(); }

	/// Return the scene's test mode
	inline ETestType getTestType() const { return m_testType; }
	/// Return the scene's test threshold
	inline Float getTestThreshold() const { return m_testThresh; }

	/**
	 * \brief Set the scene's camera. 
	 *
	 * Note that the camera is not included when this Scene instance
	 * is serialized -- the camera field will be \c NULL after 
	 * unserialization. This is intentional so that the camera can 
	 * be changed without having to re-transmit the whole scene. 
	 * Hence, the camera needs to be submitted separately
	 * and re-attached on the remote side using \ref setCamera().
	 **/
	inline void setCamera(Camera *camera) { m_camera = camera; }
	/// Return the scene's camera
	inline Camera *getCamera() { return m_camera; }
	/// Return the scene's camera (const version)
	inline const Camera *getCamera() const { return m_camera.get(); }
	
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
	 * clone of the original sampler specified in the camera.
	 */
	inline Sampler *getSampler() { return m_sampler; }
	/// Return the scene's sampler
	inline const Sampler *getSampler() const { return m_sampler.get(); }

	/// Return the scene's film
	inline Film *getFilm() { return m_camera->getFilm(); }
	/// Return the scene's film
	inline const Film *getFilm() const { return m_camera->getFilm(); }
	/// Set the scene's film
	inline void setFilm(Film *film) { m_camera->setFilm(film); }
	
	/// Return the scene's kd-tree accelerator
	inline ShapeKDTree *getKDTree() { return m_kdtree; }
	/// Return the scene's kd-tree accelerator
	inline const ShapeKDTree *getKDTree() const { return m_kdtree.get(); }

	/// Return the a list of all subsurface integrators
	inline const std::vector<Subsurface *> &getSubsurfaceIntegrators() const { return m_ssIntegrators; }

	/// Return the scene's background luminaire (if there is one)
	inline const Luminaire *getBackgroundLuminaire() const { return m_backgroundLuminaire.get(); }
	/// Does the scene have a background luminaire?
	inline bool hasBackgroundLuminaire() const { return m_backgroundLuminaire.get() != NULL; }

	/// Return the scene's triangular meshes
	inline std::vector<TriMesh *> &getMeshes() { return m_meshes; }
	/// Return the scene's triangular meshes
	inline const std::vector<TriMesh *> &getMeshes() const { return m_meshes; }
	/// Return the scene's shapes (including triangular meshes)
	inline std::vector<Shape *> &getShapes() { return m_shapes; }
	/// Return the scene's shapes (including triangular meshes)
	inline const std::vector<Shape *> &getShapes() const { return m_shapes; }
	/// Return the scene's luminaires
	inline std::vector<Luminaire *> &getLuminaires() { return m_luminaires; }
	/// Return the scene's luminaires
	inline const std::vector<Luminaire *> &getLuminaires() const { return m_luminaires; }
	/// Return the scene's participating media
	inline std::set<Medium *> &getMedia() { return m_media; }
	/// Return the scene's participating media
	inline const std::set<Medium *> &getMedia() const { return m_media; }
	/// Return referenced objects (such as textures, BSDFs)
	inline std::vector<ConfigurableObject *> &getReferencedObjects() { return m_objects; }
	/// Return referenced objects (such as textures, BSDFs)
	inline const std::vector<ConfigurableObject *> &getReferencedObjects() const { return m_objects; }

	/// Return the name of the file containing the original description of this scene
	inline const fs::path getSourceFile() const { return m_sourceFile; }
	/// Set the name of the file containing the original description of this scene
	void setSourceFile(const fs::path &name) { m_sourceFile = name; }
	/// Return the render output filename
	inline const fs::path getDestinationFile() const { return m_destinationFile; }
	/// Set the render output filename
	void setDestinationFile(const fs::path &name) { m_destinationFile = name; }
	/// Does the destination file already exist?
	inline bool destinationExists() const { return m_camera->getFilm()->destinationExists(m_destinationFile); }

	/// Set the block resolution used to split images into parallel workloads
	inline void setBlockSize(int size) { m_blockSize = size; }
	/// Return the block resolution used to split images into parallel workloads
	inline int getBlockSize() const { return m_blockSize; }

	/// Serialize the whole scene to a network/file stream
	void serialize(Stream *stream, InstanceManager *manager) const;
	/* NetworkedObject implementation */
	virtual void bindUsedResources(ParallelProcess *proc) const;
	virtual void wakeup(std::map<std::string, SerializableObject *> &params);

	/// Return a string representation
	std::string toString() const;

	//! @}
	// =============================================================

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Scene();

	/// Add a shape to the scene
	void addShape(Shape *shape);
private:
	ref<ShapeKDTree> m_kdtree;
	ref<Camera> m_camera;
	ref<Integrator> m_integrator;
	ref<Sampler> m_sampler;
	ref<Luminaire> m_backgroundLuminaire;
	std::vector<TriMesh *> m_meshes;
	std::vector<Shape *> m_shapes;
	std::vector<Luminaire *> m_luminaires;
	std::vector<Subsurface *> m_ssIntegrators;
	std::vector<ConfigurableObject *> m_objects;
	std::vector<NetworkedObject *> m_netObjects;
	std::set<Medium *> m_media;
	fs::path m_sourceFile;
	fs::path m_destinationFile;
	DiscretePDF m_luminairePDF;
	AABB m_aabb;
	BSphere m_bsphere;
	bool m_importanceSampleLuminaires;
	ETestType m_testType;
	Float m_testThresh;
	int m_blockSize;
};

MTS_NAMESPACE_END

#include <mitsuba/render/records.inl>

#endif /* __SCENE_H */

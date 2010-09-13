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

#if !defined(__SCENE_H)
#define __SCENE_H

#include <mitsuba/core/netobject.h>
#include <mitsuba/core/pdf.h>
#include <mitsuba/core/aabb.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/render/kdtree.h>
#include <mitsuba/render/camera.h>
#include <mitsuba/render/luminaire.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/medium.h>

MTS_NAMESPACE_BEGIN

/**
 * Scene data structure: holds information on surfaces, luminaires
 * and participating media and coordinates rendering jobs. Also provides
 * useful query routines mostly used by the integrator implementations.
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

	/* ==================================================================== */
    /*                     Initialization & Rendering                       */
    /* ==================================================================== */

	/// Construct a new, empty scene
	Scene(const Properties &props);

	/// Create a shallow clone of a scene
	Scene(Scene *scene);

	/// Unserialize a scene from a binary data stream
	Scene(Stream *stream, InstanceManager *manager);

	/**
	 * Initialization step - _must_ be called before using any of the
	 * methods provided by <tt>Scene</tt>.
	 */
	void initialize();

	/**
	 * Pre-process step - should be called after initialize() and 
	 * before rendering the scene. This might do a variety of things, 
	 * such as constructing photon maps or executing distributed overture 
	 * passes. Progress is tracked by sending status messages to a provided 
	 * render queue. The parameter <tt>job</tt> is required to discern 
	 * multiple render jobs occurring in parallel. The last three parameters 
	 * are resource IDs of the associated scene, camera and sample generator,
	 * which have been made available to  all local and remote workers.
	 * Returns true upon successful completion.
	 */
	void preprocess(RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID);

	/**
	 * Render the scene as seen by the scene's main camera. Progress is tracked
	 * by sending status messages to a provided render queue. The parameter
	 * <tt>job</tt> is required to discern multiple render jobs occurring in 
	 * parallel. The last three parameters are resource IDs of the associated 
	 * scene, camera and sample generator, which have been made available to 
	 * all local and remote workers. Returns true upon successful completion.
	 */
	bool render(RenderQueue *queue, const RenderJob *job,
			int sceneResID, int cameraResID, int samplerResID);
	
	/// Write out the current (partially rendered) image
	void flush();

	/**
	 * This can be called asynchronously to cancel a running render job.
	 * In this case, <tt>render()</tt> will quit with a return value of 
	 * <tt>false</tt>.
	 */
	void cancel();

	/// Post-process step after rendering. Parameters are explained above
	void postprocess(RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID);

	/// Add a child node to the scene
	void addChild(const std::string &name, ConfigurableObject *child);

	/** \brief Configure this object (called _once_ after construction
	   and addition of all child ConfigurableObjects.) */
	void configure();

	/* ==================================================================== */
    /*                         Integrator interface                         */
    /* ==================================================================== */

	/// Check for an intersection with scene objects
	inline bool rayIntersect(const Ray &pRay, Intersection &its) const {
		return m_kdtree->rayIntersect(pRay, its);
	}

	/// Cast a shadow ray
	inline bool isOccluded(const Point &p1, const Point &p2) const {
		Ray ray(p1, p2-p1);
		ray.mint = ShadowEpsilon;
		ray.maxt = 1-ShadowEpsilon;
		return m_kdtree->rayIntersect(ray);
	}

	/// Return an axis-aligned bounding box containing the whole scene
	inline const AABB &getAABB() const {
		return m_aabb;
	}
	
	/// Return a bounding sphere containing the whole scene
	inline const BSphere &getBSphere() const {
		return m_bsphere;
	}

	/**
	 * Sample a visible point on a luminaire (ideally uniform wrt. the solid angle of p)
	 * @param p
	 *     An arbitrary point in the scene
	 * @param lRec
	 *    A luminaire sampling record, which will hold information such as the
	 *    probability density, associated measure etc.
	 * @param testVisibility
	 *    If this is true, a shadow-ray will be cast to ensure that no surface
	 *    blocks the path lRec.sRec.p <-> p.
	 * @return
     *    true if sampling was successful
	 */
	bool sampleLuminaire(const Point &p,
		LuminaireSamplingRecord &lRec, const Point2 &sample, bool testVisibility = true) const;

	/**
	 * Sample a visible point on a luminaire (ideally uniform wrt. the solid angle of p). Takes
	 * a surface intersection record (some luminaires make use of this to provide improved sampling)
	 * @param its 
	 *     An arbitrary surface intersection
	 * @param lRec
	 *    A luminaire sampling record, which will hold information such as the
	 *    probability density, associated measure etc.
	 * @param testVisibility
	 *    If this is true, a shadow-ray will be cast to ensure that no surface
	 *    blocks the path lRec.sRec.p <-> p.
	 * @return
     *    true if sampling was successful
	 */
	bool sampleLuminaire(const Intersection &its,
		LuminaireSamplingRecord &lRec, const Point2 &sample, bool testVisibility = true) const;

	/**
	 * Convenience method - similar to sampleLuminaire(), but also attenuates
	 * lRec.Le by the integrated extinction coefficient on the path lRec.sRec.p <-> p.
	 */
	inline bool sampleLuminaireAttenuated(const Point &p,
		LuminaireSamplingRecord &lRec, const Point2 &sample, bool testVisibility = true) const {
		if (sampleLuminaire(p, lRec, sample, testVisibility)) {
			lRec.Le *= getAttenuation(Ray(p, lRec.sRec.p-p, 0, 1));
			return true;
		}
		return false;
	}

	/**
	 * Convenience method - similar to sampleLuminaire(), but also attenuates
	 * lRec.Le by the integrated extinction coefficient on the path lRec.sRec.p <-> p.
	 */
	inline bool sampleLuminaireAttenuated(const Intersection &its,
		LuminaireSamplingRecord &lRec, const Point2 &sample, bool testVisibility = true) const {
		if (sampleLuminaire(its, lRec, sample, testVisibility)) {
			lRec.Le *= getAttenuation(Ray(its.p, lRec.sRec.p-its.p, 0, 1));
			return true;
		}
		return false;
	}

	/**
	 * Return the probability of generating a sample on a luminaire if
	 * sampleLuminaire() is called with the point 'p'.
	 */
	Float pdfLuminaire(const Point &p,
		const LuminaireSamplingRecord &lRec) const;

	/**
	 * Return the probability of generating a sample on a luminaire if
	 * sampleLuminaire() is called with the surface interation 'its'
	 */
	Float pdfLuminaire(const Intersection &its,
		const LuminaireSamplingRecord &lRec) const;

	/**
	 * Sample a particle leaving a luminaire and return a ray describing its path
	 * as well as record containing detailed probability density information. 
	 * Two uniformly distributed 2D samples are required.
	 */
	void sampleEmission(EmissionRecord &lRec, Point2 &sample1, Point2 &sample2) const;

	/**
	 * Sample only the spatial dimension of the emission sampling strategy
	 * implemented in <tt>sampleEmission</tt>. An examplary use of this
	 * is bidirectional path tracing, where the directionally varying part
	 * is handed similarly to a BSDF, which modulates the spatially
	 * dependent radiance component. After the function call terminates,
	 * the area density as well as the spatially dependent radiance component
	 * will be stored in <tt>eRec</tt>.
	 */
	void sampleEmissionArea(EmissionRecord &lRec, Point2 &sample) const;

	/**
	 * As above, but handles only the directional part. Must be called *after*
	 * sampleEmissionArea(). The return value is to be understood as a BRDF,
	 * which modulates the emitted energy.
	 */
	Spectrum sampleEmissionDirection(EmissionRecord &lRec, Point2 &sample) const;

	/**
	 * Given an emitted particle, populate the emission record with the relevant spectra
	 * and probability densities.
	 */
	void pdfEmission(EmissionRecord &eRec) const;

	/**
	 * Return the background radiance for a ray, which did not hit anything.
	 */
	Spectrum LeBackground(const Ray &ray) const;

	/**
	 * Return the background radiance for a ray, which did not hit anything.
	 * (attenuated by any participating media along the ray)
	 */
	inline Spectrum LeBackgroundAttenuated(const Ray &ray) const {
		if (!m_backgroundLuminaire)
			return Spectrum(0.0f);
		return LeBackground(ray) * getAttenuation(ray);
	}

	/// Calculate the attenuation along the ray segment [mint, maxt]
	Spectrum getAttenuation(const Ray &ray) const;

	/// Does the scene contain participating media?
	inline bool hasMedia() const { return !m_media.empty(); }

	/**
	 * In the presence of participating media, sample a traveled distance
	 * along a ray. Should ideally importance sample with respect to the
	 * sample contribution.
	 * Returns false if the ray passes through all media (or
	 * if there were none). Assumes that none of the media overlap!
	 */
	bool sampleDistance(const Ray &ray, Float maxDist, MediumSamplingRecord &mRec,
		Sampler *sampler) const;

	/// Return the scene's test mode
	inline ETestType getTestType() const { return m_testType; }
	/// Return the scene's test threshold
	inline Float getTestThreshold() const { return m_testThresh; }

	/**
	 * Set the scene's camera. Note that the camera is not included
	 * when this Scene instance is serialized -- the camera field
	 * will be <tt>NULL</tt> after unserialization. This is intentional
	 * so that the camera can be changed without having to re-transmit
	 * the whole scene. Hence, the camera needs to be submitted separately
	 * and re-attached on the remote side using <tt>setCamera</tt>.
	 **/
	inline void setCamera(Camera *camera) { m_camera = camera; }
	/// Return the scene's camera
	inline Camera *getCamera() { return m_camera; }
	/// Return the scene's camera (const version)
	inline const Camera *getCamera() const { return m_camera.get(); }
	
	/**
	 * Set the scene's integrator. Note that the integrator is not included
	 * when this Scene instance is serialized -- the integrator field
	 * will be <tt>NULL</tt> after unserialization. This is intentional
	 * so that the integrator can be changed without having to re-transmit
	 * the whole scene. Hence, the integrator needs to be submitted separately
	 * and re-attached on the remote side using <tt>setIntegrator</tt>.
	 **/
	inline void setIntegrator(Integrator *integrator) { m_integrator = integrator; }
	/// Return the scene's integrator
	inline Integrator *getIntegrator() { return m_integrator; }
	/// Return the scene's integrator (const version)
	inline const Integrator *getIntegrator() const { return m_integrator.get(); }

	/**
	 * Set the scene's sampler. Note that the sampler is not included
	 * when this Scene instance is serialized -- the sampler field
	 * will be <tt>NULL</tt> after unserialization. This is intentional
	 * so that the sampler can be changed without having to re-transmit
	 * the whole scene. Hence, the sampler needs to be submitted separately
	 * and re-attached on the remote side using <tt>setSampler</tt>.
	 **/
	inline void setSampler(Sampler *sampler) { m_sampler = sampler; }
	/// Return the scene's sampler
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
	inline KDTree *getKDTree() { return m_kdtree; }
	/// Return the scene's kd-tree accelerator
	inline const KDTree *getKDTree() const { return m_kdtree.get(); }

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
	inline std::vector<Medium *> &getMedia() { return m_media; }
	/// Return the scene's participating media
	inline const std::vector<Medium *> &getMedia() const { return m_media; }
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

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Scene();

	/// Add a shape to the scene
	void addShape(Shape *shape);
private:
	ref<KDTree> m_kdtree;
	ref<Camera> m_camera;
	ref<Integrator> m_integrator;
	ref<Sampler> m_sampler;
	ref<Luminaire> m_backgroundLuminaire;
	std::vector<TriMesh *> m_meshes;
	std::vector<Shape *> m_shapes;
	std::vector<Luminaire *> m_luminaires;
	std::vector<Medium *> m_media;
	std::vector<Subsurface *> m_ssIntegrators;
	std::vector<ConfigurableObject *> m_objects;
	std::vector<NetworkedObject *> m_netObjects;
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

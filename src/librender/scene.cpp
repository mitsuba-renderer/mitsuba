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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/core/plugin.h>

#define DEFAULT_BLOCKSIZE 32

MTS_NAMESPACE_BEGIN

Scene::Scene() 
 : NetworkedObject(Properties()), m_blockSize(DEFAULT_BLOCKSIZE) {
	m_kdtree = new ShapeKDTree();
	m_testType = ENone;
	m_testThresh = 0.0f;
	m_importanceSampleLuminaires = true;
}

Scene::Scene(const Properties &props)
 : NetworkedObject(props), m_blockSize(DEFAULT_BLOCKSIZE) {
	m_kdtree = new ShapeKDTree();
	/* When test case mode is active (Mitsuba is started with the -t parameter), 
	  this specifies the type of test performed. Mitsuba will expect a reference 
	  solution file of the name <tt>&lt;sceneName&gt;.ref</tt>. When set to 
	  <tt>t-test</tt>, a two-sided t-test on equality to the reference will be 
	  performed at the (1 - <tt>testThresh</tt>) level (99% by default). 
	  When set to <tt>relerr</tt>, the test will fail if the relative error 
	  exceeds <tt>testThresh</tt> (default: 0.01%).
    */
	std::string testType = props.getString("testType", "none");
	if (testType == "none")
		m_testType = ENone;
	else if (testType == "t-test")
		m_testType = ETTest;
	else if (testType == "relerr")
		m_testType = ERelativeError;
	else
		Log(EError, "Unknown test mode \"%s\" specified (must be \"t-test\" or \"relerr\")",
			testType.c_str());
	/* Error threshold for use with <tt>testType</tt> */
	m_testThresh = props.getFloat("testThresh", 0.01f);
	/* By default, luminaire sampling chooses a luminaire with a probability 
	  dependent on the emitted power. Setting this parameter to false switches 
	  to uniform sampling. */
	m_importanceSampleLuminaires = props.getBoolean("importanceSampleLuminaires", true);
	/* kd-tree construction: Enable primitive clipping? Generally leads to a 
	  significant improvement of the resulting tree. */
	if (props.hasProperty("kdClip"))
		m_kdtree->setClip(props.getBoolean("kdClip"));
	/* kd-tree construction: Relative cost of a triangle intersection operation 
	   in the surface area heuristic. */
	if (props.hasProperty("kdIntersectionCost"))
		m_kdtree->setQueryCost(props.getFloat("kdIntersectionCost"));
	/* kd-tree construction: Relative cost of a kd-tree traversal operation 
	   in the surface area heuristic. */
	if (props.hasProperty("kdTraversalCost"))
		m_kdtree->setTraversalCost(props.getFloat("kdTraversalCost"));
	/* kd-tree construction: Bonus factor for cutting away regions of empty space */
	if (props.hasProperty("kdEmptySpaceBonus"))
		m_kdtree->setEmptySpaceBonus(props.getFloat("kdEmptySpaceBonus"));
	/* kd-tree construction: A kd-tree node containing this many or fewer 
	   primitives will not be split */
	if (props.hasProperty("kdStopPrims"))
		m_kdtree->setStopPrims(props.getInteger("kdStopPrims"));
	/* kd-tree construction: Maximum tree depth */
	if (props.hasProperty("kdMaxDepth"))
		m_kdtree->setMaxDepth(props.getInteger("kdMaxDepth"));
	/* kd-tree construction: Specify the number of primitives, at which the 
	   builder will switch from (approximate) Min-Max binning to the accurate 
	   O(n log n) SAH-based optimization method. */
	if (props.hasProperty("kdExactPrimitiveThreshold"))
		m_kdtree->setExactPrimitiveThreshold(props.getInteger("kdExactPrimitiveThreshold"));
	/* kd-tree construction: use multiple processors? */
	if (props.hasProperty("kdParallelBuild"))
		m_kdtree->setParallelBuild(props.getBoolean("kdParallelBuild"));
	/* kd-tree construction: specify whether or not bad splits can be "retracted". */
	if (props.hasProperty("kdRetract"))
		m_kdtree->setRetract(props.getBoolean("kdRetract"));
	/* kd-tree construction: Set the number of bad refines allowed to happen
	   in succession before a leaf node will be created.*/
	if (props.hasProperty("kdMaxBadRefines"))
		m_kdtree->setMaxBadRefines(props.getInteger("kdMaxBadRefines"));
}

Scene::Scene(Scene *scene) : NetworkedObject(Properties()) {
	m_kdtree = scene->m_kdtree;
	m_testType = scene->m_testType;
	m_testThresh = scene->m_testThresh;
	m_blockSize = scene->m_blockSize;
	m_aabb = scene->m_aabb;
	m_bsphere = scene->m_bsphere;
	m_backgroundLuminaire = scene->m_backgroundLuminaire;
	m_camera = scene->m_camera;
	m_integrator = scene->m_integrator;
	m_sourceFile = scene->m_sourceFile;
	m_destinationFile = scene->m_destinationFile;
	m_luminairePDF = scene->m_luminairePDF;
	m_importanceSampleLuminaires = scene->m_importanceSampleLuminaires;
	m_shapes = scene->m_shapes;
	for (size_t i=0; i<m_shapes.size(); ++i)
		m_shapes[i]->incRef();
	m_meshes = scene->m_meshes;
	for (size_t i=0; i<m_meshes.size(); ++i)
		m_meshes[i]->incRef();
	m_luminaires = scene->m_luminaires;
	for (size_t i=0; i<m_luminaires.size(); ++i)
		m_luminaires[i]->incRef();
	m_media = scene->m_media;
	for (std::set<Medium *>::iterator it = m_media.begin();
			it != m_media.end(); ++it)
		(*it)->incRef();
	m_ssIntegrators = scene->m_ssIntegrators;
	for (size_t i=0; i<m_ssIntegrators.size(); ++i)
		m_ssIntegrators[i]->incRef();
	m_objects = scene->m_objects;
	for (size_t i=0; i<m_objects.size(); ++i)
		m_objects[i]->incRef();
	m_netObjects = scene->m_netObjects;
	for (size_t i=0; i<m_netObjects.size(); ++i)
		m_netObjects[i]->incRef();
}


Scene::Scene(Stream *stream, InstanceManager *manager) 
 : NetworkedObject(stream, manager) {
	m_kdtree = new ShapeKDTree();
	m_kdtree->setQueryCost(stream->readFloat());
	m_kdtree->setTraversalCost(stream->readFloat());
	m_kdtree->setEmptySpaceBonus(stream->readFloat());
	m_kdtree->setStopPrims(stream->readInt());
	m_kdtree->setClip(stream->readBool());
	m_kdtree->setMaxDepth(stream->readUInt());
	m_kdtree->setExactPrimitiveThreshold(stream->readUInt());
	m_kdtree->setParallelBuild(stream->readBool());
	m_kdtree->setRetract(stream->readBool());
	m_kdtree->setMaxBadRefines(stream->readUInt());
	m_importanceSampleLuminaires = stream->readBool();
	m_testType = (ETestType) stream->readInt();
	m_testThresh = stream->readFloat();
	m_blockSize = stream->readInt();
	m_aabb = AABB(stream);
	m_bsphere = BSphere(stream);
	m_backgroundLuminaire = static_cast<Luminaire *>(manager->getInstance(stream));
	size_t count = stream->readSize();
	for (size_t i=0; i<count; ++i) {
		Shape *shape = static_cast<Shape *>(manager->getInstance(stream));
		shape->incRef();
		m_shapes.push_back(shape);
	}
	count = stream->readSize();
	for (size_t i=0; i<count; ++i) {
		TriMesh *trimesh = static_cast<TriMesh *>(manager->getInstance(stream));
		trimesh->incRef();
		m_meshes.push_back(trimesh);
	}
	count = stream->readSize();
	for (size_t i=0; i<count; ++i) {
		Luminaire *luminaire = static_cast<Luminaire *>(manager->getInstance(stream));
		luminaire->incRef();
		m_luminaires.push_back(luminaire);
	}
	count = stream->readSize();
	for (size_t i=0; i<count; ++i) {
		Medium *medium = static_cast<Medium *>(manager->getInstance(stream));
		medium->incRef();
		m_media.insert(medium);
	}
	count = stream->readSize();
	for (size_t i=0; i<count; ++i) {
		Subsurface *ssIntegrator = static_cast<Subsurface *>(manager->getInstance(stream));
		ssIntegrator->incRef();
		m_ssIntegrators.push_back(ssIntegrator);
	}
	count = stream->readSize();
	for (size_t i=0; i<count; ++i) {
		ConfigurableObject *obj = static_cast<ConfigurableObject *>(manager->getInstance(stream));
		obj->incRef();
		m_objects.push_back(obj);
	}
	count = stream->readSize();
	for (size_t i=0; i<count; ++i) {
		NetworkedObject *obj = static_cast<NetworkedObject *>(manager->getInstance(stream));
		m_netObjects.push_back(obj); // Do not increase the ref. count
	}
	initialize();
}

Scene::~Scene() {
	for (size_t i=0; i<m_shapes.size(); i++)
		m_shapes[i]->decRef();
	for (size_t i=0; i<m_meshes.size(); i++)
		m_meshes[i]->decRef();
	for (size_t i=0; i<m_objects.size(); i++)
		m_objects[i]->decRef();
	for (size_t i=0; i<m_ssIntegrators.size(); i++)
		m_ssIntegrators[i]->decRef();
	for (size_t i=0; i<m_luminaires.size(); i++)
		m_luminaires[i]->decRef();
	for (std::set<Medium *>::iterator it = m_media.begin();
			it != m_media.end(); ++it)
		(*it)->decRef();
}

void Scene::bindUsedResources(ParallelProcess *proc) const {
	for (size_t i=0; i<m_netObjects.size(); ++i)
		m_netObjects[i]->bindUsedResources(proc);
}

void Scene::wakeup(std::map<std::string, SerializableObject *> &params) {
	for (size_t i=0; i<m_netObjects.size(); ++i)
		m_netObjects[i]->wakeup(params);
}

void Scene::configure() {
	if (m_integrator == NULL) {
		/* Create a direct integrator by default */
		m_integrator = static_cast<Integrator *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Integrator), Properties("direct")));
		m_integrator->configure();
	}
	if (m_camera == NULL) {
		Properties props("perspective");
		/* Create a perspective camera with 45deg. FOV, which can see the whole scene */
		AABB aabb;
		for (size_t i=0; i<m_shapes.size(); ++i)
			aabb.expandBy(m_shapes[i]->getAABB());
		if (aabb.isValid()) {
			Log(EInfo, "No camera found! Adding a default camera.");
			Point center = aabb.getCenter();
			Vector extents = aabb.getExtents();
			Float maxExtents = std::max(extents.x, extents.y);
			Float distance = maxExtents/(2.0f * std::tan(45 * .5f * M_PI/180));

			props.setTransform("toWorld", Transform::translate(Vector(center.x, 
					center.y, aabb.min.z - distance)));
			props.setFloat("fov", 45.0f);

			m_camera = static_cast<Camera *> (PluginManager::getInstance()->
					createObject(MTS_CLASS(Camera), props));
			m_camera->configure();
			m_sampler = m_camera->getSampler();
		} else {
			m_camera = static_cast<Camera *> (PluginManager::getInstance()->
					createObject(MTS_CLASS(Camera), props));
			m_camera->configure();
			m_sampler = m_camera->getSampler();
		}
	}

	m_integrator->configureSampler(m_sampler);

	/**
	 * make it possible to serialize the integrator by 
	 * itself without all the extra baggage
	 */
	m_integrator->setParent(NULL); 
}

void Scene::initialize() {
	if (!m_kdtree->isBuilt()) {
		/* Expand all geometry */
		std::vector<Shape *> tempShapes;
		tempShapes.reserve(m_shapes.size());
		m_shapes.swap(tempShapes);
		for (size_t i=0; i<tempShapes.size(); ++i) {
			addShape(tempShapes[i]);
			tempShapes[i]->decRef();
		}

		/* Build the kd-tree */
		m_kdtree->build();

		m_aabb = m_kdtree->getAABB();
		m_bsphere = m_kdtree->getBSphere();
	}

	if (!m_luminairePDF.isReady()) {
		if (m_luminaires.size() == 0) {
			Log(EWarn, "No luminaires found -- adding a sky luminaire");
			Properties skyProps("sky");
			skyProps.setFloat("skyScale", 0.1f);
			//skyProps.setFloat("sunScale", 0.1f);
			skyProps.setBoolean("extend", true);
			ref<Luminaire> luminaire = static_cast<Luminaire *>(
				PluginManager::getInstance()->createObject(MTS_CLASS(Luminaire), skyProps));
			addChild(luminaire);
			luminaire->configure();
		}

		/* Calculate a discrete PDF to importance sample luminaires */
		for (std::vector<Luminaire *>::iterator it = m_luminaires.begin();
			it != m_luminaires.end(); ++it) {
			(*it)->preprocess(this);
			/* Add with a probability proportional to the luminaire's power */
			if (m_importanceSampleLuminaires)
				m_luminairePDF.put((*it)->getSamplingWeight());
			else
				m_luminairePDF.put(1.0f);
		}
		m_luminairePDF.build();
	} else {
		for (std::vector<Luminaire *>::iterator it = m_luminaires.begin();
			it != m_luminaires.end(); ++it) 
			(*it)->preprocess(this);
	}
}

bool Scene::preprocess(RenderQueue *queue, const RenderJob *job, 
		int sceneResID, int cameraResID, int samplerResID) {
	initialize();

	/* Pre-process step for the main scene integrator */
	if (!m_integrator->preprocess(this, queue, job,
		sceneResID, cameraResID, samplerResID))
		return false;

	/* Pre-process step for all sub-surface integrators */
	for (std::vector<Subsurface *>::iterator it = m_ssIntegrators.begin();
		it != m_ssIntegrators.end(); ++it) 
		if (!(*it)->preprocess(this, queue, job, 
			sceneResID, cameraResID, samplerResID))
			return false;
	return true;
}

bool Scene::render(RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID) {
	m_camera->getFilm()->clear();
	return m_integrator->render(this, queue, job, sceneResID, 
		cameraResID, samplerResID);
}

void Scene::cancel() {
	for (std::vector<Subsurface *>::iterator it = m_ssIntegrators.begin();
		it != m_ssIntegrators.end(); ++it) 
		(*it)->cancel();
	m_integrator->cancel();
}

void Scene::flush() {
	m_camera->getFilm()->develop(m_destinationFile);
}

void Scene::postprocess(RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID) {
	m_integrator->postprocess(this, queue, job, sceneResID, 
		cameraResID, samplerResID);
	m_camera->getFilm()->develop(m_destinationFile);
}

Float Scene::pdfLuminaire(const Point &p,
		const LuminaireSamplingRecord &lRec, bool delta) const {
	const Luminaire *luminaire = lRec.luminaire;
	Float luminance;

	if (m_importanceSampleLuminaires)
		luminance = luminaire->getSamplingWeight();
	else 
		luminance = 1.0f;

	/* Calculate the probability of importance sampling this luminaire */
	const Float fraction = luminance / m_luminairePDF.getOriginalSum();
	return luminaire->pdf(p, lRec, delta) * fraction;
}

bool Scene::sampleLuminaire(const Point &p, Float time,
		LuminaireSamplingRecord &lRec, const Point2 &s,
		bool testVisibility) const {
	Point2 sample(s);
	Float lumPdf;
	size_t index = m_luminairePDF.sampleReuse(sample.x, lumPdf);
	const Luminaire *luminaire = m_luminaires[index];
	luminaire->sample(p, lRec, sample);

	if (lRec.pdf != 0) {
		if (testVisibility) {
			Vector dir = lRec.sRec.p - p;
			Float length = dir.length();
			Ray ray(p, dir/length, Epsilon, length*(1-ShadowEpsilon), time);
			if (m_kdtree->rayIntersect(ray))
				return false;
		}
		lRec.pdf *= lumPdf;
		lRec.value /= lRec.pdf;
		lRec.luminaire = luminaire;
		return true;
	} else {
		return false;
	}
}

Spectrum Scene::getTransmittance(const Point &p1, const Point &p2,
		Float time, const Medium *medium, Sampler *sampler) const {
	if (m_media.size() == 0) {
		Vector dir = p2-p1;
		Float length = dir.length();
		Ray ray(p1, dir/length, Epsilon, length*(1-ShadowEpsilon), time);
		return Spectrum(m_kdtree->rayIntersect(ray) ? 0.0f : 1.0f);
	} else {
		Vector d = p2 - p1;
		Float remaining = d.length();
		d /= remaining;
		remaining *= 1-Epsilon;

		const Shape *shape;
		Ray ray(p1, d, time);
		Spectrum transmittance(1.0f);
		int iterations = 0;

		while (remaining > 0) {
			Normal n;
			Float t;

			bool surface = rayIntersect(ray, t, shape, n);

			if (medium)
				transmittance *= medium->getTransmittance(
					Ray(ray, 0, std::min(t, remaining)), sampler);

			if (!surface) 
				break;

			ray.o = ray(t);
			remaining -= t;

			if (remaining > 0) {
				if (shape->isOccluder())
					return Spectrum(0.0f);
				else if (shape->isMediumTransition())
					medium = dot(n, d) > 0 ? shape->getExteriorMedium()
						: shape->getInteriorMedium();
				if (++iterations > 100) { /// Just a precaution..
					Log(EWarn, "getTransmittance(): round-off error issues?");
					break;
				}
			}
		}

		return transmittance;
	}
}

bool Scene::attenuatedRayIntersect(const Ray &_ray, const Medium *medium,
		Intersection &its, bool &indexMatchedMediumTransition,
		Spectrum &transmittance, Sampler *sampler) const {
	Ray ray(_ray);
	transmittance = Spectrum(1.0f);
	int iterations = 0;

	while (true) {
		bool surface = m_kdtree->rayIntersect(ray, its);

		if (medium) 
			transmittance *= medium->getTransmittance(Ray(ray, 0, its.t), sampler);

		if (!surface)
			return false;
		else if (its.shape->isOccluder())
			return true;
		else if (its.shape->isMediumTransition()) {
			medium = dot(its.geoFrame.n, ray.d) > 0 ?
				  its.shape->getExteriorMedium()
				: its.shape->getInteriorMedium();
			indexMatchedMediumTransition = true;
		}

		ray.o = ray(its.t);
		ray.mint = Epsilon;

		if (++iterations > 100) { /// Just a precaution..
			Log(EWarn, "attenuatedRayIntersect(): round-off error issues?");
			return false;
		}
	}
}

bool Scene::sampleAttenuatedLuminaire(const Point &p, Float time, 
	const Medium *medium, LuminaireSamplingRecord &lRec, 
	const Point2 &s, Sampler *sampler) const {
	Point2 sample(s);
	Float lumPdf;
	const Luminaire *luminaire = m_luminaires[
		m_luminairePDF.sampleReuse(sample.x, lumPdf)];
	luminaire->sample(p, lRec, sample);

	if (lRec.pdf != 0) {
		lRec.value *= getTransmittance(p, lRec.sRec.p, time, medium, sampler);
		if (lRec.value.isZero())
			return false;
		lRec.pdf *= lumPdf;
		lRec.value /= lRec.pdf;
		lRec.luminaire = luminaire;
		return true;
	}
	return false;
}

bool Scene::sampleAttenuatedLuminaire(const Intersection &its, 
	const Medium *medium, LuminaireSamplingRecord &lRec, 
	const Point2 &s, Sampler *sampler) const {
	Point2 sample(s);
	Float lumPdf;
	const Luminaire *luminaire = m_luminaires[
		m_luminairePDF.sampleReuse(sample.x, lumPdf)];
	luminaire->sample(its.p, lRec, sample);

	if (lRec.pdf != 0) {
		if (its.isMediumTransition())
			medium = its.getTargetMedium(lRec.sRec.p - its.p);
		lRec.value *= getTransmittance(its.p, lRec.sRec.p, its.time, medium, sampler);
		if (lRec.value.isZero())
			return false;
		lRec.pdf *= lumPdf;
		lRec.value /= lRec.pdf;
		lRec.luminaire = luminaire;
		return true;
	}
	return false;
}


void Scene::sampleEmission(EmissionRecord &eRec, Point2 &sample1, Point2 &sample2) const {
	Float lumPdf;
	size_t index = m_luminairePDF.sampleReuse(sample1.x, lumPdf);
	const Luminaire *luminaire = m_luminaires[index];
	luminaire->sampleEmission(eRec, sample1, sample2);
	eRec.pdfArea *= lumPdf;
	eRec.luminaire = luminaire;
	Float cosTheta = (eRec.luminaire->getType() & Luminaire::EOnSurface)
		? absDot(eRec.sRec.n, eRec.d) : 1;
	eRec.value *= cosTheta / (eRec.pdfArea * eRec.pdfDir);
}

void Scene::sampleEmissionArea(EmissionRecord &eRec, Point2 &sample) const {
	Float lumPdf;
	size_t index = m_luminairePDF.sampleReuse(sample.x, lumPdf);
	const Luminaire *luminaire = m_luminaires[index];
	luminaire->sampleEmissionArea(eRec, sample);
	eRec.pdfArea *= lumPdf;
	eRec.luminaire = luminaire;
}

void Scene::pdfEmission(EmissionRecord &eRec, bool delta) const {
	const Luminaire *luminaire = eRec.luminaire;
	Float luminance;
	if (m_importanceSampleLuminaires)
		luminance = luminaire->getSamplingWeight();
	else 
		luminance = 1.0f;
	/* Calculate the probability of importance sampling this luminaire */
	const Float fraction = luminance / m_luminairePDF.getOriginalSum();

	luminaire->pdfEmission(eRec, delta);
	eRec.pdfArea *= fraction;
}

void Scene::addChild(const std::string &name, ConfigurableObject *child) {
	const Class *cClass = child->getClass();

	if (cClass->derivesFrom(MTS_CLASS(NetworkedObject)) &&
	   !cClass->derivesFrom(MTS_CLASS(Integrator))) 
		m_netObjects.push_back(static_cast<NetworkedObject *>(child));
	if (cClass->derivesFrom(MTS_CLASS(Camera))) {
		AssertEx(m_camera == NULL, "There can only be one camera per scene");
		m_camera = static_cast<Camera *>(child);
		m_sampler = m_camera->getSampler();
	} else if (cClass->derivesFrom(MTS_CLASS(Integrator))) {
		AssertEx(m_integrator == NULL, "There can only be one integrator per scene");
		m_integrator = static_cast<Integrator *>(child);
	} else if (cClass->derivesFrom(MTS_CLASS(Texture))
			|| cClass->derivesFrom(MTS_CLASS(BSDF))
			|| cClass->derivesFrom(MTS_CLASS(PhaseFunction))) {
		ConfigurableObject *obj= static_cast<ConfigurableObject *>(child);
		obj->incRef();
		m_objects.push_back(obj);
	} else if (cClass->derivesFrom(MTS_CLASS(Medium))) {
		Medium *medium = static_cast<Medium *>(child);
		if (m_media.find(medium) == m_media.end()) {
			medium->incRef();
			m_media.insert(medium);
		}
	} else if (cClass->derivesFrom(MTS_CLASS(Luminaire))) {
		Luminaire *luminaire = static_cast<Luminaire *>(child);

		if (luminaire->isCompound()) {
			int index = 0;
			do {
				ref<Luminaire> element = luminaire->getElement(index++);
				if (element == NULL)
					break;
				addChild(name, element);
			} while (true);
			return;
		}

		luminaire->incRef();
		m_luminaires.push_back(luminaire);
		if (luminaire->isBackgroundLuminaire()) {
			AssertEx(m_backgroundLuminaire.get() == NULL,
				"The scene may only contain one background luminaire");
			m_backgroundLuminaire = luminaire;
		}
	} else if (cClass->derivesFrom(MTS_CLASS(Shape))) {
		ref<Shape> shape = static_cast<Shape *>(child);
		shape->incRef();
		m_shapes.push_back(shape);
	} else if (cClass->derivesFrom(MTS_CLASS(Scene))) {
		ref<Scene> scene = static_cast<Scene *>(child);
		/* A scene from somewhere else has been included.
		   Add all of its contents */
		for (size_t i=0; i<scene->getLuminaires().size(); ++i) {
			Luminaire *lum = scene->getLuminaires()[i];
			lum->setParent(this);
			addChild("luminaire", lum);
		}
		for (size_t i=0; i<scene->getShapes().size(); ++i) {
			Shape *shape = scene->getShapes()[i];
			shape->setParent(this);
			addChild("shape", shape);
		}
		for (size_t i=0; i<scene->getReferencedObjects().size(); ++i) {
			ConfigurableObject *obj = scene->getReferencedObjects()[i];
			obj->setParent(this);
			addChild("object", obj);
		}
		for (std::set<Medium *>::iterator it = scene->getMedia().begin();
				it != scene->getMedia().end(); ++it) {
			Medium *medium = *it;
			medium->setParent(this);
			addChild("medium", medium);
		}
		if (scene->getIntegrator() != NULL) {
			scene->getIntegrator()->setParent(this);
			addChild("integrator", scene->getIntegrator());
		}
		if (scene->getCamera() != NULL) {
			scene->getCamera()->setParent(this);
			addChild("integrator", scene->getCamera());
		}
	} else {
		Log(EError, "Scene: Invalid child node! (\"%s\")",
			cClass->getName().c_str());
	}
}

void Scene::addShape(Shape *shape) {
	if (shape->isCompound()) {
		int index = 0;
		do {
			ref<Shape> element = shape->getElement(index++);
			if (element == NULL)
				break;
			addShape(element);
		} while (true);
	} else {
		if (shape->isLuminaire()) {
			if (std::find(m_luminaires.begin(), m_luminaires.end(), 
					shape->getLuminaire()) == m_luminaires.end()) {
				m_luminaires.push_back(shape->getLuminaire());
				shape->getLuminaire()->incRef();
			}
		}
		if (shape->hasSubsurface()) {
			m_netObjects.push_back(static_cast<NetworkedObject *>(shape->getSubsurface()));
			if (std::find(m_ssIntegrators.begin(), m_ssIntegrators.end(), 
					shape->getSubsurface()) == m_ssIntegrators.end()) {
				m_ssIntegrators.push_back(shape->getSubsurface());
				shape->getSubsurface()->incRef();
			}
		}

		Medium *iMedium = shape->getInteriorMedium(),
		       *eMedium = shape->getExteriorMedium();

		if (eMedium != NULL && m_media.find(eMedium) == m_media.end()) {
			m_media.insert(eMedium);
			eMedium->incRef();
		}

		if (iMedium != NULL && m_media.find(iMedium) == m_media.end()) {
			m_media.insert(iMedium);
			iMedium->incRef();
		}

		if (shape->getClass()->derivesFrom(MTS_CLASS(TriMesh))) {
			if (std::find(m_meshes.begin(), m_meshes.end(), shape)
					== m_meshes.end()) {
				m_meshes.push_back(static_cast<TriMesh *>(shape));
				shape->incRef();
			}
		}

		shape->incRef();
		m_kdtree->addShape(shape);
		m_shapes.push_back(shape);
	}
}

void Scene::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);

	stream->writeFloat(m_kdtree->getQueryCost());
	stream->writeFloat(m_kdtree->getTraversalCost());
	stream->writeFloat(m_kdtree->getEmptySpaceBonus());
	stream->writeInt(m_kdtree->getStopPrims());
	stream->writeBool(m_kdtree->getClip());
	stream->writeUInt(m_kdtree->getMaxDepth());
	stream->writeUInt(m_kdtree->getExactPrimitiveThreshold());
	stream->writeBool(m_kdtree->getParallelBuild());
	stream->writeBool(m_kdtree->getRetract());
	stream->writeUInt(m_kdtree->getMaxBadRefines());
	stream->writeBool(m_importanceSampleLuminaires);
	stream->writeInt(m_testType);
	stream->writeFloat(m_testThresh);
	stream->writeInt(m_blockSize);
	m_aabb.serialize(stream);
	m_bsphere.serialize(stream);
	manager->serialize(stream, m_backgroundLuminaire.get());
	stream->writeSize(m_shapes.size());
	for (size_t i=0; i<m_shapes.size(); ++i) 
		manager->serialize(stream, m_shapes[i]);
	stream->writeSize(m_meshes.size());
	for (size_t i=0; i<m_meshes.size(); ++i) 
		manager->serialize(stream, m_meshes[i]);
	stream->writeSize(m_luminaires.size());
	for (size_t i=0; i<m_luminaires.size(); ++i) 
		manager->serialize(stream, m_luminaires[i]);
	stream->writeSize(m_media.size());
	for (std::set<Medium *>::const_iterator it = m_media.begin();
			it != m_media.end(); ++it)
		manager->serialize(stream, *it);
	stream->writeSize(m_ssIntegrators.size());
	for (size_t i=0; i<m_ssIntegrators.size(); ++i) 
		manager->serialize(stream, m_ssIntegrators[i]);
	stream->writeSize(m_objects.size());
	for (size_t i=0; i<m_objects.size(); ++i) 
		manager->serialize(stream, m_objects[i]);
	stream->writeSize(m_netObjects.size());
	for (size_t i=0; i<m_netObjects.size(); ++i) 
		manager->serialize(stream, m_netObjects[i]);
}

std::string Scene::toString() const {
	std::ostringstream oss;

	oss << "Scene[" << endl
		<< "  testType = " << ((m_testType == ETTest) ? "t-test" : "relerr") << ", " << endl
		<< "  testThresh = " << m_testThresh << ", " << endl
		<< "  importanceSampleLuminaires = " << (int) m_importanceSampleLuminaires << ", " << endl
		<< "  camera = " << indent(m_camera.toString()) << "," << endl
		<< "  sampler = " << indent(m_sampler.toString()) << "," << endl
		<< "  integrator = " << indent(m_integrator.toString()) << "," << endl
		<< "  kdtree = " << indent(m_kdtree.toString()) << "," << endl
		<< "  backgroundLuminaire = " << indent(m_backgroundLuminaire.toString()) << "," << endl
		<< "  meshes = " << indent(containerToString(m_meshes.begin(), m_meshes.end())) << "," << endl
		<< "  shapes = " << indent(containerToString(m_shapes.begin(), m_shapes.end())) << "," << endl
		<< "  luminaires = " << indent(containerToString(m_luminaires.begin(), m_luminaires.end())) << "," << endl
		<< "  media = " << indent(containerToString(m_media.begin(), m_media.end())) << "," << endl
		<< "  ssIntegrators = " << indent(containerToString(m_ssIntegrators.begin(), m_ssIntegrators.end())) << "," << endl
		<< "  objects = " << indent(containerToString(m_objects.begin(), m_objects.end())) << endl;
	oss << "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS_S(Scene, false, ConfigurableObject)
MTS_NAMESPACE_END

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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/core/plugin.h>

MTS_NAMESPACE_BEGIN

Scene::Scene(const Properties &props)
 : NetworkedObject(props), m_blockSize(32) {
	m_kdtree = new KDTree();
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
	m_testThresh = props.getFloat("testThresh", 0.01);
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
		m_kdtree->setIntersectionCost(props.getFloat("kdIntersectionCost"));
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
	for (size_t i=0; i<m_media.size(); ++i)
		m_media[i]->incRef();
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
	m_kdtree = new KDTree();
	m_kdtree->setIntersectionCost(stream->readFloat());
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
	int count = stream->readInt();
	for (int i=0; i<count; ++i) {
		Shape *shape = static_cast<Shape *>(manager->getInstance(stream));
		shape->incRef();
		m_shapes.push_back(shape);
	}
	count = stream->readInt();
	for (int i=0; i<count; ++i) {
		TriMesh *trimesh = static_cast<TriMesh *>(manager->getInstance(stream));
		trimesh->incRef();
		m_meshes.push_back(trimesh);
	}
	count = stream->readInt();
	for (int i=0; i<count; ++i) {
		Luminaire *luminaire = static_cast<Luminaire *>(manager->getInstance(stream));
		luminaire->incRef();
		m_luminaires.push_back(luminaire);
	}
	count = stream->readInt();
	for (int i=0; i<count; ++i) {
		Medium *medium = static_cast<Medium *>(manager->getInstance(stream));
		medium->incRef();
		m_media.push_back(medium);
	}
	count = stream->readInt();
	for (int i=0; i<count; ++i) {
		Subsurface *ssIntegrator = static_cast<Subsurface *>(manager->getInstance(stream));
		ssIntegrator->incRef();
		m_ssIntegrators.push_back(ssIntegrator);
	}
	count = stream->readInt();
	for (int i=0; i<count; ++i) {
		ConfigurableObject *obj = static_cast<ConfigurableObject *>(manager->getInstance(stream));
		obj->incRef();
		m_objects.push_back(obj);
	}
	count = stream->readInt();
	for (int i=0; i<count; ++i) {
		NetworkedObject *obj = static_cast<NetworkedObject *>(manager->getInstance(stream));
		m_netObjects.push_back(obj); // Do not increase the ref. count
	}
	initialize();
}

Scene::~Scene() {
	for (unsigned int i=0; i<m_shapes.size(); i++)
		m_shapes[i]->decRef();
	for (unsigned int i=0; i<m_meshes.size(); i++)
		m_meshes[i]->decRef();
	for (unsigned int i=0; i<m_media.size(); i++) 
		m_media[i]->decRef();
	for (unsigned int i=0; i<m_objects.size(); i++)
		m_objects[i]->decRef();
	for (unsigned int i=0; i<m_ssIntegrators.size(); i++)
		m_ssIntegrators[i]->decRef();
	for (unsigned int i=0; i<m_luminaires.size(); i++)
		m_luminaires[i]->decRef();
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
				createObject(Integrator::m_theClass, Properties("direct")));
		m_integrator->configure();
	}
	if (m_camera == NULL) {
		Log(EWarn, "No camera found! Adding a default camera.");

		Properties props("perspective");
		/* Create a perspective camera with 45deg. FOV, which can see the whole scene */
		AABB aabb;
		for (size_t i=0; i<m_shapes.size(); ++i)
			aabb.expandBy(m_shapes[i]->getAABB());
		for (size_t i=0; i<m_media.size(); ++i)
			aabb.expandBy(m_media[i]->getAABB());
		if (!aabb.isValid())
			Log(EError, "Unable to set up a default camera -- does the scene contain anything at all?");
		Point center = aabb.getCenter();
		Vector extents = aabb.getExtents();
		Float maxExtents = std::max(extents.x, extents.y);
		Float distance = maxExtents/(2.0f * std::tan(45 * .5f * M_PI/180));

		props.setTransform("toWorld", Transform::translate(Vector(center.x, center.y, aabb.getMinimum().x - distance)));
		props.setFloat("fov", 45.0f);

		m_camera = static_cast<Camera *> (PluginManager::getInstance()->createObject(Camera::m_theClass, props));
		m_camera->configure();
		m_sampler = m_camera->getSamplerX();
	}

	if (m_media.size() > 1)
		Log(EError, "Scenes are currently restricted to at most one participating medium.");

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
		for (unsigned int i=0; i<m_media.size(); i++) {
			const AABB &aabb = m_media[i]->getAABB();
			m_aabb.expandBy(aabb);
			for (int j=0; j<8; ++j)
				m_bsphere.expandBy(aabb.getCorner(j));
		}
	}

	if (!m_luminairePDF.isReady()) {
		if (m_luminaires.size() == 0) {
			Log(EWarn, "No luminaires found -- adding a constant environment source");
			Properties constantProps("constant");
			constantProps.setSpectrum("intensity", Spectrum(0.8f));
			addChild("", PluginManager::getInstance()->createObject(Luminaire::m_theClass, constantProps));
		}

		/* Calculate a discrete PDF to importance sample luminaires */
		for (std::vector<Luminaire *>::iterator it = m_luminaires.begin();
			it != m_luminaires.end(); ++it) {
			(*it)->preprocess(this);
			/* Add with a probability proportional to the luminaire's power */
			if (m_importanceSampleLuminaires)
				m_luminairePDF.put((*it)->getPower().getLuminance());
			else
				m_luminairePDF.put(1.0f);
		}
		m_luminairePDF.build();
	}
}

bool Scene::preprocess(RenderQueue *queue, const RenderJob *job, 
		int sceneResID, int cameraResID, int samplerResID) {
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

	/* Pre-process step for all participating media */
	for (std::vector<Medium *>::iterator it = m_media.begin();
		it != m_media.end(); ++it) 
		(*it)->preprocess(this, queue, job,
			sceneResID, cameraResID, samplerResID);
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
		const LuminaireSamplingRecord &lRec) const {
	const Luminaire *luminaire = lRec.luminaire;
	Float luminance;

	if (m_importanceSampleLuminaires)
		luminance = luminaire->getPower().getLuminance();
	else 
		luminance = 1.0f;

	/* Calculate the probability of importance sampling this luminaire */
	const Float fraction = luminance / m_luminairePDF.getOriginalSum();
	return luminaire->pdf(p, lRec) * fraction;
}

Float Scene::pdfLuminaire(const Intersection &its,
		const LuminaireSamplingRecord &lRec) const {
	const Luminaire *luminaire = lRec.luminaire;
	Float luminance;

	if (m_importanceSampleLuminaires)
		luminance = luminaire->getPower().getLuminance();
	else 
		luminance = 1.0f;

	/* Calculate the probability of importance sampling this luminaire */
	const Float fraction = luminance / m_luminairePDF.getOriginalSum();
	return luminaire->pdf(its, lRec) * fraction;
}

bool Scene::sampleLuminaire(const Point &p,
		LuminaireSamplingRecord &lRec, const Point2 &s,
		bool testVisibility) const {
	Point2 sample(s);
	Float lumPdf;
	unsigned int index = m_luminairePDF.sampleReuse(sample.x, lumPdf);
	const Luminaire *luminaire = m_luminaires[index];
	luminaire->sample(p, lRec, sample);

	if (lRec.pdf != 0) {
		if (testVisibility && isOccluded(p, lRec.sRec.p)) 
			return false;
		lRec.pdf *= lumPdf;
		lRec.Le /= lRec.pdf;
		lRec.luminaire = luminaire;
		return true;
	} else {
		return false;
	}
}

bool Scene::sampleLuminaire(const Intersection &its,
		LuminaireSamplingRecord &lRec, const Point2 &s,
		bool testVisibility) const {
	Point2 sample(s);
	Float lumPdf;
	unsigned int index = m_luminairePDF.sampleReuse(sample.x, lumPdf);
	const Luminaire *luminaire = m_luminaires[index];
	luminaire->sample(its, lRec, sample);

	if (lRec.pdf != 0) {
		if (testVisibility && isOccluded(its.p, lRec.sRec.p)) 
			return false;
		lRec.pdf *= lumPdf;
		lRec.Le /= lRec.pdf;
		lRec.luminaire = luminaire;
		return true;
	} else {
		return false;
	}
}

void Scene::sampleEmission(EmissionRecord &eRec, Point2 &sample1, Point2 &sample2) const {
	Float lumPdf;
	unsigned int index = m_luminairePDF.sampleReuse(sample1.x, lumPdf);
	const Luminaire *luminaire = m_luminaires[index];
	luminaire->sampleEmission(eRec, sample1, sample2);
	eRec.pdfArea *= lumPdf;
	eRec.luminaire = luminaire;
	Float cosTheta = (eRec.luminaire->getType() & Luminaire::EOnSurface) ? absDot(eRec.sRec.n, eRec.d) : 1;
	eRec.P *= cosTheta / (eRec.pdfArea * eRec.pdfDir);
}

void Scene::sampleEmissionArea(EmissionRecord &eRec, Point2 &sample) const {
	Float lumPdf;
	unsigned int index = m_luminairePDF.sampleReuse(sample.x, lumPdf);
	const Luminaire *luminaire = m_luminaires[index];
	luminaire->sampleEmissionArea(eRec, sample);
	eRec.pdfArea *= lumPdf;
	eRec.luminaire = luminaire;
}

Spectrum Scene::sampleEmissionDirection(EmissionRecord &eRec, Point2 &sample) const {
	return eRec.luminaire->sampleEmissionDirection(eRec, sample);
}

void Scene::pdfEmission(EmissionRecord &eRec) const {
	const Luminaire *luminaire = eRec.luminaire;
	Float luminance;
	if (m_importanceSampleLuminaires)
		luminance = luminaire->getPower().getLuminance();
	else 
		luminance = 1.0f;
	/* Calculate the probability of importance sampling this luminaire */
	const Float fraction = luminance / m_luminairePDF.getOriginalSum();

	luminaire->pdfEmission(eRec);
	eRec.pdfArea *= fraction;
}

Spectrum Scene::LeBackground(const Ray &ray) const {
	Spectrum Le(0.0f);
	for (std::vector<Luminaire *>::const_iterator it = m_luminaires.begin();
		it != m_luminaires.end(); ++it)
		Le += (*it)->Le(ray);
	return Le;
}

Spectrum Scene::getAttenuation(const Ray &ray) const {
	Spectrum extinction(0.0f);
	for (std::vector<Medium *>::const_iterator it = 
			m_media.begin(); it != m_media.end(); ++it) {
		extinction += (*it)->tau(ray);
	}
	if (extinction.isBlack())
		return Spectrum(1.0f);
	else
		return (-extinction).exp();
}

bool Scene::sampleDistance(const Ray &ray, Float maxDist, MediumSamplingRecord &mRec, 
			Sampler *sampler) const {
	if (m_media.size() > 0) {
		return m_media[0]->sampleDistance(ray, maxDist, mRec, sampler);
	} else {
		mRec.pdf = 1.0f;
		mRec.miWeight = 1.0f;
		mRec.attenuation = Spectrum(1.0f);
		return false;
	}
}

void Scene::addChild(const std::string &name, ConfigurableObject *child) {
	const Class *cClass = child->getClass();

	if (cClass->derivesFrom(NetworkedObject::m_theClass) &&
	   !cClass->derivesFrom(Integrator::m_theClass)) 
		m_netObjects.push_back(static_cast<NetworkedObject *>(child));
	if (cClass->derivesFrom(Camera::m_theClass)) {
		AssertEx(m_camera == NULL, "There can only be one camera per scene");
		m_camera = static_cast<Camera *>(child);
		m_sampler = m_camera->getSamplerX();
	} else if (cClass->derivesFrom(Integrator::m_theClass)) {
		AssertEx(m_integrator == NULL, "There can only be one integrator per scene");
		m_integrator = static_cast<Integrator *>(child);
	} else if (cClass->derivesFrom(Texture::m_theClass)
			|| cClass->derivesFrom(BSDF::m_theClass)) {
		ConfigurableObject *obj= static_cast<ConfigurableObject *>(child);
		obj->incRef();
		m_objects.push_back(obj);
	} else if (cClass->derivesFrom(Medium::m_theClass)) {
		Medium *medium = static_cast<Medium *>(child);
		medium->incRef();
		m_media.push_back(medium);
	} else if (cClass->derivesFrom(Luminaire::m_theClass)) {
		Luminaire *luminaire = static_cast<Luminaire *>(child);
		luminaire->incRef();
		m_luminaires.push_back(luminaire);
		if (luminaire->isBackgroundLuminaire()) {
			AssertEx(m_backgroundLuminaire.get() == NULL,
				"The scene may only contain one background luminaire");
			m_backgroundLuminaire = luminaire;
		}
	} else if (cClass->derivesFrom(Shape::m_theClass)) {
		ref<Shape> shape = static_cast<Shape *>(child);
		shape->incRef();
		m_shapes.push_back(shape);
	} else if (cClass->derivesFrom(Scene::m_theClass)) {
		ref<Scene> scene = static_cast<Scene *>(child);
		/* A scene from somewhere else has been included.
		   Add all of its contents */
		for (unsigned int i=0; i<scene->getLuminaires().size(); ++i) {
			Luminaire *lum = scene->getLuminaires()[i];
			lum->setParent(this);
			addChild("luminaire", lum);
		}
		for (unsigned int i=0; i<scene->getShapes().size(); ++i) {
			Shape *shape = scene->getShapes()[i];
			shape->setParent(this);
			addChild("shape", shape);
		}
		for (unsigned int i=0; i<scene->getReferencedObjects().size(); ++i) {
			ConfigurableObject *obj = scene->getReferencedObjects()[i];
			obj->setParent(this);
			addChild("object", obj);
		}
		for (unsigned int i=0; i<scene->getMedia().size(); ++i) {
			Medium *medium = scene->getMedia()[i];
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
		if (shape->getClass()->derivesFrom(TriMesh::m_theClass)) {
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

	stream->writeFloat(m_kdtree->getIntersectionCost());
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
	stream->writeUInt((unsigned int) m_shapes.size());
	for (size_t i=0; i<m_shapes.size(); ++i) 
		manager->serialize(stream, m_shapes[i]);
	stream->writeUInt((unsigned int) m_meshes.size());
	for (size_t i=0; i<m_meshes.size(); ++i) 
		manager->serialize(stream, m_meshes[i]);
	stream->writeUInt((unsigned int) m_luminaires.size());
	for (size_t i=0; i<m_luminaires.size(); ++i) 
		manager->serialize(stream, m_luminaires[i]);
	stream->writeUInt((unsigned int) m_media.size());
	for (size_t i=0; i<m_media.size(); ++i) 
		manager->serialize(stream, m_media[i]);
	stream->writeUInt((unsigned int) m_ssIntegrators.size());
	for (size_t i=0; i<m_ssIntegrators.size(); ++i) 
		manager->serialize(stream, m_ssIntegrators[i]);
	stream->writeUInt((unsigned int) m_objects.size());
	for (size_t i=0; i<m_objects.size(); ++i) 
		manager->serialize(stream, m_objects[i]);
	stream->writeUInt((unsigned int) m_netObjects.size());
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
		<< "  meshes = " << indent(listToString(m_meshes)) << "," << endl
		<< "  shapes = " << indent(listToString(m_shapes)) << "," << endl
		<< "  luminaires = " << indent(listToString(m_luminaires)) << "," << endl
		<< "  media = " << indent(listToString(m_media)) << "," << endl
		<< "  ssIntegrators = " << indent(listToString(m_ssIntegrators)) << "," << endl
		<< "  objects = " << indent(listToString(m_objects)) << endl;
	oss << "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS_S(Scene, false, ConfigurableObject)
MTS_NAMESPACE_END

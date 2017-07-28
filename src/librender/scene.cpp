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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>

#define DEFAULT_BLOCKSIZE 32

MTS_NAMESPACE_BEGIN

// ===========================================================================
//         Constructors, destructor and serialization-related code
// ===========================================================================

Scene::Scene()
 : NetworkedObject(Properties()), m_blockSize(DEFAULT_BLOCKSIZE) {
    m_kdtree = new ShapeKDTree();
    m_sourceFile = new fs::path();
    m_destinationFile = new fs::path();
}

Scene::Scene(const Properties &props)
 : NetworkedObject(props), m_blockSize(DEFAULT_BLOCKSIZE) {
    m_kdtree = new ShapeKDTree();
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
    m_sourceFile = new fs::path();
    m_destinationFile = new fs::path();
}

Scene::Scene(Scene *scene) : NetworkedObject(Properties()) {
    m_kdtree = scene->m_kdtree;
    m_blockSize = scene->m_blockSize;
    m_aabb = scene->m_aabb;
    m_environmentEmitter = scene->m_environmentEmitter;
    m_sensor = scene->m_sensor;
    m_integrator = scene->m_integrator;
    m_sourceFile = new fs::path(*scene->m_sourceFile);
    m_destinationFile = new fs::path(*scene->m_destinationFile);
    m_emitterPDF = scene->m_emitterPDF;
    m_shapes = scene->m_shapes;
    m_sensors = scene->m_sensors;
    m_meshes = scene->m_meshes;
    m_emitters = scene->m_emitters;
    m_media = scene->m_media;
    m_ssIntegrators = scene->m_ssIntegrators;
    m_objects = scene->m_objects;
    m_netObjects = scene->m_netObjects;
    m_specialShapes = scene->m_specialShapes;
    m_degenerateSensor = scene->m_degenerateSensor;
    m_degenerateEmitters = scene->m_degenerateEmitters;
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
    m_blockSize = stream->readUInt();
    m_degenerateSensor = stream->readBool();
    m_degenerateEmitters = stream->readBool();
    m_aabb = AABB(stream);
    m_environmentEmitter = static_cast<Emitter *>(manager->getInstance(stream));
    m_sourceFile = new fs::path(stream->readString());
    m_destinationFile = new fs::path(stream->readString());

    size_t count = stream->readSize();
    m_shapes.reserve(count);
    for (size_t i=0; i<count; ++i)
        m_shapes.push_back(static_cast<Shape *>(manager->getInstance(stream)));
    count = stream->readSize();
    m_specialShapes.reserve(count);
    for (size_t i=0; i<count; ++i)
        m_specialShapes.push_back(static_cast<Shape *>(manager->getInstance(stream)));
    count = stream->readSize();
    m_meshes.reserve(count);
    for (size_t i=0; i<count; ++i)
        m_meshes.push_back(static_cast<TriMesh *>(manager->getInstance(stream)));
    count = stream->readSize();
    m_sensors.reserve(count);
    for (size_t i=0; i<count; ++i)
        m_sensors.push_back(static_cast<Sensor *>(manager->getInstance(stream)));
    count = stream->readSize();
    m_emitters.reserve(count);
    for (size_t i=0; i<count; ++i)
        m_emitters.push_back(static_cast<Emitter *>(manager->getInstance(stream)));
    count = stream->readSize();
    m_media.reserve(count);
    for (size_t i=0; i<count; ++i)
        m_media.push_back(static_cast<Medium *>(manager->getInstance(stream)));
    count = stream->readSize();
    m_ssIntegrators.reserve(count);
    for (size_t i=0; i<count; ++i)
        m_ssIntegrators.push_back(static_cast<Subsurface *>(manager->getInstance(stream)));
    count = stream->readSize();
    m_objects.reserve(count);
    for (size_t i=0; i<count; ++i)
        m_objects.push_back(static_cast<ConfigurableObject *>(manager->getInstance(stream)));
    count = stream->readSize();
    m_netObjects.reserve(count);
    for (size_t i=0; i<count; ++i)
        m_netObjects.push_back(static_cast<NetworkedObject *>(manager->getInstance(stream)));

    initialize();
}

Scene::~Scene() {
    delete m_destinationFile;
    delete m_sourceFile;
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
    stream->writeUInt(m_blockSize);
    stream->writeBool(m_degenerateSensor);
    stream->writeBool(m_degenerateEmitters);
    m_aabb.serialize(stream);
    manager->serialize(stream, m_environmentEmitter.get());
    stream->writeString(m_sourceFile->string());
    stream->writeString(m_destinationFile->string());

    stream->writeSize(m_shapes.size());
    for (size_t i=0; i<m_shapes.size(); ++i)
        manager->serialize(stream, m_shapes[i].get());

    stream->writeSize(m_specialShapes.size());
    for (size_t i=0; i<m_specialShapes.size(); ++i)
        manager->serialize(stream, m_specialShapes[i].get());

    stream->writeSize(m_meshes.size());
    for (size_t i=0; i<m_meshes.size(); ++i)
        manager->serialize(stream, m_meshes[i]);

    stream->writeSize(m_sensors.size());
    for (size_t i=0; i<m_sensors.size(); ++i)
        manager->serialize(stream, m_sensors[i].get());

    stream->writeSize(m_emitters.size());
    for (size_t i=0; i<m_emitters.size(); ++i)
        manager->serialize(stream, m_emitters[i].get());

    stream->writeSize(m_media.size());
    for (ref_vector<Medium>::const_iterator it = m_media.begin();
            it != m_media.end(); ++it)
        manager->serialize(stream, it->get());

    stream->writeSize(m_ssIntegrators.size());
    for (ref_vector<Subsurface>::const_iterator it = m_ssIntegrators.begin();
            it != m_ssIntegrators.end(); ++it)
        manager->serialize(stream, it->get());

    stream->writeSize(m_objects.size());
    for (ref_vector<ConfigurableObject>::const_iterator it = m_objects.begin();
            it != m_objects.end(); ++it)
        manager->serialize(stream, it->get());

    stream->writeSize(m_netObjects.size());
    for (ref_vector<NetworkedObject>::const_iterator it = m_netObjects.begin();
            it != m_netObjects.end(); ++it)
        manager->serialize(stream, it->get());
}

// ===========================================================================
//          Scene initialization, rendering, and miscellaneous methods
// ===========================================================================

void Scene::bindUsedResources(ParallelProcess *proc) const {
    for (ref_vector<NetworkedObject>::const_iterator it = m_netObjects.begin();
            it != m_netObjects.end(); ++it)
        it->get()->bindUsedResources(proc);
}

void Scene::wakeup(ConfigurableObject *,
        std::map<std::string, SerializableObject *> &params) {
    for (ref_vector<NetworkedObject>::iterator it = m_netObjects.begin();
            it != m_netObjects.end(); ++it)
        (*it)->wakeup(this, params);
}

void Scene::setSensor(Sensor *sensor) {
    m_sensor = sensor;
    m_degenerateSensor = sensor->getType() & Sensor::EDeltaPosition;
}

void Scene::removeSensor(Sensor *sensor) {
    if (!sensor)
        return;
    ref<Sensor> oldSensor = sensor;
    m_sensors.erase(std::remove(m_sensors.begin(),
        m_sensors.end(), oldSensor));
}

void Scene::addSensor(Sensor *sensor) {
    ref<Sensor> newSensor = sensor;
    if (!newSensor || std::find(m_sensors.begin(),
            m_sensors.end(), newSensor) != m_sensors.end())
        return;
    m_sensors.push_back(newSensor);
}

void Scene::configure() {
    if (m_integrator == NULL) {
        /* Create a direct integrator by default */
        m_integrator = static_cast<Integrator *> (PluginManager::getInstance()->
                createObject(MTS_CLASS(Integrator), Properties("direct")));
        m_integrator->configure();
    }

    if (m_sensor == NULL) {
        if (m_sensors.size() == 0) {
            Log(EInfo, "No sensors found! Adding a perspective camera..");
            Properties props("perspective");
            props.setFloat("fov", 45.0f);

            /* Create a perspective camera with a 45 deg. field of view
               and positioned so that it can see the entire scene */
            AABB aabb;
            for (ref_vector<Shape>::iterator it = m_shapes.begin();
                    it != m_shapes.end(); ++it)
                aabb.expandBy(it->get()->getAABB());
            if (aabb.isValid()) {
                Point center = aabb.getCenter();
                Vector extents = aabb.getExtents();
                Float maxExtentsXY = std::max(extents.x, extents.y);
                Float distance = maxExtentsXY/(2.0f * std::tan(45 * .5f * M_PI/180));
                Float maxExtentsXYZ = std::max(extents.z, maxExtentsXY);

                props.setFloat("farClip", maxExtentsXYZ * 5 + distance);
                props.setFloat("nearClip", distance / 100);

                props.setFloat("focusDistance", distance + extents.z/2);
                props.setTransform("toWorld", Transform::translate(Vector(center.x,
                        center.y, aabb.min.z - distance)));
            }
            Sensor *sensor = static_cast<Sensor *> (PluginManager::getInstance()->
                    createObject(MTS_CLASS(Sensor), props));
            sensor->configure();
            m_sensors.push_back(sensor);
        }
        m_sensor = m_sensors[m_sensors.size()-1];
    }
    m_sampler = m_sensor->getSampler();

    m_integrator->configureSampler(this, m_sampler);
}

void Scene::invalidate() {
    m_kdtree = new ShapeKDTree();
}

void Scene::initialize() {
    if (!m_kdtree->isBuilt()) {
        /* Expand all geometry */
        ref_vector<Shape> temp;
        temp.reserve(m_shapes.size());

        m_shapes.ensureUnique();
        m_shapes.swap(temp);
        size_t primitiveCount = 0, effPrimitiveCount = 0;

        for (size_t i=0; i<temp.size(); ++i) {
            addShape(temp[i]);
            primitiveCount += temp[i]->getPrimitiveCount();
            effPrimitiveCount += temp[i]->getEffectivePrimitiveCount();
            temp[i] = NULL;
        }
        if (primitiveCount != effPrimitiveCount) {
            Log(EDebug, "Scene contains " SIZE_T_FMT " primitives. Due to "
                "instancing or other kinds of procedural geometry, the effective number of primitives is "
                SIZE_T_FMT ".", primitiveCount, effPrimitiveCount);
        }

        /* Build the kd-tree */
        m_kdtree->build();

        m_aabb = m_kdtree->getAABB();
    }

    /* Make sure that there are no duplicates */
    m_emitters.ensureUnique();
    m_media.ensureUnique();
    m_ssIntegrators.ensureUnique();
    m_objects.ensureUnique();
    m_netObjects.ensureUnique();

    if (!m_emitterPDF.isNormalized()) {
        if (m_emitters.size() == 0) {
            Log(EWarn, "No emitters found -- adding sun & sky.");
            /* This is not a particularly realistic sky -- it extends below the
               horizon and uses an enlarged sun :). This is done to get better
               results for arbitrary input (and with a path tracer). */

            Properties skyProps("sunsky");
            skyProps.setFloat("scale", 2);
            skyProps.setTransform("toWorld", Transform::rotate(Vector(0,1,0), -120.0f));
            skyProps.setBoolean("extend", true);
            skyProps.setFloat("sunRadiusScale", 15);
            ref<Emitter> emitter = static_cast<Emitter *>(
                PluginManager::getInstance()->createObject(MTS_CLASS(Emitter), skyProps));
            addChild(emitter);
            emitter->configure();
        }

        /* Calculate a discrete PDF to importance sample emitters */
        for (ref_vector<Emitter>::iterator it = m_emitters.begin();
                it != m_emitters.end(); ++it)
            m_emitterPDF.append(it->get()->getSamplingWeight());

        m_emitterPDF.normalize();
    }

    initializeBidirectional();
}

void Scene::initializeBidirectional() {
    m_aabb = m_kdtree->getAABB();
    m_degenerateEmitters = true;
    m_specialShapes.clear();

    if (m_sensor) {
        ref<Shape> shape = m_sensor->createShape(this);
        if (shape != NULL)
            m_specialShapes.push_back(shape);
        m_aabb.expandBy(m_sensor->getAABB());

        m_degenerateSensor = m_sensor->getType() & Sensor::EDeltaPosition;
    }

    AABB aabb(m_aabb);
    for (ref_vector<Emitter>::iterator it = m_emitters.begin();
            it != m_emitters.end(); ++it) {
        Emitter *emitter = it->get();

        ref<Shape> shape = emitter->createShape(this);
        if (shape != NULL)
            m_specialShapes.push_back(shape);

        aabb.expandBy(emitter->getAABB());
        if (!(emitter->getType() & Emitter::EDeltaPosition))
            m_degenerateEmitters = false;
}
    m_aabb = aabb;
}

bool Scene::preprocess(RenderQueue *queue, const RenderJob *job,
        int sceneResID, int sensorResID, int samplerResID) {

    initialize();

    /* Pre-process step for the main scene integrator */
    if (!m_integrator->preprocess(this, queue, job,
        sceneResID, sensorResID, samplerResID))
        return false;

    /* Pre-process step for all sub-surface integrators (each one in independence) */
    for (ref_vector<Subsurface>::iterator it = m_ssIntegrators.begin();
            it != m_ssIntegrators.end(); ++it)
        (*it)->setActive(false);

    for (ref_vector<Subsurface>::iterator it = m_ssIntegrators.begin();
        it != m_ssIntegrators.end(); ++it)
        if (!(*it)->preprocess(this, queue, job,
                sceneResID, sensorResID, samplerResID))
            return false;

    for (ref_vector<Subsurface>::iterator it = m_ssIntegrators.begin();
            it != m_ssIntegrators.end(); ++it)
        (*it)->setActive(true);

    return true;
}

bool Scene::render(RenderQueue *queue, const RenderJob *job,
        int sceneResID, int sensorResID, int samplerResID) {
    m_sensor->getFilm()->clear();
    return m_integrator->render(this, queue, job, sceneResID,
        sensorResID, samplerResID);
}

void Scene::cancel() {
    for (ref_vector<Subsurface>::iterator it = m_ssIntegrators.begin();
            it != m_ssIntegrators.end(); ++it)
        (*it)->cancel();
    m_integrator->cancel();
}

void Scene::flush(RenderQueue *queue, const RenderJob *job) {
    m_sensor->getFilm()->develop(this, queue->getRenderTime(job));
}

void Scene::setDestinationFile(const fs::path &name) {
    *m_destinationFile = name;
}

void Scene::setSourceFile(const fs::path &name) {
    *m_sourceFile = name;
}

void Scene::postprocess(RenderQueue *queue, const RenderJob *job,
        int sceneResID, int sensorResID, int samplerResID) {
    m_integrator->postprocess(this, queue, job, sceneResID,
        sensorResID, samplerResID);
    m_sensor->getFilm()->develop(this, queue->getRenderTime(job));
}

void Scene::addChild(const std::string &name, ConfigurableObject *child) {
    const Class *cClass = child->getClass();

    if (cClass->derivesFrom(MTS_CLASS(NetworkedObject)) &&
       !cClass->derivesFrom(MTS_CLASS(Integrator)))
        m_netObjects.push_back(static_cast<NetworkedObject *>(child));

    if (cClass->derivesFrom(MTS_CLASS(Sensor))) {
        m_sensors.push_back(static_cast<Sensor *>(child));
    } else if (cClass->derivesFrom(MTS_CLASS(Integrator))) {
        AssertEx(m_integrator == NULL, "There can only be one integrator per scene");
        m_integrator = static_cast<Integrator *>(child);
    } else if (cClass->derivesFrom(MTS_CLASS(Texture))
            || cClass->derivesFrom(MTS_CLASS(BSDF))
            || cClass->derivesFrom(MTS_CLASS(Subsurface))
            || cClass->derivesFrom(MTS_CLASS(PhaseFunction))) {
        m_objects.push_back(static_cast<ConfigurableObject *>(child));
    } else if (cClass->derivesFrom(MTS_CLASS(Medium))) {
        m_media.push_back(static_cast<Medium *>(child));
    } else if (cClass->derivesFrom(MTS_CLASS(Emitter))) {
        Emitter *emitter = static_cast<Emitter *>(child);

        if (emitter->isCompound()) {
            size_t index = 0;
            do {
                ref<Emitter> element = emitter->getElement(index++);
                if (element == NULL)
                    break;
                addChild(name, element);
            } while (true);
            return;
        }

        if (emitter->isEnvironmentEmitter()) {
            if (m_environmentEmitter != NULL)
                Log(EError, "The scene may only contain one environment emitter");
            m_environmentEmitter = emitter;
        }

        m_emitters.push_back(emitter);
    } else if (cClass->derivesFrom(MTS_CLASS(Shape))) {
        Shape *shape = static_cast<Shape *>(child);
        if (shape->isSensor()) // determine sensors as early as possible
            addSensor(shape->getSensor());
        m_shapes.push_back(shape);
    } else if (cClass->derivesFrom(MTS_CLASS(Scene))) {
        ref<Scene> scene = static_cast<Scene *>(child);
        /* A scene from somewhere else has been included.
           Add all of its contents */
        ref_vector<Emitter> &emitters = scene->getEmitters();
        ref_vector<Sensor> &sensors = scene->getSensors();
        ref_vector<Shape> &shapes = scene->getShapes();
        ref_vector<ConfigurableObject> &refObjects = scene->getReferencedObjects();

        for (size_t i=0; i<emitters.size(); ++i)
            addChild(emitters[i]);

        for (size_t i=0; i<sensors.size(); ++i)
            addChild(sensors[i]);

        for (size_t i=0; i<shapes.size(); ++i)
            addChild(shapes[i]);

        for (ref_vector<ConfigurableObject>::iterator it = refObjects.begin();
                it != refObjects.end(); ++it)
            addChild(it->get());

        for (ref_vector<Medium>::iterator it = scene->getMedia().begin();
                it != scene->getMedia().end(); ++it)
            addChild(it->get());

        if (scene->getIntegrator() != NULL)
            addChild(scene->getIntegrator());

        if (scene->getSensor() != NULL)
            addChild(scene->getSensor());
    } else {
        ConfigurableObject::addChild(name, child);
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
        if (shape->isSensor() && !m_sensors.contains(shape->getSensor()))
            m_sensors.push_back(shape->getSensor());
        if (shape->isEmitter())
            m_emitters.push_back(shape->getEmitter());
        if (shape->hasSubsurface()) {
            m_netObjects.push_back(shape->getSubsurface());
            m_ssIntegrators.push_back(shape->getSubsurface());
        }

        Medium *iMedium = shape->getInteriorMedium(),
               *eMedium = shape->getExteriorMedium();

        if (eMedium != NULL)
            m_media.push_back(eMedium);

        if (iMedium != NULL)
            m_media.push_back(iMedium);

        if (shape->getClass()->derivesFrom(MTS_CLASS(TriMesh)))
            m_meshes.push_back(static_cast<TriMesh *>(shape));

        m_kdtree->addShape(shape);
        m_shapes.push_back(shape);
    }
}

std::string Scene::toString() const {
    std::ostringstream oss;

    oss << "Scene[" << endl
        << "  sensor = " << indent(m_sensor.toString()) << "," << endl
        << "  sampler = " << indent(m_sampler.toString()) << "," << endl
        << "  integrator = " << indent(m_integrator.toString()) << "," << endl
        << "  kdtree = " << indent(m_kdtree.toString()) << "," << endl
        << "  environmentEmitter = " << indent(m_environmentEmitter.toString()) << "," << endl
        << "  shapes = " << indent(containerToString(m_shapes.begin(), m_shapes.end())) << "," << endl
        << "  emitters = " << indent(containerToString(m_emitters.begin(), m_emitters.end())) << "," << endl
        << "  media = " << indent(containerToString(m_media.begin(), m_media.end())) << "," << endl
        << "  sensors = " << indent(containerToString(m_sensors.begin(), m_sensors.end())) << "," << endl
        << "  ssIntegrators = " << indent(containerToString(m_ssIntegrators.begin(), m_ssIntegrators.end())) << "," << endl
        << "  objects = " << indent(containerToString(m_objects.begin(), m_objects.end())) << endl
        << "]";
    return oss.str();
}

// ===========================================================================
//                 Ray tracing and transmittance computations
// ===========================================================================

static StatsCounter mediumInconsistencies("General", "Detected medium inconsistencies");

Spectrum Scene::evalTransmittance(const Point &p1, bool p1OnSurface, const Point &p2, bool p2OnSurface,
        Float time, const Medium *medium, int &interactions, Sampler *sampler) const {
    Vector d = p2 - p1;
    Float remaining = d.length();
    d /= remaining;

    Float lengthFactor = p2OnSurface ? (1-ShadowEpsilon) : 1;
    Ray ray(p1, d, p1OnSurface ? Epsilon : 0, remaining * lengthFactor, time);
    Spectrum transmittance(1.0f);
    Intersection its;
    int maxInteractions = interactions;
    interactions = 0;

    while (remaining > 0) {
        Normal n;
        bool surface = rayIntersect(ray, its.t, its.shape, its.geoFrame.n, its.uv);

        if (surface && (interactions == maxInteractions ||
            !(its.getBSDF()->getType() & BSDF::ENull))) {
            /* Encountered an occluder -- zero transmittance. */
            return Spectrum(0.0f);
        }

        if (medium)
            transmittance *= medium->evalTransmittance(
                Ray(ray, 0, std::min(its.t, remaining)), sampler);

        if (!surface || transmittance.isZero())
            break;

        const BSDF *bsdf = its.getBSDF();

        its.p = ray.o;
        its.geoFrame = Frame(its.geoFrame.n);
        its.hasUVPartials = false;
        Vector wo = its.geoFrame.toLocal(ray.d);
        BSDFSamplingRecord bRec(its, -wo, wo, ERadiance);
        bRec.typeMask = BSDF::ENull;
        transmittance *= bsdf->eval(bRec, EDiscrete);

        if (its.isMediumTransition()) {
            if (medium != its.getTargetMedium(-d)) {
                ++mediumInconsistencies;
                return Spectrum(0.0f);
            }
            medium = its.getTargetMedium(d);
        }

        if (++interactions > 100) { /// Just a precaution..
            Log(EWarn, "evalTransmittance(): round-off error issues?");
            break;
        }

        ray.o = ray(its.t);
        remaining -= its.t;
        ray.maxt = remaining * lengthFactor;
        ray.mint = Epsilon;
    }

    return transmittance;
}

// ===========================================================================
//             Ray tracing support for bidirectional algorithms
// ===========================================================================

bool Scene::rayIntersectAll(const Ray &ray) const {
    if (rayIntersect(ray))
        return true;

    Float mint = ray.mint;
    if (mint == Epsilon)
        mint *= std::max(std::max(std::max(std::abs(ray.o.x),
            std::abs(ray.o.y)), std::abs(ray.o.z)), Epsilon);

    for (size_t i=0; i<m_specialShapes.size(); ++i) {
        if (m_specialShapes[i]->rayIntersect(ray, mint, ray.maxt))
            return true;
    }

    return false;
}

bool Scene::rayIntersectAll(const Ray &ray, Float &t,
            ConstShapePtr &shapePtr, Normal &n, Point2 &uv) const {
    bool result = rayIntersect(ray, t, shapePtr, n, uv);
    if (m_specialShapes.size() == 0)
        return result;

    uint8_t buffer[MTS_KD_INTERSECTION_TEMP];
    Float tempT, maxt = result ? t : ray.maxt;

    Float mint = ray.mint;
    if (mint == Epsilon)
        mint *= std::max(std::max(std::max(std::abs(ray.o.x),
            std::abs(ray.o.y)), std::abs(ray.o.z)), Epsilon);


    for (size_t i=0; i<m_specialShapes.size(); ++i) {
        const Shape *shape = m_specialShapes[i].get();

        if (shape->rayIntersect(ray, mint, maxt, tempT, buffer)) {
            /// Uh oh... -- much unnecessary work is done here
            Intersection its;
            its.t = tempT;
            shape->fillIntersectionRecord(ray, buffer, its);
            maxt = t = tempT;
            shapePtr = shape;
            n = its.geoFrame.n;
            uv = its.uv;
            result = true;
        }
    }

    return result;
}

bool Scene::rayIntersectAll(const Ray &ray, Intersection &its) const {
    bool result = rayIntersect(ray, its);
    if (m_specialShapes.size() == 0)
        return result;

    uint8_t buffer[MTS_KD_INTERSECTION_TEMP];
    Float maxt = result ? its.t : ray.maxt;
    Float mint = ray.mint;
    if (mint == Epsilon)
        mint *= std::max(std::max(std::max(std::abs(ray.o.x),
            std::abs(ray.o.y)), std::abs(ray.o.z)), Epsilon);
    Float tempT;

    for (size_t i=0; i<m_specialShapes.size(); ++i) {
        const Shape *shape = m_specialShapes[i].get();

        if (shape->rayIntersect(ray, mint, maxt, tempT, buffer)) {
            its.t = tempT;
            shape->fillIntersectionRecord(ray, buffer, its);
            result = true;
        }
    }

    return result;
}

Spectrum Scene::evalTransmittanceAll(const Point &p1, bool p1OnSurface, const Point &p2, bool p2OnSurface,
        Float time, const Medium *medium, int &interactions, Sampler *sampler) const {
    Vector d = p2 - p1;
    Float remaining = d.length();
    d /= remaining;

    Float lengthFactor = p2OnSurface ? (1-ShadowEpsilon) : 1;
    Ray ray(p1, d, p1OnSurface ? Epsilon : 0, remaining * lengthFactor, time);
    Spectrum transmittance(1.0f);
    Intersection its;
    int maxInteractions = interactions;
    interactions = 0;

    while (remaining > 0) {
        Normal n;
        bool surface = rayIntersectAll(ray, its.t, its.shape, its.geoFrame.n, its.uv);

        if (surface && (interactions == maxInteractions ||
            !(its.getBSDF()->getType() & BSDF::ENull))) {
            /* Encountered an occluder -- zero transmittance. */
            return Spectrum(0.0f);
        }

        if (medium)
            transmittance *= medium->evalTransmittance(
                Ray(ray, 0, std::min(its.t, remaining)), sampler);

        if (!surface || transmittance.isZero())
            break;

        const BSDF *bsdf = its.getBSDF();

        its.p = ray.o;
        its.geoFrame = Frame(its.geoFrame.n);
        its.hasUVPartials = false;
        Vector wo = its.geoFrame.toLocal(ray.d);
        BSDFSamplingRecord bRec(its, -wo, wo, ERadiance);
        bRec.typeMask = BSDF::ENull;
        transmittance *= bsdf->eval(bRec, EDiscrete);

        if (its.isMediumTransition()) {
            if (medium != its.getTargetMedium(-d)) {
                ++mediumInconsistencies;
                return Spectrum(0.0f);
            }
            medium = its.getTargetMedium(d);
        }

        if (++interactions > 100) { /// Just a precaution..
            Log(EWarn, "evalTransmittanceAll(): round-off error issues?");
            break;
        }

        ray.o = ray(its.t);
        remaining -= its.t;
        ray.maxt = remaining * lengthFactor;
        ray.mint = Epsilon;
    }

    return transmittance;
}

// ===========================================================================
//                Emission and direct illumination sampling
// ===========================================================================

Spectrum Scene::sampleEmitterDirect(DirectSamplingRecord &dRec,
        const Point2 &_sample, bool testVisibility) const {
    Point2 sample(_sample);

    /* Randomly pick an emitter */
    Float emPdf;
    size_t index = m_emitterPDF.sampleReuse(sample.x, emPdf);
    const Emitter *emitter = m_emitters[index].get();
    Spectrum value = emitter->sampleDirect(dRec, sample);

    if (dRec.pdf != 0) {
        if (testVisibility) {
            Ray ray(dRec.ref, dRec.d, Epsilon,
                    dRec.dist*(1-ShadowEpsilon), dRec.time);
            if (m_kdtree->rayIntersect(ray))
                return Spectrum(0.0f);
        }
        dRec.object = emitter;
        dRec.pdf *= emPdf;
        value /= emPdf;
        return value;
    } else {
        return Spectrum(0.0f);
    }
}

Spectrum Scene::sampleAttenuatedEmitterDirect(DirectSamplingRecord &dRec,
        const Medium *medium, int &interactions, const Point2 &_sample, Sampler *sampler) const {
    Point2 sample(_sample);

    /* Randomly pick an emitter */
    Float emPdf;
    size_t index = m_emitterPDF.sampleReuse(sample.x, emPdf);
    const Emitter *emitter = m_emitters[index].get();
    Spectrum value = emitter->sampleDirect(dRec, sample);

    if (dRec.pdf != 0) {
        value *= evalTransmittance(dRec.ref, false,
            dRec.p, emitter->isOnSurface(), dRec.time, medium,
            interactions, sampler) / emPdf;
        dRec.object = emitter;
        dRec.pdf *= emPdf;
        return value;
    } else {
        return Spectrum(0.0f);
    }
}

Spectrum Scene::sampleAttenuatedEmitterDirect(DirectSamplingRecord &dRec,
        const Intersection &its, const Medium *medium, int &interactions,
        const Point2 &_sample, Sampler *sampler) const {
    Point2 sample(_sample);

    /* Randomly pick an emitter */
    Float emPdf;
    size_t index = m_emitterPDF.sampleReuse(sample.x, emPdf);
    const Emitter *emitter = m_emitters[index].get();
    Spectrum value = emitter->sampleDirect(dRec, sample);

    if (dRec.pdf != 0) {
        if (its.shape && its.isMediumTransition())
            medium = its.getTargetMedium(dRec.d);
        value *= evalTransmittance(its.p, true, dRec.p, emitter->isOnSurface(),
                dRec.time, medium, interactions, sampler) / emPdf;
        dRec.object = emitter;
        dRec.pdf *= emPdf;
        return value;
    } else {
        return Spectrum(0.0f);
    }
}

Spectrum Scene::sampleSensorDirect(DirectSamplingRecord &dRec,
        const Point2 &sample, bool testVisibility) const {
    Spectrum value = m_sensor->sampleDirect(dRec, sample);

    if (dRec.pdf != 0) {
        if (testVisibility) {
            Ray ray(dRec.ref, dRec.d, Epsilon,
                    dRec.dist*(1-ShadowEpsilon), dRec.time);
            if (m_kdtree->rayIntersect(ray))
                return Spectrum(0.0f);
        }
        dRec.object = m_sensor.get();
        return value;
    } else {
        return Spectrum(0.0f);
    }
}

Spectrum Scene::sampleAttenuatedSensorDirect(DirectSamplingRecord &dRec,
        const Medium *medium, int &interactions, const Point2 &sample, Sampler *sampler) const {
    Spectrum value = m_sensor->sampleDirect(dRec, sample);

    if (dRec.pdf != 0) {
        value *= evalTransmittance(dRec.ref, false, dRec.p, m_sensor->isOnSurface(),
            dRec.time, medium, interactions, sampler);
        dRec.object = m_sensor.get();
        return value;
    } else {
        return Spectrum(0.0f);
    }
}

Spectrum Scene::sampleAttenuatedSensorDirect(DirectSamplingRecord &dRec,
        const Intersection &its, const Medium *medium, int &interactions,
        const Point2 &sample, Sampler *sampler) const {
    Spectrum value = m_sensor->sampleDirect(dRec, sample);

    if (dRec.pdf != 0) {
        if (its.shape && its.isMediumTransition())
            medium = its.getTargetMedium(dRec.d);
        value *= evalTransmittance(its.p, true, dRec.p, m_sensor->isOnSurface(),
                dRec.time, medium, interactions, sampler);
        dRec.object = m_sensor.get();
        return value;
    } else {
        return Spectrum(0.0f);
    }
}

Float Scene::pdfEmitterDirect(const DirectSamplingRecord &dRec) const {
    const Emitter *emitter = static_cast<const Emitter *>(dRec.object);
    return emitter->pdfDirect(dRec) * pdfEmitterDiscrete(emitter);
}

Float Scene::pdfSensorDirect(const DirectSamplingRecord &dRec) const {
    return m_sensor->pdfDirect(dRec);
}

Spectrum Scene::sampleEmitterPosition(
        PositionSamplingRecord &pRec,
        const Point2 &_sample) const {
    Point2 sample(_sample);

    /* Randomly pick an emitter */
    Float emPdf;
    size_t index = m_emitterPDF.sampleReuse(sample.x, emPdf);
    const Emitter *emitter = m_emitters[index].get();

    Spectrum value = emitter->samplePosition(pRec, sample);

    pRec.object = emitter;
    pRec.pdf *= emPdf;

    return value / emPdf;
}

Float Scene::pdfEmitterPosition(const PositionSamplingRecord &pRec) const {
    const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
    return emitter->pdfPosition(pRec) * pdfEmitterDiscrete(emitter);
}

Spectrum Scene::sampleEmitterRay(Ray &ray,
        const Emitter* &emitter,
        const Point2 &spatialSample,
        const Point2 &directionalSample,
        Float time) const {

    Point2 sample(spatialSample);

    /* Randomly pick an emitter */
    Float emPdf;
    size_t index = m_emitterPDF.sampleReuse(sample.x, emPdf);
    emitter = m_emitters[index].get();

    return emitter->sampleRay(ray, sample, directionalSample, time) / emPdf;
}

MTS_IMPLEMENT_CLASS_S(Scene, false, ConfigurableObject)
MTS_NAMESPACE_END

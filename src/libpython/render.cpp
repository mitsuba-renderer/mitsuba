#include "base.h"
#include <mitsuba/render/scene.h>
#include <mitsuba/render/scenehandler.h>
#include <mitsuba/render/renderqueue.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/render/noise.h>
#include "../shapes/instance.h"

using namespace mitsuba;

static bool intersection_get_hasUVPartials(const Intersection &its) { return its.hasUVPartials;  }
static void intersection_set_hasUVPartials(Intersection &its, bool value) { its.hasUVPartials = value;  }
static uint32_t intersection_get_primIndex(const Intersection &its) { return its.primIndex;  }
static void intersection_set_primIndex(Intersection &its, uint32_t value) { its.primIndex = value;  }
static Shape * intersection_get_shape(const Intersection &its) { return const_cast<Shape *>(its.shape);  }
static void intersection_set_shape(Intersection &its, Shape * value) { its.shape = value;  }
static Shape * intersection_get_instance(const Intersection &its) { return const_cast<Shape *>(its.instance);  }
static void intersection_set_instance(Intersection &its, Shape * value) { its.instance = value;  }

static const Intersection &bsdfsamplingrecord_get_its(const BSDFSamplingRecord &bRec) {
    return bRec.its;
}

static bp::tuple texture_evalGradient(const Texture *texture, const Intersection &its) {
    Spectrum grad[2];
    texture->evalGradient(its, grad);
    return bp::make_tuple(grad[0], grad[1]);
}

static unsigned int bsdf_getType_1(const BSDF *bsdf) {
    return bsdf->getType();
}

static unsigned int bsdf_getType_2(const BSDF *bsdf, int index) {
    return bsdf->getType(index);
}

static bp::tuple bsdf_sample(const BSDF *bsdf, BSDFSamplingRecord &bRec, const Point2 &sample) {
    Float pdf;
    Spectrum result = bsdf->sample(bRec, pdf, sample);
    return bp::make_tuple(result, pdf);
}

static bp::list shapekdtree_getShapes(const ShapeKDTree *kdtree) {
    const std::vector<const Shape *> &shapes = kdtree->getShapes();
    bp::list list;
    for (size_t i=0; i<shapes.size(); ++i)
        list.append(bp::object(ref<Shape>(const_cast<Shape *>(shapes[i]))));
    return list;
}

static Shape *shape_getShapeGroup(Shape *shape) {
    if (shape->getClass()->getName() != "Instance")
        return NULL;
    return static_cast<Instance *>(shape)->getShapeGroup();
}

static bp::object shape_getElement(Shape *shape, int idx) {
    return cast(shape->getElement(idx));
}

static AABB shapekdtree_getAABB(const ShapeKDTree *kdtree) {
    return kdtree->getAABB();
}

static bool shapekdtree_isBuilt(const ShapeKDTree *kdtree) {
    return kdtree->isBuilt();
}

static bp::object shapekdtree_rayIntersect(const ShapeKDTree *kdtree, const Ray &ray) {
    Intersection its;

    if (!kdtree->rayIntersect(ray, its))
        return bp::object();

    return bp::object(its);
}


static ShapeKDTree* shape_getKDTree(const Shape *shape) {
    return const_cast<ShapeKDTree *>(static_cast<const ShapeKDTree *>(shape->getKDTree()));
}

static bp::object shape_rayIntersect(const Shape *shape, const Ray &ray, Float mint, Float maxt) {
    uint8_t temp[MTS_KD_INTERSECTION_TEMP];
    Float t;

    if (!shape->rayIntersect(ray, mint, maxt, t, temp))
        return bp::object();

    Intersection its;
    its.shape = shape;
    its.t = t;
    shape->fillIntersectionRecord(ray, temp, its);

    return bp::object(its);
}

static bp::object scene_rayIntersect(const Scene *scene, const Ray &ray) {
    Intersection its;

    if (!scene->rayIntersect(ray, its))
        return bp::object();

    return bp::object(its);
}

static bp::object scene_rayIntersectAll(const Scene *scene, const Ray &ray) {
    Intersection its;

    if (!scene->rayIntersectAll(ray, its))
        return bp::object();

    return bp::object(its);
}

static bp::tuple shape_getCurvature(const Shape *shape, const Intersection &its, bool shadingFrame) {
    Float H, K;
    shape->getCurvature(its, H, K, shadingFrame);
    return bp::make_tuple(H, K);
}

static bp::tuple shape_getNormalDerivative(const Shape *shape, const Intersection &its, bool shadingFrame) {
    Vector dpdu, dpdv;
    shape->getNormalDerivative(its, dpdu, dpdv, shadingFrame);
    return bp::make_tuple(dpdu, dpdv);
}

static ref<Scene> loadScene1(const fs::path &filename) {
    return SceneHandler::loadScene(filename);
}

static ref<Scene> loadScene2(const fs::path &filename, const StringMap &params) {
    SceneHandler::ParameterMap pmap;
    for (StringMap::const_iterator it = params.begin(); it != params.end(); ++it)
        pmap[it->first]=it->second;
    return SceneHandler::loadScene(filename, pmap);
}

static bp::list scene_getSensors(Scene *scene) {
    bp::list list;
    ref_vector<Sensor> &sensors = scene->getSensors();
    for (size_t i=0; i<sensors.size(); ++i)
        list.append(cast(sensors[i].get()));
    return list;
}

static bp::object scene_getSensor(Scene *scene) { return cast(scene->getSensor()); }
static bp::object scene_getIntegrator(Scene *scene) { return cast(scene->getIntegrator()); }

static bp::list scene_getMeshes(Scene *scene) {
    bp::list list;
    std::vector<TriMesh *> &meshes = scene->getMeshes();
    for (size_t i=0; i<meshes.size(); ++i)
        list.append(cast(meshes[i]));
    return list;
}

static bp::list scene_getShapes(Scene *scene) {
    bp::list list;
    ref_vector<Shape> &shapes = scene->getShapes();
    for (size_t i=0; i<shapes.size(); ++i)
        list.append(cast(shapes[i].get()));
    return list;
}

static bp::list scene_getEmitters(Scene *scene) {
    bp::list list;
    ref_vector<Emitter> &emitters = scene->getEmitters();
    for (size_t i=0; i<emitters.size(); ++i)
        list.append(cast(emitters[i].get()));
    return list;
}

static bp::list scene_getMedia(Scene *scene) {
    bp::list list;
    ref_vector<Medium> &media = scene->getMedia();
    for (size_t i=0; i<media.size(); ++i)
        list.append(cast(media[i].get()));
    return list;
}

typedef InternalArray<uint32_t>     InternalUInt32Array;
typedef InternalArray<Point3>       InternalPoint3Array;
typedef InternalArray<Normal>       InternalNormalArray;
typedef InternalArray<Point2>       InternalPoint2Array;
typedef InternalArray<Color3>       InternalColor3Array;
typedef InternalArray<TangentSpace> InternalTangentSpaceArray;

static InternalUInt32Array trimesh_getTriangles(TriMesh *triMesh) {
    BOOST_STATIC_ASSERT(sizeof(Triangle) == 3*sizeof(uint32_t));
    return InternalUInt32Array(triMesh, (uint32_t *) triMesh->getTriangles(), triMesh->getTriangleCount()*3);
}

static InternalPoint3Array trimesh_getVertexPositions(TriMesh *triMesh) {
    return InternalPoint3Array(triMesh, triMesh->getVertexPositions(), triMesh->getVertexCount());
}

static InternalNormalArray trimesh_getVertexNormals(TriMesh *triMesh) {
    return InternalNormalArray(triMesh, triMesh->getVertexNormals(), triMesh->getVertexCount());
}

static InternalPoint2Array trimesh_getVertexTexcoords(TriMesh *triMesh) {
    return InternalPoint2Array(triMesh, triMesh->getVertexTexcoords(), triMesh->getVertexCount());
}

static InternalColor3Array trimesh_getVertexColors(TriMesh *triMesh) {
    return InternalColor3Array(triMesh, triMesh->getVertexColors(), triMesh->getVertexCount());
}

static InternalTangentSpaceArray trimesh_getUVTangents(TriMesh *triMesh) {
    return InternalTangentSpaceArray(triMesh, triMesh->getUVTangents(), triMesh->getVertexCount());
}

static ref<TriMesh> trimesh_fromBlender(const std::string &name,
        size_t faceCount, size_t facePtr, size_t vertexCount, size_t vertexPtr, size_t uvPtr, size_t colPtr, short matID) {
    return TriMesh::fromBlender(name, faceCount, reinterpret_cast<void *>(facePtr), vertexCount,
        reinterpret_cast<void *>(vertexPtr), reinterpret_cast<void *>(uvPtr),
        reinterpret_cast<void *>(colPtr), matID);
}

class RenderListenerWrapper : public RenderListener {
public:
    RenderListenerWrapper(PyObject *self) : m_self(self), m_locked(false) { Py_INCREF(m_self); }

    void workBeginEvent(const RenderJob *job, const RectangularWorkUnit *wu, int worker) {
        CALLBACK_SYNC_GIL();
        try {
            bp::call_method<void>(m_self, "workBeginEvent", bp::ptr(job), bp::ptr(wu), worker);
        } catch (bp::error_already_set &) { check_python_exception(); }
    }

    void workEndEvent(const RenderJob *job, const ImageBlock *wr, bool cancelled) {
        CALLBACK_SYNC_GIL();
        try {
            bp::call_method<void>(m_self, "workEndEvent", bp::ptr(job), bp::ptr(wr), cancelled);
        } catch (bp::error_already_set &) { check_python_exception(); }
    }

    void workCanceledEvent(const RenderJob *job, const Point2i &offset,
            const Vector2i &size) {
        CALLBACK_SYNC_GIL();
        try {
            bp::call_method<void>(m_self, "workCanceledEvent", bp::ptr(job), offset, size);
        } catch (bp::error_already_set &) { check_python_exception(); }
    }

    void refreshEvent(const RenderJob *job) {
        CALLBACK_SYNC_GIL();
        try {
            bp::call_method<void>(m_self, "refreshEvent", bp::ptr(job));
        } catch (bp::error_already_set &) { check_python_exception(); }
    }

    void finishJobEvent(const RenderJob *job, bool cancelled) {
        CALLBACK_SYNC_GIL();
        try {
            bp::call_method<void>(m_self, "finishJobEvent", bp::ptr(job), cancelled);
        } catch (bp::error_already_set &) { check_python_exception(); }
    }

    virtual ~RenderListenerWrapper() {
        Py_DECREF(m_self);
    }
private:
    PyObject *m_self;
    bool m_locked;
};

static void renderQueue_join(RenderQueue *queue) {
    ReleaseGIL gil;
    queue->join();
}

static void renderQueue_waitLeft(RenderQueue *queue, size_t count) {
    ReleaseGIL gil;
    queue->waitLeft(count);
}

static void scene_cancel(Scene *scene) {
    ReleaseGIL gil;
    scene->cancel();
}

static void renderJob_cancel(RenderJob *job) {
    ReleaseGIL gil;
    job->cancel();
}

bp::tuple Sensor_sampleRay(Sensor *sensor, const Point2 &samplePosition, const Point2 &apertureSample, Float timeSample) {
    Ray ray;
    Spectrum result = sensor->sampleRay(ray, samplePosition, apertureSample, timeSample);
    return bp::make_tuple(result, ray);
}

bp::tuple Sensor_sampleRayDifferential(Sensor *sensor, const Point2 &samplePosition, const Point2 &apertureSample, Float timeSample) {
    RayDifferential ray;
    Spectrum result = sensor->sampleRay(ray, samplePosition, apertureSample, timeSample);
    return bp::make_tuple(result, ray);
}

bp::tuple Sensor_eval(Sensor *sensor, const Intersection &its, const Vector &d) {
    Point2 samplePos;
    Spectrum result = sensor->eval(its, d, samplePos);
    return bp::make_tuple(result, samplePos);
}

bp::tuple Sensor_getSamplePosition(Sensor *sensor, const PositionSamplingRecord &pRec, const DirectionSamplingRecord &dRec) {
    Point2 samplePos;
    bool result = sensor->getSamplePosition(pRec, dRec, samplePos);
    return bp::make_tuple(result, samplePos);
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getBitmap_overloads, getBitmap, 0, 1)

void export_render() {
    bp::object renderModule(
        bp::handle<>(bp::borrowed(PyImport_AddModule("mitsuba.render"))));
    bp::scope().attr("render") = renderModule;
    PyObject *oldScope = bp::detail::current_scope;

    BP_SETSCOPE(renderModule);
    renderModule.attr("__path__") = "mitsuba.render";

    BP_INTERNAL_ARRAY(InternalUInt32Array);
    BP_INTERNAL_ARRAY(InternalPoint3Array);
    BP_INTERNAL_ARRAY(InternalPoint2Array);
    BP_INTERNAL_ARRAY(InternalColor3Array);
    BP_INTERNAL_ARRAY(InternalNormalArray);
    BP_INTERNAL_ARRAY(InternalTangentSpaceArray);

    bp::enum_<ETransportMode>("ETransportMode")
        .value("ERadiance", ERadiance)
        .value("EImportance", EImportance)
        .export_values();

    bp::enum_<EMeasure>("EMeasure")
        .value("EInvalidMeasure", EInvalidMeasure)
        .value("ESolidAngle", ESolidAngle)
        .value("ELength", ELength)
        .value("EArea", EArea)
        .value("EDiscrete", EDiscrete)
        .export_values();

    BP_CLASS(ShapeKDTree, Object, bp::no_init)
        .def("getPrimitiveCount", &ShapeKDTree::getPrimitiveCount)
        .def("addShape", &ShapeKDTree::addShape)
        .def("build", &ShapeKDTree::build)
        .def("isBuilt", &shapekdtree_isBuilt)
        .def("getAABB", &shapekdtree_getAABB, BP_RETURN_VALUE)
        .def("getShapes", &shapekdtree_getShapes)
        .def("rayIntersect", &shapekdtree_rayIntersect);

    Sampler *(Scene::*scene_getSampler)(void) = &Scene::getSampler;
    Film *(Scene::*scene_getFilm)(void) = &Scene::getFilm;
    ShapeKDTree *(Scene::*scene_getKDTree)(void) = &Scene::getKDTree;

    BP_CLASS(Scene, NetworkedObject, bp::init<>())
        .def(bp::init<Properties>())
        .def(bp::init<Scene *>())
        .def(bp::init<Stream *, InstanceManager *>())
        .def("initialize", &Scene::initialize)
        .def("invalidate", &Scene::invalidate)
        .def("preprocess", &Scene::preprocess)
        .def("render", &Scene::render)
        .def("postprocess", &Scene::postprocess)
        .def("flush", &Scene::flush)
        .def("cancel", scene_cancel)
        .def("rayIntersect", &scene_rayIntersect)
        .def("rayIntersectAll", &scene_rayIntersectAll)
        .def("evalTransmittance", &Scene::evalTransmittance)
        .def("evalTransmittanceAll", &Scene::evalTransmittanceAll)
        .def("sampleEmitterDirect", &Scene::sampleEmitterDirect)
        .def("sampleSensorDirect", &Scene::sampleSensorDirect)
        .def("pdfEmitterDirect", &Scene::pdfEmitterDirect)
        .def("pdfSensorDirect", &Scene::pdfSensorDirect)
        .def("getAABB", &Scene::getAABB, BP_RETURN_VALUE)
        .def("getBSphere", &Scene::getBSphere, BP_RETURN_VALUE)
        .def("getBlockSize", &Scene::getBlockSize)
        .def("setBlockSize", &Scene::setBlockSize)
        .def("getSourceFile", &Scene::getSourceFile, BP_RETURN_VALUE)
        .def("setSourceFile", &Scene::setSourceFile)
        .def("getDestinationFile", &Scene::getDestinationFile, BP_RETURN_VALUE)
        .def("setDestinationFile", &Scene::setDestinationFile)
        .def("destinationExists", &Scene::destinationExists)
        .def("hasEnvironmentEmitter", &Scene::hasEnvironmentEmitter)
        .def("getEnvironmentEmitter", &Scene::getEnvironmentEmitter, BP_RETURN_VALUE)
        .def("hasDegenerateSensor", &Scene::hasDegenerateSensor)
        .def("hasDegenerateEmitters", &Scene::hasDegenerateEmitters)
        .def("hasMedia", &Scene::hasMedia)
        .def("addSensor", &Scene::addSensor)
        .def("removeSensor", &Scene::removeSensor)
        .def("getSensor", &scene_getSensor, BP_RETURN_VALUE)
        .def("setSensor", &Scene::setSensor)
        .def("getSensors", &scene_getSensors)
        .def("getIntegrator", &scene_getIntegrator, BP_RETURN_VALUE)
        .def("setIntegrator", &Scene::setIntegrator)
        .def("getSampler", scene_getSampler, BP_RETURN_VALUE)
        .def("setSampler", &Scene::setSampler)
        .def("getFilm", scene_getFilm, BP_RETURN_VALUE)
        .def("getShapes", &scene_getShapes)
        .def("getMeshes", &scene_getMeshes)
        .def("getEmitters", &scene_getEmitters)
        .def("getMedia", &scene_getMedia)
        .def("getKDTree", scene_getKDTree, BP_RETURN_VALUE);

    BP_CLASS(Sampler, ConfigurableObject, bp::no_init)
        .def("clone", &Sampler::clone, BP_RETURN_VALUE)
        .def("generate", &Sampler::generate)
        .def("advance", &Sampler::advance)
        .def("setSampleIndex", &Sampler::setSampleIndex)
        .def("getSampleIndex", &Sampler::getSampleIndex)
        .def("getSampleCount", &Sampler::getSampleCount)
        .def("getProperties", &Sampler::getProperties, BP_RETURN_VALUE)
        .def("next1D", &Sampler::next1D)
        .def("next2D", &Sampler::next2D, BP_RETURN_VALUE)
        .def("request1DArray", &Sampler::request1DArray)
        .def("request2DArray", &Sampler::request2DArray);

    bp::class_<SceneHandler, boost::noncopyable>("SceneHandler", bp::no_init)
        .def("loadScene", &loadScene1, BP_RETURN_VALUE)
        .def("loadScene", &loadScene2, BP_RETURN_VALUE)
        .staticmethod("loadScene");

    Scene *(RenderJob::*renderJob_getScene)(void) = &RenderJob::getScene;
    RenderQueue *(RenderJob::*renderJob_getRenderQueue)(void) = &RenderJob::getRenderQueue;
    BP_CLASS(RenderJob, Thread, (bp::init<const std::string &, Scene *, RenderQueue *, bp::optional<int, int, int, bool, bool> >()))
        .def("flush", &RenderJob::flush)
        .def("cancel", renderJob_cancel)
        .def("wait", &RenderJob::wait)
        .def("isInteractive", &RenderJob::isInteractive)
        .def("setInteractive", &RenderJob::setInteractive)
        .def("getScene", renderJob_getScene, BP_RETURN_VALUE)
        .def("getRenderQueue", renderJob_getRenderQueue, BP_RETURN_VALUE);

    BP_CLASS(RenderQueue, Object, bp::init<>())
        .def("getJobCount", &RenderQueue::getJobCount)
        .def("addJob", &RenderQueue::addJob)
        .def("removeJob", &RenderQueue::removeJob)
        .def("getRenderTime", &RenderQueue::getRenderTime)
        .def("registerListener", &RenderQueue::registerListener)
        .def("unregisterListener", &RenderQueue::unregisterListener)
        .def("waitLeft", renderQueue_waitLeft)
        .def("join", renderQueue_join)
        .def("flush", &RenderQueue::flush);

    BP_STRUCT(Intersection, bp::init<>())
        .def_readwrite("t", &Intersection::t)
        .def_readwrite("p", &Intersection::p)
        .def_readwrite("geoFrame", &Intersection::geoFrame)
        .def_readwrite("shFrame", &Intersection::shFrame)
        .def_readwrite("uv", &Intersection::uv)
        .def_readwrite("dpdu", &Intersection::dpdu)
        .def_readwrite("dpdv", &Intersection::dpdv)
        .def_readwrite("dudx", &Intersection::dudx)
        .def_readwrite("dudy", &Intersection::dudy)
        .def_readwrite("dvdx", &Intersection::dvdx)
        .def_readwrite("dvdy", &Intersection::dvdy)
        .def_readwrite("time", &Intersection::time)
        .def_readwrite("color", &Intersection::color)
        .def_readwrite("wi", &Intersection::wi)
        .add_property("shape", bp::make_function(&intersection_get_shape, BP_RETURN_VALUE), &intersection_set_shape)
        .add_property("instance", bp::make_function(&intersection_get_instance, BP_RETURN_VALUE), &intersection_set_instance)
        .add_property("hasUVPartials", &intersection_get_hasUVPartials, &intersection_set_hasUVPartials)
        .add_property("primIndex", &intersection_get_primIndex, &intersection_set_primIndex)
        .def("toWorld", &Intersection::toWorld, BP_RETURN_VALUE)
        .def("toLocal", &Intersection::toLocal, BP_RETURN_VALUE)
        .def("isEmitter", &Intersection::isEmitter)
        .def("isSensor", &Intersection::isSensor)
        .def("hasSubsurface", &Intersection::hasSubsurface)
        .def("isMediumTransition", &Intersection::isMediumTransition)
        .def("isValid", &Intersection::isValid)
        .def("Le", &Intersection::Le)
        .def("LoSub", &Intersection::LoSub)
        .def("__repr__", &Intersection::toString);

    BP_STRUCT(TangentSpace, bp::init<>())
        .def(bp::init<Vector, Vector>())
        .def(bp::init<Stream *>())
        .def_readwrite("dpdu", &TangentSpace::dpdu)
        .def_readwrite("dpdv", &TangentSpace::dpdv)
        .def("serialize", &TangentSpace::serialize)
        .def("__repr__", &TangentSpace::toString);

    BP_STRUCT(PositionSamplingRecord, bp::init<>())
        .def(bp::init<Float>())
        .def(bp::init<Intersection, EMeasure>())
        .def_readwrite("p", &PositionSamplingRecord::p)
        .def_readwrite("time", &PositionSamplingRecord::time)
        .def_readwrite("n", &PositionSamplingRecord::n)
        .def_readwrite("pdf", &PositionSamplingRecord::pdf)
        .def_readwrite("measure", &PositionSamplingRecord::measure)
        .def_readwrite("uv", &PositionSamplingRecord::uv)
        .def_readwrite("object", &PositionSamplingRecord::object)
        .def("__repr__", &PositionSamplingRecord::toString);

    BP_STRUCT(DirectionSamplingRecord, bp::init<>())
        .def(bp::init<Vector, EMeasure>())
        .def(bp::init<Intersection, EMeasure>())
        .def_readwrite("d", &DirectionSamplingRecord::d)
        .def_readwrite("pdf", &DirectionSamplingRecord::pdf)
        .def_readwrite("measure", &DirectionSamplingRecord::measure)
        .def("__repr__", &DirectionSamplingRecord::toString);

    BP_SUBSTRUCT(DirectSamplingRecord, PositionSamplingRecord, bp::init<>())
        .def(bp::init<Point, Float>())
        .def(bp::init<Intersection>())
        .def(bp::init<MediumSamplingRecord>())
        .def_readwrite("ref", &DirectSamplingRecord::ref)
        .def_readwrite("refN", &DirectSamplingRecord::refN)
        .def_readwrite("d", &DirectSamplingRecord::d)
        .def_readwrite("dist", &DirectSamplingRecord::dist)
        .def("setQuery", &DirectSamplingRecord::setQuery)
        .def("__repr__", &DirectSamplingRecord::toString);

    Medium *(Shape::*shape_getInteriorMedium)(void) = &Shape::getInteriorMedium;
    Medium *(Shape::*shape_getExteriorMedium)(void) = &Shape::getExteriorMedium;
    Sensor *(Shape::*shape_getSensor)(void) = &Shape::getSensor;
    Emitter *(Shape::*shape_getEmitter)(void) = &Shape::getEmitter;
    Subsurface *(Shape::*shape_getSubsurface)(void) = &Shape::getSubsurface;
    BSDF *(Shape::*shape_getBSDF)(void) = &Shape::getBSDF;

    BP_CLASS(Shape, ConfigurableObject, bp::no_init)
        .def("getName", &Shape::getName, BP_RETURN_VALUE)
        .def("isCompound", &Shape::isCompound)
        .def("getElement", shape_getElement)
        .def("getSurfaceArea", &Shape::getSurfaceArea)
        .def("getAABB", &Shape::getAABB, BP_RETURN_VALUE)
        .def("getClippedAABB", &Shape::getClippedAABB, BP_RETURN_VALUE)
        .def("createTriMesh", &Shape::createTriMesh, BP_RETURN_VALUE)
        .def("rayIntersect", &shape_rayIntersect)
        .def("getNormalDerivative", &shape_getNormalDerivative)
        .def("getCurvature", &shape_getCurvature)
        .def("samplePosition", &Shape::samplePosition)
        .def("pdfPosition", &Shape::pdfPosition)
        .def("sampleDirect", &Shape::sampleDirect)
        .def("pdfDirect", &Shape::pdfDirect)
        .def("getInteriorMedium", shape_getInteriorMedium, BP_RETURN_VALUE)
        .def("getExteriorMedium", shape_getExteriorMedium, BP_RETURN_VALUE)
        .def("isMediumTransition", &Shape::isMediumTransition)
        .def("hasSubsurface", &Shape::hasSubsurface)
        .def("getSubsurface", shape_getSubsurface, BP_RETURN_VALUE)
        .def("isEmitter", &Shape::isEmitter)
        .def("getEmitter", shape_getEmitter, BP_RETURN_VALUE)
        .def("isSensor", &Shape::isSensor)
        .def("getSensor", shape_getSensor, BP_RETURN_VALUE)
        .def("hasBSDF", &Shape::hasBSDF)
        .def("getBSDF", shape_getBSDF, BP_RETURN_VALUE)
        .def("setBSDF", &Shape::setBSDF)
        .def("getKDTree", &shape_getKDTree, BP_RETURN_VALUE)
        .def("getPrimitiveCount", &Shape::getPrimitiveCount)
        .def("getEffectivePrimitiveCount", &Shape::getEffectivePrimitiveCount)
        .def("copyAttachments", &Shape::copyAttachments)
        .def("getShapeGroup", &shape_getShapeGroup, BP_RETURN_VALUE);

    void (TriMesh::*triMesh_serialize1)(Stream *stream) const = &TriMesh::serialize;
    void (TriMesh::*triMesh_serialize2)(Stream *stream, InstanceManager *) const = &TriMesh::serialize;

    BP_CLASS(TriMesh, Shape, (bp::init<std::string, size_t, size_t, bp::optional<bool, bool, bool, bool, bool> >()))
        .def(bp::init<Stream *, InstanceManager *>())
        .def(bp::init<Stream *, int>())
        .def("getTriangleCount", &TriMesh::getTriangleCount)
        .def("getTriangles", trimesh_getTriangles)
        .def("getVertexCount", &TriMesh::getVertexCount)
        .def("getVertexPositions", trimesh_getVertexPositions, BP_RETURN_VALUE)
        .def("hasVertexNormals", &TriMesh::hasVertexNormals)
        .def("getVertexNormals", trimesh_getVertexNormals, BP_RETURN_VALUE)
        .def("hasVertexColors", &TriMesh::hasVertexColors)
        .def("getVertexColors", trimesh_getVertexColors, BP_RETURN_VALUE)
        .def("hasVertexTexcoords", &TriMesh::hasVertexTexcoords)
        .def("getVertexTexcoords", trimesh_getVertexTexcoords, BP_RETURN_VALUE)
        .def("hasUVTangents", &TriMesh::hasUVTangents)
        .def("getUVTangents", trimesh_getUVTangents, BP_RETURN_VALUE)
        .def("computeUVTangents", &TriMesh::computeUVTangents)
        .def("computeNormals", &TriMesh::computeNormals)
        .def("rebuildTopology", &TriMesh::rebuildTopology)
        .def("serialize", triMesh_serialize1)
        .def("serialize", triMesh_serialize2)
        .def("writeOBJ", &TriMesh::writeOBJ)
        .def("writePLY", &TriMesh::writePLY)
        .def("fromBlender", trimesh_fromBlender)
        .staticmethod("fromBlender");

    Shape *(AbstractEmitter::*abstractemitter_getShape)(void) = &AbstractEmitter::getShape;
    Medium *(AbstractEmitter::*abstractemitter_getMedium)(void) = &AbstractEmitter::getMedium;

    BP_CLASS(AbstractEmitter, ConfigurableObject, bp::no_init)
        .def("getType", &AbstractEmitter::getType)
        .def("setWorldTransform", &AbstractEmitter::setWorldTransform)
        .def("getWorldTransform", &AbstractEmitter::getWorldTransform, BP_RETURN_VALUE)
        .def("isOnSurface", &AbstractEmitter::isOnSurface)
        .def("needsPositionSample", &AbstractEmitter::needsPositionSample)
        .def("needsDirectionSample", &AbstractEmitter::needsDirectionSample)
        .def("needsDirectSample", &AbstractEmitter::needsDirectSample)
        .def("getDirectMeasure", &AbstractEmitter::getDirectMeasure)
        .def("isDegenerate", &AbstractEmitter::isDegenerate)
        .def("getShape", abstractemitter_getShape, BP_RETURN_VALUE)
        .def("getMedium", abstractemitter_getMedium, BP_RETURN_VALUE)
        .def("createShape", &AbstractEmitter::createShape, BP_RETURN_VALUE)
        .def("getAABB", &AbstractEmitter::getAABB, BP_RETURN_VALUE);

    BP_SETSCOPE(AbstractEmitter_class);
    bp::enum_<AbstractEmitter::EEmitterType>("EEmitterType")
        .value("EDeltaDirection", AbstractEmitter::EDeltaDirection)
        .value("EDeltaPosition", AbstractEmitter::EDeltaPosition)
        .value("EOnSurface", AbstractEmitter::EOnSurface)
        .export_values();
    BP_SETSCOPE(renderModule);

    BP_CLASS(Emitter, AbstractEmitter, bp::no_init) // incomplete
        .def("eval", &Emitter::eval, BP_RETURN_VALUE)
        .def("getSamplingWeight", &Emitter::getSamplingWeight)
        .def("isEnvironmentEmitter", &Emitter::isEnvironmentEmitter)
        .def("evalEnvironment", &Emitter::evalEnvironment, BP_RETURN_VALUE)
        .def("isCompound", &Emitter::isCompound)
        .def("getElement", &Emitter::getElement, BP_RETURN_VALUE)
        .def("getBitmap", &Emitter::getBitmap, getBitmap_overloads()[BP_RETURN_VALUE]);

    BP_SETSCOPE(Emitter_class);
    bp::enum_<Emitter::EEmitterFlags>("EEmitterFlags")
        .value("EEnvironmentEmitter", Emitter::EEnvironmentEmitter)
        .export_values();
    BP_SETSCOPE(renderModule);

    Film *(Sensor::*sensor_getFilm)(void) = &Sensor::getFilm;
    Sampler *(Sensor::*sensor_getSampler)(void) = &Sensor::getSampler;

    BP_CLASS(Sensor, AbstractEmitter, bp::no_init)
        .def("getShutterOpen", &Sensor::getShutterOpen)
        .def("setShutterOpen", &Sensor::setShutterOpen)
        .def("getShutterOpenTime", &Sensor::getShutterOpenTime)
        .def("setShutterOpenTime", &Sensor::setShutterOpenTime)
        .def("sampleTime", &Sensor::sampleTime)
        .def("sampleRay", &Sensor_sampleRay)
        .def("sampleRayDifferential", &Sensor_sampleRayDifferential)
        .def("eval", &Sensor_eval, BP_RETURN_VALUE)
        .def("getSamplePosition", &Sensor_getSamplePosition)
        .def("pdfTime", &Sensor::pdfTime)
        .def("needsTimeSample", &Sensor::needsTimeSample)
        .def("needsApertureSample", &Sensor::needsApertureSample)
        .def("getFilm", sensor_getFilm, BP_RETURN_VALUE)
        .def("getSampler", sensor_getSampler, BP_RETURN_VALUE)
        .def("getAspect", &Sensor::getAspect);

    BP_SETSCOPE(Sensor_class);
    bp::enum_<Sensor::ESensorFlags>("ESensorFlags")
        .value("EDeltaTime", Sensor::EDeltaTime)
        .value("ENeedsApertureSample", Sensor::ENeedsApertureSample)
        .value("EProjectiveCamera", Sensor::EProjectiveCamera)
        .value("EPerspectiveCamera", Sensor::EPerspectiveCamera)
        .value("EOrthographicCamera", Sensor::EOrthographicCamera)
        .value("EPositionSampleMapsToPixels", Sensor::EPositionSampleMapsToPixels)
        .value("EDirectionSampleMapsToPixels", Sensor::EDirectionSampleMapsToPixels)
        .export_values();
    BP_SETSCOPE(renderModule);

    void (Film::*film_develop1)(const Scene *scene, Float renderTime) = &Film::develop;
    bool (Film::*film_develop2)(const Point2i &offset, const Vector2i &size,
        const Point2i &targetOffset, Bitmap *target) const = &Film::develop;
    ReconstructionFilter *(Film::*film_getreconstructionfilter)() = &Film::getReconstructionFilter;

    BP_CLASS(Film, ConfigurableObject, bp::no_init)
        .def("getSize", &Film::getSize, BP_RETURN_VALUE)
        .def("getCropSize", &Film::getCropSize, BP_RETURN_VALUE)
        .def("getCropOffset", &Film::getCropOffset, BP_RETURN_VALUE)
        .def("clear", &Film::clear)
        .def("setBitmap", &Film::setBitmap)
        .def("addBitmap", &Film::addBitmap)
        .def("setDestinationFile", &Film::setDestinationFile)
        .def("develop", film_develop1)
        .def("develop", film_develop2)
        .def("destinationExists", &Film::destinationExists)
        .def("hasHighQualityEdges", &Film::hasHighQualityEdges)
        .def("hasAlpha", &Film::hasAlpha)
        .def("getReconstructionFilter", film_getreconstructionfilter, BP_RETURN_VALUE);

    void (ProjectiveCamera::*projectiveCamera_setWorldTransform1)(const Transform &) = &ProjectiveCamera::setWorldTransform;
    void (ProjectiveCamera::*projectiveCamera_setWorldTransform2)(AnimatedTransform *) = &ProjectiveCamera::setWorldTransform;
    const Transform (ProjectiveCamera::*projectiveCamera_getWorldTransform1)(Float t) const = &ProjectiveCamera::getWorldTransform;
    const AnimatedTransform *(ProjectiveCamera::*projectiveCamera_getWorldTransform2)(void) const = &ProjectiveCamera::getWorldTransform;

    BP_CLASS(ProjectiveCamera, Sensor, bp::no_init)
        .def("getViewTransform", &ProjectiveCamera::getViewTransform, BP_RETURN_VALUE)
        .def("getWorldTransform", projectiveCamera_getWorldTransform1, BP_RETURN_VALUE)
        .def("getWorldTransform", projectiveCamera_getWorldTransform2, BP_RETURN_VALUE)
        .def("setWorldTransform", projectiveCamera_setWorldTransform1)
        .def("setWorldTransform", projectiveCamera_setWorldTransform2)
        .def("getProjectionTransform", &ProjectiveCamera::getProjectionTransform, BP_RETURN_VALUE)
        .def("getNearClip", &ProjectiveCamera::getNearClip)
        .def("getFarClip", &ProjectiveCamera::getFarClip)
        .def("getFocusDistance", &ProjectiveCamera::getFocusDistance)
        .def("setNearClip", &ProjectiveCamera::setNearClip)
        .def("setFarClip", &ProjectiveCamera::setFarClip)
        .def("setFocusDistance", &ProjectiveCamera::setFocusDistance);

    BP_CLASS(PerspectiveCamera, ProjectiveCamera, bp::no_init)
        .def("getXFov", &PerspectiveCamera::getXFov)
        .def("setXFov", &PerspectiveCamera::setXFov)
        .def("getYFov", &PerspectiveCamera::getYFov)
        .def("setYFov", &PerspectiveCamera::setYFov)
        .def("getDiagonalFov", &PerspectiveCamera::getDiagonalFov)
        .def("setDiagonalFov", &PerspectiveCamera::setDiagonalFov);

    BP_CLASS(Integrator, ConfigurableObject, bp::no_init)
        .def("preprocess", &Integrator::preprocess)
        .def("render", &Integrator::render)
        .def("cancel", &Integrator::cancel)
        .def("postprocess", &Integrator::postprocess)
        .def("configureSampler", &Integrator::configureSampler)
        .def("getSubIntegrator", &Integrator::getSubIntegrator, BP_RETURN_VALUE);

    BP_STRUCT(BSDFSamplingRecord, (bp::init<const Intersection &, Sampler *, ETransportMode>()))
        .def(bp::init<const Intersection &, const Vector &, ETransportMode>())
        .def(bp::init<const Intersection &, const Vector &, const Vector &, ETransportMode>())
        .add_property("its", bp::make_function(&bsdfsamplingrecord_get_its, BP_RETURN_VALUE))
        .def_readwrite("wi", &BSDFSamplingRecord::wi)
        .def_readwrite("wo", &BSDFSamplingRecord::wo)
        .def_readwrite("mode", &BSDFSamplingRecord::mode)
        .def_readwrite("typeMask", &BSDFSamplingRecord::typeMask)
        .def_readwrite("component", &BSDFSamplingRecord::component)
        .def_readwrite("sampledType", &BSDFSamplingRecord::sampledType)
        .def_readwrite("sampledComponent", &BSDFSamplingRecord::sampledComponent)
        .def("reverse", &BSDFSamplingRecord::reverse)
        .def("__repr__", &BSDFSamplingRecord::toString);

    BP_CLASS(BSDF, ConfigurableObject, bp::no_init)
        .def("getComponentCount", &BSDF::getComponentCount)
        .def("getType", &bsdf_getType_1)
        .def("getType", &bsdf_getType_2)
        .def("getMeasure", &BSDF::getMeasure)
        .def("hasComponent", &BSDF::hasComponent)
        .def("usesRayDifferentials", &BSDF::usesRayDifferentials)
        .def("getDiffuseReflectance", &BSDF::getDiffuseReflectance, BP_RETURN_VALUE)
        .def("getEta", &BSDF::getEta)
        .def("sample", &bsdf_sample, BP_RETURN_VALUE)
        .def("eval", &BSDF::eval, BP_RETURN_VALUE)
        .def("pdf", &BSDF::pdf)
        .staticmethod("getMeasure");

    BP_SETSCOPE(BSDF_class);
    bp::enum_<BSDF::EBSDFType>("EBSDFType")
        .value("ENull", BSDF::ENull)
        .value("EDiffuseReflection", BSDF::EDiffuseReflection)
        .value("EDiffuseTransmission", BSDF::EDiffuseTransmission)
        .value("EGlossyReflection", BSDF::EGlossyReflection)
        .value("EGlossyTransmission", BSDF::EGlossyTransmission)
        .value("EDeltaReflection", BSDF::EDeltaReflection)
        .value("EDeltaTransmission", BSDF::EDeltaTransmission)
        .value("EDelta1DReflection", BSDF::EDelta1DReflection)
        .value("EDelta1DTransmission", BSDF::EDelta1DTransmission)
        .value("EAnisotropic", BSDF::EAnisotropic)
        .value("ESpatiallyVarying", BSDF::ESpatiallyVarying)
        .value("ENonSymmetric", BSDF::ENonSymmetric)
        .value("EFrontSide", BSDF::EFrontSide)
        .value("EBackSide", BSDF::EBackSide)
        .value("EUsesSampler", BSDF::EUsesSampler)
        .export_values();

    bp::enum_<BSDF::ETypeCombinations>("ETypeCombinations")
        .value("EReflection", BSDF::EReflection)
        .value("ETransmission", BSDF::ETransmission)
        .value("EDiffuse", BSDF::EDiffuse)
        .value("EGlossy", BSDF::EGlossy)
        .value("ESmooth", BSDF::ESmooth)
        .value("EDelta", BSDF::EDelta)
        .value("EDelta1D", BSDF::EDelta1D)
        .value("EAll", BSDF::EAll)
        .export_values();

    BP_SETSCOPE(renderModule);

    BP_CLASS(Texture, ConfigurableObject, bp::no_init)
        .def("eval", &Texture::eval, BP_RETURN_VALUE)
        .def("evalGradient", &texture_evalGradient)
        .def("getAverage", &Texture::getAverage, BP_RETURN_VALUE)
        .def("getMinimum", &Texture::getMinimum, BP_RETURN_VALUE)
        .def("getMaximum", &Texture::getMaximum, BP_RETURN_VALUE)
        .def("getResolution", &Texture::getResolution, BP_RETURN_VALUE)
        .def("isConstant", &Texture::isConstant)
        .def("isMonochromatic", &Texture::isMonochromatic)
        .def("expand", &Texture::expand, BP_RETURN_VALUE)
        .def("usesRayDifferentials", &Texture::usesRayDifferentials)
        .def("getBitmap", &Texture::getBitmap, getBitmap_overloads()[BP_RETURN_VALUE]);

    bp::class_<Noise>("Noise")
        .def("perlinNoise", &Noise::perlinNoise)
        .def("turbulence", &Noise::turbulence)
        .def("fbm", &Noise::fbm);

    void (ImageBlock::*imageBlock_put1)(const ImageBlock *) = &ImageBlock::put;
    bool (ImageBlock::*imageBlock_put2)(const Point2 &, const Spectrum &, Float) = &ImageBlock::put;
    Bitmap *(ImageBlock::*imageBlock_getBitmap)() = &ImageBlock::getBitmap;

    BP_CLASS(ImageBlock, WorkResult, (bp::init<Bitmap::EPixelFormat, const Vector2i &, bp::optional<const ReconstructionFilter *, int, bool> >()))
        .def("setOffset", &ImageBlock::setOffset)
        .def("getOffset", &ImageBlock::getOffset, BP_RETURN_VALUE)
        .def("setSize", &ImageBlock::setSize)
        .def("getSize", &ImageBlock::getSize, BP_RETURN_VALUE)
        .def("getWidth", &ImageBlock::getWidth)
        .def("getHeight", &ImageBlock::getHeight)
        .def("setWarn", &ImageBlock::setWarn)
        .def("getWarn", &ImageBlock::getWarn)
        .def("getBorderSize", &ImageBlock::getBorderSize)
        .def("getChannelCount", &ImageBlock::getChannelCount)
        .def("getBitmap", imageBlock_getBitmap, BP_RETURN_VALUE)
        .def("clear", &ImageBlock::clear)
        .def("put", imageBlock_put1)
        .def("put", imageBlock_put2)
        .def("clone", &ImageBlock::clone, BP_RETURN_VALUE)
        .def("copyTo", &ImageBlock::copyTo);

    BP_CLASS(RectangularWorkUnit, WorkUnit, bp::init<>())
        .def("getOffset", &RectangularWorkUnit::getOffset, BP_RETURN_VALUE)
        .def("setOffset", &RectangularWorkUnit::setOffset)
        .def("getSize", &RectangularWorkUnit::getSize, BP_RETURN_VALUE)
        .def("setSize", &RectangularWorkUnit::setSize);

    bp::class_<RenderListener, ref<RenderListenerWrapper>, boost::noncopyable>
            RenderListener_class("RenderListener", bp::init<>());
    bp::register_ptr_to_python<RenderListener*>();
    bp::implicitly_convertible<ref<RenderListenerWrapper>, ref<RenderListener> >();

    RenderListener_class.def("workBeginEvent", &RenderListener::workBeginEvent)
        .def("workEndEvent", &RenderListener::workEndEvent)
        .def("workCanceledEvent", &RenderListener::workCanceledEvent)
        .def("refreshEvent", &RenderListener::refreshEvent)
        .def("finishJobEvent", &RenderListener::finishJobEvent);

    bp::detail::current_scope = oldScope;
}

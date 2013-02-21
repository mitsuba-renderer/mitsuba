#include "base.h"
#include <mitsuba/render/scene.h>
#include <mitsuba/render/scenehandler.h>
#include <mitsuba/render/renderqueue.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/render/noise.h>

using namespace mitsuba;

bool intersection_get_hasUVPartials(const Intersection &its) { return its.hasUVPartials;  }
void intersection_set_hasUVPartials(Intersection &its, bool value) { its.hasUVPartials = value;  }
uint32_t intersection_get_primIndex(const Intersection &its) { return its.primIndex;  }
void intersection_set_primIndex(Intersection &its, uint32_t value) { its.primIndex = value;  }

const Intersection &bsdfsamplingrecord_get_its(const BSDFSamplingRecord &bRec) {
	return bRec.its;
}

unsigned int bsdf_getType_1(const BSDF *bsdf) {
	return bsdf->getType();
}

unsigned int bsdf_getType_2(const BSDF *bsdf, int index) {
	return bsdf->getType(index);
}

bp::tuple bsdf_sample(const BSDF *bsdf, BSDFSamplingRecord &bRec, const Point2 &sample) {
	Float pdf;
	Spectrum result = bsdf->sample(bRec, pdf, sample);
	return bp::make_tuple(result, pdf);
}

bp::object shape_rayIntersect(const Shape *shape, const Ray &ray, Float mint, Float maxt) {
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

bp::object scene_rayIntersect(const Scene *scene, const Ray &ray) {
	Intersection its;

	if (!scene->rayIntersect(ray, its))
		return bp::object();

	return bp::object(its);
}

bp::object scene_rayIntersectAll(const Scene *scene, const Ray &ray) {
	Intersection its;

	if (!scene->rayIntersectAll(ray, its))
		return bp::object();

	return bp::object(its);
}

bp::tuple shape_getCurvature(const Shape *shape, const Intersection &its, bool shadingFrame) {
	Float H, K;
	shape->getCurvature(its, H, K, shadingFrame);
	return bp::make_tuple(H, K);
}

bp::tuple shape_getNormalDerivative(const Shape *shape, const Intersection &its, bool shadingFrame) {
	Vector dpdu, dpdv;
	shape->getNormalDerivative(its, dpdu, dpdv, shadingFrame);
	return bp::make_tuple(dpdu, dpdv);
}

ref<Scene> loadScene(const fs::path &filename, const StringMap &params) {
	SceneHandler::ParameterMap pmap;
	for (StringMap::const_iterator it = params.begin(); it != params.end(); ++it)
		pmap[it->first]=it->second;
	return SceneHandler::loadScene(filename, pmap);
}

bp::list scene_getSensors(Scene *scene) {
	bp::list list;
	ref_vector<Sensor> &sensors = scene->getSensors();
	for (size_t i=0; i<sensors.size(); ++i)
		list.append(cast(sensors[i].get()));
	return list;
}

bp::object scene_getSensor(Scene *scene) { return cast(scene->getSensor()); }
bp::object scene_getIntegrator(Scene *scene) { return cast(scene->getIntegrator()); }

bp::list scene_getMeshes(Scene *scene) {
	bp::list list;
	std::vector<TriMesh *> &meshes = scene->getMeshes();
	for (size_t i=0; i<meshes.size(); ++i)
		list.append(cast(meshes[i]));
	return list;
}

bp::list scene_getShapes(Scene *scene) {
	bp::list list;
	ref_vector<Shape> &shapes = scene->getShapes();
	for (size_t i=0; i<shapes.size(); ++i)
		list.append(cast(shapes[i].get()));
	return list;
}

bp::list scene_getEmitters(Scene *scene) {
	bp::list list;
	ref_vector<Emitter> &emitters = scene->getEmitters();
	for (size_t i=0; i<emitters.size(); ++i)
		list.append(cast(emitters[i].get()));
	return list;
}

bp::list scene_getMedia(Scene *scene) {
	bp::list list;
	ref_vector<Medium> &media = scene->getMedia();
	for (size_t i=0; i<media.size(); ++i)
		list.append(cast(media[i].get()));
	return list;
}

void export_render() {
	bp::object renderModule(
		bp::handle<>(bp::borrowed(PyImport_AddModule("mitsuba.render"))));
	bp::scope().attr("render") = renderModule;
	PyObject *oldScope = bp::detail::current_scope;

	BP_SETSCOPE(renderModule);

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

	Sampler *(Scene::*scene_getSampler)(void) = &Scene::getSampler;
	Film *(Scene::*scene_getFilm)(void) = &Scene::getFilm;

	BP_CLASS(Scene, NetworkedObject, bp::init<>())
		.def(bp::init<Properties>())
		.def(bp::init<Scene *>())
		.def(bp::init<Stream *, InstanceManager *>())
		.def("initialize", &Scene::initialize)
		.def("preprocess", &Scene::preprocess)
		.def("render", &Scene::render)
		.def("postprocess", &Scene::postprocess)
		.def("flush", &Scene::flush)
		.def("cancel", &Scene::cancel)
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
		.def("getMedia", &scene_getMedia);

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

	BP_CLASS(ReconstructionFilter, ConfigurableObject, bp::no_init)
		.def("eval", &ReconstructionFilter::eval)
		.def("evalDiscretized", &ReconstructionFilter::evalDiscretized)
		.def("getRadius", &ReconstructionFilter::getRadius)
		.def("getBorderSize", &ReconstructionFilter::getBorderSize);

	bp::class_<SceneHandler, boost::noncopyable>("SceneHandler", bp::no_init)
		.def("loadScene", &loadScene, BP_RETURN_VALUE)
		.staticmethod("loadScene");

	BP_CLASS(RenderJob, Thread, (bp::init<const std::string &, Scene *, RenderQueue *>()))
		.def(bp::init<const std::string &, Scene *, RenderQueue *, int>())
		.def(bp::init<const std::string &, Scene *, RenderQueue *, int, int>())
		.def(bp::init<const std::string &, Scene *, RenderQueue *, int, int, int>())
		.def("flush", &RenderJob::flush)
		.def("cancel", &RenderJob::cancel)
		.def("wait", &RenderJob::wait);

	BP_CLASS(RenderQueue, Object, bp::init<>())
		.def("getJobCount", &RenderQueue::getJobCount)
		.def("addJob", &RenderQueue::addJob)
		.def("removeJob", &RenderQueue::removeJob)
		.def("waitLeft", &RenderQueue::waitLeft)
		.def("join", &RenderQueue::join)
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
		.def_readwrite("shape", &Intersection::shape)
		.def_readwrite("instance", &Intersection::instance)
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
		.def("getElement", &Shape::getElement, BP_RETURN_VALUE)
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
		.def("getPrimitiveCount", &Shape::getPrimitiveCount)
		.def("getEffectivePrimitiveCount", &Shape::getEffectivePrimitiveCount);

	void (TriMesh::*triMesh_serialize1)(Stream *stream) const = &TriMesh::serialize;
	void (TriMesh::*triMesh_serialize2)(Stream *stream, InstanceManager *) const = &TriMesh::serialize;

	BP_CLASS(TriMesh, Shape, (bp::init<std::string, size_t, size_t, bool, bool, bool, bool, bool>()))
		.def(bp::init<Stream *, InstanceManager *>())
		.def(bp::init<Stream *, int>())
		.def("getTriangleCount", &TriMesh::getTriangleCount)
		.def("getVertexCount", &TriMesh::getVertexCount)
		.def("hasVertexNormals", &TriMesh::hasVertexNormals)
		.def("hasVertexColors", &TriMesh::hasVertexColors)
		.def("hasVertexTexcoords", &TriMesh::hasVertexTexcoords)
		.def("hasUVTangents", &TriMesh::hasUVTangents)
		.def("computeUVTangents", &TriMesh::computeUVTangents)
		.def("computeNormals", &TriMesh::computeNormals)
		.def("rebuildTopology", &TriMesh::rebuildTopology)
		.def("serialize", triMesh_serialize1)
		.def("serialize", triMesh_serialize2)
		.def("writeOBJ", &TriMesh::writeOBJ);

	BP_CLASS(Sensor, ConfigurableObject, bp::no_init) // incomplete
		.def("getShutterOpen", &Sensor::getShutterOpen)
		.def("setShutterOpen", &Sensor::setShutterOpen)
		.def("getShutterOpenTime", &Sensor::getShutterOpenTime)
		.def("setShutterOpenTime", &Sensor::setShutterOpenTime);

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
		.value("ENull,", BSDF::ENull)
		.value("EDiffuseReflection,", BSDF::EDiffuseReflection)
		.value("EDiffuseTransmission,", BSDF::EDiffuseTransmission)
		.value("EGlossyReflection", BSDF::EGlossyReflection)
		.value("EGlossyTransmission,", BSDF::EGlossyTransmission)
		.value("EDeltaReflection,", BSDF::EDeltaReflection)
		.value("EDeltaTransmission,", BSDF::EDeltaTransmission)
		.value("EDelta1DReflection,", BSDF::EDelta1DReflection)
		.value("EDelta1DTransmission,", BSDF::EDelta1DTransmission)
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
		.def("getAverage", &Texture::getAverage, BP_RETURN_VALUE)
		.def("getMinimum", &Texture::getMinimum, BP_RETURN_VALUE)
		.def("getMaximum", &Texture::getMaximum, BP_RETURN_VALUE)
		.def("getResolution", &Texture::getResolution, BP_RETURN_VALUE)
		.def("isConstant", &Texture::isConstant)
		.def("usesRayDifferentials", &Texture::usesRayDifferentials);

	bp::class_<Noise>("Noise")
		.def("perlinNoise", &Noise::perlinNoise)
		.def("turbulence", &Noise::turbulence)
		.def("fbm", &Noise::fbm);

	bp::detail::current_scope = oldScope;
}

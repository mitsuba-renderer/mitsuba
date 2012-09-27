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

ref<Scene> loadScene(const fs::path &filename, const StringMap &params) {
	SceneHandler::ParameterMap pmap;
	for (StringMap::const_iterator it = params.begin(); it != params.end(); ++it)
		pmap[it->first]=it->second;
	return SceneHandler::loadScene(filename, pmap);
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
		.def("getAABB", &Scene::getAABB, BP_RETURN_VALUE)
		.def("getBSphere", &Scene::getBSphere, BP_RETURN_VALUE)
		.def("getBlockSize", &Scene::getBlockSize)
		.def("setBlockSize", &Scene::setBlockSize)
		.def("getSourceFile", &Scene::getSourceFile, BP_RETURN_VALUE)
		.def("setSourceFile", &Scene::setSourceFile)
		.def("getDestinationFile", &Scene::getDestinationFile, BP_RETURN_VALUE)
		.def("setDestinationFile", &Scene::setDestinationFile)
		.def("destinationExists", &Scene::destinationExists);

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
		.def("fbm", &Noise::fbm)
		.def("turbulence", &Noise::turbulence);

	bp::detail::current_scope = oldScope;
}

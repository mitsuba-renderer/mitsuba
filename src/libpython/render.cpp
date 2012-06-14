#include "base.h"
#include <mitsuba/render/scene.h>
#include <mitsuba/render/scenehandler.h>
#include <mitsuba/render/renderqueue.h>
#include <mitsuba/render/renderjob.h>

using namespace mitsuba;

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

	bp::detail::current_scope = oldScope;
}

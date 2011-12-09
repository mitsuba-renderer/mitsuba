#include "base.h"
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/appender.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/aabb.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/sched_remote.h>
#include <mitsuba/core/netobject.h>
#include <mitsuba/core/sstream.h>
#include <mitsuba/core/sshstream.h>
#include <mitsuba/render/scenehandler.h>

using namespace mitsuba;

void initializeFramework() {
	/* Initialize the core framework */
	Class::staticInitialization();
	PluginManager::staticInitialization();
	Statistics::staticInitialization();
	Thread::staticInitialization();
	Logger::staticInitialization();
	Spectrum::staticInitialization();
	Scheduler::staticInitialization();
	SHVector::staticInitialization();
	SceneHandler::staticInitialization();
}

void shutdownFramework() {
	/* Shutdown the core framework */
	SceneHandler::staticShutdown();
	SHVector::staticShutdown();
	Scheduler::staticShutdown();
	Spectrum::staticShutdown();
	Logger::staticShutdown();
	Thread::staticShutdown();
	Statistics::staticShutdown();
	PluginManager::staticShutdown();
	Class::staticShutdown();
}

class spectrum_wrapper {
public:
	static Float get(const Spectrum &spec, int i) {
		if (i < 0 || i >= SPECTRUM_SAMPLES) {
			SLog(EError, "Index %i is out of range!", i);
			return 0.0f;
		}
		return spec[i];
	}
	
	static void set(Spectrum &spec, int i, Float value) {
		if (i < 0 || i >= SPECTRUM_SAMPLES) 
			SLog(EError, "Index %i is out of range!", i);
		else
			spec[i] = value;
	}

	static int len(Spectrum &) {
		return SPECTRUM_SAMPLES;
	}
};

class properties_wrapper {
public:
	static bp::object get(const Properties &props, std::string name) {
		if (!props.hasProperty(name))
			SLog(EError, "Properties: keyword \"%s\" not found!", name.c_str());

		switch (props.getType(name)) {
			case Properties::EBoolean:
				return bp::object(props.getBoolean(name));
			case Properties::EString:
				return bp::object(props.getString(name));
			case Properties::EInteger:
				return bp::object(props.getInteger(name));
			case Properties::EFloat:
				return bp::object(props.getFloat(name));
			case Properties::EVector:
				return bp::object(props.getVector(name));
			case Properties::EPoint:
				return bp::object(props.getPoint(name));
			case Properties::ETransform:
				return bp::object(props.getTransform(name));
			case Properties::ESpectrum:
				return bp::object(props.getSpectrum(name));
			default:
				SLog(EError, "Properties: type of keyword \"%s\" is not supported!", name.c_str());
				return bp::object();
		}
	}

	static void set(Properties &props, const std::string name, bp::object value) {
		bp::extract<std::string> extractString(value);
		bp::extract<bool> extractBoolean(value);
		bp::extract<int> extractInteger(value);
		bp::extract<Float> extractFloat(value);
		bp::extract<Vector> extractVector(value);
		bp::extract<Point> extractPoint(value);
		bp::extract<Transform> extractTransform(value);
		bp::extract<Spectrum> extractSpectrum(value);

		if (extractString.check()) {
			props.setString(name, extractString());
		} else if (extractBoolean.check() && PyObject_IsInstance(value.ptr(), (PyObject *) &PyBool_Type)) {
			props.setBoolean(name, extractBoolean());
		} else if (extractInteger.check()) {
			props.setInteger(name, extractInteger());
		} else if (extractFloat.check()) {
			props.setFloat(name, extractFloat());
		} else if (extractPoint.check()) {
			props.setPoint(name, extractPoint());
		} else if (extractVector.check()) {
			props.setVector(name, extractVector());
		} else if (extractTransform.check()) {
			props.setTransform(name, extractTransform());
		} else if (extractSpectrum.check()) {
			props.setSpectrum(name, extractSpectrum());
		} else {
			SLog(EError, "Properties: type of keyword \"%s\" is not supported!", name.c_str());
		}
	}
};

struct path_to_python_str {
	static PyObject* convert(fs::path const& path) {
		return boost::python::incref(
			boost::python::object(path.file_string()).ptr());
	}
};


static void Matrix4x4_setItem(Matrix4x4 *matrix, bp::tuple tuple, Float value) {
	if (bp::len(tuple) != 2)
		SLog(EError, "Invalid matrix indexing operation, required a tuple of length 2");
	int i = bp::extract<int>(tuple[0]);
	int j = bp::extract<int>(tuple[1]);

	if (i < 0 || j < 0 || i >= 4 || j >= 4)
		SLog(EError, "Index (%i, %i) is out of bounds!", i, j);

	matrix->operator()(i, j) = value;
}

static Float Matrix4x4_getItem(Matrix4x4 *matrix, bp::tuple tuple) {
	if (bp::len(tuple) != 2)
		SLog(EError, "Invalid matrix indexing operation, required a tuple of length 2");
	int i = bp::extract<int>(tuple[0]);
	int j = bp::extract<int>(tuple[1]);

	if (i < 0 || j < 0 || i >= 4 || j >= 4)
		SLog(EError, "Index (%i, %i) is out of bounds!", i, j);

	return matrix->operator()(i, j);
}

Class *object_getClass(Object *object) {
	return const_cast<Class *>(object->getClass());
}

Class *class_forName(const char *name) {
	return const_cast<Class *>(Class::forName(name));
}

Class *class_getSuperClass(Class *theClass) {
	return const_cast<Class *>(theClass->getSuperClass());
}

static ref<SerializableObject> instance_manager_getinstance(InstanceManager *manager, Stream *stream) {
	return manager->getInstance(stream);
}

void appender_logProgress(Appender *appender, Float progress, const std::string &name,
	const std::string &formatted, const std::string &eta) {
	appender->logProgress(progress, name, formatted, eta, NULL);
}

void logger_logProgress(Logger *logger, Float progress, const std::string &name,
	const std::string &formatted, const std::string &eta) {
	logger->logProgress(progress, name, formatted, eta, NULL);
}
	
void mts_log(ELogLevel level, const std::string &msg) {
	bp::object traceback(bp::import("traceback"));
	bp::object extract_stack(traceback.attr("extract_stack"));
	bp::object stack = extract_stack();
	bp::object top(stack[bp::len(stack)-1]);
	std::string module = bp::extract<std::string>(top[2]);
	Thread::getThread()->getLogger()->log(level, 
		NULL, bp::extract<const char *>(top[0]), 
		bp::extract<int>(top[1]), "%s%s: %s", 
		module.c_str(), module[0] == '<' ? "" : "()",
		msg.c_str());
}

class FormatterWrapper : public Formatter {
public:
	FormatterWrapper(PyObject *self) : m_self(self) { Py_INCREF(m_self); }

	std::string format(ELogLevel logLevel, const Class *theClass,
			const Thread *thread, const std::string &text, 
			const char *file, int line) {
		return bp::call_method<std::string>(m_self, "format", logLevel, 
				bp::ptr(const_cast<Class *>(theClass)),
				bp::ptr(const_cast<Thread *>(thread)), text, file, line);
	}

	virtual ~FormatterWrapper() {
		Py_DECREF(m_self);
	}
private:
	PyObject *m_self;
};

class AppenderWrapper : public Appender {
public:
	AppenderWrapper(PyObject *self) : m_self(self) { Py_INCREF(m_self); }

	void append(ELogLevel level, const std::string &text) {
		bp::call_method<void>(m_self, "append", level, text);
	}

	void logProgress(Float progress, const std::string &name,
		const std::string &formatted, const std::string &eta, 
		const void *ptr) {
		bp::call_method<void>(m_self, "logProgress", name, formatted, eta);
	}

	virtual ~AppenderWrapper() {
		Py_DECREF(m_self);
	}
private:
	PyObject *m_self;
};

static Spectrum *spectrum_array_constructor(bp::list list) {
	Float spec[SPECTRUM_SAMPLES];
	if (bp::len(list) != SPECTRUM_SAMPLES)
		SLog(EError, "Spectrum: expected %i arguments", SPECTRUM_SAMPLES);

	for (int i=0; i<bp::len(list); ++i)
		spec[i] = bp::extract<Float>(list[i]);

	return new Spectrum(spec);
}

static Point ray_eval(const Ray &ray, Float t) {
	return ray(t);
}

void aabb_expandby_aabb(AABB *aabb, const AABB &aabb2) { aabb->expandBy(aabb2); }
void aabb_expandby_point(AABB *aabb, const Point &p) { aabb->expandBy(p); }
Float aabb_distanceto_aabb(AABB *aabb, const AABB &aabb2) { return aabb->distanceTo(aabb2); }
Float aabb_distanceto_point(AABB *aabb, const Point &p) { return aabb->distanceTo(p); }
Float aabb_sqrdistanceto_aabb(AABB *aabb, const AABB &aabb2) { return aabb->squaredDistanceTo(aabb2); }
Float aabb_sqrdistanceto_point(AABB *aabb, const Point &p) { return aabb->squaredDistanceTo(p); }
bool aabb_contains_aabb(AABB *aabb, const AABB &aabb2) { return aabb->contains(aabb2); }
bool aabb_contains_point(AABB *aabb, const Point &p) { return aabb->contains(p); }
bp::object bsphere_rayIntersect(BSphere *bsphere, const Ray &ray) {
	Float nearT, farT;
	if (bsphere->rayIntersect(ray, nearT, farT))
		return bp::make_tuple(nearT, farT);
	else
		return bp::object();

}

bp::object aabb_rayIntersect(AABB *aabb, const Ray &ray) {
	Float nearT, farT;
	if (aabb->rayIntersect(ray, nearT, farT))
		return bp::make_tuple(nearT, farT);
	else
		return bp::object();

}

Vector transform_mul_vector(Transform *transform, const Vector &vector) { return transform->operator()(vector); }
Vector4 transform_mul_vector4(Transform *transform, const Vector4 &vector) { return transform->operator()(vector); }
Normal transform_mul_normal(Transform *transform, const Normal &normal) { return transform->operator()(normal); }
Point transform_mul_point(Transform *transform, const Point &point) { return transform->operator()(point); }
Ray transform_mul_ray(Transform *transform, const Ray &ray) { return transform->operator()(ray); }
Transform transform_mul_transform(Transform *transform, const Transform &other) { return *transform * other; }

ConfigurableObject *pluginmgr_create(PluginManager *manager, bp::dict dict) {
	Properties properties;
	bp::list list = dict.items();
	std::map<std::string, ConfigurableObject *> children;

	for (int i=0; i<bp::len(list); ++i) {
		bp::tuple tuple = bp::extract<bp::tuple>(list[i]);
		std::string name = bp::extract<std::string>(tuple[0]);
		bp::extract<bp::dict> extractDict(tuple[1]);
		bp::extract<std::string> extractString(tuple[1]);
		bp::extract<ConfigurableObject *> extractConfigurableObject(tuple[1]);

		if (name == "type") {
			if (!extractString.check())
				SLog(EError, "'type' property must map to a string!");
			else
				properties.setPluginName(extractString());
		} else if (extractDict.check()) {
			children[name] = pluginmgr_create(manager, extractDict());
		} else if (extractConfigurableObject.check()) {
			children[name] = extractConfigurableObject();
		} else {
			properties_wrapper::set(properties, name, tuple[1]);
		}
	}

	ConfigurableObject *object = manager->createObject(properties);
	for (std::map<std::string, ConfigurableObject *>::iterator it = children.begin();
		it != children.end(); ++it) 
		object->addChild(it->first, it->second);
	object->configure();
	return object;
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(fromLinearRGB_overloads, fromLinearRGB, 3, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(fromXYZ_overloads, fromXYZ, 3, 4)

void export_core() {
	bp::to_python_converter<fs::path, path_to_python_str>();
	bp::implicitly_convertible<std::string, fs::path>();

	bp::object coreModule(
		bp::handle<>(bp::borrowed(PyImport_AddModule("mitsuba.core"))));
	bp::scope().attr("core") = coreModule;
	PyObject *oldScope = bp::detail::current_scope;

	BP_SETSCOPE(coreModule);

	/* Basic STL containers */
	bp::class_<StringVector>("StringVector")
		.def(bp::vector_indexing_suite<StringVector>());
	bp::class_<StringMap>("StringMap")
		.def(bp::map_indexing_suite<StringMap>());

	bp::enum_<ELogLevel>("ELogLevel")
		.value("ETrace", ETrace)
		.value("EDebug", EDebug)
		.value("EInfo", EInfo)
		.value("EWarn", EWarn)
		.value("EError", EError)
		.export_values();

	bp::def("Log", &mts_log);

	bp::class_<Class, boost::noncopyable>("Class", bp::no_init)
		.def("getName", &Class::getName, BP_RETURN_CONSTREF)
		.def("isAbstract", &Class::isAbstract)
		.def("isInstantiable", &Class::isInstantiable)
		.def("isSerializable", &Class::isSerializable)
		.def("derivesFrom", &Class::derivesFrom)
		.def("getSuperClass", &class_getSuperClass, BP_RETURN_INTREF)
		.def("forName", &class_forName, BP_RETURN_INTREF)
		.def("unserialize", &Class::unserialize, BP_RETURN_VALUE)
		.def("instantiate", &Class::instantiate, BP_RETURN_VALUE)
		.staticmethod("forName");
	bp::register_ptr_to_python<Class*>();

	bp::class_<Object, ref<Object>, boost::noncopyable>("Object", bp::no_init)
		.def("getRefCount", &Object::getRefCount)
		.def("getClass", &object_getClass, BP_RETURN_INTREF)
		.def("__str__", &Object::toString);
	bp::register_ptr_to_python<Object*>();

	BP_CLASS(Stream, Object, bp::no_init)
		.def("setByteOrder", &Stream::setByteOrder)
		.def("getByteOrder", &Stream::getByteOrder)
		.def("getHostByteOrder", &Stream::getHostByteOrder)
		.def("truncate", &Stream::truncate)
		.def("setPos", &Stream::setPos)
		.def("getPos", &Stream::getPos)
		.def("getSize", &Stream::getSize)
		.def("flush", &Stream::flush)
		.def("canWrite", &Stream::canWrite)
		.def("canRead", &Stream::canRead)
		.def("skip", &Stream::skip)
		.def("copyTo", &Stream::copyTo)
		.def("writeString", &Stream::writeString)
		.def("readString", &Stream::readString)
		.def("writeLine", &Stream::writeLine)
		.def("readLine", &Stream::readLine)
		.def("writeChar", &Stream::writeChar)
		.def("readChar", &Stream::readChar)
		.def("writeUChar", &Stream::writeUChar)
		.def("readUChar", &Stream::readUChar)
		.def("writeShort", &Stream::writeShort)
		.def("readShort", &Stream::readShort)
		.def("writeUShort", &Stream::writeUShort)
		.def("readUShort", &Stream::readUShort)
		.def("writeInt", &Stream::writeInt)
		.def("readInt", &Stream::readInt)
		.def("writeUInt", &Stream::writeUInt)
		.def("readUInt", &Stream::readUInt)
		.def("writeLong", &Stream::writeLong)
		.def("readLong", &Stream::readLong)
		.def("writeULong", &Stream::writeULong)
		.def("readULong", &Stream::readULong)
		.def("writeFloat", &Stream::writeFloat)
		.def("readFloat", &Stream::readFloat)
		.def("writeSingle", &Stream::writeSingle)
		.def("readSingle", &Stream::readSingle)
		.def("writeDouble", &Stream::writeDouble)
		.def("readDouble", &Stream::readDouble);

	BP_SETSCOPE(Stream_class);
	bp::enum_<Stream::EByteOrder>("EByteOrder")
		.value("EBigEndian", Stream::EBigEndian)
		.value("ELittleEndian", Stream::ELittleEndian)
		.value("ENetworkByteOrder", Stream::ENetworkByteOrder)
		.export_values();
	BP_SETSCOPE(coreModule);

	BP_CLASS(FileStream, Stream, bp::init<>())
		.def(bp::init<std::string, FileStream::EFileMode>())
		.def("getPath", &FileStream::getPath, 
				BP_RETURN_CONSTREF)
		.def("open", &FileStream::open)
		.def("close", &FileStream::close)
		.def("remove", &FileStream::remove);

	BP_CLASS(SocketStream, Stream, (bp::init<std::string, int>()))
		.def("getPeer", &SocketStream::getPeer, BP_RETURN_CONSTREF)
		.def("getReceivedBytes", &SocketStream::getReceivedBytes)
		.def("getSentBytes", &SocketStream::getSentBytes);

	BP_CLASS(SSHStream, Stream, (bp::init<std::string, std::string, const StringVector &>()))
		.def(bp::init<std::string, std::string, const StringVector &, int>())
		.def(bp::init<std::string, std::string, const StringVector &, int, int>())
		.def("getUserName", &SSHStream::getUserName, BP_RETURN_CONSTREF)
		.def("getHostName", &SSHStream::getHostName, BP_RETURN_CONSTREF)
		.def("getReceivedBytes", &SSHStream::getReceivedBytes)
		.def("getSentBytes", &SSHStream::getSentBytes);

	BP_SETSCOPE(FileStream_class);
	bp::enum_<FileStream::EFileMode>("EFileMode")
		.value("EReadOnly", FileStream::EReadOnly)
		.value("EReadWrite", FileStream::EReadWrite)
		.value("ETruncWrite", FileStream::ETruncWrite)
		.value("ETruncReadWrite", FileStream::ETruncReadWrite)
		.value("EAppendWrite", FileStream::EAppendWrite)
		.value("EAppendReadWrite", FileStream::EAppendReadWrite)
		.export_values();
	BP_SETSCOPE(coreModule);

	BP_CLASS(SerializableObject, Object, bp::no_init)
		.def("serialize", &SerializableObject::serialize);
	
	ConfigurableObject *(ConfigurableObject::*cobject_get_parent)() = &ConfigurableObject::getParent;
	void (ConfigurableObject::*cobject_add_child_1)(ConfigurableObject *) = &ConfigurableObject::addChild;
	void (ConfigurableObject::*cobject_add_child_2)(const std::string &, ConfigurableObject *) = &ConfigurableObject::addChild;

	BP_CLASS(ConfigurableObject, SerializableObject, bp::no_init)
		.def("getParent", cobject_get_parent, BP_RETURN_VALUE)
		.def("setParent", &ConfigurableObject::setParent)
		.def("addChild", cobject_add_child_1)
		.def("addChild", cobject_add_child_2)
		.def("configure", &ConfigurableObject::configure);

	BP_CLASS(NetworkedObject, ConfigurableObject, bp::no_init)
		.def("bindUsedResources", &NetworkedObject::bindUsedResources);

	Thread *(Thread::*thread_get_parent)() = &Thread::getParent;
	BP_CLASS(Thread, Object, bp::no_init)
		.def("getID", &Thread::getID)
		.def("getStackSize", &Thread::getStackSize)
		.def("setPriority", &Thread::setPriority)
		.def("getPriority", &Thread::getPriority)
		.def("setCritical", &Thread::setCritical)
		.def("getCritical", &Thread::getCritical)
		.def("setName", &Thread::setName)
		.def("getName", &Thread::getName, BP_RETURN_CONSTREF)
		.def("getParent", thread_get_parent, BP_RETURN_VALUE)
		.def("setLogger", &Thread::setLogger)
		.def("getLogger", &Thread::getLogger, BP_RETURN_VALUE)
		.def("setFileResolver", &Thread::setFileResolver)
		.def("getFileResolver", &Thread::getFileResolver, BP_RETURN_VALUE)
		.def("getThread", &Thread::getThread, BP_RETURN_VALUE)
		.def("isRunning", &Thread::isRunning)
		.def("registerUnmanagedThread", &Thread::registerUnmanagedThread, BP_RETURN_VALUE)
		.def("sleep", &Thread::sleep)
		.def("detach", &Thread::detach)
		.def("join", &Thread::join)
		.def("start", &Thread::start)
		.staticmethod("sleep")
		.staticmethod("getThread")
		.staticmethod("registerUnmanagedThread");

	BP_SETSCOPE(Thread_class);
	bp::enum_<Thread::EThreadPriority>("EThreadPriority")
		.value("EIdlePriority", Thread::EIdlePriority)
		.value("ELowestPriority", Thread::ELowestPriority)
		.value("ELowPriority", Thread::ELowPriority)
		.value("ENormalPriority", Thread::ENormalPriority)
		.value("EHighPriority", Thread::EHighPriority)
		.value("EHighestPriority", Thread::EHighestPriority)
		.value("ERealtimePriority", Thread::ERealtimePriority)
		.export_values();
	BP_SETSCOPE(coreModule);

	BP_WRAPPED_CLASS(Appender, AppenderWrapper, Object, bp::init<>())
		.def("append", &Appender::append)
		.def("logProgress", &appender_logProgress);

	BP_WRAPPED_CLASS(Formatter, FormatterWrapper, Object, bp::init<>())
		.def("format", &Formatter::format);
	
	Appender *(Logger::*logger_get_appender)(size_t) = &Logger::getAppender;
	BP_CLASS(Logger, Object, bp::init<ELogLevel>())
		.def("logProgress", logger_logProgress)
		.def("setLogLevel", &Logger::setLogLevel)
		.def("getLogLevel", &Logger::getLogLevel)
		.def("setErrorLevel", &Logger::setErrorLevel)
		.def("getErrorLevel", &Logger::getErrorLevel)
		.def("addAppender", &Logger::addAppender)
		.def("removeAppender", &Logger::removeAppender)
		.def("clearAppenders", &Logger::clearAppenders)
		.def("getAppenderCount", &Logger::getAppenderCount)
		.def("getAppender", logger_get_appender, BP_RETURN_VALUE)
		.def("getFormatter", &Logger::getFormatter, BP_RETURN_VALUE)
		.def("setFormatter", &Logger::setFormatter)
		.def("getWarningCount", &Logger::getWarningCount);

	BP_CLASS(InstanceManager, Object, bp::init<>())
		.def("serialize", &InstanceManager::serialize)
		.def("getInstance", &instance_manager_getinstance, BP_RETURN_VALUE);

	bp::class_<ContinuousSpectrum, boost::noncopyable>("ContinuousSpectrum", bp::no_init)
		.def("eval", &ContinuousSpectrum::eval)
		.def("average", &ContinuousSpectrum::average)
		.def("__str__", &ContinuousSpectrum::toString);

	bp::class_<InterpolatedSpectrum, bp::bases<ContinuousSpectrum>, boost::noncopyable>
			("InterpolatedSpectrum", bp::init<>())
		.def(bp::init<size_t>())
		.def(bp::init<fs::path>())
		.def("append", &InterpolatedSpectrum::append)
		.def("clear", &InterpolatedSpectrum::clear)
		.def("zeroExtend", &InterpolatedSpectrum::zeroExtend);

	BP_CLASS(Bitmap, Object, (bp::init<int, int, int>()))
		.def(bp::init<Bitmap::EFileFormat, Stream *>())
		.def("clone", &Bitmap::clone, BP_RETURN_VALUE)
		.def("clear", &Bitmap::clear)
		.def("save", &Bitmap::save)
		.def("setTitle", &Bitmap::setTitle)
		.def("getTitle", &Bitmap::getTitle, BP_RETURN_INTREF)
		.def("setAuthor", &Bitmap::setAuthor)
		.def("getAuthor", &Bitmap::getAuthor, BP_RETURN_INTREF)
		.def("setComment", &Bitmap::setComment)
		.def("getComment", &Bitmap::getComment, BP_RETURN_INTREF)
		.def("setGamma", &Bitmap::setGamma)
		.def("getGamma", &Bitmap::getGamma)
		.def("getWidth", &Bitmap::getWidth)
		.def("getHeight", &Bitmap::getHeight)
		.def("getBitsPerPixel", &Bitmap::getBitsPerPixel)
		.def("getSize", &Bitmap::getSize);

	BP_SETSCOPE(Bitmap_class);
	bp::enum_<Bitmap::EFileFormat>("EFileFormat")
		.value("EPNG", Bitmap::EPNG)
		.value("EEXR", Bitmap::EEXR)
		.value("ETGA", Bitmap::ETGA)
		.value("EBMP", Bitmap::EBMP)
		.value("EJPEG", Bitmap::EJPEG)
		.export_values();
	BP_SETSCOPE(coreModule);

	BP_CLASS(FileResolver, Object, bp::init<>())
		.def("resolve", &FileResolver::resolve, BP_RETURN_VALUE)
		.def("resolveAbsolute", &FileResolver::resolveAbsolute, BP_RETURN_VALUE)
		.def("clone", &FileResolver::clone, BP_RETURN_VALUE)
		.def("addPath", &FileResolver::addPath)
		.def("clear", &FileResolver::clear);

	void (Random::*random_seed_random)(Random *) = &Random::seed;
	void (Random::*random_seed_uint64_t)(uint64_t) = &Random::seed;
	BP_CLASS(Random, Object, bp::init<>())
		.def(bp::init<uint64_t>())
		.def(bp::init<Random *>())
		.def(bp::init<Stream *, InstanceManager *>())
		.def("set", &Random::set)
		.def("seed", random_seed_random)
		.def("seed", random_seed_uint64_t)
		.def("nextULong", &Random::nextULong)
		.def("nextUInt", &Random::nextUInt)
		.def("nextSize", &Random::nextSize)
		.def("nextFloat", &Random::nextFloat)
		.def("serialize", &Random::serialize);

	ConfigurableObject *(PluginManager::*pluginmgr_createobject_1)(const Properties &) = &PluginManager::createObject;
	ConfigurableObject *(PluginManager::*pluginmgr_createobject_2)(const Class *, const Properties &) = &PluginManager::createObject;

	BP_CLASS(PluginManager, Object, bp::no_init)
		.def("ensurePluginLoaded", &PluginManager::ensurePluginLoaded)
		.def("getLoadedPlugins", &PluginManager::getLoadedPlugins)
		.def("create", pluginmgr_create, BP_RETURN_VALUE)
		.def("createObject", pluginmgr_createobject_1, BP_RETURN_VALUE)
		.def("createObject", pluginmgr_createobject_2, BP_RETURN_VALUE)
		.def("getInstance", &PluginManager::getInstance, BP_RETURN_VALUE)
		.staticmethod("getInstance");

	BP_CLASS(Statistics, Object, bp::no_init)
		.def("getStats", &Statistics::getStats, BP_RETURN_VALUE)
		.def("printStats", &Statistics::printStats)
		.def("getInstance", &Statistics::getInstance, BP_RETURN_VALUE)
		.staticmethod("getInstance");

	BP_CLASS(WorkUnit, Object, bp::no_init)
		.def("set", &WorkUnit::set)
		.def("load", &WorkUnit::load)
		.def("save", &WorkUnit::save);
	
	BP_CLASS(WorkResult, Object, bp::no_init)
		.def("load", &WorkResult::load)
		.def("save", &WorkResult::save);
	
	BP_CLASS(WorkProcessor, SerializableObject, bp::no_init)
		.def("createWorkUnit", &WorkProcessor::createWorkUnit, BP_RETURN_VALUE)
		.def("createWorkResult", &WorkProcessor::createWorkResult, BP_RETURN_VALUE)
		.def("clone", &WorkProcessor::clone, BP_RETURN_VALUE)
		.def("prepare", &WorkProcessor::prepare)
		.def("process", &WorkProcessor::process);

	BP_CLASS(ParallelProcess, Object, bp::no_init)
		.def("generateWork", &ParallelProcess::generateWork)
		.def("processResult", &ParallelProcess::processResult)
		.def("handleCancellation", &ParallelProcess::handleCancellation)
		.def("getReturnStatus", &ParallelProcess::getReturnStatus)
		.def("createWorkProcessor", &ParallelProcess::createWorkProcessor)
		.def("bindResource", &ParallelProcess::bindResource)
		.def("isLocal", &ParallelProcess::isLocal)
		.def("getLogLevel", &ParallelProcess::getLogLevel)
		.def("getRequiredPlugins", &ParallelProcess::getRequiredPlugins, BP_RETURN_VALUE);
	
	BP_SETSCOPE(ParallelProcess_class);
	bp::enum_<ParallelProcess::EStatus>("EStatus")
		.value("EUnknown", ParallelProcess::EUnknown)
		.value("EPause", ParallelProcess::EPause)
		.value("ESuccess", ParallelProcess::ESuccess)
		.value("EFailure", ParallelProcess::EFailure)
		.export_values();
	BP_SETSCOPE(coreModule);
	
	BP_CLASS(Worker, Thread, bp::no_init)
		.def("getCoreCount", &Worker::getCoreCount)
		.def("isRemoteWorker", &Worker::isRemoteWorker);

	BP_CLASS(LocalWorker, Worker, bp::init<const std::string>()) 
		.def(bp::init<const std::string, Thread::EThreadPriority>());

	BP_CLASS(RemoteWorker, Worker, (bp::init<const std::string, Stream *>()))
		.def("getNodeName", &RemoteWorker::getNodeName, BP_RETURN_VALUE);

	bp::class_<SerializableObjectVector>("SerializableObjectVector")
		.def(bp::vector_indexing_suite<SerializableObjectVector>());

	bool (Scheduler::*scheduler_cancel)(ParallelProcess *)= &Scheduler::cancel;
	BP_CLASS(Scheduler, Object, bp::no_init)
		.def("schedule", &Scheduler::schedule)
		.def("wait", &Scheduler::wait)
		.def("cancel", scheduler_cancel)
		.def("registerResource", &Scheduler::registerResource)
		.def("registerManifoldResource", &Scheduler::registerManifoldResource)
		.def("retainResource", &Scheduler::retainResource)
		.def("unregisterResource", &Scheduler::unregisterResource)
		.def("getResourceID", &Scheduler::getResourceID)
		.def("registerWorker", &Scheduler::registerWorker)
		.def("unregisterWorker", &Scheduler::unregisterWorker)
		.def("getWorkerCount", &Scheduler::getWorkerCount)
		.def("getLocalWorkerCount", &Scheduler::getLocalWorkerCount)
		.def("getWorker", &Scheduler::getWorker, BP_RETURN_VALUE)
		.def("start", &Scheduler::start)
		.def("pause", &Scheduler::pause)
		.def("stop", &Scheduler::stop)
		.def("getCoreCount", &Scheduler::getCoreCount)
		.def("hasLocalWorkers", &Scheduler::hasLocalWorkers)
		.def("hasRemoteWorkers", &Scheduler::hasRemoteWorkers)
		.def("getInstance", &Scheduler::getInstance, BP_RETURN_VALUE)
		.def("isRunning", &Scheduler::isRunning)
		.def("isBusy", &Scheduler::isBusy)
		.staticmethod("getInstance");

	BP_STRUCT(Spectrum, bp::init<>())
		.def("__init__", bp::make_constructor(spectrum_array_constructor))
		.def(bp::init<Float>())
		.def(bp::init<Stream *>())
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self -= bp::self)
		.def(bp::self *= Float())
		.def(bp::self * Float())
		.def(bp::self *= bp::self)
		.def(bp::self * bp::self)
		.def(bp::self / Float())
		.def(bp::self /= Float())
		.def(bp::self /= bp::self)
		.def(bp::self / bp::self)
		.def("isValid", &Spectrum::isValid)
		.def("isNaN", &Spectrum::isNaN)
		.def("average", &Spectrum::average)
		.def("sqrt", &Spectrum::sqrt)
		.def("exp", &Spectrum::exp)
		.def("pow", &Spectrum::pow)
		.def("clampNegative", &Spectrum::clampNegative)
		.def("min", &Spectrum::min)
		.def("max", &Spectrum::max)
		.def("isZero", &Spectrum::isZero)
		.def("eval", &Spectrum::eval)
		.def("getLuminance", &Spectrum::getLuminance)
		.def("fromXYZ", &Spectrum::fromXYZ, fromXYZ_overloads())
		.def("toXYZ", &Spectrum::toXYZ)
		.def("fromLinearRGB", &Spectrum::fromLinearRGB, fromLinearRGB_overloads())
		.def("toLinearRGB", &Spectrum::toLinearRGB)
		.def("fromSRGB", &Spectrum::fromSRGB)
		.def("toSRGB", &Spectrum::toSRGB)
		.def("fromContinuousSpectrum", &Spectrum::fromContinuousSpectrum)
		.def("serialize", &Spectrum::serialize)
		.def("__str__", &Spectrum::toString)
		.def("__len__", &spectrum_wrapper::len)
		.def("__getitem__", &spectrum_wrapper::get)
		.def("__setitem__", &spectrum_wrapper::set);

	BP_SETSCOPE(Spectrum_struct);
	bp::enum_<Spectrum::EConversionIntent>("EConversionIntent")
		.value("EReflectance", Spectrum::EReflectance)
		.value("EIlluminant", Spectrum::EIlluminant)
		.export_values();
	BP_SETSCOPE(coreModule);

	bp::class_<Properties> properties("Properties");
	properties
		.def(bp::init<std::string>())
		.def("getPluginName", &Properties::getPluginName, BP_RETURN_CONSTREF)
		.def("setPluginName", &Properties::setPluginName)
		.def("getID", &Properties::getPluginName, BP_RETURN_CONSTREF)
		.def("setID", &Properties::setPluginName)
		.def("getType", &Properties::getType)
		.def("getNames", &Properties::getNames)
		.def("hasProperty", &Properties::hasProperty)
		.def("wasQueried", &Properties::wasQueried)
		.def("markQueried", &Properties::markQueried)
		.def("__getitem__", &properties_wrapper::get)
		.def("__setitem__", &properties_wrapper::set)
		.def("__contains__", &Properties::hasProperty)
		.def("__str__", &Properties::toString);

	BP_SETSCOPE(properties);
	bp::enum_<Properties::EPropertyType>("EPropertyType")
		.value("EBoolean", Properties::EBoolean)
		.value("EInteger", Properties::EInteger)
		.value("EFloat", Properties::EFloat)
		.value("EPoint", Properties::EPoint)
		.value("ETransform", Properties::ETransform)
		.value("ESpectrum", Properties::ESpectrum)
		.value("EString", Properties::EString)
		.value("EData", Properties::EData)
		.export_values();
	BP_SETSCOPE(coreModule);

	BP_STRUCT(Vector2, bp::init<>())
		.def(bp::init<Float, Float>())
		.def(bp::init<Point2>())
		.def_readwrite("x", &Vector2::x)
		.def_readwrite("y", &Vector2::y);

	BP_STRUCT(Vector2i, bp::init<>())
		.def(bp::init<int, int>())
		.def(bp::init<Point2i>())
		.def_readwrite("x", &Vector2i::x)
		.def_readwrite("y", &Vector2i::y);

	BP_STRUCT(Vector3, bp::init<>())
		.def(bp::init<Float, Float, Float>())
		.def(bp::init<Point3>())
		.def(bp::init<Normal>())
		.def_readwrite("x", &Vector3::x)
		.def_readwrite("y", &Vector3::y)
		.def_readwrite("z", &Vector3::z);

	BP_STRUCT(Normal, bp::init<>())
		.def(bp::init<Float, Float, Float>())
		.def(bp::init<Vector>())
		.def_readwrite("x", &Normal::x)
		.def_readwrite("y", &Normal::y)
		.def_readwrite("z", &Normal::z);

	BP_STRUCT(Vector3i, bp::init<>())
		.def(bp::init<int, int, int>())
		.def(bp::init<Point3i>())
		.def_readwrite("x", &Vector3i::x)
		.def_readwrite("y", &Vector3i::y)
		.def_readwrite("z", &Vector3i::z);

	BP_STRUCT(Vector4, bp::init<>())
		.def(bp::init<Float, Float, Float, Float>())
		.def(bp::init<Point4>())
		.def_readwrite("x", &Vector4::x)
		.def_readwrite("y", &Vector4::y)
		.def_readwrite("z", &Vector4::z)
		.def_readwrite("w", &Vector4::w);

	BP_STRUCT(Vector4i, bp::init<>())
		.def(bp::init<int, int, int, int>())
		.def(bp::init<Point4i>())
		.def_readwrite("x", &Vector4i::x)
		.def_readwrite("y", &Vector4i::y)
		.def_readwrite("z", &Vector4i::z)
		.def_readwrite("w", &Vector4i::w);

	BP_STRUCT(Point2, bp::init<>())
		.def(bp::init<Float, Float>())
		.def(bp::init<Vector2>())
		.def_readwrite("x", &Point2::x)
		.def_readwrite("y", &Point2::y);

	BP_STRUCT(Point2i, bp::init<>())
		.def(bp::init<int, int>())
		.def(bp::init<Vector2i>())
		.def_readwrite("x", &Point2i::x)
		.def_readwrite("y", &Point2i::y);

	BP_STRUCT(Point3, bp::init<>())
		.def(bp::init<Float, Float, Float>())
		.def(bp::init<Vector3>())
		.def(bp::init<Normal>())
		.def_readwrite("x", &Point3::x)
		.def_readwrite("y", &Point3::y)
		.def_readwrite("z", &Point3::z);
	
	BP_STRUCT(Point3i, bp::init<>())
		.def(bp::init<int, int, int>())
		.def(bp::init<Vector3i>())
		.def_readwrite("x", &Point3i::x)
		.def_readwrite("y", &Point3i::y)
		.def_readwrite("z", &Point3i::z);

	BP_STRUCT(Point4, bp::init<>())
		.def(bp::init<Float, Float, Float, Float>())
		.def(bp::init<Vector4>())
		.def_readwrite("x", &Point4::x)
		.def_readwrite("y", &Point4::y)
		.def_readwrite("z", &Point4::z)
		.def_readwrite("w", &Point4::w);

	BP_STRUCT(Point4i, bp::init<>())
		.def(bp::init<int, int, int, int>())
		.def(bp::init<Vector4i>())
		.def_readwrite("x", &Point4i::x)
		.def_readwrite("y", &Point4i::y)
		.def_readwrite("z", &Point4i::z)
		.def_readwrite("w", &Point4i::w);

	BP_IMPLEMENT_VECTOR_OPS(Normal, Float, 3);
	BP_IMPLEMENT_VECTOR_OPS(Vector2i, int, 2);
	BP_IMPLEMENT_VECTOR_OPS(Vector3i, int, 3);
	BP_IMPLEMENT_VECTOR_OPS(Vector4i, int, 3);
	BP_IMPLEMENT_VECTOR_OPS(Vector2, Float, 2);
	BP_IMPLEMENT_VECTOR_OPS(Vector3, Float, 3);
	BP_IMPLEMENT_VECTOR_OPS(Vector4, Float, 3);
	BP_IMPLEMENT_POINT_OPS(Point2i, int, 2);
	BP_IMPLEMENT_POINT_OPS(Point3i, int, 3);
	BP_IMPLEMENT_POINT_OPS(Point4i, int, 3);
	BP_IMPLEMENT_POINT_OPS(Point2, Float, 2);
	BP_IMPLEMENT_POINT_OPS(Point3, Float, 3);
	BP_IMPLEMENT_POINT_OPS(Point4, Float, 3);

	Float (*dot2)(const Vector2 &, const Vector2 &) = &dot;
	Float (*dot3)(const Vector3 &, const Vector3 &) = &dot;
	Float (*dot4)(const Vector4 &, const Vector4 &) = &dot;
	Float (*absDot2)(const Vector2 &, const Vector2 &) = &absDot;
	Float (*absDot3)(const Vector3 &, const Vector3 &) = &absDot;
	Float (*absDot4)(const Vector4 &, const Vector4 &) = &absDot;
	Vector2 (*normalize2)(const Vector2 &) = &normalize;
	Vector3 (*normalize3)(const Vector3 &) = &normalize;
	Vector4 (*normalize4)(const Vector4 &) = &normalize;
	Vector3 (*cross3)(const Vector3 &, const Vector3 &) = &cross;
	bp::def("dot", dot2);
	bp::def("dot", dot3);
	bp::def("dot", dot4);
	bp::def("absDot", absDot2);
	bp::def("absDot", absDot3);
	bp::def("absDot", absDot4);
	bp::def("normalize", normalize2);
	bp::def("normalize", normalize3);
	bp::def("normalize", normalize4);
	bp::def("cross", cross3);

	bp::scope().attr("Vector") = bp::scope().attr("Vector3");
	bp::scope().attr("Point") = bp::scope().attr("Point3");

	bp::class_<Matrix4x4>("Matrix4x4", bp::init<Float>())
		.def(bp::init<Stream *>())
		.def("__setitem__", &Matrix4x4_setItem)
		.def("__getitem__", &Matrix4x4_getItem)
		.def("setIdentity", &Matrix4x4::setIdentity)
		.def("isZero", &Matrix4x4::isZero)
		.def("trace", &Matrix4x4::trace)
		.def("det", &Matrix4x4::det)
		.def("serialize", &Matrix4x4::serialize)
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self + Float())
		.def(bp::self += Float())
		.def(bp::self - bp::self)
		.def(bp::self -= bp::self)
		.def(bp::self - Float())
		.def(bp::self -= Float())
		.def(bp::self * Float())
		.def(Float() * bp::self)
		.def(bp::self *= Float())
		.def(bp::self * bp::self)
		.def(bp::self *= bp::self)
		.def(bp::self / Float())
		.def(bp::self /= Float())
		.def("__str__", &Matrix4x4::toString);

	bp::class_<Ray>("Ray", bp::init<>())
		.def(bp::init<Ray &>())
		.def(bp::init<Ray &, Float, Float>())
		.def(bp::init<Point, Vector, Float>())
		.def(bp::init<Point, Vector, Float, Float, Float>())
		.def_readwrite("o", &Ray::o)
		.def_readwrite("d", &Ray::d)
		.def_readwrite("dRcp", &Ray::dRcp)
		.def_readwrite("mint", &Ray::mint)
		.def_readwrite("maxt", &Ray::maxt)
		.def_readwrite("time", &Ray::time)
		.def("setOrigin", &Ray::setOrigin)
		.def("setDirection", &Ray::setDirection)
		.def("setTime", &Ray::setTime)
		.def("eval", &ray_eval, BP_RETURN_VALUE)
		.def("__str__", &Ray::toString);

	bp::class_<BSphere>("BSphere", bp::init<>())
		.def(bp::init<BSphere>())
		.def(bp::init<Point, Float>())
		.def(bp::init<Stream *>())
		.def_readwrite("center", &BSphere::center)
		.def_readwrite("radius", &BSphere::radius)
		.def("isEmpty", &BSphere::isEmpty)
		.def("expandBy", &BSphere::expandBy)
		.def("contains", &BSphere::contains)
		.def(bp::self == bp::self)
		.def(bp::self != bp::self)
		.def("rayIntersect", &bsphere_rayIntersect)
		.def("serialize", &BSphere::serialize)
		.def("__str__", &BSphere::toString);

	bp::class_<AABB>("AABB", bp::init<>())
		.def(bp::init<AABB>())
		.def(bp::init<Point>())
		.def(bp::init<Point, Point>())
		.def(bp::init<Stream *>())
		.def_readwrite("min", &AABB::min)
		.def_readwrite("max", &AABB::max)
		.def("getSurfaceArea", &AABB::getSurfaceArea)
		.def("getVolume", &AABB::getVolume)
		.def("getCorner", &AABB::getCorner, BP_RETURN_VALUE)
		.def("overlaps", &AABB::overlaps)
		.def("getCenter", &AABB::getCenter, BP_RETURN_VALUE)
		.def("reset", &AABB::reset)
		.def("clip", &AABB::clip)
		.def("isValid", &AABB::isValid)
		.def("expandBy", &aabb_expandby_aabb)
		.def("expandBy", &aabb_expandby_point)
		.def("distanceTo", &aabb_distanceto_aabb)
		.def("distanceTo", &aabb_distanceto_point)
		.def("squaredDistanceTo", &aabb_sqrdistanceto_aabb)
		.def("squaredDistanceTo", &aabb_sqrdistanceto_point)
		.def("contains", &aabb_contains_aabb)
		.def("contains", &aabb_contains_point)
		.def("getLargestAxis", &AABB::getLargestAxis)
		.def("getShortestAxis", &AABB::getShortestAxis)
		.def("getExtents", &AABB::getExtents, BP_RETURN_VALUE)
		.def(bp::self == bp::self)
		.def(bp::self != bp::self)
		.def("rayIntersect", &aabb_rayIntersect)
		.def("getBSphere", &AABB::getBSphere, BP_RETURN_VALUE)
		.def("serialize", &AABB::serialize)
		.def("__str__", &AABB::toString);
	
	bp::class_<Frame>("Frame", bp::init<>())
		.def(bp::init<Vector, Vector, Normal>())
		.def(bp::init<Vector, Vector, Vector>())
		.def(bp::init<Normal>())
		.def(bp::init<Stream *>())
		.def_readwrite("s", &Frame::s)
		.def_readwrite("t", &Frame::t)
		.def_readwrite("n", &Frame::n)
		.def("serialize", &Frame::serialize)
		.def("toLocal", &Frame::toLocal, BP_RETURN_VALUE)
		.def("toWorld", &Frame::toWorld, BP_RETURN_VALUE)
		.def("__str__", &Frame::toString)
		.def("cosTheta", &Frame::cosTheta)
		.def("sinTheta", &Frame::sinTheta)
		.def("sinTheta2", &Frame::sinTheta2)
		.def("tanTheta", &Frame::tanTheta)
		.def("sinPhi", &Frame::sinPhi)
		.def("cosPhi", &Frame::cosPhi)
		.def("sinPhi2", &Frame::sinPhi2)
		.def("cosPhi2", &Frame::cosPhi2)
		.staticmethod("cosTheta")
		.staticmethod("sinTheta")
		.staticmethod("sinTheta2")
		.staticmethod("tanTheta")
		.staticmethod("sinPhi")
		.staticmethod("cosPhi")
		.staticmethod("sinPhi2")
		.staticmethod("cosPhi2");
	
	bp::class_<Transform>("Transform", bp::init<>())
		.def(bp::init<Stream *>())
		.def(bp::init<const Matrix4x4 &>())
		.def(bp::init<const Matrix4x4 &, const Matrix4x4 &>())
		.def("inverse", &Transform::inverse, BP_RETURN_VALUE)
		.def("getMatrix", &Transform::getMatrix, BP_RETURN_VALUE)
		.def("getInverseMatrix", &Transform::getInverseMatrix, BP_RETURN_VALUE)
		.def("det3x3", &Transform::det3x3)
		.def("hasScale", &Transform::hasScale)
		.def("isIdentity", &Transform::isIdentity)
		.def("serialize", &Transform::serialize)
		.def("__mul__", &transform_mul_transform, BP_RETURN_VALUE)
		.def("__mul__", &transform_mul_point, BP_RETURN_VALUE)
		.def("__mul__", &transform_mul_vector, BP_RETURN_VALUE)
		.def("__mul__", &transform_mul_vector4, BP_RETURN_VALUE)
		.def("__mul__", &transform_mul_normal, BP_RETURN_VALUE)
		.def("__mul__", &transform_mul_ray, BP_RETURN_VALUE)
		.def("__str__", &Transform::toString)
		.def("translate", &Transform::translate, BP_RETURN_VALUE)
		.def("rotate", &Transform::rotate, BP_RETURN_VALUE)
		.def("scale", &Transform::scale, BP_RETURN_VALUE)
		.def("lookAt", &Transform::lookAt, BP_RETURN_VALUE)
		.def("perspective", &Transform::perspective, BP_RETURN_VALUE)
		.def("orthographic", &Transform::orthographic, BP_RETURN_VALUE)
		.def("glPerspective", &Transform::glPerspective, BP_RETURN_VALUE)
		.def("glFrustum", &Transform::glFrustum, BP_RETURN_VALUE)
		.def("glOrthographic", &Transform::glOrthographic, BP_RETURN_VALUE)
		.def("fromFrame", &Transform::fromFrame, BP_RETURN_VALUE)
		.staticmethod("translate")
		.staticmethod("rotate")
		.staticmethod("scale")
		.staticmethod("lookAt")
		.staticmethod("perspective")
		.staticmethod("orthographic")
		.staticmethod("glPerspective")
		.staticmethod("glFrustum")
		.staticmethod("glOrthographic")
		.staticmethod("fromFrame");

	/* Functions from utility.h */
	bp::def("fresnel", &fresnel);
	bp::def("fresnelDielectric", &fresnelDielectric);
	bp::def("fresnelConductor", &fresnelConductor, BP_RETURN_VALUE);
	bp::def("fresnelDiffuseReflectance", &fresnelDiffuseReflectance);

	bp::detail::current_scope = oldScope;
}

BOOST_PYTHON_MODULE(mitsuba) {
	bp::object package = bp::scope();
	package.attr("__path__") = "mitsuba";

	/* Automatically take care of the framework
	   initialization / shutdown */
	initializeFramework();
	atexit(shutdownFramework);

	export_core();
	export_render();
}

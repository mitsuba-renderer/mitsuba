#include "core.h"
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/appender.h>

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
}

void shutdownFramework() {
	/* Shutdown the core framework */
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
				return bp::object(props.getFloat(name));
			case Properties::EPoint:
				return bp::object(props.getFloat(name));
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
			props.setPoint(name, Point(extractVector()));
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
	
static ref<SerializableObject> instance_manager_getinstance(InstanceManager *manager, Stream *stream) {
	return manager->getInstance(stream);
}

void appender_logProgress(Appender *appender, Float progress, const std::string &name,
	const std::string &formatted, const std::string &eta) {
	appender->logProgress(progress, name, formatted, eta, NULL);
}
	
void mts_log(ELogLevel level, const std::string &msg) {
	bp::object traceback(bp::import("traceback"));
	bp::object extract_stack(traceback.attr("extract_stack"));
	bp::object top(extract_stack()[1]);
	Thread::getThread()->getLogger()->log(level, 
		NULL, bp::extract<const char *>(top[0]), 
		bp::extract<int>(top[1]),
		"%s: %s", (const char *) bp::extract<const char *>(top[2]), msg.c_str());
}

class AppenderWrapper : public Appender, public bp::wrapper<Appender> {
public:
	AppenderWrapper(PyObject *self) : m_self(self) { }

	void append(ELogLevel level, const std::string &text) {
		get_override("append")(level, text);
	}

	void logProgress(Float progress, const std::string &name,
		const std::string &formatted, const std::string &eta, 
		const void *ptr) {
		get_override("logProgress")(progress, name, formatted, eta);
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

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(fromLinearRGB_overloads, fromLinearRGB, 3, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(fromXYZ_overloads, fromXYZ, 3, 4)

void export_core() {
	boost::python::to_python_converter<
		fs::path, path_to_python_str>();

	bp::object coreModule(
		bp::handle<>(bp::borrowed(PyImport_AddModule("mitsuba.core"))));
	bp::scope().attr("core") = coreModule;
	PyObject *oldScope = bp::detail::current_scope;

	BP_SETSCOPE(coreModule);

	bp::enum_<ELogLevel>("ELogLevel")
		.value("ETrace", ETrace)
		.value("EDebug", EDebug)
		.value("EInfo", EInfo)
		.value("EWarn", EWarn)
		.value("EError", EError)
		.export_values();

	bp::def("Log", &mts_log);

	bp::class_<Object, ref<Object>, boost::noncopyable>("Object", bp::no_init)
		.def("getRefCount", &Object::getRefCount)
		.def("__str__", &Object::toString);

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

	BP_CLASS(SerializableObject, Stream, bp::no_init)
		.def("serialize", &SerializableObject::serialize);
	
	typedef Object* (Thread::*get_parent_type)();

	BP_CLASS(Thread, Object, bp::no_init)
		.def("getID", &Thread::getID)
		.def("getStackSize", &Thread::getStackSize)
		.def("setPriority", &Thread::setPriority)
		.def("getPriority", &Thread::getPriority)
		.def("setCritical", &Thread::setCritical)
		.def("getCritical", &Thread::getCritical)
		.def("setName", &Thread::setName)
		.def("getName", &Thread::getName, BP_RETURN_CONSTREF)
//		.def("getParent", (get_parent_type) &Thread::getParent)
//		.def("setLogger", &Thread::setLogger)
//		.def("getLogger", &Thread::getLogger)
//		.def("setFileResolver", &Thread::setFileResolver)
//		.def("getFileResolver", &Thread::getFileResolver)
		.def("getThread", &Thread::getThread, bp::return_internal_reference<>())
		.def("isRunning", &Thread::isRunning)
		.def("sleep", &Thread::sleep)
		.def("detach", &Thread::detach)
		.def("join", &Thread::join)
		.def("start", &Thread::start)
		.staticmethod("sleep")
		.staticmethod("getThread");

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

	BP_CLASS(InstanceManager, Object, bp::init<>())
		.def("serialize", &InstanceManager::serialize)
		.def("getInstance", &instance_manager_getinstance);

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

	bp::detail::current_scope = oldScope;
}

BOOST_PYTHON_MODULE(mitsuba) {
	bp::object package = bp::scope();
	package.attr("__path__") = "mitsuba";

	/* Basic STL containers */
	bp::class_<StringVector>("StringVector")
		.def(bp::vector_indexing_suite<StringVector>());
	bp::class_<StringMap>("StringMap")
		.def(bp::map_indexing_suite<StringMap>());

	/* Automatically take care of the framework
	   initialization / shutdown */
	initializeFramework();
	atexit(shutdownFramework);

	export_core();
}

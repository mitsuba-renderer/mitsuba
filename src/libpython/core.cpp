#include "core.h"
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/properties.h>

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

template <typename T> class fixedsize_wrapper {
public:
	typedef typename T::value_type value_type;

	static value_type get(const T &vec, int i) {
		if (i < 0 || i >= T::dim) {
			SLog(EError, "Index %i is out of range!", i);
			return 0.0f;
		}
		return vec[i];
	}
	
	static void set(T &vec, int i, value_type value) {
		if (i < 0 || i >= T::dim) 
			SLog(EError, "Index %i is out of range!", i);
		else
			vec[i] = value;
	}

	static int len(T &) {
		return T::dim;
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


void Matrix4x4_setItem(Matrix4x4 *matrix, bp::tuple tuple, Float value) {
	if (bp::len(tuple) != 2)
		SLog(EError, "Invalid matrix indexing operation, required a tuple of length 2");
	int i = bp::extract<int>(tuple[0]);
	int j = bp::extract<int>(tuple[1]);

	if (i < 0 || j < 0 || i >= 4 || j >= 4)
		SLog(EError, "Index (%i, %i) is out of bounds!", i, j);

	matrix->operator()(i, j) = value;
}

Float Matrix4x4_getItem(Matrix4x4 *matrix, bp::tuple tuple) {
	if (bp::len(tuple) != 2)
		SLog(EError, "Invalid matrix indexing operation, required a tuple of length 2");
	int i = bp::extract<int>(tuple[0]);
	int j = bp::extract<int>(tuple[1]);

	if (i < 0 || j < 0 || i >= 4 || j >= 4)
		SLog(EError, "Index (%i, %i) is out of bounds!", i, j);

	return matrix->operator()(i, j);
}
	
ref<SerializableObject> instance_manager_getinstance(InstanceManager *manager, Stream *stream) {
	return manager->getInstance(stream);
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(setString_overloads, setString, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getString_overloads, getString, 1, 2)

void export_core() {
	boost::python::to_python_converter<
		fs::path, path_to_python_str>();

	bp::object coreModule(
		bp::handle<>(bp::borrowed(PyImport_AddModule("mitsuba.core"))));
	bp::scope().attr("core") = coreModule;
	bp::scope scope = coreModule;

	bp::class_<Object, ref<Object>, boost::noncopyable>("Object", bp::no_init)
		.def("getRefCount", &Object::getRefCount)
		.def("__str__", &Object::toString);

	BP_DECLARE_CLASS(Stream, Object) stream("Stream", bp::no_init);

	bp::scope streamScope = stream;
	bp::enum_<Stream::EByteOrder>("EByteOrder")
		.value("EBigEndian", Stream::EBigEndian)
		.value("ELittleEndian", Stream::ELittleEndian)
		.value("ENetworkByteOrder", Stream::ENetworkByteOrder)
		.export_values();

	stream.def("setByteOrder", &Stream::setByteOrder)
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

	bp::scope scope2 = coreModule;

	BP_DECLARE_CLASS(FileStream, Stream) fstream("FileStream");
	bp::scope fstreamScope = fstream;
	bp::enum_<FileStream::EFileMode>("EFileMode")
		.value("EReadOnly", FileStream::EReadOnly)
		.value("EReadWrite", FileStream::EReadWrite)
		.value("ETruncWrite", FileStream::ETruncWrite)
		.value("ETruncReadWrite", FileStream::ETruncReadWrite)
		.value("EAppendWrite", FileStream::EAppendWrite)
		.value("EAppendReadWrite", FileStream::EAppendReadWrite)
		.export_values();

	fstream
		.def(bp::init<std::string, FileStream::EFileMode>())
		.def("getPath", &FileStream::getPath, 
				bp::return_value_policy<bp::copy_const_reference>())
		.def("open", &FileStream::open)
		.def("close", &FileStream::close)
		.def("remove", &FileStream::remove);

	bp::scope scope3 = coreModule;

	BP_DECLARE_CLASS(SerializableObject, Stream)
		.def("serialize", &SerializableObject::serialize);

	BP_DECLARE_CLASS(InstanceManager, Object)
		.def("serialize", &InstanceManager::serialize)
		.def("getInstance", &instance_manager_getinstance);

	BP_DECLARE_CLASS(ContinuousSpectrum, SerializableObject)
		.def("eval", &ContinuousSpectrum::eval)
		.def("average", &ContinuousSpectrum::average);
	
	BP_DECLARE_CLASS(InterpolatedSpectrum, ContinuousSpectrum)
		.def("append", &InterpolatedSpectrum::append)
		.def("clear", &InterpolatedSpectrum::clear)
		.def("zeroExtend", &ContinuousSpectrum::zeroExtend);

	bp::class_<Properties> properties("Properties");
	bp::scope propertiesScope = properties;
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

	properties
		.def(bp::init<std::string>())
		.def("getPluginName", &Properties::getPluginName,
			bp::return_value_policy<bp::copy_const_reference>())
		.def("setPluginName", &Properties::setPluginName)
		.def("getID", &Properties::getPluginName,
			bp::return_value_policy<bp::copy_const_reference>())
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

	bp::scope scope4 = coreModule;
	bp::class_<Vector2>("Vector2", bp::init<Float, Float>())
		.def(bp::init<Float>())
		.def(bp::init<Point2>())
		.def(bp::init<Stream *>())
		.def_readwrite("x", &Vector2::x)
		.def_readwrite("y", &Vector2::y)
		.def("length", &Vector2::length)
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self -= bp::self)
		.def(bp::self *= Float())
		.def(bp::self * Float())
		.def(bp::self / Float())
		.def(bp::self /= Float())
		.def("serialize", &Vector2::serialize)
		.def("__str__", &Vector2::toString)
		.def("__len__", &fixedsize_wrapper<Vector2>::len)
		.def("__getitem__", &fixedsize_wrapper<Vector2>::get)
		.def("__setitem__", &fixedsize_wrapper<Vector2>::set);

	bp::class_<Point2>("Point2", bp::init<Float, Float>())
		.def(bp::init<Float>())
		.def(bp::init<Vector2>())
		.def(bp::init<Stream *>())
		.def_readwrite("x", &Point2::x)
		.def_readwrite("y", &Point2::y)
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self *= Float())
		.def(bp::self * Float())
		.def(bp::self / Float())
		.def(bp::self /= Float())
		.def("serialize", &Point2::serialize)
		.def("__str__", &Point2::toString)
		.def("__len__", &fixedsize_wrapper<Point2>::len)
		.def("__getitem__", &fixedsize_wrapper<Point2>::get)
		.def("__setitem__", &fixedsize_wrapper<Point2>::set);

	bp::class_<Vector>("Vector", bp::init<Float, Float, Float>())
		.def(bp::init<Float>())
		.def(bp::init<Point>())
		.def(bp::init<Normal>())
		.def(bp::init<Stream *>())
		.def_readwrite("x", &Vector::x)
		.def_readwrite("y", &Vector::y)
		.def_readwrite("z", &Vector::z)
		.def("length", &Vector::length)
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self -= bp::self)
		.def(bp::self *= Float())
		.def(bp::self * Float())
		.def(bp::self / Float())
		.def(bp::self /= Float())
		.def("serialize", &Vector::serialize)
		.def("__str__", &Vector::toString)
		.def("__len__", &fixedsize_wrapper<Vector>::len)
		.def("__getitem__", &fixedsize_wrapper<Vector>::get)
		.def("__setitem__", &fixedsize_wrapper<Vector>::set);

	bp::class_<Normal>("Normal", bp::init<Float, Float, Float>())
		.def(bp::init<Float>())
		.def(bp::init<Vector>())
		.def(bp::init<Stream *>())
		.def_readwrite("x", &Normal::x)
		.def_readwrite("y", &Normal::y)
		.def_readwrite("z", &Normal::z)
		.def("length", &Normal::length)
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self -= bp::self)
		.def(bp::self *= Float())
		.def(bp::self * Float())
		.def(bp::self / Float())
		.def(bp::self /= Float())
		.def("serialize", &Normal::serialize)
		.def("__str__", &Normal::toString)
		.def("__len__", &fixedsize_wrapper<Normal>::len)
		.def("__getitem__", &fixedsize_wrapper<Normal>::get)
		.def("__setitem__", &fixedsize_wrapper<Normal>::set);

	bp::class_<Point>("Point", bp::init<Float, Float, Float>())
		.def(bp::init<Float>())
		.def(bp::init<Vector>())
		.def(bp::init<Stream *>())
		.def_readwrite("x", &Point::x)
		.def_readwrite("y", &Point::y)
		.def_readwrite("z", &Point::z)
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self *= Float())
		.def(bp::self * Float())
		.def(bp::self / Float())
		.def(bp::self /= Float())
		.def("serialize", &Point::serialize)
		.def("__str__", &Point::toString)
		.def("__len__", &fixedsize_wrapper<Point>::len)
		.def("__getitem__", &fixedsize_wrapper<Point>::get)
		.def("__setitem__", &fixedsize_wrapper<Point>::set);

	bp::class_<Vector4>("Vector4", bp::init<Float, Float, Float, Float>())
		.def(bp::init<Float>())
		.def(bp::init<Point4>())
		.def(bp::init<Stream *>())
		.def_readwrite("x", &Vector4::x)
		.def_readwrite("y", &Vector4::y)
		.def_readwrite("z", &Vector4::z)
		.def_readwrite("w", &Vector4::w)
		.def("length", &Vector4::length)
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self -= bp::self)
		.def(bp::self *= Float())
		.def(bp::self * Float())
		.def(bp::self / Float())
		.def(bp::self /= Float())
		.def("serialize", &Vector4::serialize)
		.def("__str__", &Vector4::toString)
		.def("__len__", &fixedsize_wrapper<Vector4>::len)
		.def("__getitem__", &fixedsize_wrapper<Vector4>::get)
		.def("__setitem__", &fixedsize_wrapper<Vector4>::set);

	bp::class_<Point4>("Point4", bp::init<Float, Float, Float, Float>())
		.def(bp::init<Float>())
		.def(bp::init<Vector4>())
		.def(bp::init<Stream *>())
		.def_readwrite("x", &Point4::x)
		.def_readwrite("y", &Point4::y)
		.def_readwrite("z", &Point4::z)
		.def_readwrite("w", &Point4::w)
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self *= Float())
		.def(bp::self * Float())
		.def(bp::self / Float())
		.def(bp::self /= Float())
		.def("serialize", &Point4::serialize)
		.def("__str__", &Point4::toString)
		.def("__len__", &fixedsize_wrapper<Point4>::len)
		.def("__getitem__", &fixedsize_wrapper<Point4>::get)
		.def("__setitem__", &fixedsize_wrapper<Point4>::set);

	bp::class_<Vector2i>("Vector2i", bp::init<int, int>())
		.def(bp::init<int>())
		.def(bp::init<Point2i>())
		.def(bp::init<Stream *>())
		.def_readwrite("x", &Vector2i::x)
		.def_readwrite("y", &Vector2i::y)
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self -= bp::self)
		.def(bp::self *= int())
		.def(bp::self * int())
		.def(bp::self / int())
		.def(bp::self /= int())
		.def("serialize", &Vector2i::serialize)
		.def("__str__", &Vector2i::toString)
		.def("__len__", &fixedsize_wrapper<Vector2i>::len)
		.def("__getitem__", &fixedsize_wrapper<Vector2i>::get)
		.def("__setitem__", &fixedsize_wrapper<Vector2i>::set);

	bp::class_<Point2i>("Point2i", bp::init<int, int>())
		.def(bp::init<int>())
		.def(bp::init<Vector2i>())
		.def(bp::init<Stream *>())
		.def_readwrite("x", &Point2i::x)
		.def_readwrite("y", &Point2i::y)
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self *= int())
		.def(bp::self * int())
		.def(bp::self / int())
		.def(bp::self /= int())
		.def("serialize", &Point2i::serialize)
		.def("__str__", &Point2i::toString)
		.def("__len__", &fixedsize_wrapper<Point2i>::len)
		.def("__getitem__", &fixedsize_wrapper<Point2i>::get)
		.def("__setitem__", &fixedsize_wrapper<Point2i>::set);

	bp::class_<Vector3i>("Vector3i", bp::init<int, int, int>())
		.def(bp::init<int>())
		.def(bp::init<Point3i>())
		.def(bp::init<Stream *>())
		.def_readwrite("x", &Vector3i::x)
		.def_readwrite("y", &Vector3i::y)
		.def_readwrite("z", &Vector3i::z)
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self -= bp::self)
		.def(bp::self *= int())
		.def(bp::self * int())
		.def(bp::self / int())
		.def(bp::self /= int())
		.def("serialize", &Vector3i::serialize)
		.def("__str__", &Vector3i::toString)
		.def("__len__", &fixedsize_wrapper<Vector3i>::len)
		.def("__getitem__", &fixedsize_wrapper<Vector3i>::get)
		.def("__setitem__", &fixedsize_wrapper<Vector3i>::set);

	bp::class_<Point3i>("Point3i", bp::init<int, int, int>())
		.def(bp::init<int>())
		.def(bp::init<Vector3i>())
		.def(bp::init<Stream *>())
		.def_readwrite("x", &Point3i::x)
		.def_readwrite("y", &Point3i::y)
		.def_readwrite("z", &Point3i::z)
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self *= int())
		.def(bp::self * int())
		.def(bp::self / int())
		.def(bp::self /= int())
		.def("serialize", &Point3i::serialize)
		.def("__str__", &Point3i::toString)
		.def("__len__", &fixedsize_wrapper<Point3i>::len)
		.def("__getitem__", &fixedsize_wrapper<Point3i>::get)
		.def("__setitem__", &fixedsize_wrapper<Point3i>::set);

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
		.def(bp::self *= Float())
		.def(bp::self * bp::self)
		.def(bp::self *= bp::self)
		.def(bp::self / Float())
		.def(bp::self /= Float())
		.def("__str__", &Matrix4x4::toString);
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

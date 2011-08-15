#include "mtspy.h"
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/transform.h>

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


void export_core() {
	bp::object coreModule(
		bp::handle<>(bp::borrowed(PyImport_AddModule("mtspy.core"))));
	bp::scope().attr("core") = coreModule;
	bp::scope coreScope = coreModule;

	bp::class_<Object, ref<Object>, boost::noncopyable>("Object", bp::no_init)
		.def("getRefCount", &Object::getRefCount)
		.def("__str__", &Object::toString);
	
	bp::class_<Stream, ref<Stream>, bp::bases<Object>, boost::noncopyable>("Stream", bp::no_init);

	bp::enum_<Stream::EByteOrder>("EByteOrder")
		.value("EBigEndian", Stream::EBigEndian)
		.value("ELittleEndian", Stream::ELittleEndian)
		.value("ENetworkByteOrder", Stream::ENetworkByteOrder)
		.export_values();

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

BOOST_PYTHON_MODULE(mtspy) {
	bp::object package = bp::scope();
	package.attr("__path__") = "mtspy";

	bp::def("initializeFramework", initializeFramework);
	bp::def("shutdownFramework", shutdownFramework);

	export_core();
}

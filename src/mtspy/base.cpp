#include <mitsuba/mitsuba.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/sched.h>

namespace boost{
	namespace python{
		template <typename T> T* get_pointer(mitsuba::ref<T> & p) {
			return p.get();
		}
		template <typename T> const T* get_pointer(const mitsuba::ref<T> & p) {
			return p.get();
		}
	}
}

#include <boost/python.hpp>

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

template <typename T> boost::python::str as_str(const T &t) {
	std::string str = t.toString();
	return boost::python::str(str.c_str(), str.length());
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

/* Support ref<..> smart pointers */
namespace boost {
	namespace python {
		template <typename T> struct pointee<mitsuba::ref<T> > {
			typedef T type;
		};	
	}
}

class A: public Object {
public:
	virtual void test() { cout << "A::test()" << endl; }
};

class B: public A {
public:
	virtual void test() { cout << "B::test()" << endl; }

	ref<A> getter() { return this; }
};


void export_core() {
	namespace bp = boost::python;

	bp::object coreModule(
		bp::handle<>(bp::borrowed(PyImport_AddModule("mtspy.core"))));
	bp::scope().attr("core") = coreModule;

	bp::scope coreScope = coreModule;

	bp::class_<Object, boost::noncopyable>("Object", bp::no_init)
		.def("getRefCount", &Object::getRefCount)
		.def("__str__", &Object::toString);


	bp::class_<A, bp::bases<Object>, boost::noncopyable >("A")
		.def("test", &A::test);
	
	bp::class_<B, bp::bases<A>, boost::noncopyable >("B")
		.def("test", &B::test).def("getter", &B::getter);
	
	bp::register_ptr_to_python< ref<Object> >();
	bp::register_ptr_to_python< ref<A> >();
	bp::register_ptr_to_python< ref<B> >();

	bp::class_<Vector2>("Vector2", bp::init<Float, Float>())
		.def(bp::init<Float>())
		.def(bp::init<Point2>())
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
		.def("__str__", as_str<Vector2>)
		.def("__len__", &fixedsize_wrapper<Vector2>::len)
		.def("__getitem__", &fixedsize_wrapper<Vector2>::get)
		.def("__setitem__", &fixedsize_wrapper<Vector2>::set);

	bp::class_<Point2>("Point2", bp::init<Float, Float>())
		.def(bp::init<Float>())
		.def(bp::init<Vector2>())
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
		.def("__str__", as_str<Point2>)
		.def("__len__", &fixedsize_wrapper<Point2>::len)
		.def("__getitem__", &fixedsize_wrapper<Point2>::get)
		.def("__setitem__", &fixedsize_wrapper<Point2>::set);

	bp::class_<Vector>("Vector", bp::init<Float, Float, Float>())
		.def(bp::init<Float>())
		.def(bp::init<Point>())
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
		.def("__str__", as_str<Vector>)
		.def("__len__", &fixedsize_wrapper<Vector>::len)
		.def("__getitem__", &fixedsize_wrapper<Vector>::get)
		.def("__setitem__", &fixedsize_wrapper<Vector>::set);

	bp::class_<Point>("Point", bp::init<Float, Float, Float>())
		.def(bp::init<Float>())
		.def(bp::init<Vector>())
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
		.def("__str__", as_str<Point>)
		.def("__len__", &fixedsize_wrapper<Point>::len)
		.def("__getitem__", &fixedsize_wrapper<Point>::get)
		.def("__setitem__", &fixedsize_wrapper<Point>::set);

	bp::class_<Vector4>("Vector4", bp::init<Float, Float, Float, Float>())
		.def(bp::init<Float>())
		.def(bp::init<Point4>())
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
		.def("__str__", as_str<Vector4>)
		.def("__len__", &fixedsize_wrapper<Vector4>::len)
		.def("__getitem__", &fixedsize_wrapper<Vector4>::get)
		.def("__setitem__", &fixedsize_wrapper<Vector4>::set);

	bp::class_<Point4>("Point4", bp::init<Float, Float, Float, Float>())
		.def(bp::init<Float>())
		.def(bp::init<Vector4>())
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
		.def("__str__", as_str<Point4>)
		.def("__len__", &fixedsize_wrapper<Point4>::len)
		.def("__getitem__", &fixedsize_wrapper<Point4>::get)
		.def("__setitem__", &fixedsize_wrapper<Point4>::set);

	bp::class_<Vector2i>("Vector2i", bp::init<int, int>())
		.def(bp::init<int>())
		.def(bp::init<Point2i>())
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
		.def("__str__", as_str<Vector2i>)
		.def("__len__", &fixedsize_wrapper<Vector2i>::len)
		.def("__getitem__", &fixedsize_wrapper<Vector2i>::get)
		.def("__setitem__", &fixedsize_wrapper<Vector2i>::set);

	bp::class_<Point2i>("Point2i", bp::init<int, int>())
		.def(bp::init<int>())
		.def(bp::init<Vector2i>())
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
		.def("__str__", as_str<Point2i>)
		.def("__len__", &fixedsize_wrapper<Point2i>::len)
		.def("__getitem__", &fixedsize_wrapper<Point2i>::get)
		.def("__setitem__", &fixedsize_wrapper<Point2i>::set);

	bp::class_<Vector3i>("Vector3i", bp::init<int, int, int>())
		.def(bp::init<int>())
		.def(bp::init<Point3i>())
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
		.def("__str__", as_str<Vector3i>)
		.def("__len__", &fixedsize_wrapper<Vector3i>::len)
		.def("__getitem__", &fixedsize_wrapper<Vector3i>::get)
		.def("__setitem__", &fixedsize_wrapper<Vector3i>::set);

	bp::class_<Point3i>("Point3i", bp::init<int, int, int>())
		.def(bp::init<int>())
		.def(bp::init<Vector3i>())
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
		.def("__str__", as_str<Point3i>)
		.def("__len__", &fixedsize_wrapper<Point3i>::len)
		.def("__getitem__", &fixedsize_wrapper<Point3i>::get)
		.def("__setitem__", &fixedsize_wrapper<Point3i>::set);
}

BOOST_PYTHON_MODULE(mtspy) {
	using namespace boost::python;

	object package = scope();
	package.attr("__path__") = "mtspy";

	def("initializeFramework", initializeFramework);
	def("shutdownFramework", shutdownFramework);

	export_core();
}

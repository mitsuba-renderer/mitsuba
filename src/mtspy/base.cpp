#include <mitsuba/mitsuba.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/sched.h>
#include <boost/python.hpp>

using namespace boost::python;
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

BOOST_PYTHON_MODULE(mtspy) {
	def("initializeFramework", initializeFramework);
	def("shutdownFramework", shutdownFramework);

	class_<ref<Object>, boost::noncopyable>("Object");
//		.def("getRefCount", &Object::getRefCount)
//		.def("toString", &Object::toString);

	class_<Vector2>("Vector2", init<Float, Float>())
		.def(init<Float>())
		.def(init<Point2>())
		.def_readwrite("x", &Vector2::x)
		.def_readwrite("y", &Vector2::y)
		.def("length", &Vector2::length)
		.def(self != self)
		.def(self == self)
		.def(-self)
		.def(self + self)
		.def(self += self)
		.def(self - self)
		.def(self -= self)
		.def(self *= Float())
		.def(self * Float())
		.def(self / Float())
		.def(self /= Float())
		.def("__str__", as_str<Vector2>)
		.def("__len__", &fixedsize_wrapper<Vector2>::len)
		.def("__getitem__", &fixedsize_wrapper<Vector2>::get)
		.def("__setitem__", &fixedsize_wrapper<Vector2>::set);

	class_<Point2>("Point2", init<Float, Float>())
		.def(init<Float>())
		.def(init<Vector2>())
		.def_readwrite("x", &Point2::x)
		.def_readwrite("y", &Point2::y)
		.def(self != self)
		.def(self == self)
		.def(-self)
		.def(self + self)
		.def(self += self)
		.def(self - self)
		.def(self *= Float())
		.def(self * Float())
		.def(self / Float())
		.def(self /= Float())
		.def("__str__", as_str<Point2>)
		.def("__len__", &fixedsize_wrapper<Point2>::len)
		.def("__getitem__", &fixedsize_wrapper<Point2>::get)
		.def("__setitem__", &fixedsize_wrapper<Point2>::set);

	class_<Vector>("Vector", init<Float, Float, Float>())
		.def(init<Float>())
		.def(init<Point>())
		.def_readwrite("x", &Vector::x)
		.def_readwrite("y", &Vector::y)
		.def_readwrite("z", &Vector::z)
		.def("length", &Vector::length)
		.def(self != self)
		.def(self == self)
		.def(-self)
		.def(self + self)
		.def(self += self)
		.def(self - self)
		.def(self -= self)
		.def(self *= Float())
		.def(self * Float())
		.def(self / Float())
		.def(self /= Float())
		.def("__str__", as_str<Vector>)
		.def("__len__", &fixedsize_wrapper<Vector>::len)
		.def("__getitem__", &fixedsize_wrapper<Vector>::get)
		.def("__setitem__", &fixedsize_wrapper<Vector>::set);

	class_<Point>("Point", init<Float, Float, Float>())
		.def(init<Float>())
		.def(init<Vector>())
		.def_readwrite("x", &Point::x)
		.def_readwrite("y", &Point::y)
		.def_readwrite("z", &Point::z)
		.def(self != self)
		.def(self == self)
		.def(-self)
		.def(self + self)
		.def(self += self)
		.def(self - self)
		.def(self *= Float())
		.def(self * Float())
		.def(self / Float())
		.def(self /= Float())
		.def("__str__", as_str<Point>)
		.def("__len__", &fixedsize_wrapper<Point>::len)
		.def("__getitem__", &fixedsize_wrapper<Point>::get)
		.def("__setitem__", &fixedsize_wrapper<Point>::set);

	class_<Vector4>("Vector4", init<Float, Float, Float, Float>())
		.def(init<Float>())
		.def(init<Point4>())
		.def_readwrite("x", &Vector4::x)
		.def_readwrite("y", &Vector4::y)
		.def_readwrite("z", &Vector4::z)
		.def_readwrite("w", &Vector4::w)
		.def("length", &Vector4::length)
		.def(self != self)
		.def(self == self)
		.def(-self)
		.def(self + self)
		.def(self += self)
		.def(self - self)
		.def(self -= self)
		.def(self *= Float())
		.def(self * Float())
		.def(self / Float())
		.def(self /= Float())
		.def("__str__", as_str<Vector4>)
		.def("__len__", &fixedsize_wrapper<Vector4>::len)
		.def("__getitem__", &fixedsize_wrapper<Vector4>::get)
		.def("__setitem__", &fixedsize_wrapper<Vector4>::set);

	class_<Point4>("Point4", init<Float, Float, Float, Float>())
		.def(init<Float>())
		.def(init<Vector4>())
		.def_readwrite("x", &Point4::x)
		.def_readwrite("y", &Point4::y)
		.def_readwrite("z", &Point4::z)
		.def_readwrite("w", &Point4::w)
		.def(self != self)
		.def(self == self)
		.def(-self)
		.def(self + self)
		.def(self += self)
		.def(self - self)
		.def(self *= Float())
		.def(self * Float())
		.def(self / Float())
		.def(self /= Float())
		.def("__str__", as_str<Point4>)
		.def("__len__", &fixedsize_wrapper<Point4>::len)
		.def("__getitem__", &fixedsize_wrapper<Point4>::get)
		.def("__setitem__", &fixedsize_wrapper<Point4>::set);

	class_<Vector2i>("Vector2i", init<int, int>())
		.def(init<int>())
		.def(init<Point2i>())
		.def_readwrite("x", &Vector2i::x)
		.def_readwrite("y", &Vector2i::y)
		.def(self != self)
		.def(self == self)
		.def(-self)
		.def(self + self)
		.def(self += self)
		.def(self - self)
		.def(self -= self)
		.def(self *= int())
		.def(self * int())
		.def(self / int())
		.def(self /= int())
		.def("__str__", as_str<Vector2i>)
		.def("__len__", &fixedsize_wrapper<Vector2i>::len)
		.def("__getitem__", &fixedsize_wrapper<Vector2i>::get)
		.def("__setitem__", &fixedsize_wrapper<Vector2i>::set);

	class_<Point2i>("Point2i", init<int, int>())
		.def(init<int>())
		.def(init<Vector2i>())
		.def_readwrite("x", &Point2i::x)
		.def_readwrite("y", &Point2i::y)
		.def(self != self)
		.def(self == self)
		.def(-self)
		.def(self + self)
		.def(self += self)
		.def(self - self)
		.def(self *= int())
		.def(self * int())
		.def(self / int())
		.def(self /= int())
		.def("__str__", as_str<Point2i>)
		.def("__len__", &fixedsize_wrapper<Point2i>::len)
		.def("__getitem__", &fixedsize_wrapper<Point2i>::get)
		.def("__setitem__", &fixedsize_wrapper<Point2i>::set);

	class_<Vector3i>("Vector3i", init<int, int, int>())
		.def(init<int>())
		.def(init<Point3i>())
		.def_readwrite("x", &Vector3i::x)
		.def_readwrite("y", &Vector3i::y)
		.def_readwrite("z", &Vector3i::z)
		.def(self != self)
		.def(self == self)
		.def(-self)
		.def(self + self)
		.def(self += self)
		.def(self - self)
		.def(self -= self)
		.def(self *= int())
		.def(self * int())
		.def(self / int())
		.def(self /= int())
		.def("__str__", as_str<Vector3i>)
		.def("__len__", &fixedsize_wrapper<Vector3i>::len)
		.def("__getitem__", &fixedsize_wrapper<Vector3i>::get)
		.def("__setitem__", &fixedsize_wrapper<Vector3i>::set);

	class_<Point3i>("Point3i", init<int, int, int>())
		.def(init<int>())
		.def(init<Vector3i>())
		.def_readwrite("x", &Point3i::x)
		.def_readwrite("y", &Point3i::y)
		.def_readwrite("z", &Point3i::z)
		.def(self != self)
		.def(self == self)
		.def(-self)
		.def(self + self)
		.def(self += self)
		.def(self - self)
		.def(self *= int())
		.def(self * int())
		.def(self / int())
		.def(self /= int())
		.def("__str__", as_str<Point3i>)
		.def("__len__", &fixedsize_wrapper<Point3i>::len)
		.def("__getitem__", &fixedsize_wrapper<Point3i>::get)
		.def("__setitem__", &fixedsize_wrapper<Point3i>::set);
}

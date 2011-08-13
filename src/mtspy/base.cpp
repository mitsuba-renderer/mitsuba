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
	class_<Vector>("Vector", init<Float, Float, Float>())
		.def(init<Float>())
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
}

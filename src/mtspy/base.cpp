#include <mitsuba/mitsuba.h>
#include <boost/python.hpp>

using namespace boost::python;
using namespace mitsuba;

std::string hello() {
	return "Hello world!";
}

BOOST_PYTHON_MODULE(mtspy) {
	def("hello", hello, "Return hello world...");
}

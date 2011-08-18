/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__MTSPY_H)
#define __MTSPY_H

#include <mitsuba/mitsuba.h>

#define BP_STRUCT(Name, Init) \
	bp::class_<Name> Name ##_struct(#Name, Init); \
	Name ##_struct

#define BP_CLASS(Name, Base, Init) \
	bp::class_<Name, ref<Name>, bp::bases<Base>, boost::noncopyable> Name ##_class(#Name, Init); \
	Name ##_class

#define BP_WRAPPED_CLASS(Name, Wrapper, Base, Init) \
	bp::class_<Name, ref<Wrapper>, bp::bases<Base>, boost::noncopyable> Name ##_class(#Name, Init); \
	Name ##_class

#define BP_IMPLEMENT_VECTOR_OPS(Name, Scalar, Size) \
	Name ##_struct \
		.def(bp::init<Stream *>()) \
		.def(bp::self != bp::self) \
		.def(bp::self == bp::self) \
		.def(-bp::self) \
		.def(bp::self + bp::self) \
		.def(bp::self += bp::self) \
		.def(bp::self - bp::self) \
		.def(bp::self -= bp::self) \
		.def(bp::self *= Scalar()) \
		.def(bp::self * Scalar()) \
		.def(Scalar() * bp::self) \
		.def(bp::self / Scalar()) \
		.def(bp::self /= Scalar()) \
		.def("serialize", &Name::serialize) \
		.def("__str__", &Name::toString) \
		.def("__len__", &FixedSizeSupport<Name, Scalar, Size>::len) \
		.def("__getitem__", &FixedSizeSupport<Name, Scalar, Size>::get) \
		.def("__setitem__", &FixedSizeSupport<Name, Scalar, Size>::set)

#define BP_IMPLEMENT_POINT_OPS(Name, Scalar, Size) \
	Name ##_struct \
		.def(bp::init<Stream *>()) \
		.def(bp::self != bp::self) \
		.def(bp::self == bp::self) \
		.def(Scalar() * bp::self) \
		.def(-bp::self) \
		.def(bp::self + Name::vector_type()) \
		.def(bp::self += Name::vector_type()) \
		.def(bp::self - Name::vector_type()) \
		.def(bp::self -= Name::vector_type()) \
		.def(bp::self - bp::self) \
		.def(bp::self *= Scalar()) \
		.def(bp::self * Scalar()) \
		.def(bp::self / Scalar()) \
		.def(bp::self /= Scalar()) \
		.def("serialize", &Name::serialize) \
		.def("__str__", &Name::toString) \
		.def("__len__", &FixedSizeSupport<Name, Scalar, Size>::len) \
		.def("__getitem__", &FixedSizeSupport<Name, Scalar, Size>::get) \
		.def("__setitem__", &FixedSizeSupport<Name, Scalar, Size>::set)
			
#define BP_SETSCOPE(value) do { \
		bp::detail::current_scope = value.ptr(); \
	} while (0);

#define BP_RETURN_CONSTREF bp::return_value_policy<bp::copy_const_reference>()
#define BP_RETURN_INTREF bp::return_internal_reference<>

namespace boost {
	namespace python {
		template <typename T> T* get_pointer(mitsuba::ref<T> & p) {
			return p.get();
		}
		template <typename T> const T* get_pointer(const mitsuba::ref<T> & p) {
			return p.get();
		}
	}
}

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>

namespace bp = boost::python;

/* Support ref<..> smart pointers */
namespace boost {
	namespace python {
		template <typename T> struct pointee< mitsuba::ref<T> > {
			typedef T type;
		};
	}
}

template <typename T, typename Scalar, int Size> class FixedSizeSupport {
public:
	static Scalar get(const T &value, int i) {
		using namespace mitsuba;

		if (i < 0 || i >= Size) {
			SLog(EError, "Index %i is out of range!", i);
			return 0.0f;
		}
		return value[i];
	}

	static void set(T &value, int i, Scalar arg) {
		using namespace mitsuba;

		if (i < 0 || i >= Size) 
			SLog(EError, "Index %i is out of range!", i);
		else
			value[i] = arg;
	}

	static int len(const T &) {
		return Size;
	}
};

typedef std::vector<std::string> StringVector;
typedef std::map<std::string, std::string> StringMap;

#endif /* __MTSPY_H */


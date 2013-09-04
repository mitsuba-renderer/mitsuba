/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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

#if !defined(__PYTHON_BASE_H)
#define __PYTHON_BASE_H

#include <mitsuba/mitsuba.h>

#if defined(_MSC_VER)
#pragma warning(disable : 4244) // 'return' : conversion from 'Py_ssize_t' to 'unsigned int', possible loss of data
#pragma warning(disable : 4267) // 'return' : conversion from 'size_t' to 'long', possible loss of data
#endif

#define BP_STRUCT(Name, Init) \
	bp::class_<Name> Name ##_struct(#Name, Init); \
	bp::register_ptr_to_python<Name*>(); \
	Name ##_struct

#define BP_SUBSTRUCT(Name, Base, Init) \
	bp::class_<Name, bp::bases<Base> > Name ##_struct(#Name, Init); \
	bp::register_ptr_to_python<Name*>(); \
	Name ##_struct

#define BP_CLASS(Name, Base, Init) \
	bp::class_<Name, ref<Name>, bp::bases<Base>, boost::noncopyable> Name ##_class(#Name, Init); \
	bp::register_ptr_to_python<Name*>(); \
	Name ##_class

#define BP_WRAPPED_CLASS(Name, Wrapper, Base, Init) \
	bp::class_<Name, ref<Wrapper>, bp::bases<Base>, boost::noncopyable> Name ##_class(#Name, Init); \
	bp::register_ptr_to_python<Name*>(); \
	bp::implicitly_convertible<ref<Wrapper>, ref<Name> >(); \
	Name ##_class

#define BP_IMPLEMENT_VECTOR_OPS(Name, Scalar, Size) \
	Name ##_struct \
		.def(bp::init<Stream *>()) \
		.def(bp::init<Name>()) \
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
		.def("length", &Name::length) \
		.def("lengthSquared", &Name::lengthSquared) \
		.def("__repr__", &Name::toString) \
		.def("__len__", &FixedSizeSupport<Name, Scalar, Size>::len) \
		.def("__getitem__", &FixedSizeSupport<Name, Scalar, Size>::get) \
		.def("__setitem__", &FixedSizeSupport<Name, Scalar, Size>::set)

#define BP_IMPLEMENT_POINT_OPS(Name, Scalar, Size) \
	Name ##_struct \
		.def(bp::init<Stream *>()) \
		.def(bp::init<Name>()) \
		.def(bp::self != bp::self) \
		.def(bp::self == bp::self) \
		.def(Scalar() * bp::self) \
		.def(-bp::self) \
		.def(bp::self + Name::VectorType()) \
		.def(bp::self += Name::VectorType()) \
		.def(bp::self - Name::VectorType()) \
		.def(bp::self -= Name::VectorType()) \
		.def(bp::self - bp::self) \
		.def(bp::self *= Scalar()) \
		.def(bp::self * Scalar()) \
		.def(bp::self / Scalar()) \
		.def(bp::self /= Scalar()) \
		.def("serialize", &Name::serialize) \
		.def("__repr__", &Name::toString) \
		.def("__len__", &FixedSizeSupport<Name, Scalar, Size>::len) \
		.def("__getitem__", &FixedSizeSupport<Name, Scalar, Size>::get) \
		.def("__setitem__", &FixedSizeSupport<Name, Scalar, Size>::set)

#define BP_SETSCOPE(value) do { \
		bp::detail::current_scope = value.ptr(); \
	} while (0);

#define BP_RETURN_CONSTREF bp::return_value_policy<bp::copy_const_reference>()
#define BP_RETURN_NONCONSTREF bp::return_value_policy<bp::copy_non_const_reference>()
#define BP_RETURN_VALUE bp::return_value_policy<bp::return_by_value>()
#define BP_RETURN_INTREF bp::return_internal_reference<>()

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
			return (Scalar) 0;
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

namespace mitsuba {
	class SerializableObject;
	class ConfigurableObject;
};

typedef std::vector<std::string> StringVector;
typedef std::vector<mitsuba::SerializableObject *> SerializableObjectVector;
typedef std::map<std::string, std::string> StringMap;

extern void export_core();
extern void export_render();
extern bp::object cast(mitsuba::ConfigurableObject *obj);

#endif /* __PYTHON_BASE_H */


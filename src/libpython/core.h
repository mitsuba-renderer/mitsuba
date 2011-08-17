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

#define BP_DECLARE_CLASS(Name, Base) \
	bp::class_<Name, ref<Name>, bp::bases<Base>, boost::noncopyable>

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

template <typename T> void registerClass() {
	boost::python::register_ptr_to_python< mitsuba::ref<T> >();
}

typedef std::vector<std::string> StringVector;
typedef std::map<std::string, std::string> StringMap;

#endif /* __MTSPY_H */


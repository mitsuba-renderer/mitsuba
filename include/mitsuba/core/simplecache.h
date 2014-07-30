/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

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

#pragma once
#if !defined(__MITSUBA_CORE_CACHE_H_)
#define __MITSUBA_CORE_CACHE_H_

#include <mitsuba/core/tls.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Generic thread-local storage for caching evaluations of
 * expensive function calls
 *
 * This class implements a simple and efficient one-entry storage for caching
 * the result of a call to an expensive-to-evaluate function that has no side
 * effects.
 *
 * The target application is a situation where multiple threads are causing
 * calls to such a function and each individual thread may invoke it a few
 * times in a row using the same input argument.
 *
 * This storage class provides the means to avoid re-evaluating the function
 * in this case. By isolating threads from one another, it does not suffer
 * heavy costs for inter-thread synchronizations.
 *
 * Here is an example snippet:
 *
 * \code
 * struct MyFunctor {
 *     inline void operator()(const Point &input, Float &output) {
 *          // .... Perform expensive function call / computation .....
 *     }
 * };
 *
 * void test() {
 *     SimpleCache<Point, Float> myCache;
 *     MyFunctor functor;
 *
 *     // First call -- evaluate the functor for the given input
 *     Float result = myCache.get(functor, Point(1,2,3));
 *
 *     // Now, the evaluation can uses the cached value
 *     Float result2 = myCache.get(functor, Point(1,2,3));
 * }
 *
 * \endcode
 *
 * \tparam ArgType
 *     Argument type of the function whose return values should be cached
 *
 * \tparam ReturnType
 *     Return type of the function whose return values should be cached
 */
template <typename ArgType, typename ReturnType> class SimpleCache
	: protected PrimitiveThreadLocal< std::pair<ArgType, ReturnType> > {
protected:
	typedef std::pair<ArgType, ReturnType>     ValueType;
	typedef PrimitiveThreadLocal<ValueType>    ParentType;
public:
	SimpleCache() : ParentType() { }

	/**
	 * \brief Return the cache entry for the argument \c argument
	 * or run \c UpdateFunctor to compute it
	 */
	template <typename UpdateFunctor> inline ReturnType &get(const UpdateFunctor &functor, const ArgType &argument) {
		bool existed;
		ValueType *value = (ValueType *) this->m_base.get(existed);

		if (EXPECT_NOT_TAKEN(!existed || value->first != argument)) {
			value->first = argument;
			functor(value->first, value->second);
		}

		return value->second;
	}
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_CACHE_H_ */

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
#if !defined(__MITSUBA_CORE_BRENT_H_)
#define __MITSUBA_CORE_BRENT_H_

#include <mitsuba/mitsuba.h>
#include <boost/function.hpp>

MTS_NAMESPACE_BEGIN

/**
 * \brief Brent's method nonlinear zero finder
 *
 * The implementation is transcribed from the Apache Commons
 * Java implementation. The supplied function is required to be
 * continuous, but not necessarily smooth.
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE BrentSolver {
public:
	/// Return value of \ref BrentSolver::solve()
	struct Result {
		bool success;
		size_t iterations;
		Float x;
		Float y;

		/// Create a new result instance
		inline Result(bool success, size_t iterations, Float x, Float y)
			: success(success), iterations(iterations), x(x), y(y) { }

		/// Return a string representation of the result
		inline std::string toString() const {
			std::ostringstream oss;
			oss << "BrentSolver::Result["
				<< "success=" << success << ", "
				<< "iterations=" << iterations << ", "
				<< "x=" << x << ", "
				<< "y=" << y << "]";
			return oss.str();
		}
	};

	/**
	 * \brief Create a new Brent-style solver with the
	 * specified accuracy requirements
	 *
	 * \param maxIterations
	 *      Max. number of successive iterations (default: 100)
	 * \param absAccuracy
	 *      Absolute accuracy requirement -- the iterations will stop
	 *      when |f(x)| < absAccuracy.
	 * \param absAccuracyPos
	 *      Absolute accuracy requirement of the position -- the
	 *      iterations will stop when |minX-maxX| < absAccuracyPos.
	 * \param absAccuracyPos
	 *      Absolute accuracy requirement of the position -- the
	 *      iterations will stop when |minX-maxX|/minX < relAccuracyPos.
	 */
	inline BrentSolver(size_t maxIterations = 100,
			Float absAccuracy = 1e-6f,
			Float absAccuracyPos = 1e-6f,
			Float relAccuracyPos = 1e-6f)
		: m_maxIterations(maxIterations),
		  m_absAccuracy(absAccuracy),
		  m_absAccuracyPos(absAccuracyPos),
		  m_relAccuracyPos(relAccuracyPos) { }

	/**
	 * \brief Find a zero in the given interval.
	 *
	 * Requires that the values of the function at the endpoints
	 * have opposite signs.
	 *
	 * \param min the lower bound for the interval.
	 * \param max the upper bound for the interval.
	 * \return the value where the function is zero
	 */
	Result solve(const boost::function<Float (Float)> &func,
			Float min, Float max) const;

	/**
	 * \brief Find a zero in the given interval with an initial guess
	 *
	 * Requires that the values of the function at the endpoints
	 * have opposite signs (note that it is allowed to have endpoints
	 * with the same sign if the initial point has opposite sign
	 * function-wise).
	 *
	 * \param min the lower bound for the interval.
	 * \param max the upper bound for the interval.
	 * \param initial the start value to use (must be set to min
	 *    if no initial point is known)
	 * \return the value where the function is zero
	 */
	Result solve(const boost::function<Float (Float)> &func,
			Float min, Float max, Float initial) const;

	/**
	 * Find a zero starting search according to the three provided points.
	 *
	 * \param x0 old approximation for the root
	 * \param y0 function value at the approximation for the root
	 * \param x1 last calculated approximation for the root
	 * \param y1 function value at the last calculated approximation
	 *    for the root
	 * \param x2 bracket point (must be set to x0 if no bracket point is
	 *   known, this will force starting with linear interpolation)
	 * \param y2 function value at the bracket point.
	 * \return the value where the function is zero
	 */
	Result solve(const boost::function<Float (Float)> &func,
				 Float x0, Float y0,
				 Float x1, Float y1,
				 Float x2, Float y2) const;
protected:
	size_t m_maxIterations;
	Float m_absAccuracy;
	Float m_absAccuracyPos;
	Float m_relAccuracyPos;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_BRENT_H_ */

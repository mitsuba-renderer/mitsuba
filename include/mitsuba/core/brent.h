/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__BRENT_H)
#define __BRENT_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Brent's method nonlinear zero finder
 *
 * The implementation is transcribed from the Apache Commons 
 * Java implementation. The function \a Func is required to be
 * continuous, but not necessarily smooth.
 */
template <typename T, typename Func> class BrentSolver {
public:
	typedef T value_type;

	/// Return value of \ref BrentSolver::solve()
	struct Result {
		bool success;
		int iterations;
		value_type x;
		value_type y;

		/// Create a new result instance
		inline Result(bool success, int iterations, value_type x, value_type y)
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
	inline BrentSolver(int maxIterations = 100, 
			value_type absAccuracy = 1e-6f,
			value_type absAccuracyPos = 1e-6f,
			value_type relAccuracyPos = 1e-6f) 
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
	Result solve(Func &f, value_type min, value_type max) const {
		// return the first endpoint if it is good enough
		value_type yMin = f(min);
		if (std::abs(yMin) <= m_absAccuracy)
			return Result(true, 0, min, yMin);

		// return the second endpoint if it is good enough
		value_type yMax = f(max);
		if (std::abs(yMax) <= m_absAccuracy)
			return Result(true, 0, max, yMax);

		value_type sign = yMin * yMax;
		if (sign > 0) {
			SLog(EWarn, "BrentSolver: Function values at the endpoints do not have different signs -- "
				"endpoints: [%f, %f], values: [%f, %f]", min, max, yMin, yMax);
			return Result(false, 0, 0, 0);
		} else {
			// solve using only the first endpoint as initial guess
			return solve(f, min, yMin, max, yMax, min, yMin);
		}
	}

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
	Result solve(Func &f, value_type min, value_type max, value_type initial) const {
		if (initial < min || initial > max) {
			SLog(EWarn, "BrentSolver: Invalid interval: lower=%f, initial=%f, upper=%f",
				min, max, initial);
			return Result(false, 0, 0, 0);
		}

		// return the initial guess if it is good enough
		value_type yInitial = f(initial);
		if (std::abs(yInitial) <= m_absAccuracy) 
			return Result(true, 0, initial, yInitial);

		// return the first endpoint if it is good enough
		value_type yMin = f(min);
		if (std::abs(yMin) <= m_absAccuracy)
			return Result(true, 0, min, yMin);

		// reduce interval if min and initial bracket the root
        if (yInitial * yMin < 0)
            return solve(f, min, yMin, initial, yInitial, min, yMin);

		// return the second endpoint if it is good enough
		value_type yMax = f(max);
		if (std::abs(yMax) <= m_absAccuracy)
			return Result(true, 0, max, yMax);

		// reduce interval if initial and max bracket the root
        if (yInitial * yMax < 0) 
            return solve(f, initial, yInitial, max, yMax, initial, yInitial);

		SLog(EWarn, "BrentSolver: Function values at the endpoints do not have different signs -- "
			"endpoints: [%f, %f], values: [%f, %f]", min, max, yMin, yMax);

		return Result(false, 0, 0, 0);
	}

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
	Result solve(Func &f,
				 value_type x0, value_type y0,
				 value_type x1, value_type y1,
				 value_type x2, value_type y2) const {
		value_type delta = x1 - x0;
		value_type oldDelta = delta;

		int i = 0;
		while (i < m_maxIterations) {
			if (std::abs(y2) < std::abs(y1)) {
				// use the bracket point if is better than last approximation
				x0 = x1;
				x1 = x2;
				x2 = x0;
				y0 = y1;
				y1 = y2;
				y2 = y0;
			}
			if (std::abs(y1) <= m_absAccuracy) {
				// Avoid division by very small values. Assume
				// the iteration has converged (the problem may
				// still be ill conditioned)
				return Result(true, i, x1, y1);
			}
			value_type dx = x2 - x1;
			value_type tolerance =
				std::max(m_relAccuracyPos * std::abs(x1), m_absAccuracyPos);

			if (std::abs(dx) <= tolerance) 
				return Result(true, i, x1, y1);
			if ((std::abs(oldDelta) < tolerance) ||
					(std::abs(y0) <= std::abs(y1))) {
				// Force bisection.
				delta = (value_type) 0.5f * dx;
				oldDelta = delta;
			} else {
				value_type r3 = y1 / y0;
				value_type p;
				value_type p1;
				// the equality test (x0 == x2) is intentional,
				// it is part of the original Brent's method,
				// it should NOT be replaced by proximity test
				if (x0 == x2) {
					// Linear interpolation.
					p = dx * r3;
					p1 = 1 - r3;
				} else {
					// Inverse quadratic interpolation.
					value_type r1 = y0 / y2;
					value_type r2 = y1 / y2;
					p = r3 * (dx * r1 * (r1 - r2) - (x1 - x0) * (r2 - 1));
					p1 = (r1 - 1) * (r2 - 1) * (r3 - 1);
				}
				if (p > 0) {
					p1 = -p1;
				} else {
					p = -p;
				}
				if (2 * p >= (value_type) 1.5f * dx * p1 - std::abs(tolerance * p1) ||
						p >= std::abs((value_type) 0.5f * oldDelta * p1)) {
					// Inverse quadratic interpolation gives a value
					// in the wrong direction, or progress is slow.
					// Fall back to bisection.
					delta = (value_type) 0.5f * dx;
					oldDelta = delta;
				} else {
					oldDelta = delta;
					delta = p / p1;
				}
			}
			// Save old X1, Y1
			x0 = x1;
			y0 = y1;
			// Compute new X1, Y1
			if (std::abs(delta) > tolerance) {
				x1 = x1 + delta;
			} else if (dx > 0) {
				x1 = x1 + (value_type) 0.5f * tolerance;
			} else if (dx <= 0) {
				x1 = x1 - (value_type) 0.5f * tolerance;
			}
			y1 = f(x1);
			if ((y1 > 0) == (y2 > 0)) {
				x2 = x0;
				y2 = y0;
				delta = x1 - x0;
				oldDelta = delta;
			}
			i++;
		}
		SLog(EWarn, "BrentSolver: Maximum number of iterations (%i) exceeded!",
			m_maxIterations);
		return Result(false, i, x1, y1);
	}
private:
	int m_maxIterations;
	value_type m_absAccuracy;
	value_type m_absAccuracyPos;
	value_type m_relAccuracyPos;
}; 

MTS_NAMESPACE_END

#endif /* __BRENT_H */

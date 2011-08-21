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

#if !defined(__SPLINE_H)
#define __SPLINE_H

#include <mitsuba/core/serialization.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Simple natural cubic spline interpolation
 *
 * \ingroup libcore
 */
class MTS_EXPORT_CORE CubicSpline : public SerializableObject {
public:
	/**
	 * \brief Initialize a cubic spline with the given set of 
	 * points (must be in ascending order, with no duplicates)
	 */
	inline CubicSpline(std::vector<Float> &x, std::vector<Float> &y)
		: m_x(x), m_y(y) { }

	/**
	 * \brief Initialize an empty cubic spline and reserve memory
	 * for the specified number of points
	 */
	inline CubicSpline(size_t nPoints) {
		m_x.reserve(nPoints);
		m_y.reserve(nPoints);
		m_deriv.reserve(nPoints);
	}

	/**
	 * \brief Unserialize a cubic spline from a
	 * binary data stream
	 */
	CubicSpline(Stream *stream, InstanceManager *manager);

	/**
	 * \brief Append a point at the end -- must be called with 
	 * increasing values of \a x
	 */
	inline void append(Float x, Float y) {
		m_x.push_back(x);
		m_y.push_back(y);
	}

	/// Return the number of control points
	inline size_t getSize() const { return m_x.size(); }

	/// Clear the internal representation
	void clear();

	/// Compute the derivatives -- must be called prior to \ref eval
	void build();

	/// Evaluate the spline at \c x
	Float eval(Float x) const;

	/// Serialize this spline to a binary data stream
	void serialize(Stream *stream, InstanceManager *manager) const;

	/// Return a human-readable representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
private:
	std::vector<Float> m_x, m_y;
	std::vector<Float> m_deriv; /// 2nd derivatives
};

MTS_NAMESPACE_END

#endif /* __SPLINE_H */

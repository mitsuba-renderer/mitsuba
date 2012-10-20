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

#pragma once
#if !defined(__MITSUBA_CORE_VMF_H_)
#define __MITSUBA_CORE_VMF_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Von Mises-Fisher distribution on the 2-sphere
 *
 * This is a basic implementation, which assumes that the
 * distribution is centered around the Z-axis.
 *
 * \author Wenzel Jakob
 * \ingroup libcore
 */
struct MTS_EXPORT_CORE VonMisesFisherDistr {
public:
	/**
	 * \brief Create a new von Mises-Fisher distribution
	 * with the given concentration parameter
	 */
	VonMisesFisherDistr(Float kappa = 0);

	/// Return the concentration parameter kappa
	inline Float getKappa() const {
		return m_kappa;
	}

	/// Evaluate the distribution for a given value of cos(theta)
	Float eval(Float cosTheta) const;

	/**
	 * \brief Generate a sample from this distribution
	 *
	 * \param sample
	 *     A uniformly distributed point on <tt>[0,1]^2</tt>
	 */
	Vector sample(const Point2 &sample) const;

	/**
	 * \brief Compute an appropriate concentration parameter so that
	 * the associated vMF distribution takes on the value \c x at its peak
	 */
	static Float forPeakValue(Float x);

	/**
	 * \brief Compute an concentration parameter that approximately
	 * corresponds to the spherical convolution of two vMF distributions.
	 *
	 * For details, see "Directional Statistics" by Mardia and Jupp, p.44
	 */
	static Float convolve(Float kappa1, Float kappa2);
private:
	Float m_kappa;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_VMF_H_ */

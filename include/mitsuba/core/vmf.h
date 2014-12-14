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
#if !defined(__MITSUBA_CORE_VMF_H_)
#define __MITSUBA_CORE_VMF_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Von Mises-Fisher distribution on the 2-sphere
 *
 * This is a basic implementation, which assumes that the
 * distribution is centered around the Z-axis. All provided
 * functions are implemented in such a way that they avoid
 * issues with numerical overflow.
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
	explicit inline VonMisesFisherDistr(Float kappa = 0) : m_kappa(kappa) { }

	/// Return the concentration parameter kappa
	inline void setKappa(Float kappa) {
		m_kappa = kappa;
	}

	/// Return the concentration parameter kappa
	inline Float getKappa() const {
		return m_kappa;
	}

	/// Return the mean cosine of the distribution
	Float getMeanCosine() const;

	/// Evaluate the distribution for a given value of cos(theta)
	Float eval(Float cosTheta) const;

	/**
	 * \brief Generate a sample from this distribution
	 *
	 * \param sample
	 *     A uniformly distributed point on <tt>[0,1]^2</tt>
	 */
	Vector sample(const Point2 &sample) const;

	/// Return a string representation
	std::string toString() const;

	/**
	 * \brief Compute an appropriate concentration parameter so that
	 * the associated vMF distribution takes on the value \c x at its peak
	 */
	static Float forPeakValue(Float x);

	/**
	 * \brief Estimate the vMF concentration parameter
	 * based on the length of the mean vector that is produced
	 * by simply averaging a set of sampled directions
	 *
	 * This is an unbiased estimator [Banerjee et al. 05]
	 */
	static Float forMeanLength(Float length);

	/**
	 * \brief Compute an appropriate concentration parameter so that
	 * the associated vMF distribution has the mean cosine \c g.
	 */
	static Float forMeanCosine(Float g);

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

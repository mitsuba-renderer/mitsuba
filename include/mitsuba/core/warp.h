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
#if !defined(__MITSUBA_CORE_WARP_H_)
#define __MITSUBA_CORE_WARP_H_

#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Implements common warping techniques that map from the unit
 * square to other domains, such as spheres, hemispheres, etc.
 *
 * The main application of this class is to generate uniformly
 * distributed or weighted point sets in certain common target domains.
 */
namespace warp {
	// =============================================================
	//! @{ \name Warping techniques related to spheres and subsets
	// =============================================================

	/// Uniformly sample a vector on the unit sphere with respect to solid angles
	extern MTS_EXPORT_CORE Vector squareToUniformSphere(const Point2 &sample);

	/// Density of \ref squareToUniformSphere() with respect to solid angles
	extern MTS_EXPORT_CORE inline Float squareToUniformSpherePdf() { return INV_FOURPI; }

	/// Uniformly sample a vector on the unit hemisphere with respect to solid angles
	extern MTS_EXPORT_CORE Vector squareToUniformHemisphere(const Point2 &sample);

	/// Density of \ref squareToUniformHemisphere() with respect to solid angles
	extern MTS_EXPORT_CORE inline Float squareToUniformHemispherePdf() { return INV_TWOPI; }

	/// Sample a cosine-weighted vector on the unit hemisphere with respect to solid angles
	extern MTS_EXPORT_CORE Vector squareToCosineHemisphere(const Point2 &sample);

	/// Density of \ref squareToCosineHemisphere() with respect to solid angles
	extern MTS_EXPORT_CORE inline Float squareToCosineHemispherePdf(const Vector &d)
		{ return INV_PI * Frame::cosTheta(d); }

	/**
	 * \brief Uniformly sample a vector that lies within a given
	 * cone of angles around the Z axis
	 *
	 * \param cosCutoff Cosine of the cutoff angle
	 * \param sample A uniformly distributed sample on \f$[0,1]^2\f$
	 */
	extern MTS_EXPORT_CORE Vector squareToUniformCone(Float cosCutoff, const Point2 &sample);

	/**
	 * \brief Uniformly sample a vector that lies within a given
	 * cone of angles around the Z axis
	 *
	 * \param cosCutoff Cosine of the cutoff angle
	 * \param sample A uniformly distributed sample on \f$[0,1]^2\f$
	 */
	extern MTS_EXPORT_CORE inline Float squareToUniformConePdf(Float cosCutoff) {
		return INV_TWOPI / (1-cosCutoff);
	}

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Warping techniques that operate in the plane
	// =============================================================

	/// Uniformly sample a vector on a 2D disk
	extern MTS_EXPORT_CORE Point2 squareToUniformDisk(const Point2 &sample);

	/// Density of \ref squareToUniformDisk per unit area
	extern MTS_EXPORT_CORE inline Float squareToUniformDiskPdf() { return INV_PI; }

	/// Low-distortion concentric square to disk mapping by Peter Shirley (PDF: 1/PI)
	extern MTS_EXPORT_CORE Point2 squareToUniformDiskConcentric(const Point2 &sample);

	/// Inverse of the mapping \ref squareToUniformDiskConcentric
	extern MTS_EXPORT_CORE Point2 uniformDiskToSquareConcentric(const Point2 &p);

	/// Density of \ref squareToUniformDisk per unit area
	extern MTS_EXPORT_CORE inline Float squareToUniformDiskConcentricPdf() { return INV_PI; }

	/// Convert an uniformly distributed square sample into barycentric coordinates
	extern MTS_EXPORT_CORE Point2 squareToUniformTriangle(const Point2 &sample);

	/**
	 * \brief Sample a point on a 2D standard normal distribution
	 *
	 * Internally uses the Box-Muller transformation
	 */
	extern MTS_EXPORT_CORE Point2 squareToStdNormal(const Point2 &sample);

	/// Density of \ref squareToStdNormal per unit area
	extern MTS_EXPORT_CORE Float squareToStdNormalPdf(const Point2 &pos);

	/// Warp a uniformly distributed square sample to a 2D tent distribution
	extern MTS_EXPORT_CORE Point2 squareToTent(const Point2 &sample);

	/**
	 * \brief Warp a uniformly distributed sample on [0, 1] to a nonuniform
	 * tent distribution with nodes <tt>{a, b, c}</tt>
	 */
	extern MTS_EXPORT_CORE Float intervalToNonuniformTent(Float a, Float b, Float c, Float sample);

	//! @}
	// =============================================================
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_WARP_H_ */

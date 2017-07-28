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
#if !defined(__MITSUBA_CORE_SPLINE_H_)
#define __MITSUBA_CORE_SPLINE_H_

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/*! \addtogroup libcore */
/*! @{ */

// -----------------------------------------------------------------------
//! @{ \name Functions for evaluating and sampling Catmull-Rom splines
// -----------------------------------------------------------------------

/**
 * \brief Evaluate a cubic spline interpolant of a uniformly sampled 1D function
 *
 * This implementation relies on Catmull-Rom splines, i.e. it uses finite
 * differences to approximate the derivatives at the endpoints of each spline
 * segment.
 *
 * \param x
 *      Evaluation point
 * \param values
 *      Floating point array containing \c size regularly spaced evaluations
 *      in the range [\c min,\c max] of the function to be approximated.
 * \param size
 *      Denotes the size of the \c values array
 * \param min
 *      Position of the first knot
 * \param max
 *      Position of the last knot
 * \param extrapolate
 *      Extrapolate values when \c x is out of range? (default: \c false)
 * \return
 *      The interpolated value or zero when <tt>extrapolate=false</tt>tt>
 *      and \c x lies outside of [\c min, \c max]
 */
extern MTS_EXPORT_CORE Float evalCubicInterp1D(Float x, const Float *values,
        size_t size, Float min, Float max, bool extrapolate = false);

/**
 * \brief Evaluate a cubic spline interpolant of a \a nonuniformly sampled 1D function
 *
 * This implementation relies on Catmull-Rom splines, i.e. it uses finite
 * differences to approximate the derivatives at the endpoints of each spline
 * segment.
 *
 * \param x
 *      Evaluation point
 * \param nodes
 *      Floating point array containing \c size nonuniformly spaced values
 *      denoting positions the where the function to be interpolated was evaluated.
 *      They must be provided in \a increasing order.
 * \param values
 *      Floating point array containing function evaluations matched to
 *      the entries of \c nodes.
 * \param size
 *      Denotes the size of the \c values array
 * \param extrapolate
 *      Extrapolate values when \c x is out of range? (default: \c false)
 * \return
 *      The interpolated value or zero when <tt>extrapolate=false</tt>tt>
 *      and \c x lies outside of \a [\c min, \c max]
 */
extern MTS_EXPORT_CORE Float evalCubicInterp1DN(Float x, const Float *nodes,
        const Float *values, size_t size, bool extrapolate = false);

/**
 * \brief Computes the definite integral over a segment of a uniformly
 * sampled 1D Catmull-Rom spline interpolant
 *
 * This is useful for sampling spline segments as part of an importance
 * sampling scheme (in conjunction with \ref sampleCubicInterp1D)
 *
 * \param idx
 *      Denotes the desires spline segment (must be between 0 and size-2)
  * \param values
 *      Floating point array containing \c size regularly spaced evaluations
 *      in the range [\c min,\c max] of the function to be approximated.
 * \param size
 *      Denotes the size of the \c values array
 * \param min
 *      Position of the first knot
 * \param max
 *      Position of the last knot
 * \return
 *      The definite integral over the specified segment
 */
extern MTS_EXPORT_CORE Float integrateCubicInterp1D(size_t idx,
        const Float *values, size_t size, Float min, Float max);

/**
 * \brief Computes the definite integral over a segment of a \a nonuniformly
 * sampled 1D Catmull-Rom spline interpolant
 *
 * This is useful for sampling spline segments as part of an importance
 * sampling scheme (in conjunction with \ref sampleCubicInterp1D)
 *
 * \param idx
 *      Denotes the desires spline segment (must be between 0 and size-2)
 * \param nodes
 *      Floating point array containing \c size nonuniformly spaced values
 *      denoting positions the where the function to be interpolated was evaluated.
 *      They must be provided in \a increasing order.
 * \param values
 *      Floating point array containing function evaluations matched to
 *      the entries of \c nodes.
 * \param size
 *      Denotes the size of the \c values array
 * \return
 *      The definite integral over the specified segment
 */
extern MTS_EXPORT_CORE Float integrateCubicInterp1DN(size_t idx,
        const Float *nodes, const Float *values, size_t size);

/**
 * \brief Importance sample a segment of a uniformly sampled 1D Catmull-Rom
 * spline interpolant
 *
 * \param idx
 *      Denotes the desires spline segment (must be between 0 and size-2)
  * \param values
 *      Floating point array containing \c size regularly spaced evaluations
 *      in the range [\c min,\c max] of the function to be approximated.
 * \param size
 *      Denotes the size of the \c values array
 * \param min
 *      Position of the first knot
 * \param max
 *      Position of the last knot
 * \param sample
 *      A uniformly distributed random sample in the interval <tt>[0,1]</tt>
 * \param fval
 *      If set to a non-NULL pointer, this argument will be used to return
 *      the value of the spline at the sampled position
 * \return
 *      The sampled position
 */
extern MTS_EXPORT_CORE Float sampleCubicInterp1D(size_t idx, const Float *values,
        size_t size, Float min, Float max, Float sample, Float *fval = NULL);

/**
 * \brief Importance sample a segment of a \a nonuniformly sampled 1D Catmull-Rom
 * spline interpolant
 *
 * \param idx
 *      Denotes the desires spline segment (must be between 0 and size-2)
 * \param nodes
 *      Floating point array containing \c size nonuniformly spaced values
 *      denoting positions the where the function to be interpolated was evaluated.
 *      They must be provided in \a increasing order.
 * \param values
 *      Floating point array containing function evaluations matched to
 *      the entries of \c nodes.
 * \param size
 *      Denotes the size of the \c values array
 * \param sample
 *      A uniformly distributed random sample in the interval <tt>[0,1]</tt>
 * \param fval
 *      If set to a non-NULL pointer, this argument will be used to return
 *      the value of the spline at the sampled position
 * \return
 *      The sampled position
 */
extern MTS_EXPORT_CORE Float sampleCubicInterp1DN(size_t idx, const Float *nodes,
        const Float *values, size_t size, Float sample, Float *fval = NULL);

/**
 * \brief Evaluate a cubic spline interpolant of a uniformly sampled 2D function
 *
 * This implementation relies on a tensor product of Catmull-Rom splines, i.e. it uses
 * finite differences to approximate the derivatives at the endpoints of each spline
 * patch.
 *
 * \param p
 *      Evaluation point
 * \param values
 *      A 2D floating point array of <tt>size.x*size.y</tt> cells containing regularly
 *      spaced evaluations of the function to be interpolated on the domain <tt>[min, max]</tt>.
 *      Consecutive entries of this array correspond to increments in the 'x' coordinate.
 * \param size
 *      Denotes the size of the \c values array (along each dimension)
 * \param min
 *      Position of the first knot on each dimension
 * \param max
 *      Position of the last knot on each dimension
 * \param extrapolate
 *      Extrapolate values when \c p is out of range? (default: \c false)
 * \return
 *      The interpolated value or zero when <tt>extrapolate=false</tt>tt> and
 *      \c p lies outside of the knot range
 */
extern MTS_EXPORT_CORE Float evalCubicInterp2D(const Point2 &p, const Float *values,
        const Size2 &size, const Point2 &min, const Point2 &max, bool extrapolate = false);

/**
 * \brief Evaluate a cubic spline interpolant of a \a nonuniformly sampled 2D function
 *
 * This implementation relies on a tensor product of Catmull-Rom splines, i.e. it uses
 * finite differences to approximate the derivatives at the endpoints of each spline
 * region.
 *
 * When the underlying function is sampled on a regular grid, \ref evalCubicInterp2D()
 * should be preferred, since value lookups will be considerably faster.
 *
 * \param p
 *      Evaluation point
 * \param nodes
 *      Pointer to a list for each dimension denoting the positions where the function
 *      to be interpolated was evaluated. The <tt>i</tt>-th array must have
 *      size <tt>size[i]</tt> and contain position values in \a increasing order.
 * \param values
 *      A 2D floating point array of <tt>size.x*size.y</tt> cells containing nonuniformly
 *      spaced evaluations of the function to be interpolated on the domain <tt>[min, max]</tt>.
 *      Consecutive entries of this array correspond to increments in the 'x' coordinate.
 * \param size
 *      Denotes the size of the \c values array (along each dimension)
 * \param extrapolate
 *      Extrapolate values when \c p is out of range? (default: \c false)
 * \return
 *      The interpolated value or zero when <tt>extrapolate=false</tt>tt> and
 *      \c p lies outside of the knot range
 */
extern MTS_EXPORT_CORE Float evalCubicInterp2DN(const Point2 &p, const Float **nodes,
        const Float *values, const Size2 &size, bool extrapolate = false);

/**
 * \brief Evaluate a cubic spline interpolant of a uniformly sampled 3D function
 *
 * This implementation relies on a tensor product of Catmull-Rom splines, i.e. it uses
 * finite differences to approximate the derivatives at the endpoints of each spline
 * region.
 *
 * \param p
 *      Evaluation point of the interpolant
 * \param values
 *      A 3D floating point array of <tt>size.x*size.y*size.z</tt> cells containing regularly
 *      spaced evaluations of the function to be interpolated on the domain <tt>[min, max]</tt>.
 *      Consecutive entries of this array correspond to increments in the 'x' coordinate,
 *      then 'y', and finally 'z' increments.
 * \param size
 *      Denotes the size of the \c values array (along each dimension)
 * \param min
 *      Position of the first knot on each dimension
 * \param max
 *      Position of the last knot on each dimension
 * \param extrapolate
 *      Extrapolate values when \c p is out of range? (default: \c false)
 * \return
 *      The interpolated value or zero when <tt>extrapolate=false</tt>tt> and
 *      \c p lies outside of the knot range
 */
extern MTS_EXPORT_CORE Float evalCubicInterp3D(const Point3 &p, const Float *values,
        const Size3 &size, const Point3 &min, const Point3 &max, bool extrapolate = false);

/**
 * \brief Evaluate a cubic spline interpolant of a \a nonuniformly sampled 3D function
 *
 * This implementation relies on a tensor product of Catmull-Rom splines, i.e. it uses
 * finite differences to approximate the derivatives at the endpoints of each spline
 * region.
 *
 * When the underlying function is sampled on a regular grid, \ref evalCubicInterp3D()
 * should be preferred, since value lookups will be considerably faster.
 *
 * \param p
 *      Evaluation point
 * \param nodes
 *      Pointer to a list for each dimension denoting the positions where the function
 *      to be interpolated was evaluated. The <tt>i</tt>-th array must have
 *      size <tt>size[i]</tt> and contain position values in \a increasing order.
 * \param values
 *      A 2D floating point array of <tt>size.x*size.y</tt> cells containing nonuniformly
 *      spaced evaluations of the function to be interpolated on the domain <tt>[min, max]</tt>.
 *      Consecutive entries of this array correspond to increments in the 'x' coordinate,
 *      then 'y', and finally 'z' increments.
 * \param size
 *      Denotes the size of the \c values array (along each dimension)
 * \param extrapolate
 *      Extrapolate values when \c p is out of range? (default: \c false)
 * \return
 *      The interpolated value or zero when <tt>extrapolate=false</tt>tt> and
 *      \c p lies outside of the knot range
 */
extern MTS_EXPORT_CORE Float evalCubicInterp3DN(const Point3 &p, const Float **nodes,
        const Float *values, const Size3 &size, bool extrapolate = false);

//! @}
// -----------------------------------------------------------------------

/*! @} */

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_SPLINE_H_ */

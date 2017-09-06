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
#if !defined(__MITSUBA_RENDER_COMMON_H_)
#define __MITSUBA_RENDER_COMMON_H_

#include <mitsuba/core/ray.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Specifies the transport mode when sampling or
 * evaluating a scattering function
 *
 * \ingroup librender
 */
enum ETransportMode {
    /* Note to self: do not change these enumeration
       values, some code depends on them. */

    /// Radiance transport
    ERadiance = 0,
    /// Importance transport
    EImportance = 1,
    /// Specifies the number of supported transport modes
    ETransportModes = 2
};

/**
 * \brief A list of measures that are associated with various sampling
 * methods in Mitsuba.
 *
 * Every now and then, sampling densities consist of a sum of several terms
 * defined on spaces that have different associated measures. In this case,
 * one of the constants in \ref EMeasure can be specified to clarify the
 * component in question when performing query operations.
 *
 * \ingroup librender
 */
enum EMeasure {
    /// Invalid measure
    EInvalidMeasure = 0,
    /// Solid angle measure
    ESolidAngle = 1,
    /// Length measure
    ELength = 2,
    /// Area measure
    EArea = 3,
    /// Discrete measure
    EDiscrete = 4
};

/**
 * \brief Generic sampling record for positions
 *
 * This sampling record is used to implement techniques that draw a position
 * from a point, line, surface, or volume domain in 3D and furthermore provide
 * auxilary information about the sample.
 *
 * Apart from returning the position and (optionally) the surface normal, the
 * responsible sampling method must annotate the record with the
 * associated probability density and measure.
 *
 * \ingroup librender
 */
struct MTS_EXPORT_RENDER PositionSamplingRecord {
public:
    /// Sampled position
    Point p;

    /// Associated time value
    Float time;

    /// Sampled surface normal (if applicable)
    Normal n;

    /// Probability density at the sample
    Float pdf;

    /**
     * \brief Denotes the measure associated with the sample.
     *
     * This is necessary to deal with quantities that are defined on
     * unusual spaces, e.g. areas that have collapsed to a point
     * or a line.
     */
    EMeasure measure;

    /**
     * \brief Optional: 2D sample position associated with the record
     *
     * In some uses of this record, a sampled position may be associated
     * with an important 2D quantity, such as the texture coordinates on
     * a triangle mesh or a position on the aperture of a sensor. When
     * applicable, such positions are stored in the \c uv attribute.
     */
    Point2 uv;

    /**
     * \brief Optional: Pointer to an associated object
     *
     * In some uses of this record, sampling a position also involves
     * choosing one of several objects (shapes, emitters, ..) on which
     * the position lies. In that case, the \c object attribute stores
     * a pointer to this object.
     */
    const ConfigurableObject *object;
public:
    /// Create an invalid position sampling record
    inline PositionSamplingRecord() { }

    /**
     * \brief Create a new position sampling record that can be
     * passed e.g. to \ref Shape::samplePosition
     *
     * \param time
     *    Specifies the time that should be associated with the
     *    position sample. This only matters when things are in motion
     */
    inline PositionSamplingRecord(Float time) : time(time),
        uv(0.0f), object(NULL) { }

    /**
     * \brief Create a position sampling record
     * from a surface intersection
     *
     * This is useful to determine the hypothetical sampling
     * density on a surface after hitting it using standard
     * ray tracing. This happens for instance in path tracing
     * with multiple importance sampling.
     */
    inline PositionSamplingRecord(const Intersection &its,
        EMeasure measure = EArea);

    /// Return a human-readable description of the record
    std::string toString() const;
};

/**
 * \brief Generic sampling record for directions
 *
 * This sampling record is used to implement techniques that randomly draw a
 * unit vector from a subset of the sphere and furthermore provide
 * auxilary information about the sample.
 *
 * Apart from returning the sampled direction, the responsible sampling method
 * must annotate the record with the associated probability density
 * and measure.
 *
 * \ingroup librender
 */
struct MTS_EXPORT_RENDER DirectionSamplingRecord {
public:
    /// Sampled direction
    Vector d;

    /// Probability density at the sample
    Float pdf;

    /// Measure associated with the density function
    EMeasure measure;

public:
    /// Return a human-readable description of the record
    std::string toString() const;

    /**
     * \brief Create an uninitialized position sampling record
     *
     * The resulting data structure is meant to be used
     * to generate a new direction sample.
     *
     * \sa Emitter::sampleDirection
     */
    inline DirectionSamplingRecord() { }

    /**
     * \brief Create a direction sampling record filled with a
     * specified direction.
     *
     * The resulting data structure is meant to be used to
     * query the density of a direction sampling technique
     *
     * \sa Emitter::pdfDirection
     */
    inline DirectionSamplingRecord(const Vector &d,
            EMeasure measure = ESolidAngle)
        : d(d), measure(measure) { }

    /**
     * \brief Create a direction sampling record
     * from a surface intersection
     *
     * This is useful to determine the hypothetical sampling
     * density of a direction after hitting it using standard
     * ray tracing. This happens for instance when hitting
     * the camera aperture in bidirectional rendering
     * techniques.
     */
    inline DirectionSamplingRecord(const Intersection &its,
        EMeasure measure = ESolidAngle);
};

/**
 * \brief Record for solid-angle based area sampling techniques
 *
 * This sampling record is used to implement techniques that randomly pick
 * a position on the surface of an object with the goal of importance sampling
 * a quantity that is defined over the sphere seen from a given reference point.
 *
 * This general approach for sampling positions is named "direct" sampling
 * throughout Mitsuba motivated by direct illumination rendering techniques,
 * which represent the most important application.
 *
 * This record inherits all fields from \ref PositionSamplingRecord and
 * extends it with two useful quantities that are cached so that they don't
 * need to be recomputed many times: the unit direction and length from the
 * reference position to the sampled point.
 *
 * \ingroup librender
 */
struct MTS_EXPORT_RENDER DirectSamplingRecord : public PositionSamplingRecord {
public:
    /// Reference point for direct sampling
    Point ref;

    /**
     * \brief Optional: normal vector associated with the reference point
     *
     * When nonzero, the direct sampling method can use the normal vector
     * to sample according to the projected solid angle at \c ref.
     */
    Normal refN;

    /// Unit direction from the reference point to the target direction
    Vector d;

    /// Distance from the reference point to the target direction
    Float dist;

public:
    /// Create an invalid direct sampling record
    inline DirectSamplingRecord() { }

    /**
     * \brief Create an new direct sampling record for a reference point
     * \c ref located somewhere in space (i.e. \a not on a surface)
     *
     * \param ref
     *     The reference point
     * \param time
     *     An associated time value
     */
    inline DirectSamplingRecord(const Point &ref, Float time)
        : PositionSamplingRecord(time), ref(ref), refN(0.0f) { }

    /**
     * \brief Create an new direct sampling record for a reference point
     * \c ref located on a surface.
     *
     * \param its
     *     The reference point specified using an intersection record
     */
    inline DirectSamplingRecord(const Intersection &refIts);

    /**
     * \brief Create an new direct sampling record for a reference point
     * \c ref located in a medium
     *
     * \param mRec
     *     The reference point specified using an medium sampling record
     */
    inline DirectSamplingRecord(const MediumSamplingRecord &mRec);

    /**
     * \brief Create a direct sampling record, which can be used to \a query
     * the density of a surface position (where there reference point lies on
     * a \a surface)
     *
     * \param ray
     *     Reference to the ray that generated the intersection \c its.
     *     The ray origin must be located at \c refIts.p
     *
     * \param its
     *     A surface intersection record (usually on an emitter)
     */

    inline void setQuery(
        const Ray &ray,
        const Intersection &its,
        EMeasure measure = ESolidAngle);

    /// Return a human-readable description of the record
    std::string toString() const;
};

/// \cond
extern MTS_EXPORT_RENDER std::ostream &operator<<(std::ostream &os, const ETransportMode &mode);
extern MTS_EXPORT_RENDER std::ostream &operator<<(std::ostream &os, const EMeasure &measure);
/// \endcond

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_COMMON_H_ */

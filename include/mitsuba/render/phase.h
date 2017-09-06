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
#if !defined(__MITSUBA_RENDER_PHASE_H_)
#define __MITSUBA_RENDER_PHASE_H_

#include <mitsuba/core/netobject.h>
#include <mitsuba/render/common.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Data structure, which contains information
 * required to sample or query a phase function.
 *
 * \ingroup librender
 */
struct MTS_EXPORT_RENDER PhaseFunctionSamplingRecord {
    /**
     * \brief Reference to a Medium sampling record created
     * by \ref Medium::sampleDistance()
     */
    const MediumSamplingRecord &mRec;

    /**
     * \brief Normalized incident direction vector, which points away
     * from the scattering event.
     *
     * In Mitsuba, the direction convention for phase functions is the
     * same as for BSDFs, as opposed to much of the literature, where
     * \c wi points inwards.
     */
    Vector wi;

    /// Normalized outgoing direction vector
    Vector wo;

    /* Transported mode (radiance or importance) -- required for
       rendering with non-reciprocal phase functions */
    ETransportMode mode;

    /**
     * \brief Given a medium interaction and an incident direction,
     * construct a query record which can be used to sample an outgoing
     * direction.
     *
     * \param mRec
     *      An reference to the underlying medium sampling record
     * \param wi
     *      An incident direction in world coordinates. This should
     *      be a normalized direction vector that points \a away from
     *      the scattering event.
     * \param mode
     *      The transported mode (\ref ERadiance or \ref EImportance)
     */

    inline PhaseFunctionSamplingRecord(const MediumSamplingRecord &mRec,
        const Vector &wi, ETransportMode mode = ERadiance)
        : mRec(mRec), wi(wi), mode(mode) { }

    /*
     * \brief Given a medium interaction an an incident/exitant direction
     * pair (wi, wo), create a query record to evaluate the phase function
     * or its sampling density.
     *
     * \param mRec
     *      An reference to the underlying medium sampling record
     * \param wi
     *      An incident direction in world coordinates. This should
     *      be a normalized direction vector that points \a away from
     *      the scattering event.
     * \param wo
     *      An outgoing direction in world coordinates. This should
     *      be a normalized direction vector that points \a away from
     *      the scattering event.
     * \param mode
     *      The transported mode (\ref ERadiance or \ref EImportance)
     */
    inline PhaseFunctionSamplingRecord(const MediumSamplingRecord &mRec,
        const Vector &wi, const Vector &wo, ETransportMode mode = ERadiance)
        : mRec(mRec), wi(wi),  wo(wo), mode(mode) { }

    /**
     * \brief Reverse the direction of light transport in the record
     *
     * This function essentially swaps \c wi and \c wo and adjusts
     * \c mode appropriately, so that non-symmetric scattering
     * models can be queried in the reverse direction.
     */
    inline void reverse() {
        std::swap(wo, wi);
        mode = (ETransportMode) (1-mode);
    }

    std::string toString() const;
};

/** \brief Abstract phase function.
 * \ingroup librender
 */
class MTS_EXPORT_RENDER PhaseFunction : public ConfigurableObject {
public:
    enum EPhaseFunctionType {
        /// Completely isotropic 1/(4 pi) phase function
        EIsotropic       = 0x01,
        /// The phase function only depends on \c dot(wi,wo)
        EAngleDependence = 0x04,
        /// The opposite of \ref EAngleDependence (there is an arbitrary dependence)
        EAnisotropic     = 0x02,
        /// The phase function is non symmetric, i.e. eval(wi,wo) != eval(wo, wi)
        ENonSymmetric    = 0x08
    };

    /**
     * \brief Return information flags of this phase function,
     * combined binary OR.
     * \sa EPhaseFunctionType
     */
    inline unsigned int getType() const {
        return m_type;
    }

    /// Configure the material (called after construction by the XML parser)
    virtual void configure();

    /**
     * \brief Evaluate the phase function for an outward-pointing
     * pair of directions (wi, wo)
     */
    virtual Float eval(const PhaseFunctionSamplingRecord &pRec) const = 0;

    /**
     * \brief Sample the phase function and return the importance weight (i.e. the
     * value of the phase function divided by the probability density of the sample).
     *
     * When the probability density is not explicitly required, this function
     * should be preferred, since it is potentially faster by making use of
     * cancellations during the division.
     *
     * \param pRec    A phase function query record
     * \param sampler A sample generator
     *
     * \return The phase function value divided by the probability
     *         density of the sample
     */
    virtual Float sample(PhaseFunctionSamplingRecord &pRec,
        Sampler *sampler) const = 0;

    /**
     * \brief Sample the phase function and return the probability density \a and the
     * importance weight of the sample (i.e. the value of the phase function divided
     * by the probability density)
     *
     * \param pRec    A phase function query record
     * \param sampler A sample generator
     * \param pdf     Will record the probability with respect to solid angles
     *
     * \return The phase function value divided by the probability
     *         density of the sample
     */
    virtual Float sample(PhaseFunctionSamplingRecord &pRec,
        Float &pdf, Sampler *sampler) const = 0;

    /**
     * \brief Calculate the probability of sampling wo (given wi).
     *
     * Assuming that the phase function can be sampled exactly,
     * the default implementation just evaluates \ref eval()
     */
    virtual Float pdf(const PhaseFunctionSamplingRecord &pRec) const;

    /**
     * \brief Does this phase function require directionally varying scattering
     * and extinction coefficients?
     *
     * This is used to implement rendering of media that have an anisotropic
     * structure (cf. "A radiative transfer framework for rendering materials with
     *   anisotropic structure" by Wenzel Jakob, Adam Arbree, Jonathan T. Moon,
     *   Kavita Bala, and Steve Marschner, SIGGRAPH 2010)
     */
    virtual bool needsDirectionallyVaryingCoefficients() const;

    /**
     * \brief For anisotropic media: evaluate the directionally varying component
     * of the scattering and absorption coefficients.
     *
     * \param cosTheta
     *    Angle between the axis of rotational symmetry and the
     *    direction of propagation
     */
    virtual Float sigmaDir(Float cosTheta) const;

    /**
     * \brief Returns the maximum value take on on by \ref sigmaDirMax().
     * This is useful when implementing Woodcock tracking.
     */
    virtual Float sigmaDirMax() const;

    /**
     * \brief Returns the mean cosine (often referred to by
     * the constant "g") of this phase function
     *
     * The default implementation throws an exception
     */
    virtual Float getMeanCosine() const;

    /// Return a string representation
    virtual std::string toString() const = 0;

    MTS_DECLARE_CLASS()
protected:
    /// Create a new phase function instance
    inline PhaseFunction(const Properties &props) :
        ConfigurableObject(props) { }

    /// Unserialize a phase function
    inline PhaseFunction(Stream *stream, InstanceManager *manager) :
        ConfigurableObject(stream, manager) { }

    /// Virtual destructor
    virtual ~PhaseFunction() { }
protected:
    unsigned int m_type;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_PHASE_H_ */

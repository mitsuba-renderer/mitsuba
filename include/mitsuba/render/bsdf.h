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
#if !defined(__MITSUBA_RENDER_BSDF_H_)
#define __MITSUBA_RENDER_BSDF_H_

#include <mitsuba/core/cobject.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/common.h>
#include <mitsuba/render/shader.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief This data structured contains all information
 * that is required to sample or query a BSDF.
 *
 * \sa BSDF::f()
 * \sa BSDF::sample()
 * \ingroup librender
 * \ingroup libpython
 */
struct MTS_EXPORT_RENDER BSDFSamplingRecord {
public:
    /**
     * \brief Given a surface interaction and an incident direction,
     * construct a query record which can be used to sample an outgoing
     * direction.
     *
     * By default, all components will be sampled irregardless of
     * what measure they live on. For convenience, this function
     * uses the local incident direction vector contained in the
     * supplied intersection record. The mode of transport is
     * set to \ref ERadiance -- the \ref mode fie
     *
     * \param its
     *      An reference to the underlying intersection record
     *
     * \param sampler
     *      A source of (pseudo-) random numbers. Note that this sampler
     *      is only used when the scattering model for some reason needs
     *      more than the two unformly distributed numbers supplied in
     *      the \ref BSDF::sample() methods
     *
     * \param mode
     *      The transported mode (\ref ERadiance or \ref EImportance)
     */
    explicit inline BSDFSamplingRecord(const Intersection &its, Sampler *sampler,
            ETransportMode mode = ERadiance);

    /**
     * \brief Given a surface interaction an an incident/exitant direction
     * pair (wi, wo), create a query record to evaluate the BSDF or its
     * sampling density.
     *
     * For convenience, this function uses the local incident direction
     * vector contained in the supplied intersection record.
     *
     * \param its
     *      An reference to the underlying intersection record
     * \param wo
     *      An outgoing direction in local coordinates. This should
     *      be a normalized direction vector that points \a away from
     *      the scattering event.
     * \param mode
     *      The transported mode (\ref ERadiance or \ref EImportance)
     */
    inline BSDFSamplingRecord(const Intersection &its, const Vector &wo,
        ETransportMode mode = ERadiance);

    /**
     * \brief Given a surface interaction an an incident/exitant direction
     * pair (wi, wo), create a query record to evaluate the BSDF or its
     * sampling density.
     *
     * \param its
     *      An reference to the underlying intersection record
     * \param wi
     *      An incident direction in local coordinates. This should
     *      be a normalized direction vector that points \a away from
     *      the scattering event.
     * \param wo
     *      An outgoing direction in local coordinates. This should
     *      be a normalized direction vector that points \a away from
     *      the scattering event.
     * \param mode
     *      The transported mode (\ref ERadiance or \ref EImportance)
     *
     */
    inline BSDFSamplingRecord(const Intersection &its,
        const Vector &wi, const Vector &wo,
        ETransportMode mode = ERadiance);

    /**
     * \brief Reverse the direction of light transport in the record
     *
     * This function essentially swaps \c wi and \c wo and adjusts
     * \c mode appropriately, so that non-symmetric scattering
     * models can be queried in the reverse direction.
     */
    inline void reverse();

    /// Return a string representation
    std::string toString() const;
public:
    /// Reference to the underlying surface interaction
    const Intersection &its;

    /**
     * \brief Pointer to a \ref Sampler instance (optional).
     *
     * Some BSDF implementations can significantly improve
     * the quality of their importance sampling routines
     * when having access to extra random numbers. This
     * attribute provides a means of providing this
     * capability to the BSDF.
     */
    Sampler *sampler;

    /**
     * \brief Normalized incident direction in local coordinates
     *
     * Mitsuba uses the convention that \c wi and \c wo
     * point away from the scattering event
     */
    Vector wi;

    /**
     * \brief Normalized outgoing direction in local coordinates
     *
     * Mitsuba uses the convention that \c wi and \c wo
     * point away from the scattering event
     */
    Vector wo;

    /// Relative index of refraction in the sampled direction
    Float eta;

    /** \brief Transported mode (radiance or importance)
     *
     * This information is required for rendering with non-reciprocal
     * BSDFs such as transmission through a dielectric material
     */
    ETransportMode mode;

    /**
     * \brief Bit mask containing the requested BSDF component types that
     * should be sampled/evaluated.
     *
     * Set to \c BSDF::EAll by default. After sampling has been performed,
     * the component type is stored inside \ref sampledType.
     *
     * \sa BSDF::EBSDFType
     */
    unsigned int typeMask;

    /**
     * \brief Integer value specifying the requested BSDF component index that
     * should be sampled/evaluated (for multi-lobed BSDFs).
     *
     * After sampling has been performed, the component index is stored
     * inside \ref sampledComponent.
     */
    int component;

    /**
     * \brief Stores the component type that was sampled by \ref BSDF::sample()
     * \sa BSDF::EBSDFType
     */
    unsigned int sampledType;

    /**
     * \brief Stores the component index that was sampled by \ref BSDF::sample()
     */
    int sampledComponent;
};


/**
 * \brief Abstract %BSDF base-class.
 *
 * This class implements an abstract interface to all BSDF plugins in Mitsuba.
 * It exposes functions for evaluating and sampling the model, and it allows
 * querying the probability density of the sampling method. Smooth
 * two-dimensional density functions, as well as degenerate one-dimensional
 * and discrete densities are all handled within the same framework.
 *
 * For improved flexibility with respect to the various rendering algorithms,
 * this class can sample and evaluate a complete BSDF, but it also allows to
 * pick and choose individual components of multi-lobed BSDFs based on their
 * properties and component indices. This selection is specified using a
 * special record that is provided along with every query.
 *
 * \ref BSDFSamplingRecord.
 * \ingroup librender
 * \ingroup libpython
 */
class MTS_EXPORT_RENDER BSDF : public ConfigurableObject, public HWResource {
public:
    /**
     * \brief This list of flags is used to classify the different
     * types of lobes that are implemented in a BSDF instance.
     *
     * They are also useful for picking out individual components
     * by setting combinations in \ref BSDFSamplingRecord::typeMask.
     */
    enum EBSDFType {
        // =============================================================
        //                      BSDF lobe types
        // =============================================================

        /// 'null' scattering event, i.e. particles do not undergo deflection
        ENull                 = 0x00001,
        /// Ideally diffuse reflection
        EDiffuseReflection    = 0x00002,
        /// Ideally diffuse transmission
        EDiffuseTransmission  = 0x00004,
        /// Glossy reflection
        EGlossyReflection     = 0x00008,
        /// Glossy transmission
        EGlossyTransmission   = 0x00010,
        /// Reflection into a discrete set of directions
        EDeltaReflection      = 0x00020,
        /// Transmission into a discrete set of directions
        EDeltaTransmission    = 0x00040,
        /// Reflection into a 1D space of directions
        EDelta1DReflection    = 0x00080,
        /// Transmission into a 1D space of directions
        EDelta1DTransmission  = 0x00100,

        // =============================================================
        //!                   Other lobe attributes
        // =============================================================
        /// The lobe is not invariant to rotation around the normal
        EAnisotropic          = 0x01000,
        /// The BSDF depends on the UV coordinates
        ESpatiallyVarying     = 0x02000,
        /// Flags non-symmetry (e.g. transmission in dielectric materials)
        ENonSymmetric         = 0x04000,
        /// Supports interactions on the front-facing side
        EFrontSide            = 0x08000,
        /// Supports interactions on the back-facing side
        EBackSide             = 0x10000,
        /// Uses extra random numbers from the supplied sampler instance
        EUsesSampler          = 0x20000
    };

    /// Convenient combinations of flags from \ref EBSDFType
    enum ETypeCombinations {
        /// Any reflection component (scattering into discrete, 1D, or 2D set of directions)
        EReflection   = EDiffuseReflection | EDeltaReflection
            | EDelta1DReflection | EGlossyReflection,
        /// Any transmission component (scattering into discrete, 1D, or 2D set of directions)
        ETransmission = EDiffuseTransmission | EDeltaTransmission
            | EDelta1DTransmission | EGlossyTransmission | ENull,
        /// Diffuse scattering into a 2D set of directions
        EDiffuse      = EDiffuseReflection | EDiffuseTransmission,
        /// Non-diffuse scattering into a 2D set of directions
        EGlossy       = EGlossyReflection | EGlossyTransmission,
        /// Scattering into a 2D set of directions
        ESmooth       = EDiffuse | EGlossy,
        /// Scattering into a discrete set of directions
        EDelta        = ENull | EDeltaReflection | EDeltaTransmission,
        /// Scattering into a 1D space of directions
        EDelta1D      = EDelta1DReflection | EDelta1DTransmission,
        /// Any kind of scattering
        EAll          = EDiffuse | EGlossy | EDelta | EDelta1D
    };

    /// Return the number of components of this BSDF
    inline int getComponentCount() const {
        return (int) m_components.size();
    }

    /**
     * \brief Return a listing of this BSDF's component types and
     * properties, combined using binary OR.
     * \sa EBSDFType
     */
    inline unsigned int getType() const {
        return m_combinedType;
    }

    /**
     * Returns the classification flags of a specific component
     * \sa EBSDFType
     */
    inline unsigned int getType(int component) const {
        return m_components[component];
    }

    /**
     * \brief Return the measure corresponding to a particular
     * component type
     */
    inline static EMeasure getMeasure(unsigned int componentType) {
        if (componentType & ESmooth) {
            return ESolidAngle;
        } else if (componentType & EDelta) {
            return EDiscrete;
        } else if (componentType & EDelta1D) {
            return ELength;
        } else {
            Log(EError, "getMeasure(): Invalid component type!");
            return ESolidAngle; // will never be reached
        }
    }

    /// Test whether this BSDF contains a certain type of component
    inline bool hasComponent(unsigned int type) const {
        return (type & m_combinedType) != 0;
    }

    /// Return whether this BSDF makes use of ray differentials
    inline bool usesRayDifferentials() const {
        return m_usesRayDifferentials;
    }

    /// Return the diffuse reflectance value (if any)
    virtual Spectrum getDiffuseReflectance(const Intersection &its) const;

    /// Return the specular reflectance value (if any)
    virtual Spectrum getSpecularReflectance(const Intersection &its) const {
        return Spectrum(0.0f);
    }

    /**
     * \brief Sample the BSDF and return the importance weight (i.e. the
     * value of the BSDF divided by the probability density of the sample).
     *
     * When the probability density is not explicitly required, this function
     * should be preferred, since it is potentially faster by making use of
     * cancellations during the division.
     *
     * If a component mask or a specific component index is specified, the
     * sample is drawn from the matching component, if it exists. Depending
     * on the provided transport type, either the BSDF or its adjoint version
     * is used.
     *
     * \param bRec    A BSDF query record
     * \param sample  A uniformly distributed sample on \f$[0,1]^2\f$
     *
     * \return The BSDF value divided by the probability density of the sample
     *         sample (multiplied by the cosine foreshortening factor when a
     *         non-delta component is sampled) A zero spectrum means that
     *         sampling failed.
     *
     * \remark This function is not exposed by the Python API. See the other
     *         sample function instead.
     *
     */
    virtual Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const = 0;

    /**
     * \brief Sample the BSDF and return the probability density \a and the
     * importance weight of the sample (i.e. the value of the BSDF divided
     * by the probability density)
     *
     * If a component mask or a specific component index is specified, the
     * sample is drawn from the matching component, if it exists. Depending
     * on the provided transport type, either the BSDF or its adjoint version
     * is used.
     *
     * When sampling a continuous/non-delta component, this method also
     * multiplies by the cosine foreshorening factor with respect to the
     * sampled direction.
     *
     * \param bRec    A BSDF query record
     * \param sample  A uniformly distributed sample on \f$[0,1]^2\f$
     * \param pdf     Will record the probability with respect to solid angles
     *                (or the discrete probability when a delta component is sampled)
     *
     * \return The BSDF value (multiplied by the cosine foreshortening
     *         factor when a non-delta component is sampled). A zero spectrum
     *         means that sampling failed.
     *
     * \remark From Python, this function is is called using the syntax
     *         <tt>value, pdf = bsdf.sample(bRec, sample)</tt>
     */
    virtual Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf,
        const Point2 &sample) const = 0;

    /**
     * \brief Evaluate the BSDF f(wi, wo) or its adjoint version f^{*}(wi, wo)
     *
     * This method allows to query the BSDF as a whole or pick out
     * individual components. When querying a smooth (i.e. non-degenerate)
     * component, it already multiplies the result by the cosine
     * foreshortening factor with respect to the outgoing direction.
     *
     * \param bRec
     *     A record with detailed information on the BSDF query
     *
     * \param measure
     *     Specifies the measure of the component. This is necessary
     *     to handle BSDFs, whose components live on spaces with
     *     different measures. (E.g. a diffuse material with an
     *     ideally smooth dielectric coating).
     */
    virtual Spectrum eval(const BSDFSamplingRecord &bRec,
        EMeasure measure = ESolidAngle) const = 0;

    /**
     * \brief Compute the probability of sampling \c bRec.wo (given
     * \c bRec.wi).
     *
     * This method provides access to the probability density that
     * would result when supplying the same BSDF query record to the
     * \ref sample() method. It correctly handles changes in probability
     * when only a subset of the components is chosen for sampling
     * (this can be done using the \ref BSDFSamplingRecord::component and
     * \ref BSDFSamplingRecord::typeMask fields).
     *
     * \param bRec
     *     A record with detailed information on the BSDF query
     *
     * \param measure
     *     Specifies the measure of the component. This is necessary
     *     to handle BSDFs, whose components live on spaces with
     *     different measures. (E.g. a diffuse material with an
     *     ideally smooth dielectric coating).
     */
    virtual Float pdf(const BSDFSamplingRecord &bRec,
        EMeasure measure = ESolidAngle) const = 0;

    /**
     * \brief For transmissive BSDFs: return the material's
     * relative index of refraction
     *
     * The default implementation returns <tt>1.0</tt>.
     *
     * \return interior IOR / exteriorIOR
     */
    virtual Float getEta() const;

    /**
     * \brief For rough BSDFs: return the root mean square
     * surface roughness of the given BSDF component
     *
     * An infinite value indicates a component that is ideally diffuse
     */
    virtual Float getRoughness(const Intersection &its, int index) const;

    /**
     * \brief Sometimes, BSDF models make use of a perturbed frame for
     * internal shading computations (e.g. bump maps). This function
     * exposes this internal frame.
     *
     * By default, it returns <tt>its.shFrame</tt>
     */
    virtual Frame getFrame(const Intersection &its) const;

    /**
     * \brief Sometimes, BSDF models make use of a perturbed frame for
     * internal shading computations (e.g. bump maps). This function
     * computes the derivative of this frame with respect to the UV
     * parameterization of the underlying shape.
     */
    virtual void getFrameDerivative(const Intersection &its,
            Frame &du, Frame &dv) const;

    // =============================================================
    //! @{ \name ConfigurableObject interface
    // =============================================================

    /// Configure the material (called after construction by the XML parser)
    virtual void configure();

    /// Add a child object
    virtual void addChild(const std::string &string, ConfigurableObject *obj);

    /// Add an unnamed child
    inline void addChild(ConfigurableObject *child) { addChild("", child); }

    /// Serialize this object to a stream
    virtual void serialize(Stream *stream, InstanceManager *manager) const;

    /// Set the parent object
    virtual void setParent(ConfigurableObject *parent);

    //! @}
    // =============================================================

    MTS_DECLARE_CLASS()
protected:
    /// Create a new BSDF instance
    BSDF(const Properties &props);

    /// Unserialize a BSDF instance
    BSDF(Stream *stream, InstanceManager *manager);

    /**
     * \brief Convenience function to ensure energy conservation
     *
     * This function determines the component-wise maximum of the
     * texture \c tex and checks if it is below \c max. If yes,
     * it returns the texture unmodified. Otherwise, it wraps
     * the texture into a \ref ScaleTexture instance (with a
     * scaling factor chosen so that the desired maximum \c max
     * is abided) and prints a warning.
     */
    Texture *ensureEnergyConservation(Texture *tex,
        const std::string &paramName, Float max) const;

    /**
     * \brief Convenience function to ensure energy conservation
     *
     * This function determines the component-wise maximum of the
     * sum \c tex1 + \c tex2 and checks if it is below \c max. If yes,
     * it returns the texture unmodified. Otherwise, it wraps
     * each the texture into a \ref ScaleTexture instance (with a
     * scaling factor chosen so that the desired maximum \c max
     * is abided) and prints a warning.
     */
    std::pair<Texture *, Texture *> ensureEnergyConservation(
        Texture *tex1, Texture *tex2, const std::string &paramName1,
        const std::string &paramName2, Float max) const;

    /// Virtual destructor
    virtual ~BSDF();
protected:
    std::vector<unsigned int> m_components;
    unsigned int m_combinedType;
    bool m_usesRayDifferentials;
    bool m_ensureEnergyConservation;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_BSDF_H_ */

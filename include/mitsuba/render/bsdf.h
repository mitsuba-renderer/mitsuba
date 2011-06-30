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

#if !defined(__BSDF_H)
#define __BSDF_H

#include <mitsuba/core/cobject.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/common.h>
#include <mitsuba/render/shader.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Data structure, which contains information 
 * required to sample or query a BSDF. 
 */
struct MTS_EXPORT_RENDER BSDFQueryRecord {
public:
	/**
	 * \brief Given a surface interaction and an incident direction, 
	 * construct a query record which can be used to sample
	 * an outgoing direction.
	 *
	 * For convenience, this function uses the local incident direction 
	 * vector contained in the supplied intersection record.
	 */
	inline BSDFQueryRecord(const Intersection &its);

	/**
	 * \brief Given a surface interaction an an incident/exitant direction 
	 * pair (wi, wo), create a BSDF query record to evaluate f(wi, wo). 
	 *
	 * For convenience, this function uses the local incident direction 
	 * vector contained in the supplied intersection record.
	 */
	inline BSDFQueryRecord(const Intersection &its, const Vector &wo);

	/**
	 * \brief Given a surface interaction an an incident/exitant direction 
	 * pair (wi, wo), create a BSDF query record to evaluate f(wi, wo).
	 */
	inline BSDFQueryRecord(const Intersection &its, 
		const Vector &wi, const Vector &wo); 

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
	 * attribute provides a means of providing this capability
	 * to the BSDF.
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

	/** \brief Transported quantity (radiance or importance)
	 * 
	 * This information is required for rendering with non-reciprocal 
	 * BSDFs such as transmission through a dielectric material
	 */
	ETransportQuantity quantity;

	/**
	 * \brief Bit mask containing the requested BSDF component types that 
	 * should be sampled/evaluated.
	 *
	 * Set to \c 0xFFFFFFFF by default. After sampling has been performed, 
	 * the component type is stored inside \ref sampledType.
	 * 
	 * \sa BSDF::EBSDFType
	 */
	unsigned int typeMask;

	/**
	 * \brief Stores the component type that was sampled by \ref BSDF::sample()
	 * \sa BSDF::EBSDFType
	 */
	unsigned int sampledType;

	/**
	 * \brief Integer value specifying the requested BSDF component index that 
	 * should be sampled/evaluated (for multi-lobed BSDFs).
	 *
	 * After sampling has been performed, the component index is stored
	 * inside \ref sampledComponent.
	 */
	int component;

	/**
	 * \brief Stores the component index that was sampled by \ref BSDF::sample()
	 */
	int sampledComponent;
};

/** 
 * \brief Abstract BSDF base-class.
 *
 * This class implements an abstract interface to all BSDF plugins in Mitsuba.
 * It exposes functions for evaluating and sampling the model, and it allows
 * querying the probability density of the sampling method. Both smooth
 * and Dirac delta densities are handled within the same framework.
 *
 * For improved flexibility with respect to the various rendering algorithms, 
 * this class can sample and evaluate a complete BRDF, but it also allows to
 * pick and choose individual components of multi-lobed BRDFs based on their 
 * properties and component indices.
 *
 * \sa BSDFQueryRecord
 */
class MTS_EXPORT_RENDER BSDF : public ConfigurableObject, public HWResource {
public:
	/**
	 * BSDF classification types, can be combined using binary OR.
	 */
	enum EBSDFType {
		EUnknown              = 0x00000,
		/// Perfect diffuse reflection 
		EDiffuseReflection    = 0x00001, 
		/// Perfect diffuse transmission
		EDiffuseTransmission  = 0x00002,
		/// Reflection using a delta function
		EDeltaReflection      = 0x00004,
		/// Transmission using a delta function
		EDeltaTransmission    = 0x00008,
		/// Glossy reflection
		EGlossyReflection     = 0x00010,
		/// Glossy transmission
		EGlossyTransmission   = 0x00020,
		/// Reflection is not invariant to rotation
		EAnisotropic          = 0x01000,
		/// Supports interactions on the front-facing side
		EFrontSide            = 0x02000,
		/// Supports interactions on the back-facing side
		EBackSide             = 0x04000,
		/// Can use a sampler instance to improve the sampling method
		ECanUseSampler        = 0x10000  
	};

	/// Type combinations
	enum ETypeCombinations {
		EDiffuse      = EDiffuseReflection | EDiffuseTransmission,
		EGlossy       = EGlossyReflection | EGlossyTransmission,
		EDelta        = EDeltaReflection | EDeltaTransmission,
		EReflection   = EDiffuseReflection | EDeltaReflection | EGlossyReflection,
		ETransmission = EDiffuseTransmission | EDeltaTransmission | EGlossyTransmission,
		ENonDelta     = EDiffuse | EGlossy,
		EAll          = EDiffuse | EGlossy | EDelta
	};

	/// Return the number of components of this BSDF
	inline int getComponentCount() const {
		return m_componentCount;
	}

	/**
	 * \brief Return a listing this BSDF's component types and
	 * properties, combined using binary OR.
	 * \sa EBSDFType
	 */
	inline unsigned int getType() const {
		return m_combinedType;
	}

	/**
	 * Returns the BSDF type for a specific component
	 * \sa EBSDFType
	 */
	inline unsigned int getType(int component) const {
		return m_type[component];
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
	virtual Spectrum getDiffuseReflectance(const Intersection &its) const = 0;

	/**
	 * \brief Sample the BSDF and divide by the probability of the sample. 
	 *
	 * When the probability density is not explicitly required, this function
	 * should be preferred, since it is potentially faster by making use of
	 * cancellations during the division.
	 * 
	 * If a component mask or a specific component index is given, the sample
	 * is drawn only from the matching component. Depending on the transport type
	 * either the BSDF or its adjoint version is used. 
	 *
	 * \return The BSDF value divided by the sample probability
	 */
	virtual Spectrum sample(BSDFQueryRecord &bRec, const Point2 &sample) const = 0;

	/**
	 * \brief Convenience method - similar to sample(), but also multiplies
	 * by the cosine foreshorening factor with respect to the sampled direction.
	 *
	 * \return The BSDF value multiplied by \c absDot(N, bRec.wo) and divided by
	 * the sample probability
	 */
	inline Spectrum sampleCos(BSDFQueryRecord &bRec, const Point2 &_sample) const {
		Spectrum bsdfVal = sample(bRec, _sample);
		if (bsdfVal.isZero())
			return bsdfVal; // bRec.wo is undefined, play safe
		return bsdfVal * std::abs(Frame::cosTheta(bRec.wo));
	}

	/**
	 * \c Sample the BSDF and explicitly provide the probability density
	 * of the sampled direction. 
	 *
	 * If a component mask or a specific component index is given, the 
	 * sample is drawn from the matching component. Depending on the 
	 * transport type, either the BSDF or its adjoint version is used. 
	 *
	 * \return The BSDF value (\a not divided by the probability)
	 */
	virtual Spectrum sample(BSDFQueryRecord &bRec, Float &pdf, 
		const Point2 &sample) const;

	/**
	 * \brief Convenience method - similar to sample(), but also multiplies
	 * by the cosine foreshorening factor with respect to the sampled direction.
	 *
	 * \return The BSDF value multiplied by \c absDot(N, bRec.wo) (and
	 * \a not divided by the sample probability)
	 */
	inline Spectrum sampleCos(BSDFQueryRecord &bRec, Float &pdf,
			const Point2 &_sample) const {
		Spectrum bsdfVal(sample(bRec, pdf, _sample));
		if (bsdfVal.isZero())
			return bsdfVal; // bRec.wo is undefined, play safe
		return bsdfVal * std::abs(Frame::cosTheta(bRec.wo));
	}

	/// Evaluate the BSDF f(wi, wo) or its adjoint version f^{*}(wi, wo)
	virtual Spectrum f(const BSDFQueryRecord &bRec) const = 0;

	/**
	 * \brief Evaluate the BSDF f(wi, wo) or its adjoint version f^{*}(wi, wo).
	 *
	 * Also multiplies by the cosine foreshorening factor with respect
	 * to the outgoing direction.
	 */
	inline Spectrum fCos(const BSDFQueryRecord &bRec) const  {
		return f(bRec) * std::abs(Frame::cosTheta(bRec.wo));
	}

	/// Calculate the probability of sampling wi (given wo) -- continuous version
	virtual Float pdf(const BSDFQueryRecord &bRec) const = 0;

	/// Calculate the probability of sampling wi (given wo) -- degenerate 0D (Dirac delta) version
	virtual Float pdfDelta(const BSDFQueryRecord &bRec) const;

	/// Evaluate the transfport from wi to wo -- degenerate 0D (Dirac delta) component
	virtual Spectrum fDelta(const BSDFQueryRecord &bRec) const;

	/// Return the name of this BSDF
	inline const std::string &getName() const { return m_name; }

	/// Set the name of this BSDF
	inline void setName(const std::string &name) { m_name = name; }

	/// Add a child object
	virtual void addChild(const std::string &string, ConfigurableObject *obj);

	/// Serialize this object to a stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Set the parent object
	virtual void setParent(ConfigurableObject *parent);

	MTS_DECLARE_CLASS()
protected:
	/// Create a new BSDF instance
	BSDF(const Properties &props);
	
	/// Unserialize a BSDF instance
	BSDF(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~BSDF();
protected:
	unsigned int *m_type;
	unsigned int m_combinedType;
	int m_componentCount;
	bool m_usesRayDifferentials;
	std::string m_name;
};

MTS_NAMESPACE_END

#endif /* __BSDF_H */

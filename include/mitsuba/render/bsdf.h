#if !defined(__BSDF_H)
#define __BSDF_H

#include <mitsuba/render/records.h>
#include <mitsuba/render/shader.h>

MTS_NAMESPACE_BEGIN

/** \brief Abstract BxDF base-class, which supports reflection and
 * transmission using a common interface for sampling and evaluation.
 * A BSDF can consist of several components, which can optionally be
 * independently sampled and/or evaluated.
 */
class MTS_EXPORT_RENDER BSDF : public ConfigurableObject, public HWResource {
public:
	/**
	 * BSDF classification types, can be combined using binary OR.
	 */
	enum EBSDFType {
		EUnknown              = 0x00,
		EDiffuseReflection    = 0x01, /* Perfect diffuse reflection */
		EDiffuseTransmission  = 0x02, /* Perfect diffuse transmission */
		EDeltaReflection      = 0x04, /* Reflection using a delta function */
		EDeltaTransmission    = 0x08, /* Transmission using a delta function */
		EGlossyReflection     = 0x10, /* Glossy reflection */
		EGlossyTransmission   = 0x20  /* Glossy transmission */
	};

	/// Type combinations
	enum ETypeCombinations {
		EDiffuse      = EDiffuseReflection | EDiffuseTransmission,
		EGlossy       = EGlossyReflection | EGlossyTransmission,
		EDelta        = EDeltaReflection | EDeltaTransmission,
		EReflection   = EDiffuseReflection | EDeltaReflection | EGlossyReflection,
		ETransmission = EDiffuseTransmission | EDeltaTransmission | EGlossyTransmission,
		ENonDelta     = EDiffuse | EGlossy,
		EAll          = EDiffuse | EDelta | EGlossy
	};

	/// Return the number of components of this BxDF
	inline int getComponentCount() const {
		return m_componentCount;
	}

	/// Return an (OR-ed) integer listing this BxDF's components
	inline unsigned int getType() const {
		return m_combinedType;
	}

	/// Returns the BxDF type for a specific component
	inline unsigned int getType(int component) const {
		return m_type[component];
	}

	/// Test whether this BSDF contains a certain type of component
	inline bool hasComponent(unsigned int type) const {
		return (type & m_combinedType) != 0;
	}

	/// Return whether this BxDF makes use of ray differentials
	inline bool usesRayDifferentials() const {
		return m_usesRayDifferentials;
	}

	/// Return the diffuse reflectance value (if any)
	virtual Spectrum getDiffuseReflectance(const Intersection &its) const = 0;

	/**
	 * Sample the BxDF and divide by the probability of the sample. If a
	 * component mask or a specific component index is given, the sample
	 * is drawn from the matching component. Depending on the transport type
	 * the BSDF or its adjoint version is used. 
	 * @return The BSDF value divided by the sample probability
	 */
	virtual Spectrum sample(BSDFQueryRecord &bRec) const = 0;
	
	/**
	 * Convenience method - similar to sample(), but also multiplies
	 * by the cosine factor of the sampled direction.
	 * @return The BSDF value divided by the sample probability
	 */
	inline Spectrum sampleCos(BSDFQueryRecord &bRec) const {
		Spectrum bsdfVal = sample(bRec);
		return bsdfVal * std::abs(Frame::cosTheta(bRec.wo));
	}

	/**
	 * Sample the BxDF and explicitly provided the probability density
	 * of the sampled direction. If a component mask or a specific 
	 * component index is given, the sample is drawn from the matching 
	 * component. Depending on the transport type the BSDF or its adjoint 
	 * version is used. 
	 * @return The BSDF value (not divided by the probability)
	 */
	virtual Spectrum sample(BSDFQueryRecord &bRec, Float &pdf) const;

	/**
	 * Convenience method - similar to sample(), but also multiplies
	 * by the cosine factor of the sampled direction.
	 * @return The BSDF value (not divided by the probability)
	 */
	inline Spectrum sampleCos(BSDFQueryRecord &bRec, Float &pdf) const {
		Spectrum bsdfVal(sample(bRec, pdf));
		return bsdfVal * std::abs(Frame::cosTheta(bRec.wo));
	}

	/// Evaluate the BxDF f(wi, wo) or its adjoint version f^{*}(wi, wo)
	virtual Spectrum f(const BSDFQueryRecord &bRec) const = 0;

	/**
	 * Evaluate the BxDF f(wi, wo) or its adjoint version f^{*}(wi, wo).
	 * Also multiplies by the cosine factor of the sampled direction.
	 */
	inline Spectrum fCos(const BSDFQueryRecord &bRec) const  {
		return f(bRec) * std::abs(Frame::cosTheta(bRec.wo));
	}

	/// Calculate the probability of sampling wi (given wo) -- continuous component
	virtual Float pdf(const BSDFQueryRecord &bRec) const = 0;

	/// Calculate the probability of sampling wi (given wo) -- 0D (discrete) component
	virtual Float pdfDelta(const BSDFQueryRecord &bRec) const;

	/// Evaluate the transfport from wi to wo -- 0D (discrete) component
	virtual Spectrum fDelta(const BSDFQueryRecord &bRec) const;

	/// Return the name of this BSDF
	inline const std::string &getName() const { return m_name; }

	/// Set the name of this BSDF
	inline void setName(const std::string &name) { m_name = name; }

	/// Add a child object
	virtual void addChild(const std::string &string, ConfigurableObject *obj);

	/// Serialize this object to a stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

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

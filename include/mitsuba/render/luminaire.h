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

#if !defined(__LUMINAIRE_H)
#define __LUMINAIRE_H

#include <mitsuba/render/shape.h>
#include <mitsuba/render/shader.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Data structure used by the direct illumination / shadow ray
 * sampling methods in the class \ref Luminaire.
 * \ingroup librender
 */
struct MTS_EXPORT_RENDER LuminaireSamplingRecord {
public:
	/// Create an invalid shadow ray sampling record
	inline LuminaireSamplingRecord() : luminaire(NULL) { }

	/// Create a shadow ray sampling record based on a surface intersection
	inline LuminaireSamplingRecord(const Intersection &its, const Vector &direction);

	/// Return a string representation
	std::string toString() const;
public:
	/// Associated luminaire
	const Luminaire *luminaire;

	/// Data record of the associated shape sample
	ShapeSamplingRecord sRec;

	/// Direction vector pointing away from the light source
	Vector d;

	/**
	 * \brief Probability density (wrt. solid angle) of the sampled 
	 * point on the luminaire
	 */
	Float pdf;

	/**
	 * \brief Emitted radiance at \c p into direction \c d divided by
	 * the associated probability. Already contains the geometric term
	 * and optionally also transmittance when generated via 
	 * \ref Scene::sampleAttenuatedLuminaire.
	 */
	Spectrum value;
};

/**
 * \brief Data structure used to record information associated with
 * emission sampling in the class \ref Luminaire.
 * \ingroup librender
 */
struct MTS_EXPORT_RENDER EmissionRecord {
public:
	/**
	 * \brief This class supports a special \a preview mode when
	 * sampling emissions for use in a VPL-style rendering algorithm
	 */
	enum ESamplingType {
		ENormal,
		EPreview
	};

	/// Construct a luminaire sampling record that can be used to query a luminaire
	inline EmissionRecord(const Luminaire *luminaire, 
			const ShapeSamplingRecord &sRec, const Vector &d) 
		: luminaire(luminaire), type(ENormal), sRec(sRec), d(d) { }

	inline EmissionRecord() : luminaire(NULL), type(ENormal) { }

	/// Return a string representation
	std::string toString() const;
public:
	/// Associated luminaire
	const Luminaire *luminaire;

	ESamplingType type;

	/// Data record of the associated shape sample
	ShapeSamplingRecord sRec;

	/// Direction vector pointing away from the light source
	Vector d;

	/**
	 * \brief Stores the spatial component of the radiant
	 * emittance
	 *
	 * When the record was populated using Scene::sampleEmission(),
	 * \c P will also be modulated by the directional scattering
	 * distribution and divided by the associated sampling densities.
	 */
	Spectrum value;

	/// Area probability density
	Float pdfArea;

	/// Directional probability density (wrt. solid angle)
	Float pdfDir;
};

/**
 * \brief Abstract implementation of a luminaire. Supports emission and
 * direct illumination sampling strategies, and computes related probabilities.
 * \ingroup librender
 */
class MTS_EXPORT_RENDER Luminaire : public ConfigurableObject, public HWResource {
public:
	/**
	 * \brief Non-exhaustive list of flags that can be used to characterize 
	 * light sources in \ref getType()
	 */
	enum EType {
		/// The light source has a degenerate directional density
		EDeltaDirection = 0x1,

		/// The light source is perfectly diffuse with respect to direction.
		EDiffuseDirection = 0x02,

		/// The light source is associated with a surface
		EOnSurface = 0x04,

		/// The light source has a degenerate spatial density
		EDeltaPosition = 0x8,

		/// The light source has a degenerate spatial \a or directional density
		EDelta = EDeltaDirection | EDeltaPosition
	};
	
	// =============================================================
	//! @{ \name General information
	// =============================================================

	/// Return the name of this luminaire
	inline const std::string &getName() const { return m_name; }

	/// Return the luminaire type (a combination of the properties in \ref EType)
	inline int getType() const { return m_type; }

	/**
	 * \brief Return an estimate of the total amount of power emitted 
	 * by this luminaire.
	 */
	virtual Spectrum getPower() const;

	/// Is this luminaire intersectable (e.g. can it be encountered by a tracing a ray)?
	inline bool isIntersectable() const { return m_intersectable; }

	/// Specify the medium that surrounds the luminaire
	inline void setMedium(Medium *medium) { m_medium = medium; }

	/// Return a pointer to the medium that surrounds the luminaire
	inline Medium *getMedium() { return m_medium.get(); }

	/// Return a pointer to the medium that surrounds the luminaire (const version)
	inline const Medium *getMedium() const { return m_medium.get(); }

	/**
	 * \brief Return the luminaire's sampling weight
	 *
	 * This is used by the luminaire importance sampling
	 * routines in \ref Scene.
	 */
	inline Float getSamplingWeight() const { return m_samplingWeight; }

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Direct illumination sampling strategies
	// =============================================================

	/**
	 * \brief Shadow ray sampling routine: Given an arbitrary 3D position,
	 * generate a sample point on the luminaire and fill the supplied sampling
	 * record with relevant information.
	 *
	 * Sampling is ideally done with respect to solid angle at \c p.
	 */
	virtual void sample(const Point &p, 
		LuminaireSamplingRecord &lRec, const Point2 &sample) const;

	/**
	 * \brief Calculate the solid angle density for generating this sample
	 * using the luminaire sampling strategy implemented by this class.
	 *
	 * When \c delta is set to true, only components with a Dirac delta density
	 * are considered in the query. Otherwise, they are left out.
	 */
	virtual Float pdf(const Point &p, 
		const LuminaireSamplingRecord &lRec, bool delta) const;
	
	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Emission sampling strategies
	// =============================================================

	/**
	 * \brief Sample a particle leaving this luminaire and return a ray 
	 * describing its path as well as record containing detailed probability
	 * density information. 
	 *
	 * Two uniformly distributed 2D samples are required. This method does
	 * exactly the same as calling \c sampleEmissionArea and 
	 * \c sampleEmissionDirection in sequence, modulating \c eRec.Le
	 * by the return value of the latter and dividing by the product
	 * of the spatial and directional sampling densities.
	 */
	virtual void sampleEmission(EmissionRecord &eRec,
		const Point2& areaSample, const Point2 &dirSample) const;

	/**
	 * \brief Sample only the spatial part of the emission sampling strategy
	 * implemented in \c sampleEmission.
	 *
	 * An examplary use of this method is bidirectional path tracing or MLT, 
	 * where the area and direction sampling steps take place in different 
	 * vertices.
	 *
	 * After the function call terminates, the area density as well as the
	 * spatially dependent emittance component will be stored in \c eRec.
	 */
	virtual void sampleEmissionArea(EmissionRecord &lRec,
			const Point2 &sample) const;

	/**
	 * \brief Sample only the directional part of the emission sampling strategy
	 * implemented in \c sampleEmission.
	 *
	 * Can only be called \a after a preceding invocation of
	 * \ref sampleEmissionArea() with the same emission sampling record. 
	 *
	 * The return value of this function should be used to modulate the spatial 
	 * component of the radiant emittance obtained in \ref sampleEmissionArea.
	 */
	virtual Spectrum sampleEmissionDirection(EmissionRecord &lRec, 
			const Point2 &sample) const;

	/**
	 * \brief Given an emitted particle, populate the emission record with the 
	 * relevant probability densities.
	 *
	 * When \c delta is set to true, only components with a Dirac delta density
	 * are considered in the query. Otherwise, they are left out.
	 */
	virtual void pdfEmission(EmissionRecord &eRec, bool delta) const;

	/**
	 * \brief Evaluate the spatial component of the radiant emittance at a
	 * point on the luminaire (ignoring any directional variations).
	 */
	virtual Spectrum evalArea(const EmissionRecord &eRec) const;

	/**
	 * \brief Evaluate the directional emission distribution of this light source
	 * at a given point (ignoring the spatial component).
	 *
	 * This function is normalized so that it integrates to one.
	 */
	virtual Spectrum evalDirection(const EmissionRecord &eRec) const;

	//! @}
	// =============================================================
	
	// =============================================================
	//! @{ \name Area luminaire support 
	// =============================================================

	/**
	 * \brief This function is specific to area luminaires (i.e. luminaires
	 * that can be intersected by a \ref Scene::rayIntersect). It returns the
	 * radiant emittance into a given direction.
	 *
	 * This is function is used when an area light source has been hit by a 
	 * ray in a path tracing-style integrator, and it subsequently needs to
	 * be queried for the emitted radiance along the negative ray direction.
	 
	 * The default implementation throws an exception, which states that
	 * the method is not implemented.
	 */
	virtual Spectrum Le(const ShapeSamplingRecord &sRec,
			const Vector &d) const;

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Background luminaire support 
	// =============================================================

	/// Is this a background luminaire (e.g. an environment map?)
	virtual bool isBackgroundLuminaire() const;

	/**
	 * \brief Return the radiant emittance along a ray which does not
	 * intersect any scene objects.
	 *
	 * The default implementation throws an exception, which states that
	 * the method is not implemented.
	 */
	virtual Spectrum Le(const Ray &ray) const;

	/**
	 * \brief This function fills an emission sampling record with relevant
	 * information for the supplied ray (which doesn't intersect \a any 
	 * scene objects).
	 * 
	 * The record will be populated with the radiant emittance, as well as
	 * the spatial and directional densities for sampling an emitted ray along
	 * the opposite ray direction.
	 *
	 * This function is only relevant to background luminaires. The 
	 * default implementation throws an exception, which states that
	 * the method is not implemented.
	 *
	 * \return \c true upon success
	 */
	virtual bool createEmissionRecord(EmissionRecord &eRec, const Ray &ray) const;

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Miscellaneous
	// =============================================================

	/// Is this a compound luminaire consisting of several sub-objects?
	virtual bool isCompound() const;

	/**
	 * \brief Return a sub-element of a compound luminaire. 
	 *
	 * When expanding luminaires, the scene will repeatedly call this
	 * function with increasing indices. Returning \a NULL indicates
	 * that no more are available.
	 */
	virtual Luminaire *getElement(int i);

	/// Serialize this luminaire to a binary data stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Optional pre-process step before rendering starts
	virtual void preprocess(const Scene *scene);

	/// Add a child (e.g. a medium reference) to this luminaire
	void addChild(const std::string &name, ConfigurableObject *child);
	/// Add an unnamed child
	inline void addChild(ConfigurableObject *child) { addChild("", child); }

	//! @}
	// =============================================================

	MTS_DECLARE_CLASS()
protected:
    /// Create a new luminaire
    Luminaire(const Properties &props);
    
	/// Unserialize a luminaire
    Luminaire(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~Luminaire();
protected:
	Transform m_worldToLuminaire, m_luminaireToWorld;
	Float m_samplingWeight;
	ref<Medium> m_medium;
	int m_type;
	bool m_intersectable;
	std::string m_name;
};

MTS_NAMESPACE_END

#endif /* __LUMINAIRE_H */

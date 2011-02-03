/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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
 * Data structure used to record information associated with
 * sampled shadow rays
 */
struct MTS_EXPORT_RENDER LuminaireSamplingRecord {
public:
	/// Create an invalid record
	inline LuminaireSamplingRecord() : luminaire(NULL) { }
	
	/**
	 * When a ray strikes a luminaire that is part of the scene,
	 * the associated intersection record can be converted into
	 * a luminaire sampling record in order to query the luminaire
	 * for emitted radiance. (defined in records.inl)
	 */
	LuminaireSamplingRecord(const Intersection &its, const Vector &direction);

	/// Return a string representation
	std::string toString() const;
public:
	/// Associated luminaire
	const Luminaire *luminaire;

	/// Data record of the associated shape sample
	ShapeSamplingRecord sRec;

	/// Direction vector pointing away from the light source
	Vector d;

	/// Probability density of the sampled point on the luminaire
	Float pdf;

	/**
	 * Emitted radiance at 'p' into direction 'd' divided by the associated
	 * probability. Already contains the geometric term and optionally 
	 * attenuation when generated via Scene::sampleLuminaireAttenuated.
	 */
	Spectrum Le;
};

/**
 * Data structure used to record information associated with
 * luminaire emission sampling
 */
struct MTS_EXPORT_RENDER EmissionRecord {
public:
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
	 * Radiant emittance at the sampled point. When this 
	 * record was populated using Scene::sampleEmission(), 'P'
	 * has already been multiplied by the directional 
	 * scattering distribution and divided by the associated 
	 * sampling densities.
	 */
	Spectrum P;

	/// Area probability density
	Float pdfArea;

	/// Directional probability density (wrt. solid angle)
	Float pdfDir;
};

/**
 * Abstract implementation of a luminaire. Supports emission and
 * direct illumination sampling strategies and computes related probabilities.
 */
class MTS_EXPORT_RENDER Luminaire : public ConfigurableObject, public HWResource {
public:
	enum EType {
		EDeltaDirection = 0x1,
		EDiffuseDirection = 0x02,
		EOnSurface = 0x04,
		EDeltaPosition = 0x8,
		EDelta = EDeltaDirection | EDeltaPosition
	};

	/// ================= Direct illumination sampling ================= 
	/**
	 * Return the radiant emittance into a given direction. This is
	 * primarily used when an area light source has been hit by a ray, 
	 * and it subsequently needs to be queried for the emitted radiance
	 * passing into the opposite direction.
	 */
	virtual Spectrum Le(const LuminaireSamplingRecord &lRec) const = 0;

	/**
	 * Return the radiant emittance along a ray which does not
	 * intersect any scene objects. The default implementation
	 * returns zero - this can be used to implement background light sources.
	 */
	virtual Spectrum Le(const Ray &ray) const;

	/**
	 * Shadow ray sampling routine: Given an arbitrary 3D position
	 * generate a sample point on the luminaire and fill the 
	 * sampling record with all associated information.
	 * Sampling is ideally wrt. solid angle at 'p'.
	 */
	virtual void sample(const Point &p, 
		LuminaireSamplingRecord &lRec, const Point2 &sample) const = 0;

	/**
	 * Shadow ray sampling routine: Given a surface intersection, 
	 * generate a sample point on the luminaire and fill the 
	 * sampling record with all associated information.
	 * Sampling is ideally wrt. solid angle at 'its'.
	 */
	virtual void sample(const Intersection &its, 
		LuminaireSamplingRecord &lRec, const Point2 &sample) const = 0;

	/**
	 * Calculate the probability of generating this sample using
	 * the luminaire sampling strategy implemented by this class.
	 */
	virtual Float pdf(const Point &p, 
		const LuminaireSamplingRecord &lRec, bool delta) const = 0;

	/**
	 * Calculate the probability of generating this sample using
	 * the luminaire sampling strategy implemented by this class
	 * (as above, but for surface interactions)
	 */
	virtual Float pdf(const Intersection &its, 
		const LuminaireSamplingRecord &lRec, bool delta) const = 0;

	/// ================= Emission sampling ================= 

	/**
	 * Sample a particle leaving this luminaire and return a ray describing its path
	 * as well as record containing detailed probability density information. 
	 * Two uniformly distributed 2D samples are required.
	 */
	virtual void sampleEmission(EmissionRecord &eRec,
		const Point2& areaSample, const Point2 &dirSample) const = 0;

	/**
	 * Sample only the spatial dimension of the emission sampling strategy
	 * implemented in <tt>sampleEmission</tt>. An examplary use of this is in
	 * bidirectional path tracing or MLT, where the area and direction sampling
	 * steps take place in different vertices (essentially, the directional
	 * variation is similar to a BSDF that modulates the spatially dependent
	 * radiance component). After the function call terminates, the area density 
	 * as well as the spatially dependent emittance will be stored in <tt>eRec</tt>.
	 */
	virtual void sampleEmissionArea(EmissionRecord &lRec, const Point2 &sample) const = 0;

	/**
	 * As above, but handles only the directional part. Must be called *after*
	 * sampleEmissionArea(). The return value is to be understood as a BRDF,
	 * which modulates the radiant emittance.
	 */
	virtual Spectrum sampleEmissionDirection(EmissionRecord &lRec, const Point2 &sample) const = 0;

	/**
	 * Given an emitted particle, populate the emission record with the relevant 
	 * probability densities. When \a delta is set to true, only components
	 * with a Dirac delta density are queried.
	 */
	virtual void pdfEmission(EmissionRecord &eRec, bool delta) const = 0;

	/**
	 * Evaluate the directional scattering distribution of this light source
	 * at a given point. Similar to Le(), except that this function is 
	 * normalized like a BSDF.
	 */
	virtual Spectrum f(const EmissionRecord &eRec) const = 0;

	/**
	 * Evaluate the radiant emittance at a point on the luminaire 
	 * (ignoring any directional variations).
	 */
	virtual Spectrum fArea(const EmissionRecord &eRec) const = 0;

	/// ================= Misc. ================= 

	/**
	 * \brief Fill the supplied emission record with information matching 
	 * the associated ray.
	 *
	 * This function is only relevant to background luminaires. The default
	 * implementation throws an exception, which states that the method is
	 * unimplemented.
	 *
	 * \return \a true upon success
	 */
	virtual bool createEmissionRecord(EmissionRecord &eRec, const Ray &ray) const;

	/**
	 * Return an estimate of the total amount of power emitted 
	 * by this luminaire.
	 */
	virtual Spectrum getPower() const = 0;

	/// Is this a background luminaire (e.g. an environment map?)
	virtual bool isBackgroundLuminaire() const;

	/// Is this luminaire intersetable (e.g. can it be found by a tracing a ray)?
	inline bool isIntersectable() const { return m_intersectable; }

	/// Serialize this luminaire to disk
	virtual void serialize(Stream *stream, InstanceManager *manager) const;
	
	/// Return the name of this luminaire
	inline const std::string &getName() const { return m_name; }

	/// Return the luminaire type
	inline int getType() const { return m_type; }

	/// Optional pre-process step before rendering starts
	virtual void preprocess(const Scene *scene);

	/**
	 * \brief Return the luminaire's sampling weight. This is used by
	 * the luminaire importance sampling routines in \ref Scene.
	 */
	inline Float getSamplingWeight() const { return m_samplingWeight; }

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
	int m_type;
	bool m_intersectable;
	std::string m_name;
};

MTS_NAMESPACE_END

#endif /* __LUMINAIRE_H */

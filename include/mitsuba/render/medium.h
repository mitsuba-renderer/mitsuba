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

#if !defined(__MEDIUM_H)
#define __MEDIUM_H

#include <mitsuba/core/netobject.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN
/**
 * Data record associated with the sampling procedure responsible for
 * choosing a point on the in-scattering line integral (while solving 
 * the radiative transfer equation using Monte Carlo methods).
 */
struct MTS_EXPORT_RENDER MediumSamplingRecord {
public:
	inline MediumSamplingRecord() : medium(NULL) { }

	/// Return a string representation
	std::string toString() const;
public:
	/* Traveled distance */
	Float t;

	/* Interaction point */
	Point p;

	/* Local particle orientation */
	Vector orientation;

	/* Reference to the associated medium */
	const Medium *medium;

	/* Specifies the attenuation along the segment [mint, t].
	   When sampling a distance fails, this contains the 
	   attenuation along the whole ray.
	*/
	Spectrum attenuation;

	/* The medium's absorption coefficient at that point */
	Spectrum sigmaA;

	/* The medium's scattering coefficient at that point */
	Spectrum sigmaS;

	/**
	 * Can contain two things:
	 * If a medium interaction occurred, this records the probability 
	 * of having sampled the point p. Otherwise, it contains the
	 * probability of moving through the medium without an interaction.
	 */
	Float pdf;

	/// Max. albedo over all spectral samples
	Float albedo;

	/// Multiple importance sampling weight
	Float miWeight;
};

/** \brief Abstract phase function.
 *
 * The convention used here is that the incident and exitant
 * direction arguments point away from the scattering event (similar to BSDFs).
 */
class MTS_EXPORT_RENDER PhaseFunction : public ConfigurableObject {
public:
	enum ESampledType {
		ENormal = 0,
		EDelta
	};

	/// Evaluate the phase function for an outward-pointing pair of directions (wi, wo)
	virtual Spectrum f(const MediumSamplingRecord &mRec, const Vector &wi, const Vector &wo) const = 0;

	/**
	 * Importance sample the phase function. Division by the
	 * PDF is included in the returned value.
	 */
	virtual Spectrum sample(const MediumSamplingRecord &mRec, const Vector &wi, Vector &wo, 
		ESampledType &sampledType, const Point2 &sample) const = 0;

	/**
	 * Importance sample the phase function - do not divide
	 * by the PDF and return it explicitly.
	 */
	virtual Spectrum sample(const MediumSamplingRecord &mRec, const Vector &wi, Vector &wo, 
		ESampledType &sampledType, Float &pdf, const Point2 &sample) const = 0;

	/**
	 * Calculate the probability of sampling wo (given wi). Assuming
	 * that the phase function can be sampled exactly, the default 
	 * implementation just evaluates f()
	 */
	virtual Float pdf(const MediumSamplingRecord &mRec, const Vector &wi, const Vector &wo) const;

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
};

/** \brief Abstract participating medium 
 */
class MTS_EXPORT_RENDER Medium : public NetworkedObject {
public:
	/**
	 * Possibly perform a pre-process task. The last three parameters are
	 * resource IDs of the associated scene, camera and sample generator,
	 * which have been made available to all local and remote workers.
	 */
	virtual void preprocess(const Scene *scene, RenderQueue *queue, 
		const RenderJob *job, int sceneResID, int cameraResID, int samplerResID);

	/** 
	 * Optical thickness computation - integrates 
	 * the extinction coefficient along a ray segment [mint, maxt]
	 */
	virtual Spectrum tau(const Ray &ray) const = 0;
	
	/**
	 * Sample a distance along a ray - should ideally importance 
	 * sample with respect to the attenuation.
	 * Returns false if the maximum distance was exceeded.
	 */
	virtual bool sampleDistance(const Ray &ray, Float maxDist,
		MediumSamplingRecord &mRec, Sampler *sampler) const = 0;

	/// Return the phase function of this medium
	inline const PhaseFunction *getPhaseFunction() const { return m_phaseFunction.get(); }

	/** \brief Configure the object (called _once_ after construction
	   and addition of all child ConfigurableObjects. */
	virtual void configure();

	/// Serialize this medium to a stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Add a child ConfigurableObject
	virtual void addChild(const std::string &name, ConfigurableObject *child);

	/// Return a bounding volume
	inline const AABB &getAABB() const { return m_aabb; }

	/// Return a string representation
	virtual std::string toString() const = 0;

	MTS_DECLARE_CLASS()
protected:
	/// Create a new participating medium instance
	Medium(const Properties &props);
	
	/// Unserialize a participating medium
	Medium(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~Medium() { }
protected:	
	Spectrum m_sigmaA;
	Spectrum m_sigmaS;
	Spectrum m_sigmaT;
	Float m_albedo;
	AABB m_aabb;
	ref<PhaseFunction> m_phaseFunction;
	Float m_sizeMultiplier;
};

MTS_NAMESPACE_END

#endif /* __MEDIUM_H */

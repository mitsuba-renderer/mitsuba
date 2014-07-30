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
#if !defined(__MITSUBA_RENDER_MEDIUM_H_)
#define __MITSUBA_RENDER_MEDIUM_H_

#include <mitsuba/core/netobject.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN
/**
 * \brief Data record for sampling a point on the in-scattering
 * integral of the RTE
 *
 * \sa Medium::sampleDistance()
 * \ingroup librender
 */
struct MTS_EXPORT_RENDER MediumSamplingRecord {
public:

	/// Traveled distance
	Float t;

	/// Location of the scattering interaction
	Point p;

	/// Time value associated with the medium scattering event
	Float time;

	/// Local particle orientation at \ref p
	Vector orientation;

	/**
	 * \brief Specifies the transmittance along the segment [mint, t]
	 *
	 * When sampling a distance fails, this contains the
	 * transmittance along the whole ray segment [mint, maxDist].
	 */
	Spectrum transmittance;

	/// The medium's absorption coefficient at \ref p
	Spectrum sigmaA;

	/// The medium's scattering coefficient at \ref p
	Spectrum sigmaS;

	/// Records the probability density of sampling a medium interaction at p
	Float pdfSuccess;

	/**
	 * \brief Records the probability density of sampling a medium
	 * interaction in the reverse direction
	 *
	 * This is essentially the density of obtained by calling \ref sampleDistance,
	 * but starting at \c p and stopping at \c ray.o. These probabilities
	 * are important for bidirectional methods.
	 */
	Float pdfSuccessRev;

	/**
	 * When the \ref Medium::sampleDistance() is successful, this function
	 * returns the probability of \a not having generated a medium interaction
	 * until \ref t. Otherwise, it records the probability of
	 * not generating any interactions in the whole interval [mint, maxt].
	 * This probability is assumed to be symmetric with respect to
	 * sampling from the other direction, which is why there is no
	 * \c pdfFailureRev field.
	 */
	Float pdfFailure;

	/// Pointer to the associated medium
	const Medium *medium;

public:
	inline MediumSamplingRecord() : medium(NULL) { }

	/// Return a pointer to the phase function
	inline const PhaseFunction *getPhaseFunction() const;

	/// Return a string representation
	std::string toString() const;
};

/** \brief Abstract participating medium
 * \ingroup librender
 */
class MTS_EXPORT_RENDER Medium : public NetworkedObject {
public:
	// =============================================================
	//! @{ \name Medium sampling strategy
	// =============================================================

	/**
	 * \brief Sample a distance along the ray segment [mint, maxt]
	 *
	 * Should ideally importance sample with respect to the transmittance.
	 * It is assumed that the ray has a normalized direction value.
	 *
	 * \param ray      Ray, along which a distance should be sampled
	 * \param mRec     Medium sampling record to be filled with the result
	 * \return         \c false if the maximum distance was exceeded, or if
	 *                 no interaction inside the medium could be sampled.
	 */
	virtual bool sampleDistance(const Ray &ray,
		MediumSamplingRecord &mRec, Sampler *sampler) const = 0;

	/**
	 * \brief Compute the 1D density of sampling distance \a ray.maxt
	 * along the ray using the sampling strategy implemented by
	 * \a sampleDistance.
	 *
	 * The function computes the continuous densities in the case of
	 * a successful \ref sampleDistance() invocation (in both directions),
	 * as well as the Dirac delta density associated with a failure.
	 * For convenience, it also stores the transmittance along the
	 * supplied ray segment within \a mRec.
	 */
	virtual void eval(const Ray &ray,
		MediumSamplingRecord &mRec) const = 0;

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Functions for querying the medium
	// =============================================================

	/**
	 * \brief Compute the transmittance along a ray segment
	 *
	 * Computes the transmittance along a ray segment
	 * [mint, maxt] associated with the ray. It is assumed
	 * that the ray has a normalized direction value.
	 */
	virtual Spectrum evalTransmittance(const Ray &ray,
		Sampler *sampler = NULL) const = 0;

	/// Return the phase function of this medium
	inline const PhaseFunction *getPhaseFunction() const { return m_phaseFunction.get(); }

	/// Determine whether the medium is homogeneous
	virtual bool isHomogeneous() const = 0;

	/// For homogeneous media: return the absorption coefficient
	inline const Spectrum &getSigmaA() const { return m_sigmaA; }

	/// For homogeneous media: return the scattering coefficient
	inline const Spectrum &getSigmaS() const { return m_sigmaS; }

	/// For homogeneous media: return the extinction coefficient
	inline const Spectrum &getSigmaT() const { return m_sigmaT; }

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Miscellaneous
	// =============================================================

	/** \brief Configure the object (called \a once after construction
	   and addition of all child \ref ConfigurableObject instances). */
	virtual void configure();

	/// Serialize this medium to a stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Add a child ConfigurableObject
	virtual void addChild(const std::string &name, ConfigurableObject *child);
	/// Add an unnamed child
	inline void addChild(ConfigurableObject *child) { addChild("", child); }

	/// Return a string representation
	virtual std::string toString() const = 0;

	//! @}
	// =============================================================

	MTS_DECLARE_CLASS()
protected:
	/// Create a new participating medium instance
	Medium(const Properties &props);

	/// Unserialize a participating medium
	Medium(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~Medium() { }
protected:
	ref<PhaseFunction> m_phaseFunction;
	Spectrum m_sigmaA;
	Spectrum m_sigmaS;
	Spectrum m_sigmaT;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_MEDIUM_H_ */

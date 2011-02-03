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

#if !defined(__PHASE_H)
#define __PHASE_H

#include <mitsuba/core/netobject.h>
#include <mitsuba/render/common.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Data structure, which contains information 
 * required to sample or query a phase function. 
 */
struct MTS_EXPORT_RENDER PhaseFunctionQueryRecord {
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
	 * \a wi points inwards.
	 */
	Vector wi;

	/// Normalized outgoing direction vector
	Vector wo;

	/* Transported quantity (radiance or importance) -- required for 
	   rendering with non-reciprocal phase functions */
	ETransportQuantity quantity;

	inline PhaseFunctionQueryRecord(const MediumSamplingRecord &mRec,
		const Vector &wi) : mRec(mRec), wi(wi), quantity(ERadiance) {
	}

	inline PhaseFunctionQueryRecord(const MediumSamplingRecord &mRec,
		const Vector &wi, const Vector &wo) : mRec(mRec), wi(wi), 
		wo(wo), quantity(ERadiance) {
	}

	std::string toString() const;
};

/** \brief Abstract phase function.
 */
class MTS_EXPORT_RENDER PhaseFunction : public ConfigurableObject {
public:
	/**
	 * \brief Evaluate the phase function for an outward-pointing 
	 * pair of directions (wi, wo)
	 */
	virtual Spectrum f(const PhaseFunctionQueryRecord &pRec) const = 0;

	/**
	 * \brief Importance sample the phase function. 
	 *
	 * \param sampler
	 *     Sample generator
	 * \return
	 *     Weight value equal to the throughput divided by 
	 *     the probability of the sampled direction.
	 */
	virtual Spectrum sample(PhaseFunctionQueryRecord &pRec, 
		Sampler *sampler) const = 0;

	/**
	 * \brief Importance sample the phase function, but don't
	 *    divide by the computed probability.
	 * \param pdf
	 *     The probability of sampling \a pRec.wo will be returned using
	 *     this argument.
	 * \return
	 *     Phase function value for the direction pair (wi, wo)
	 */
	virtual Spectrum sample(PhaseFunctionQueryRecord &pRec,
		Float &pdf, Sampler *sampler) const = 0;

	/**
	 * \brief Calculate the probability of sampling wo (given wi).
	 *
	 * Assuming that the phase function can be sampled exactly, 
	 * the default implementation just evaluates \ref f()
	 */
	virtual Float pdf(const PhaseFunctionQueryRecord &pRec) const;

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

MTS_NAMESPACE_END

#endif /* __PHASE_H */

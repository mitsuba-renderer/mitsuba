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

#include <mitsuba/bidir/common.h>
#include <mitsuba/bidir/mutator.h>

#define MTS_BD_MEDIUM_PERTURBATION_MONOCHROMATIC 1

MTS_NAMESPACE_BEGIN

std::string EndpointRecord::toString() const {
	std::ostringstream oss;
	oss << "EndpointRecord[time=" << time << "]";
	return oss.str();
}

std::ostream &operator<<(std::ostream &os, const Mutator::EMutationType &type) {
	switch (type) {
		case Mutator::EBidirectionalMutation: os << "bidir"; break;
		case Mutator::ELensPerturbation: os << "lens"; break;
		case Mutator::ELensSubpathMutation: os << "lensSubpath"; break;
		case Mutator::ECausticPerturbation: os << "caustic"; break;
		case Mutator::EIndependentMutation: os << "indep"; break;
		case Mutator::EMultiChainPerturbation: os << "multiChain"; break;
		case Mutator::EManifoldPerturbation: os << "manifold"; break;
		default: os << "invalid"; break;
	};
	return os;
}

std::string MutationRecord::toString() const {
	std::ostringstream oss;
	oss << "MutationRecord["
		<< "type=" << type
		<< ", l=" << l
		<< ", m=" << m
		<< ", kd=" << m-l
		<< ", ka=" << ka
		<< ", weight=" << weight.toString()
		<< "]";
	return oss.str();
}

MutatorBase::MutatorBase() {
	m_mediumDensityMultiplier = 100.0f;
}

Float MutatorBase::perturbMediumDistance(Sampler *sampler, const PathVertex *vertex) {
	if (vertex->isMediumInteraction()) {
#if MTS_BD_MEDIUM_PERTURBATION_MONOCHROMATIC == 1
		/* Monochromatic version */
		const MediumSamplingRecord &mRec = vertex->getMediumSamplingRecord();
		Float sigma = (mRec.sigmaA + mRec.sigmaS).average() * m_mediumDensityMultiplier;
#else
		const MediumSamplingRecord &mRec = vertex->getMediumSamplingRecord();
		Spectrum sigmaT = (mRec.sigmaA + mRec.sigmaS) * m_mediumDensityMultiplier;
		Float sigma = sigmaT[
			std::min((int) (sampler->next1D() * SPECTRUM_SAMPLES), SPECTRUM_SAMPLES-1)];
#endif
		return (sampler->next1D() > .5 ? -1.0f : 1.0f) *
			math::fastlog(1-sampler->next1D()) / sigma;
	} else {
		return 0.0f;
	}
}

Float MutatorBase::pdfMediumPerturbation(const PathVertex *oldVertex,
		const PathEdge *oldEdge, const PathEdge *newEdge) const {
	BDAssert(oldEdge->medium && newEdge->medium);
	const MediumSamplingRecord &mRec = oldVertex->getMediumSamplingRecord();
#if MTS_BD_MEDIUM_PERTURBATION_MONOCHROMATIC == 1
	Float sigmaT = (mRec.sigmaA + mRec.sigmaS).average() * m_mediumDensityMultiplier;
	Float diff = std::abs(oldEdge->length - newEdge->length);
	return 0.5f * sigmaT*math::fastexp(-sigmaT*diff);
#else
	Spectrum sigmaT = (mRec.sigmaA + mRec.sigmaS) * m_mediumDensityMultiplier;
	Float diff = std::abs(oldEdge->length - newEdge->length);
	Float sum = 0.0f;
	for (int i=0; i<SPECTRUM_SAMPLES; ++i)
		sum += sigmaT[i]*math::fastexp(-sigmaT[i]*diff);
	return sum * (0.5f / SPECTRUM_SAMPLES);
#endif
}

MTS_IMPLEMENT_CLASS(Mutator, true, Object)
MTS_IMPLEMENT_CLASS(MutatorBase, true, Mutator)
MTS_NAMESPACE_END

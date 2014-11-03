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

#include <mitsuba/core/statistics.h>
#include <mitsuba/bidir/mut_bidir.h>
#include <mitsuba/bidir/geodist2.h>

MTS_NAMESPACE_BEGIN

static StatsCounter statsAccepted("Bidirectional mutation",
		"Acceptance rate", EPercentage);
static StatsCounter statsGenerated("Bidirectional mutation",
		"Successful generation rate", EPercentage);

BidirectionalMutator::BidirectionalMutator(const Scene *scene,
	Sampler *sampler, MemoryPool &pool, int kmin, int kmax) :
	m_scene(scene), m_sampler(sampler), m_pool(pool),
	m_kmin(kmin), m_kmax(kmax) {
}

BidirectionalMutator::~BidirectionalMutator() { }

Mutator::EMutationType BidirectionalMutator::getType() const {
	return EBidirectionalMutation;
}

Float BidirectionalMutator::suitability(const Path &path) const {
	return 1.0f;
}

bool BidirectionalMutator::sampleMutation(
		Path &source, Path &proposal, MutationRecord &muRec, const MutationRecord& sourceMuRec) {
	TwoTailedGeoDistr desiredLength(2), deletionLength(2);
	int k = source.length();

	/* Sample the desired path length of the proposal. This
	   is done using a two-tailed geometric distribution that is
	   centered around the current path length, and which respects
	   the specified minimum and maximum length constraints. */

	desiredLength.configure(k, m_kmin, m_kmax);
	int kPrime = desiredLength.sample(m_sampler->next1D());

	/* Sample the length of the deletion (in # of edges, 1 means
	   no vertices are removed). When kPrime is smaller than k,
	   we must delete at least k-kPrime+1 edges to be able to
	   achieve the desired path length.

	   When k==kPrime, we must delete *something*, or the mutation
	   is trivial, hence the conditional below expression. */

	int minDeletion = std::max((k == kPrime) ? 2 : 1, k-kPrime+1);
	deletionLength.configure(2, minDeletion, k);
	int kd = deletionLength.sample(m_sampler->next1D());

	/* Based on the desired length, this tells us how many
	   edges need to be added (k' = k - kd + ka) */
	int ka = kPrime-k+kd;

	/* Sample the left endpoint of the deleted range */
	int lMin = 0, lMax = k - kd;
	if (kd == 1 || ka == 1) {
		/* This will help to avoid certain path changes that would otherwise
		   always be rejected. Specifically, we don't want to remove the sensor
		   or emitter sample vertex, and we don't want to insert
		   vertices between a sensor/emitter sample and its supernode */
		lMin++; lMax--;
	}
	m_temp.clear();
	for (int l=lMin; l<=lMax; ++l) {
		int m = l+kd;
		if (!source.vertex(l)->isDegenerate() &&
			!source.vertex(m)->isDegenerate())
			m_temp.push_back(l);
	}
	if (m_temp.size() == 0)
		return false;

	int l = m_temp[std::min((int) (m_temp.size() *
			m_sampler->next1D()), (int) m_temp.size()-1)];
	int m = l+kd;

	/* Don't try to hit the emitter or sensor if they are degenerate */
	int sMin = 0, sMax = ka-1;
	if (l == 0 && m_scene->hasDegenerateEmitters())
		++sMin;
	else if (m == k && m_scene->hasDegenerateSensor())
		--sMax;

	/* Sample the number of SIS-type steps to take from the emitter direction */
	int s = std::min(sMin + (int) ((sMax-sMin+1) * m_sampler->next1D()), sMax);
	int t = ka - s - 1;

	/* Check a few assumptions */
	BDAssert(ka >= 1 && kd >= 1 && kd <= k
			&& l >= lMin && l <= lMax
			&& kPrime == k - kd + ka
			&& kPrime >= m_kmin
			&& kPrime <= m_kmax);

	/* Construct a mutation record */
	muRec = MutationRecord(EBidirectionalMutation, l, m, ka,
		source.getPrefixSuffixWeight(l, m));

	/* Keep some statistics */
	statsGenerated.incrementBase();
	statsAccepted.incrementBase();

	proposal.clear();
	proposal.append(source, 0, l+1);
	proposal.vertex(l) = proposal.vertex(l)->clone(m_pool);

	/* Perform a random walk from the emitter direction */
	if (proposal.randomWalk(m_scene, m_sampler, s, -1, EImportance, m_pool) != s) {
		proposal.release(l, proposal.vertexCount(), m_pool);
		return false;
	}

	/* Perform a random walk from the sensor direction */
	m_tempPath.clear();
	m_tempPath.append(source, m, k+1, true);
	m_tempPath.vertex(k-m) = m_tempPath.vertex(k-m)->clone(m_pool);

	if (m_tempPath.randomWalk(m_scene, m_sampler, t, -1, ERadiance, m_pool) != t) {
		proposal.release(l, proposal.vertexCount(), m_pool);
		m_tempPath.release(k-m, m_tempPath.vertexCount(), m_pool);
		return false;
	}

	PathEdge *connectionEdge = m_pool.allocEdge();
	proposal.append(connectionEdge);
	proposal.append(m_tempPath, 0, m_tempPath.vertexCount(), true);

	BDAssert(proposal.length() == kPrime &&
			 proposal.vertexCount() == proposal.edgeCount() + 1);

	const PathVertex
		*vsPred = l+s > 0 ? proposal.vertex(l+s-1) : NULL,
		*vtPred = l+s+2 <= kPrime ? proposal.vertex(l+s+2) : NULL;
	const PathEdge
		*vsEdge = l+s > 0 ? proposal.edge(l+s-1) : NULL,
		*vtEdge = l+s+1 < kPrime ? proposal.edge(l+s+1) : NULL;

	/* Now try to connect the two subpaths and reject
	   the proposal if there is no throughput */
	PathVertex *vs = proposal.vertex(l+s),
			   *vt = proposal.vertex(l+s+1);

	if (!PathVertex::connect(m_scene, vsPred,
			vsEdge, vs, connectionEdge, vt, vtEdge, vtPred)) {
		proposal.release(l, l+ka+1, m_pool);
		return false;
	}

	if (m >= k-1)
		proposal.vertex(kPrime-1)->updateSamplePosition(
			proposal.vertex(kPrime-2));

	++statsGenerated;
	return true;
}

Float BidirectionalMutator::pmfMutation(const Path &source, const MutationRecord &muRec) const {
	TwoTailedGeoDistr desiredLength(2), deletionLength(2);
	const int k = source.length(), m = muRec.m, l = muRec.l,
		kd = m - l, ka = muRec.ka, kPrime = k - kd + ka;
	int minDeletion = std::max((k == kPrime) ? 2 : 1, k-kPrime+1);

	/* See the sampleMutation() function for a detailed
	   description of this construction */
	int sMin = 0, sMax = ka-1;
	if (l == 0 && m_scene->hasDegenerateEmitters())
		++sMin;
	else if (m == k && m_scene->hasDegenerateSensor())
		--sMax;

	int lMin = 0, lMax = k - kd, ctr = 0;
	if (kd == 1 || ka == 1) {
		lMin++; lMax--;
	}

	for (int l=lMin; l<=lMax; ++l) {
		int m = l+kd;
		if (!source.vertex(l)->isDegenerate() &&
			!source.vertex(m)->isDegenerate())
			++ctr;
	}
	if (ctr == 0)
		return 0.0f;

	desiredLength.configure(k, m_kmin, m_kmax);
	deletionLength.configure(2, minDeletion, k);

	Float factor1 = desiredLength.pmf(kPrime);
	Float factor2 = deletionLength.pmf(kd);
	Float factor3 = 1 / (Float) ctr;
	Float factor4 = (Float) 1 / (Float) (sMax-sMin+1);

	return factor1 * factor2 * factor3 * factor4;
}

Float BidirectionalMutator::Q(const Path &source, const Path &proposal,
		const MutationRecord &muRec) const {
	const int k = source.length(), l = muRec.l,
		      m = muRec.m, ka = muRec.ka, mPrime = l+ka;

	Spectrum *importanceWeights = (Spectrum *) alloca(ka * sizeof(Spectrum)),
			 *radianceWeights  = (Spectrum *) alloca(ka * sizeof(Spectrum));

	/* Compute importance transport weights along the subpath */
	importanceWeights[0] = Spectrum(1.0f);
	for (int s=1; s<ka; ++s)
		importanceWeights[s] = importanceWeights[s-1] *
			proposal.vertex(l+s-1)->weight[EImportance] *
			proposal.edge(l+s-1)->weight[EImportance];

	/* Compute radiance transport weights along the subpath */
	radianceWeights[0] = Spectrum(1.0f);
	for (int t=1; t<ka; ++t)
		radianceWeights[t] = radianceWeights[t-1] *
			proposal.vertex(mPrime-t+1)->weight[ERadiance] *
			proposal.edge(mPrime-t)->weight[ERadiance];

	int sMin = 0, sMax = ka-1;
	if (l == 0 && m_scene->hasDegenerateEmitters())
		++sMin;
	else if (m == k && m_scene->hasDegenerateSensor())
		--sMax;

	Float result = 0.0f;
	for (int s = sMin; s <= sMax; ++s) {
		const PathEdge *edge = proposal.edge(l+s);
		const PathVertex *vs = proposal.vertex(l+s),
			  *vt = proposal.vertex(l+s+1);
		int t = ka - s - 1;

		/* Cannot connect endpoints with degenerate distributions */
		if (!vs->isConnectable() || !vt->isConnectable())
			continue;

		Spectrum weight = importanceWeights[s]
			* radianceWeights[t]
			* edge->evalCached(vs, vt, PathEdge::EEverything)
			* muRec.weight;

		Float luminance = weight.getLuminance();

		if (luminance <= RCPOVERFLOW || !std::isfinite(luminance)) {
			Log(EWarn, "Internal error: luminance = %f!", luminance);
			continue;
		}

		result += 1 / luminance;
	}

	return result * pmfMutation(source, muRec);
}

void BidirectionalMutator::accept(const MutationRecord &) {
	++statsAccepted;
}

MTS_IMPLEMENT_CLASS(BidirectionalMutator, false, Mutator)
MTS_NAMESPACE_END

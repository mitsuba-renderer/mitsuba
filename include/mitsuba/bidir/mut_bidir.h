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
#if !defined(__MITSUBA_BIDIR_MUT_BIDIR_H_)
#define __MITSUBA_BIDIR_MUT_BIDIR_H_

#include <mitsuba/bidir/mutator.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Bidirectional mutation strategy
 *
 * This class implements a slightly extended version of the bidirectional
 * mutation proposed by Veach. The main change is that it builds on top of
 * a two-tailed geometric distribution that is used to sample path
 * configuration proposals in a more flexible manner.
 *
 * \author Wenzel Jakob
 * \ingroup libbidir
 */
class MTS_EXPORT_BIDIR BidirectionalMutator : public Mutator {
public:
	/**
	 * \brief Construct a new bidirectional mutator
	 *
	 * \param scene
	 *     A pointer to the underlying scene
	 *
	 * \param sampler
	 *     A sample generator
	 *
	 * \param pool
	 *     A memory pool used to allocate new path vertices and edges
	 *
	 * \param kmin
	 *     Minimum number of edges in newly proposed paths. This can
	 *     be used to exclude direct illumination.
	 *
	 * \param kmax
	 *     Minimum number of edges in newly proposed paths.
	 */
	BidirectionalMutator(const Scene *scene, Sampler *sampler,
		MemoryPool &pool, int kmin, int kmax);

	// =============================================================
	//! @{ \name Implementation of the Mutator interface

	EMutationType getType() const;
	Float suitability(const Path &path) const;
	bool sampleMutation(Path &source, Path &proposal, MutationRecord &muRec, const MutationRecord& sourceMuRec);
	Float Q(const Path &source, const Path &proposal,
			const MutationRecord &muRec) const;
	void accept(const MutationRecord &muRec);

	//! @}
	// =============================================================

	MTS_DECLARE_CLASS()
protected:
	/**
	 * \brief Compute the probability mass associated with one
	 * of the internally implemented mutation strategies
	 */
	Float pmfMutation(const Path &source, const MutationRecord &muRec) const;

	/// Virtual destructor
	virtual ~BidirectionalMutator();
protected:
	ref<const Scene> m_scene;
	ref<Sampler> m_sampler;
	std::vector<int> m_temp;
	MemoryPool &m_pool;
	int m_kmin, m_kmax;
	Path m_tempPath;
};

MTS_NAMESPACE_END

#endif /*__MITSUBA_BIDIR_MUT_BIDIR_H_ */

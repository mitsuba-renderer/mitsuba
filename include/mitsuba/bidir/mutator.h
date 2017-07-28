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
#if !defined(__MITSUBA_BIDIR_MUTATOR_H_)
#define __MITSUBA_BIDIR_MUTATOR_H_

#include <mitsuba/bidir/path.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Generic interface to path-space mutators
 *
 * This is the superclass of all path-space mutators, such as
 * the bidirectional mutation or lens perturbation in Veach-MLT.
 */
class MTS_EXPORT_BIDIR Mutator : public Object {
public:
    /// Specifies the type of mutation implemented by the mutator
    enum EMutationType {
        EBidirectionalMutation = 0,
        ELensPerturbation,
        ELensSubpathMutation,
        EIndependentMutation,
        ECausticPerturbation,
        EMultiChainPerturbation,
        EManifoldPerturbation,
        EMutationTypeCount
    };

    /// What kind of mutations does this mutator perform?
    virtual EMutationType getType() const = 0;

    /// Determine the general "suitability" of this mutator for a given kind of path
    virtual Float suitability(const Path &path) const = 0;

    /**
     * \brief Given a path, this function produces a new proposal
     * according to the internally implemented mutation strategy
     *
     * \param source
     *     The sampling strategy implemented by the mutator
     *     will condition on this path.
     *
     * \param proposal
     *     Path data structure to be filled with the proposed mutated path
     *
     * \param muRec
     *     Data record that describes the sampled mutation strategy
     *
     * \param sourceMuRec
     *     Data record that describes the last successful mutation strategy
     *     (for the source path)
     *
     * \return \a true upon success. When the sampling step is
     *     unsuccessful (this could happen due to various
     *     reasons), the function returns <tt>false</tt>.

     */
    virtual bool sampleMutation(Path &source, Path &proposal,
            MutationRecord &muRec, const MutationRecord& sourceMuRec) = 0;

    /**
     * \brief For a pair of paths, this function computes the inverse
     * transition probability (matching the Q term in [Veach 97])
     *
     * \param source
     *     A path data structure containing the original path
     *
     * \param proposal
     *     A path data structure containing the proposed mutated path
     *
     * \param muRec
     *     Data record that describes the mutation strategy, which
     *     transformed \c source to \c proposal.
     */
    virtual Float Q(const Path &source, const Path &proposal,
        const MutationRecord &muRec) const = 0;

    /**
     * \brief Record an accepted mutation
     *
     * This function exists to allow mutators to track their
     * acceptance rate and other statistics.
     */
    virtual void accept(const MutationRecord &muRec) = 0;

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~Mutator() { }
};

/**
 * \brief Stores supplemental information about an executed mutation strategy
 *
 * These records are filled in by \ref Mutator::sampleMutation().
 */
struct MTS_EXPORT_BIDIR MutationRecord {
    Mutator::EMutationType type; ///< Type of executed mutation
    int l;                       ///< Left vertex of the affected range
    int m;                       ///< Right vertex of the affected range
    int ka;                      ///< Size of the insertion
    Spectrum weight;             ///< Spectral weight of the unchanged portion
    int extra[5];

    inline MutationRecord() { }
    inline MutationRecord(Mutator::EMutationType type, int l,
            int m, int ka, const Spectrum &weight)
     : type(type), l(l), m(m), ka(ka), weight(weight) { }

    MutationRecord reverse() const {
        MutationRecord result(type, l, l+ka, m-l, weight);
        memcpy(result.extra, extra, sizeof(extra));
        return result;
    }

    std::string toString() const;
};

extern MTS_EXPORT_BIDIR std::ostream &operator<<(
        std::ostream &os, const Mutator::EMutationType &type);

/**
 * \brief Medium-aware mutator base class
 */
class MTS_EXPORT_BIDIR MutatorBase : public Mutator {
public:
    MTS_DECLARE_CLASS()
protected:
    /// Protected constructor
    MutatorBase();

    /// Virtual destructor
    virtual ~MutatorBase() { }

    /// Perturb a distance within a medium
    Float perturbMediumDistance(Sampler *sampler,
        const PathVertex *vertex);

    /// Density function of \ref perturbMediumDistance
    Float pdfMediumPerturbation(const PathVertex *oldVertex,
        const PathEdge *oldEdge, const PathEdge *newEdge) const;

protected:
    Float m_mediumDensityMultiplier;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_BIDIR_MUTATOR_H_ */

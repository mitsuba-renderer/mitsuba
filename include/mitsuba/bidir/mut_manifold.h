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
#if !defined(__MITSUBA_BIDIR_MUT_MANIFOLD_H_)
#define __MITSUBA_BIDIR_MUT_MANIFOLD_H_

#include <mitsuba/bidir/mutator.h>

MTS_NAMESPACE_BEGIN

#define MTS_MANIFOLD_QUANTILE_SURFACE 0.9f
#define MTS_MANIFOLD_QUANTILE_MEDIUM  0.5f

/**
 * \brief Specular manifold perturbation strategy
 *
 * \author Wenzel Jakob
 * \ingroup libbidir
 */
class MTS_EXPORT_BIDIR ManifoldPerturbation : public MutatorBase {
public:
    /**
     * \brief Construct a new specular manifold perturbation strategy
     *
     * \param scene
     *     A pointer to the underlying scene
     *
     * \param sampler
     *     A sample generator
     *
     * \param pool
     *     A memory pool used to allocate new path vertices and edges
     */
    ManifoldPerturbation(const Scene *scene, Sampler *sampler,
        MemoryPool &pool,
        Float probFactor,
        bool enableOffsetManifolds,
        bool enableSpecularMedia,
        Float avgAngleChangeSurface = 0,
        Float avgAngleChangeMedium = 0);

    // =============================================================
    //! @{ \name Implementation of the Mutator interface

    EMutationType getType() const;
    Float suitability(const Path &path) const;
    bool sampleMutation(Path &source, Path &proposal,
            MutationRecord &muRec, const MutationRecord& sourceMuRec);
    Float Q(const Path &source, const Path &proposal,
            const MutationRecord &muRec) const;
    void accept(const MutationRecord &muRec);

    //! @}
    // =============================================================

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~ManifoldPerturbation();

    /// Helper function for choosing mutation strategies
    bool sampleMutationRecord(const Path &source,
        int &a, int &b, int &c, int &step);

    Float nonspecularProbSurface(Float alpha) const;
    Float nonspecularProbMedium(Float g) const;

    Float nonspecularProb(const PathVertex *vertex) const;
    inline Float specularProb(const PathVertex *vertex) const {
        return 1 - nonspecularProb(vertex);
    }

    int getSpecularChainEnd(const Path &path, int pos, int step);
protected:
    ref<const Scene> m_scene;
    ref<Sampler> m_sampler;
    mutable ref<SpecularManifold> m_manifold;
    MemoryPool &m_pool;
    Float m_probFactor, m_probFactor2;
    bool m_enableOffsetManifolds;
    bool m_enableSpecularMedia;
    static Float m_thetaDiffSurface;
    static Float m_thetaDiffMedium;
    static int m_thetaDiffSurfaceSamples;
    static int m_thetaDiffMediumSamples;
    static Mutex *m_thetaDiffMutex;
};

MTS_NAMESPACE_END

#endif /*__MITSUBA_BIDIR_MUT_MANIFOLD_H_ */

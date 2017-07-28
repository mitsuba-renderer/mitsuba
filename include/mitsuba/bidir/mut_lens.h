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
#if !defined(__MITSUBA_BIDIR_MUT_LENS_H_)
#define __MITSUBA_BIDIR_MUT_LENS_H_

#include <mitsuba/bidir/mutator.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Veach-style Lens subpath perturbation strategy
 *
 * This class implements a simple lens perturbation strategy
 * as described in Eric Veach's PhD thesis.
 *
 * \author Wenzel Jakob
 * \ingroup libbidir
 */
class MTS_EXPORT_BIDIR LensPerturbation : public MutatorBase {
public:
    /**
     * \brief Construct a new lens perturbation strategy
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
     * \param minJump
     *     Minimum jump distance in fractional pixel coordinates
     *
     * \param coveredArea
     *     Approximate fractional image plane area that is
     *     reachable using the lens perturbation
     */
    LensPerturbation(const Scene *scene, Sampler *sampler,
        MemoryPool &pool, Float minJump, Float coveredArea);

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
    virtual ~LensPerturbation();
protected:
    ref<const Scene> m_scene;
    ref<Sampler> m_sampler;
    MemoryPool &m_pool;
    Vector2 m_filmRes;
    Float m_r1, m_r2;
    Float m_logRatio;
    Float m_imagePlaneArea;
};

MTS_NAMESPACE_END

#endif /*__MITSUBA_BIDIR_MUT_LENS_H_ */

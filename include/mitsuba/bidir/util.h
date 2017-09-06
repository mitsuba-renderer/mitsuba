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
#if !defined(__MITSUBA_BIDIR_UTIL_H_)
#define __MITSUBA_BIDIR_UTIL_H_

#include <mitsuba/bidir/path.h>
#include <mitsuba/bidir/rsampler.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief A collection of powerful support functions that can be used to construct
 * bidirectional rendering algorithms
 */
class MTS_EXPORT_BIDIR BidirectionalUtils {
public:
    /**
     * \brief Render the direct illumination component of a scene
     *
     * The function is made available here, since it is used
     * by both Keleman-style and Veach-style MLT.
     */
    static ref<Bitmap> renderDirectComponent(Scene *scene, int sceneResID,
        int sensorResID, RenderQueue *queue, const RenderJob *job,
        size_t directSamples);

    /**
     * \brief Execute the first pass of a 2-pass MLT scheme.
     *
     * The function is made available here, since it is used
     * by both Keleman-style and Veach-style MLT.
     *
     * \param scene
     *     Pointer to the underlying scene
     *
     * \param sceneResID
     *     Resource ID of the scene (used for executing the first
     *     stage in parallel over multiple machines)
     *
     * \param queue
     *     Pointer to the render queue associated with the original job
     *
     * \param sizeFactor
     *     Size reduction factor to use when rendering the
     *     luminance image
     *
     * \param nestedJob
     *     Reference to a nested render job. Can be used to terminate
     *     the process from another thread
     */
    static ref<Bitmap> mltLuminancePass(Scene *scene, int sceneResID,
            RenderQueue *queue, int sizeFactor, ref<RenderJob> &nestedJob);
};

/// Restores the measure of a path vertex after going out of scope
struct RestoreMeasureHelper {
    RestoreMeasureHelper(PathVertex *vertex)
        : vertex(vertex), measure(vertex->measure) { }
    ~RestoreMeasureHelper() { vertex->measure = measure; }
    PathVertex *vertex;
    EMeasure measure;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_BIDIR_UTIL_H_ */

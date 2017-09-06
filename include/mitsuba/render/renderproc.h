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
#if !defined(__MITSUBA_RENDER_RENDERPROC_H_)
#define __MITSUBA_RENDER_RENDERPROC_H_

#include <mitsuba/render/scene.h>
#include <mitsuba/render/imageproc.h>
#include <mitsuba/render/renderqueue.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Parallel process for rendering with sampling-based integrators.
 *
 * Splits an image into independent rectangular pixel regions, which are
 * then rendered in parallel.
 *
 * \sa SamplingIntegrator
 * \ingroup librender
 */
class MTS_EXPORT_RENDER BlockedRenderProcess : public BlockedImageProcess {
public:
    BlockedRenderProcess(const RenderJob *parent, RenderQueue *queue,
        int blockSize);

    /**
     * \brief Set the pixel format associated with the rendering process
     *
     * By default, this class works with image blocks that store a set
     * of spectral samples along with an alpha value and an accumulated
     * reconstruction filter weight for each pixel. This method can be
     * used to change this, which is useful when additional information
     * should be returned (e.g. normals or depth values).
     *
     * \param pixelFormat
     *    Desired pixel format (e.g. \ref Bitmap::EMultiChannel)
     * \param channelCount
     *    Number of image channels. Only needs to be specified when
     *    setting <tt>pixelFormat=Bitmap::EMultiChannel</tt>
     * \param warnInvalid
     *    By default, the rendering process issues a warning when writing
     *    negative, infinite or NaN-valued pixels. This flag can be used
     *    to turn off the warnings.
     */
    void setPixelFormat(Bitmap::EPixelFormat pixelFormat,
        int channelCount = -1, bool warnInvalid = false);

    // ======================================================================
    //! @{ \name Implementation of the ParallelProcess interface
    // ======================================================================

    ref<WorkProcessor> createWorkProcessor() const;
    void processResult(const WorkResult *result, bool cancelled);
    void bindResource(const std::string &name, int id);
    EStatus generateWork(WorkUnit *unit, int worker);

    //! @}
    // ======================================================================

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~BlockedRenderProcess();
protected:
    ref<RenderQueue> m_queue;
    ref<Scene> m_scene;
    ref<Film> m_film;
    const RenderJob *m_parent;
    int m_resultCount;
    ref<Mutex> m_resultMutex;
    ProgressReporter *m_progress;
    int m_borderSize;
    Bitmap::EPixelFormat m_pixelFormat;
    int m_channelCount;
    bool m_warnInvalid;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_RENDERPROC_H_ */

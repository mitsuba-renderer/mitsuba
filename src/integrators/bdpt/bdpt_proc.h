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

#if !defined(__BDPT_PROC_H)
#define __BDPT_PROC_H

#include <mitsuba/render/renderproc.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/core/bitmap.h>
#include "bdpt_wr.h"

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                           Parallel process                           */
/* ==================================================================== */

/**
 * \brief Renders work units (rectangular image regions) using
 * bidirectional path tracing
 */
class BDPTProcess : public BlockedRenderProcess {
public:
    BDPTProcess(const RenderJob *parent, RenderQueue *queue,
        const BDPTConfiguration &config);

    inline const BDPTWorkResult *getResult() const { return m_result.get(); }

    /// Develop the image
    void develop();

    /* ParallelProcess impl. */
    void processResult(const WorkResult *wr, bool cancelled);
    ref<WorkProcessor> createWorkProcessor() const;
    void bindResource(const std::string &name, int id);

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~BDPTProcess() { }
private:
    ref<BDPTWorkResult> m_result;
    ref<Timer> m_refreshTimer;
    BDPTConfiguration m_config;
};

MTS_NAMESPACE_END

#endif /* __BDPT_PROC */

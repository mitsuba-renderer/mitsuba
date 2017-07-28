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

#if !defined(__IRRCACHE_PROC_H)
#define __IRRCACHE_PROC_H

#include <mitsuba/render/imageproc.h>
#include <mitsuba/render/rectwu.h>
#include <mitsuba/render/irrcache.h>

MTS_NAMESPACE_BEGIN

/**
 * This stores a number of irradiance samples, which can be sent
 * over the wire as needed. Used to implement parallel overture
 * passes
 */
class IrradianceRecordVector : public WorkResult {
public:
    IrradianceRecordVector() { }

    inline void put(const IrradianceCache::Record *rec) {
        if (rec)
            m_samples.push_back(new IrradianceCache::Record(rec));
    }

    inline size_t size() const {
        return m_samples.size();
    }

    inline void clear() {
        for (size_t i=0; i<m_samples.size(); ++i)
            delete m_samples[i];
        m_samples.clear();
    }

    inline const IrradianceCache::Record *operator[](size_t index) const {
        return m_samples[index];
    }

    /* WorkUnit interface */
    void load(Stream *stream);
    void save(Stream *stream) const;
    std::string toString() const;

    MTS_DECLARE_CLASS()
protected:
    // Virtual destructor
    virtual ~IrradianceRecordVector();
private:
    std::vector<IrradianceCache::Record *> m_samples;
};

/**
 * Parallel process for performing a distributed overture pass
 */
class OvertureProcess : public BlockedImageProcess {
public:
    OvertureProcess(const RenderJob *job, int resolution, bool gradients,
        bool clampNeighbor, bool clampScreen, Float quality);

    inline const IrradianceRecordVector *getSamples() const {
        return m_samples.get();
    }

    ref<WorkProcessor> createWorkProcessor() const;
    void processResult(const WorkResult *wr, bool cancelled);
    void bindResource(const std::string &name, int id);

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~OvertureProcess();
private:
    const RenderJob *m_job;
    ref<Scene> m_scene;
    int m_resultCount;
    ref<Mutex> m_resultMutex;
    ref<IrradianceRecordVector> m_samples;
    int m_resolution;
    bool m_gradients, m_clampNeighbor, m_clampScreen;
    Float m_quality;
    ProgressReporter *m_progress;
};

MTS_NAMESPACE_END

#endif /* __IRRCACHE_PROC_H */

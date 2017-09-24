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

#if !defined(__VCM_BASICS_H)
#define __VCM_BASICS_H

#include <mitsuba/core/plugin.h>
#include <mitsuba/render/renderproc.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/sfcurve.h>
#include <mitsuba/bidir/path.h>
#include "../pssmlt/pssmlt_sampler_base.h"

#if defined(MTS_OPENMP)
#define NANOFLANN_USE_OMP
#endif
#include <mitsuba/core/nanoflann.h>

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                           Parallel process                           */
/* ==================================================================== */

/**
 * \brief Renders work units (rectangular image regions) using
 * bidirectional path tracing
 */

struct VCMPhoton {
    size_t blockID;
    size_t pointID;
    size_t vertexID;
    Point3 pos;
    float radius;
};

struct VCMPhotonMap {
    std::vector<VCMPhoton> photons;

    inline size_t kdtree_get_point_count() const {
        return photons.size();
    }

    inline Float kdtree_distance(const Float* p1, const size_t idx_p2, size_t) const {
        const Float d0 = p1[0] - photons[idx_p2].pos[0];
        const Float d1 = p1[1] - photons[idx_p2].pos[1];
        const Float d2 = p1[2] - photons[idx_p2].pos[2];
        Float sqr_dist = d0 * d0 + d1 * d1 + d2 * d2;
        return sqr_dist;
    }

    inline Float kdtree_get_pt(const size_t idx, int dim) const {
        return photons[idx].pos[dim];
    }

    inline VCMPhoton& kdtree_get_data(const size_t idx) {
        return photons[idx];
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const {
        return false;
    }
};

inline size_t ceil_div(size_t a, size_t b) {
    return (a + b - 1) / b;
}

struct VCMStat {
    Float value;
    unsigned count;

    VCMStat() {
        value = 0.f;
        count = 0;
    }

    void accumulate(Float v, unsigned c) {
        if(c == 0) return;
        value *= count / Float(c + count);
        value += v * c / Float(c + count);
        count += c;
    }
};

class VCMWorkResultBase : public WorkResult {
    using WorkResult::WorkResult;
public:

    inline void clearPhotons() {
        m_photonMap.clear();
    }

    inline void putPhoton(const VCMPhoton& photon) {
        m_photonMap.push_back(photon);
    }

    inline const std::vector<VCMPhoton>& getPhotons() const {
        return m_photonMap;
    }
    virtual void setSize(const Vector2i &size) = 0;

    virtual void setOffset(const Point2i &offset) = 0;

    virtual void clear() = 0;

    void mergeStats(VCMWorkResultBase* result) {
        for (int i = 0; i < stats.size(); i++)
            stats[i].accumulate(result->stats[i].value, result->stats[i].count);
    }

    std::vector<VCMStat> stats;
protected:
    std::vector<VCMPhoton> m_photonMap;
};

class VCMProcessBase : public BlockedRenderProcess {
public:
    using BlockedRenderProcess::BlockedRenderProcess;

    void bindResource(const std::string &name, int id) {
        BlockedRenderProcess::bindResource(name, id);
        if (name == "sensor") {
            if (!m_sensorPathPool.size()) {
                nbx = ceil_div(m_film->getCropSize().x, m_blockSize);
                nby = ceil_div(m_film->getCropSize().y, m_blockSize);
                m_sensorPathPool.resize(nbx * nby);
                m_emitterPathPool.resize(nbx * nby);
            }
        }
    }

    typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<Float, VCMPhotonMap>,
    VCMPhotonMap, 3 /* dim */ > kd_tree_t;

    enum Phase {
        SAMPLE, EVAL
    } phase;
    size_t nbx, nby; // number of blocks to render
    std::vector<CompactBlockPathPool> m_sensorPathPool;
    std::vector<CompactBlockPathPool> m_emitterPathPool;
    VCMPhotonMap m_photonMap;
    std::shared_ptr<kd_tree_t> m_photonKDTree;
    Float m_mergeRadius;

    void buildPhotonLookupStructure() {
        m_photonKDTree.reset(new kd_tree_t(3, m_photonMap, nanoflann::KDTreeSingleIndexAdaptorParams()));
        m_photonKDTree->buildIndex();
    }

    virtual void updateRadius(int n) {
    }

    void clearPhotons() {
        m_photonMap.photons.clear();
    }

    std::vector<VCMPhoton> lookupPhotons(const PathVertex* vertex, float radius) {
        std::vector<VCMPhoton> photons;
        if (vertex->isDegenerate()) return photons;
        const Point &position = vertex->getPosition();
        std::vector<std::pair<size_t, Float> > indices_dists;
        nanoflann::RadiusResultSet<Float> resultSet(radius*radius, indices_dists); // we square here because we used squared distance metric
        m_photonKDTree->findNeighbors(resultSet, (Float*) & position, nanoflann::SearchParams());
        for (const std::pair<size_t, Float>& index_dist : indices_dists) {
            photons.push_back(m_photonMap.photons[index_dist.first]);
        }
        return photons;
    }

    void extractPhotonPath(const VCMPhoton& photon, Path& path, MemoryPool* pool = NULL, bool metropolis = false) {
        if (!metropolis) {
            if (!pool) {
                m_emitterPathPool[photon.blockID].extractPathItem(path, photon.pointID);
                return;
            }
            Path ep;
            m_emitterPathPool[photon.blockID].extractPathItem(ep, photon.pointID);
            ep.clone(path, *pool);
        } else {
            if (!pool) {
                m_sensorPathPool[photon.blockID].extractPathItem(path, photon.pointID);
                return;
            }
            Path sp;
            m_sensorPathPool[photon.blockID].extractPathItem(sp, photon.pointID);
            sp.clone(path, *pool);
        }
    }

    void processResultSample(const VCMWorkResultBase *result) {
        const std::vector<VCMPhoton>& photons = result->getPhotons();
        m_photonMap.photons.insert(m_photonMap.photons.end(), photons.begin(), photons.end());
    }
};

struct VCMConfigBase {
    int maxDepth;
    int rrDepth;
    Float phExponent;
    Float initialRadius;
    bool mergeOnly;
    bool metropolis;
};

class VCMRendererBase : public WorkProcessor {
    using WorkProcessor::WorkProcessor;
protected:

    void extractPathPair(VCMProcessBase* process, Path& emitterSubpath, Path& sensorSubpath, const RectangularWorkUnit *rect, int pointID, bool clone = false) {
        size_t blockSize = m_scene->getBlockSize();
        size_t nbx = process->nbx;
        size_t blockID = (rect->getOffset().y / blockSize) * nbx +
                rect->getOffset().x / blockSize;
        CompactBlockPathPool& sensorPathPool = process->m_sensorPathPool[blockID];
        CompactBlockPathPool& emitterPathPool = process->m_emitterPathPool[blockID];

        if (!clone) {
            sensorPathPool.extractPathItem(sensorSubpath, pointID);
            emitterPathPool.extractPathItem(emitterSubpath, pointID);
            return;
        }

        Path s, e;
        sensorPathPool.extractPathItem(s, pointID);
        emitterPathPool.extractPathItem(e, pointID);
        s.clone(sensorSubpath, m_pool);
        e.clone(emitterSubpath, m_pool);
    }

    void processSampling(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop, VCMProcessBase* m_process, VCMConfigBase* config) {
        VCMConfigBase& m_config = *config;
        const RectangularWorkUnit *rect = static_cast<const RectangularWorkUnit *> (workUnit);
        VCMWorkResultBase *result = static_cast<VCMWorkResultBase *> (workResult);
        bool needsTimeSample = m_sensor->needsTimeSample();
        Float time = m_sensor->getShutterOpen();
        size_t block_size = m_scene->getBlockSize();

        result->setOffset(rect->getOffset());
        result->setSize(rect->getSize());
        result->clear();
        m_hilbertCurve.initialize(TVector2<uint8_t>(rect->getSize()));

        size_t nbx = m_process->nbx;
        size_t blockID = (rect->getOffset().y / block_size) * nbx +
                rect->getOffset().x / block_size;

#if defined(MTS_DEBUG_FP)
        enableFPExceptions();
#endif

        Path emitterSubpath;
        Path sensorSubpath;

        /* Determine the necessary random walk depths based on properties of
           the endpoints */
        int emitterDepth = m_config.maxDepth,
                sensorDepth = m_config.maxDepth;

        /* Go one extra step if the sensor can be intersected */
        if (!m_scene->hasDegenerateSensor() && emitterDepth != -1)
            ++emitterDepth;

        /* Go one extra step if there are emitters that can be intersected */
        if (!m_scene->hasDegenerateEmitters() && sensorDepth != -1)
            ++sensorDepth;

        CompactBlockPathPool& sensorPathPool = m_process->m_sensorPathPool[blockID];
        CompactBlockPathPool& emitterPathPool = m_process->m_emitterPathPool[blockID];

        sensorPathPool.clear();
        emitterPathPool.clear();
        result->clearPhotons();

        for (size_t i = 0; i < m_hilbertCurve.getPointCount(); ++i) {
            Point2i offset = Point2i(m_hilbertCurve[i]) + Vector2i(rect->getOffset());
            m_sampler->generate(offset);
            if (stop)
                break;
            if (needsTimeSample)
                time = m_sensor->sampleTime(m_sampler->next1D());


            if (m_config.metropolis) {
                /* Perform two random walks from the sensor and emitter side */
                sensorSubpath.initialize(m_scene, time, ERadiance, m_pool);
                sensorSubpath.randomWalkFromPixel(m_scene, m_sampler, sensorDepth,
                        offset, m_config.rrDepth, m_pool);

                const Vector2i& image_size = m_sensor->getFilm()->getCropSize();
                size_t nEmitterPaths = image_size.x * image_size.y;

                Float radius = Path::estimateSensorMergingRadius(m_scene, emitterSubpath, sensorSubpath, 0, 2, nEmitterPaths,
                        m_process->m_mergeRadius);

                for (size_t k = 2; k < sensorSubpath.vertexCount(); k++) {
                    if (k > 2) Path::adjustRadius(sensorSubpath.vertexOrNull(k - 1), radius, m_config.mergeOnly);

                    if (radius == 0.f) break;

                    PathVertex* vertex = sensorSubpath.vertex(k);
                    if (!vertex->isSurfaceInteraction()) continue;
                    if (vertex->isDegenerate()) continue;
                    // add a new photon
                    VCMPhoton photon;
                    photon.pos = vertex->getPosition();
                    photon.blockID = blockID;
                    photon.pointID = i;
                    photon.vertexID = k;
                    photon.radius = radius;
                    result->putPhoton(photon);
                }

                sensorPathPool.addPathItem(sensorSubpath); // cache path into sensor Path Pool
                sensorSubpath.release(m_pool);
            } else {
                /* Start new emitter and sensor subpaths */
                emitterSubpath.initialize(m_scene, time, EImportance, m_pool);
                sensorSubpath.initialize(m_scene, time, ERadiance, m_pool);
                /* Perform a random walk using alternating steps on each path */
                Path::alternatingRandomWalkFromPixel(m_scene, m_sampler,
                        emitterSubpath, emitterDepth, sensorSubpath,
                        sensorDepth, offset, m_config.rrDepth, m_pool);

                for (size_t k = 2; k < emitterSubpath.vertexCount(); k++) {
                    PathVertex* vertex = emitterSubpath.vertex(k);
                    if (!vertex->isSurfaceInteraction()) continue;
                    if (vertex->isDegenerate()) continue;
                    // add a new photon
                    VCMPhoton photon;
                    photon.pos = vertex->getPosition();
                    photon.blockID = blockID;
                    photon.pointID = i;
                    photon.vertexID = k;
                    result->putPhoton(photon);
                }

                sensorPathPool.addPathItem(sensorSubpath); // cache path into sensor Path Pool
                emitterPathPool.addPathItem(emitterSubpath); // cache path into emitter Path Pool
                emitterSubpath.release(m_pool);
                sensorSubpath.release(m_pool);
            }
        }
        sensorPathPool.buildIndex();
        emitterPathPool.buildIndex();
    }

    ref<Scene> m_scene;
    ref<Sensor> m_sensor;
    ref<Sampler> m_sampler;
    ref<PSSMLTSamplerBase> m_pssmltSampler;
    HilbertCurve2D<uint8_t> m_hilbertCurve;
    MemoryPool m_pool;
};

class VCMIntegratorBase : public Integrator {
protected:
    bool m_cancelled;
    using Integrator::Integrator;
public:

    bool iterateVCM(VCMProcessBase* process, int sensorResID, int iter) {
        ref<Scheduler> scheduler = Scheduler::getInstance();
        process->updateRadius(iter + 1);
        process->clearPhotons();
        if (m_cancelled) return false;
        // connection phase
        // we also collect photons in this phase

        process->bindResource("sensor", sensorResID);
        process->phase = VCMProcessBase::SAMPLE;
        scheduler->schedule(process);
        scheduler->wait(process);
        if (m_cancelled) return false;

        // build photon look up structure
        process->buildPhotonLookupStructure();

        // merging phase
        process->bindResource("sensor", sensorResID);
        process->phase = VCMProcessBase::EVAL;
        scheduler->schedule(process);
        scheduler->wait(process);
        if (m_cancelled) return false;
        return true;
    }
};



MTS_NAMESPACE_END

#endif /* __VCM_PROC */

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

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>
#include "irrcache_proc.h"

MTS_NAMESPACE_BEGIN

/* Parallel overture pass implementation (worker) */
class OvertureWorker : public WorkProcessor {
public:
    OvertureWorker(int resolution, bool gradients, bool clampNeighbor,
        bool clampScreen, Float quality) : m_resolution(resolution), m_gradients(gradients),
        m_clampNeighbor(clampNeighbor), m_clampScreen(clampScreen),
        m_quality(quality) {
    }

    OvertureWorker(Stream *stream, InstanceManager *manager) {
        m_resolution = stream->readInt();
        m_gradients = stream->readBool();
        m_clampNeighbor = stream->readBool();
        m_clampScreen = stream->readBool();
        m_quality = stream->readFloat();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        stream->writeInt(m_resolution);
        stream->writeBool(m_gradients);
        stream->writeBool(m_clampNeighbor);
        stream->writeBool(m_clampScreen);
        stream->writeFloat(m_quality);
    }

    ref<WorkUnit> createWorkUnit() const {
        return new RectangularWorkUnit();
    }

    ref<WorkResult> createWorkResult() const {
        return new IrradianceRecordVector();
    }

    void prepare() {
        m_scene = static_cast<Scene *>(getResource("scene"));
        m_sensor = static_cast<Sensor *>(getResource("sensor"));
        m_subIntegrator = static_cast<SamplingIntegrator *>(getResource("subIntegrator"));
        Properties props("independent");
        props.setInteger("sampleCount", m_resolution * 3 * m_resolution);
        m_sampler = static_cast<Sampler *> (PluginManager::getInstance()->
            createObject(MTS_CLASS(Sampler), props));
        m_subIntegrator->wakeup(NULL, m_resources);

        m_irrCache = new IrradianceCache(m_scene->getAABB());
        m_irrCache->clampNeighbor(m_clampNeighbor);
        m_irrCache->clampScreen(m_clampScreen);
        m_irrCache->useGradients(m_gradients);
        m_irrCache->setQuality(m_quality);
        m_hs = new HemisphereSampler(m_resolution, 3*m_resolution);
    }

    void process(const WorkUnit *workUnit, WorkResult *workResult,
        const bool &stop) {
        const RectangularWorkUnit *rect = static_cast<const RectangularWorkUnit *>(workUnit);
        IrradianceRecordVector *result = static_cast<IrradianceRecordVector *>(workResult);
        const SamplingIntegrator *integrator = m_subIntegrator.get();

        Intersection its;
        RadianceQueryRecord rRec(m_scene, m_sampler);

        const int sx = rect->getOffset().x,
                sy = rect->getOffset().y,
                ex = sx + rect->getSize().x,
                ey = sy + rect->getSize().y;
        result->clear();

        for (int y = sy; y < ey; y++) {
            for (int x = sx; x < ex; x++) {
                if (stop)
                    break;
                Point2 pixelSample(x + .5f, y + .5f);
                RayDifferential ray;
                if (m_sensor->sampleRayDifferential(ray, pixelSample, Point2(0.5f), 0.5f).isZero())
                    continue;
                if (m_scene->rayIntersect(ray, its)) {
                    const BSDF *bsdf = its.getBSDF();
                    if ((bsdf->getType() & BSDF::EAll) != BSDF::EDiffuseReflection)
                        continue;

                    Spectrum E;
                    if (m_irrCache->get(its, E))
                        continue;

                    /* Irradiance cache miss - create a record. The following
                       generates stratified cosine-weighted samples and computes
                       rotational + translational gradients */
                    m_hs->generateDirections(its);
                    m_sampler->generate(Point2i(0));

                    for (unsigned int j=0; j<m_hs->getM(); j++) {
                        for (unsigned int k=0; k<m_hs->getN(); k++) {
                            HemisphereSampler::SampleEntry &entry = (*m_hs)(j, k);
                            entry.dist = std::numeric_limits<Float>::infinity();
                            rRec.newQuery(RadianceQueryRecord::ERadianceNoEmission
                                | RadianceQueryRecord::EDistance, m_sensor->getMedium());
                            rRec.extra = RadianceQueryRecord::ECacheQuery;
                            rRec.depth = 2;
                            entry.L = integrator->Li(RayDifferential(its.p, entry.d, 0.0f), rRec);
                            entry.dist = rRec.dist;
                            m_sampler->advance();
                        }
                    }

                    m_hs->process(its);
                    result->put(m_irrCache->put(ray, its, *m_hs));
                }
            }
        }
    }

    ref<WorkProcessor> clone() const {
        return new OvertureWorker(m_resolution, m_gradients, m_clampNeighbor,
            m_clampScreen, m_quality);
    }

    MTS_DECLARE_CLASS()
protected:
    virtual ~OvertureWorker() { }
private:
    ref<Scene> m_scene;
    ref<Sensor> m_sensor;
    ref<Sampler> m_sampler;
    ref<HemisphereSampler> m_hs;
    ref<SamplingIntegrator> m_subIntegrator;
    ref<IrradianceCache> m_irrCache;
    int m_resolution;
    bool m_gradients, m_clampNeighbor, m_clampScreen;
    Float m_quality;
};

void IrradianceRecordVector::load(Stream *stream) {
    clear();
    size_t count = stream->readUInt();
    m_samples.resize(count);
    for (size_t i=0; i<count; ++i)
        m_samples[i] = new IrradianceCache::Record(stream);
}

IrradianceRecordVector::~IrradianceRecordVector() {
    clear();
}

void IrradianceRecordVector::save(Stream *stream) const {
    stream->writeUInt((unsigned int) m_samples.size());
    for (size_t i=0; i<m_samples.size(); ++i)
        m_samples[i]->serialize(stream);
}

std::string IrradianceRecordVector::toString() const {
    std::ostringstream oss;
    oss << "IrradianceRecordVector[size="
        << m_samples.size() << "]";
    return oss.str();
}

OvertureProcess::OvertureProcess(const RenderJob *job, int resolution, bool gradients,
    bool clampNeighbor, bool clampScreen, Float quality) : m_job(job), m_resolution(resolution),
    m_gradients(gradients), m_clampNeighbor(clampNeighbor),
    m_clampScreen(clampScreen), m_quality(quality), m_progress(NULL) {
    m_resultCount = 0;
    m_resultMutex = new Mutex();
    m_samples = new IrradianceRecordVector();
}

OvertureProcess::~OvertureProcess() {
    if (m_progress)
        delete m_progress;
}

ref<WorkProcessor> OvertureProcess::createWorkProcessor() const {
    return new OvertureWorker(m_resolution, m_gradients, m_clampNeighbor,
        m_clampScreen, m_quality);
}

void OvertureProcess::processResult(const WorkResult *wr, bool cancelled) {
    const IrradianceRecordVector *result = static_cast<const IrradianceRecordVector *>(wr);
    LockGuard lock(m_resultMutex);
    for (size_t i=0; i<result->size(); ++i)
        m_samples->put((*result)[i]);
    m_progress->update(++m_resultCount);
}

void OvertureProcess::bindResource(const std::string &name, int id) {
    if (name == "scene") {
        m_scene = static_cast<Scene *>(Scheduler::getInstance()->getResource(id));
        const Film *film = m_scene->getFilm();
        Point2i offset = film->getCropOffset();
        Vector2i size = film->getCropSize();
        BlockedImageProcess::init(offset, size, 64);
        if (m_progress)
            delete m_progress;
        m_progress = new ProgressReporter("Overture pass", m_numBlocksTotal, m_job);
    }
    BlockedImageProcess::bindResource(name, id);
}

MTS_IMPLEMENT_CLASS(IrradianceRecordVector, false, WorkResult);
MTS_IMPLEMENT_CLASS_S(OvertureWorker, false, WorkProcessor);
MTS_IMPLEMENT_CLASS(OvertureProcess, false, BlockedImageProcess);
MTS_NAMESPACE_END

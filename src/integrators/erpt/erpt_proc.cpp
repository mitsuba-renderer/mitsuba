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

#include <mitsuba/bidir/mut_bidir.h>
#include <mitsuba/bidir/mut_lens.h>
#include <mitsuba/bidir/mut_caustic.h>
#include <mitsuba/bidir/mut_mchain.h>
#include <mitsuba/bidir/mut_manifold.h>
#include <mitsuba/bidir/pathsampler.h>
#include <mitsuba/bidir/util.h>
#include <mitsuba/core/sfcurve.h>
#include <boost/bind.hpp>
#include "erpt_proc.h"

//#define MTS_BD_DEBUG_HEAVY

MTS_NAMESPACE_BEGIN

static StatsCounter statsAccepted("Energy redistribution path tracing",
        "Accepted mutations", EPercentage);
static StatsCounter statsChainsPerPixel("Energy redistribution path tracing",
        "Chains started per pixel", EAverage);

/* ==================================================================== */
/*                      Worker result implementation                    */
/* ==================================================================== */

class ERPTWorkResult : public ImageBlock {
public:
    Point2i origOffset;
    Vector2i origSize;

    ERPTWorkResult(const Vector2i &size, const ReconstructionFilter *filter)
        : ImageBlock(Bitmap::ESpectrum, size, filter) { }

    void load(Stream *stream) {
        ImageBlock::load(stream);
        origOffset = Point2i(stream);
        origSize = Vector2i(stream);
    }

    void save(Stream *stream) const {
        ImageBlock::save(stream);
        origOffset.serialize(stream);
        origSize.serialize(stream);
    }
};

/* ==================================================================== */
/*                         Worker implementation                        */
/* ==================================================================== */

class ERPTRenderer : public WorkProcessor {
public:
    ERPTRenderer(const ERPTConfiguration &conf)
        : m_config(conf) {
    }

    ERPTRenderer(Stream *stream, InstanceManager *manager)
        : WorkProcessor(stream, manager) {
        m_config = ERPTConfiguration(stream);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        m_config.serialize(stream);
    }

    ref<WorkUnit> createWorkUnit() const {
        return new RectangularWorkUnit();
    }

    ref<WorkResult> createWorkResult() const {
        return new ERPTWorkResult(
            m_sensor->getFilm()->getCropSize(),
            m_sensor->getFilm()->getReconstructionFilter());
    }

    void prepare() {
        Scene *scene = static_cast<Scene *>(getResource("scene"));
        m_scene = new Scene(scene);
        m_sampler = static_cast<Sampler *>(getResource("sampler"));
        m_indepSampler = static_cast<Sampler *>(getResource("indepSampler"));
        m_sensor = static_cast<Sensor *>(getResource("sensor"));
        m_scene->removeSensor(scene->getSensor());
        m_scene->addSensor(m_sensor);
        m_scene->setSensor(m_sensor);
        m_scene->setSampler(m_sampler);
        m_scene->wakeup(NULL, m_resources);
        m_scene->initializeBidirectional();

        m_pathSampler = new PathSampler(PathSampler::EBidirectional, m_scene,
            m_sampler, m_sampler, m_sampler, m_config.maxDepth, 10,
            m_config.separateDirect, true, true);

        m_pool = &m_pathSampler->getMemoryPool();

        /* Jump sizes recommended by Eric Veach */
        Float minJump = 0.1f, coveredArea = 0.05f;

        /* Register all available mutators */
        if (m_config.bidirectionalMutation)
            m_mutators.push_back(new BidirectionalMutator(m_scene, m_indepSampler,
                *m_pool, 3, m_config.maxDepth == -1 ? INT_MAX : m_config.maxDepth + 2));

        if (m_config.lensPerturbation)
            m_mutators.push_back(new LensPerturbation(m_scene, m_indepSampler, *m_pool,
                minJump, coveredArea));

        if (m_config.multiChainPerturbation)
            m_mutators.push_back(new MultiChainPerturbation(m_scene, m_indepSampler, *m_pool,
                minJump, coveredArea));

        if (m_config.causticPerturbation)
            m_mutators.push_back(new CausticPerturbation(m_scene, m_indepSampler, *m_pool,
                minJump, coveredArea));

        if (m_config.manifoldPerturbation)
            m_mutators.push_back(new ManifoldPerturbation(m_scene, m_indepSampler, *m_pool,
                m_config.probFactor, true, true,
                m_config.avgAngleChangeSurface,
                m_config.avgAngleChangeMedium));

        if (m_mutators.size() == 0)
            Log(EError, "There must be at least one mutator!");
    }

    void pathCallback(int s, int t, Float weight, Path &path, const bool *stop) {
        if (std::isnan(weight) || std::isinf(weight) || weight < 0)
            Log(EWarn, "Invalid path weight: %f, ignoring path!", weight);

#if 0
        /* Don't run ERPT on paths that start with two diffuse vertices. It's
           usually safe to assume that these are handled well enough by BDPT */
        int k = path.length();
        if (path.vertex(k-2)->isDiffuseInteraction() && path.vertex(k-3)->isDiffuseInteraction()) {
            Spectrum value = path.getRelativeWeight() * weight / m_sampler->getSampleCount();
            m_result->put(path.getSamplePosition(), &value[0]);
            return;
        }
#endif

        Float meanChains = m_config.numChains * weight
            / (m_config.luminance * m_sampler->getSampleCount());

        /* Optional: do not launch too many chains if this is desired by the user */
        if (m_config.maxChains > 0 && meanChains > m_config.maxChains)
            meanChains = std::min(meanChains, (Float) m_config.maxChains);

        /* Decide the actual number of chains that will be launched, as well
           as their deposition energy */
        int numChains = (int) std::floor(m_indepSampler->next1D() + meanChains);
        if (numChains == 0)
            return;

        Float depositionEnergy = weight / (m_sampler->getSampleCount()
                * meanChains * m_config.chainLength);

        DiscreteDistribution suitabilities(m_mutators.size());
        std::ostringstream oss;
        Spectrum relWeight(0.0f);
        Float accumulatedWeight = 0;
        MutationRecord muRec, currentMuRec(Mutator::EMutationTypeCount, 0, 0, 0, Spectrum(0.f));
        Path *current = new Path(),
             *proposed = new Path();
        size_t mutations = 0;

        #if defined(MTS_BD_DEBUG_HEAVY)
            if (!path.verify(m_scene, EImportance, oss))
                Log(EError, "Started ERPT with an invalid path: %s", oss.str().c_str());
        #endif

        for (int chain=0; chain<numChains && !*stop; ++chain) {
            relWeight = path.getRelativeWeight();
            path.clone(*current, *m_pool);
            accumulatedWeight = 0;
            ++statsChainsPerPixel;

            for (size_t it=0; it<m_config.chainLength; ++it) {
                /* Query all mutators for their suitability */
                suitabilities.clear();
                for (size_t j=0; j<m_mutators.size(); ++j)
                    suitabilities.append(m_mutators[j]->suitability(*current));

                /* Pick a mutator according to the suitabilities */
                int mutatorIdx = -1;
                bool success = false;
                Mutator *mutator = NULL;

                if (suitabilities.normalize() == 0) {
                    /* No mutator can handle this path -- give up */
                    accumulatedWeight += m_config.chainLength - it;
                    break;
                }

                mutatorIdx = (int) suitabilities.sample(m_indepSampler->next1D());
                mutator = m_mutators[mutatorIdx].get();

                /* Sample a mutated path */
                success = mutator->sampleMutation(*current, *proposed, muRec, currentMuRec);

                statsAccepted.incrementBase(1);
                if (success) {
                    Float Qxy = mutator->Q(*current, *proposed, muRec) * suitabilities[mutatorIdx];
                    suitabilities.clear();
                    for (size_t j=0; j<m_mutators.size(); ++j)
                        suitabilities.append(m_mutators[j]->suitability(*proposed));
                    suitabilities.normalize();
                    Float Qyx = mutator->Q(*proposed, *current, muRec.reverse()) * suitabilities[mutatorIdx];
                    Float a = std::min((Float) 1, Qyx / Qxy);

                    #if defined(MTS_BD_DEBUG_HEAVY)
                        if (!proposed->verify(m_scene, EImportance, oss)) {
                            Log(EWarn, "%s proposed as %s, Qxy=%f, Qyx=%f", oss.str().c_str(),
                                    muRec.toString().c_str(), Qxy, Qyx);
                            Log(EWarn, "Original path: %s", current->toString().c_str());
                            proposed->release(muRec.l, muRec.l + muRec.ka + 1, *m_pool);
                            oss.str("");
                            continue;
                        }
                    #endif

                    if (Qxy == 0) { // be tolerant of this (can occasionally happen due to floating point inaccuracies)
                        a = 0;
                    } else if (Qxy < 0 || Qyx < 0 || std::isnan(Qxy) || std::isnan(Qyx)) {
                        #if defined(MTS_BD_DEBUG)
                            Log(EDebug, "Source path: %s", current->toString().c_str());
                            Log(EDebug, "Proposal path: %s", proposed->toString().c_str());
                            Log(EWarn, "Internal error while computing acceptance probabilities: "
                                "Qxy=%f, Qyx=%f, muRec=%s", Qxy, Qyx, muRec.toString().c_str());
                        #endif
                        a = 0;
                    }

                    accumulatedWeight += 1-a;

                    /* Accept with probability 'a' */
                    if (a == 1 || m_indepSampler->next1D() < a) {
                        Spectrum value = relWeight * (accumulatedWeight * depositionEnergy);
                        m_result->put(current->getSamplePosition(), &value[0]);

                        /* The mutation was accepted */
                        current->release(muRec.l, muRec.m+1, *m_pool);
                        std::swap(current, proposed);
                        relWeight = current->getRelativeWeight();
                        mutator->accept(muRec);
                        currentMuRec = muRec;
                        accumulatedWeight = a;
                        ++statsAccepted;
                        ++mutations;
                    } else {
                        if (a > 0) {
                            Spectrum value = proposed->getRelativeWeight() * (a * depositionEnergy);
                            m_result->put(proposed->getSamplePosition(), &value[0]);
                        }
                        /* The mutation was rejected */
                        proposed->release(muRec.l, muRec.l + muRec.ka + 1, *m_pool);
                    }
                } else {
                    accumulatedWeight += 1;
                }
            }
            if (accumulatedWeight > 0) {
                Spectrum value = relWeight * (accumulatedWeight * depositionEnergy);
                m_result->put(current->getSamplePosition(), &value[0]);
            }
            current->release(*m_pool);
        }

        /*if (mutations == 0) {
            cout << "Path was never mutated: " << path.summarize() << endl;
            cout << path.toString() << endl;
        }*/

        delete current;
        delete proposed;
    }

    void process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop) {
        const RectangularWorkUnit *rect = static_cast<const RectangularWorkUnit *>(workUnit);
        m_result = static_cast<ERPTWorkResult *>(workResult);
        m_result->origOffset = rect->getOffset();
        m_result->origSize = rect->getSize();

        m_hilbertCurve.initialize(TVector2<uint8_t>(rect->getSize()));
        m_result->clear();
        boost::function<void (int, int, Float, Path &)> callback
            = boost::bind(&ERPTRenderer::pathCallback, this, _1, _2, _3, _4, &stop);

        for (size_t i=0; i<m_hilbertCurve.getPointCount(); ++i) {
            if (stop)
                break;

            statsChainsPerPixel.incrementBase();

            Point2i offset = Point2i(m_hilbertCurve[i]) + Vector2i(rect->getOffset());
            m_sampler->generate(offset);

            for (size_t j = 0; j<m_sampler->getSampleCount(); j++) {
                m_pathSampler->samplePaths(offset, callback);
                m_sampler->advance();
            }
        }

        if (!m_pool->unused())
            Log(EError, "Internal error: detected a memory pool leak!");
        m_result = NULL;
    }

    ref<WorkProcessor> clone() const {
        return new ERPTRenderer(m_config);
    }

    MTS_DECLARE_CLASS()
private:
    ERPTConfiguration m_config;
    ref<Sensor> m_sensor;
    ref<Scene> m_scene;
    ref<Sampler> m_sampler, m_indepSampler;
    ref<PathSampler> m_pathSampler;
    ref_vector<Mutator> m_mutators;
    HilbertCurve2D<uint8_t> m_hilbertCurve;
    ERPTWorkResult *m_result;
    MemoryPool *m_pool;
};

/* ==================================================================== */
/*                           Parallel process                           */
/* ==================================================================== */

ERPTProcess::ERPTProcess(const RenderJob *job, RenderQueue *queue,
    const ERPTConfiguration &conf, const Bitmap *directImage)
    : BlockedRenderProcess(job, queue, conf.blockSize), m_job(job), m_config(conf) {
    m_directImage = directImage;
}

ref<WorkProcessor> ERPTProcess::createWorkProcessor() const {
    return new ERPTRenderer(m_config);
}

void ERPTProcess::develop() {
    LockGuard lock(m_resultMutex);
    m_film->setBitmap(m_accum->getBitmap());
    if (m_directImage)
        m_film->addBitmap(m_directImage);
    m_queue->signalRefresh(m_job);
}

void ERPTProcess::processResult(const WorkResult *wr, bool cancelled) {
    const ERPTWorkResult *result = static_cast<const ERPTWorkResult *>(wr);
    UniqueLock lock(m_resultMutex);
    m_progress->update(++m_resultCount);
    m_accum->put(result);
    develop();
    lock.unlock();
    m_queue->signalWorkCanceled(m_parent, result->origOffset, result->origSize);
}

void ERPTProcess::bindResource(const std::string &name, int id) {
    BlockedRenderProcess::bindResource(name, id);
    if (name == "sensor") {
        Film *film = static_cast<Sensor *>(Scheduler::getInstance()->getResource(id))->getFilm();

        m_accum = new ImageBlock(Bitmap::ESpectrum, film->getCropSize());
        m_accum->clear();
    }
}

MTS_IMPLEMENT_CLASS_S(ERPTRenderer, false, WorkProcessor)
MTS_IMPLEMENT_CLASS(ERPTProcess, false, BlockedRenderProcess)

MTS_NAMESPACE_END

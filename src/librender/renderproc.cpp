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

#include <mitsuba/core/statistics.h>
#include <mitsuba/core/sfcurve.h>
#include <mitsuba/render/renderproc.h>
#include <mitsuba/render/rectwu.h>

MTS_NAMESPACE_BEGIN

class BlockRenderer : public WorkProcessor {
public:
    BlockRenderer(Bitmap::EPixelFormat pixelFormat, int channelCount, int blockSize,
        int borderSize, bool warnInvalid) : m_pixelFormat(pixelFormat),
        m_channelCount(channelCount), m_blockSize(blockSize),
        m_borderSize(borderSize), m_warnInvalid(warnInvalid) { }

    BlockRenderer(Stream *stream, InstanceManager *manager) {
        m_pixelFormat = (Bitmap::EPixelFormat) stream->readInt();
        m_channelCount = stream->readInt();
        m_blockSize = stream->readInt();
        m_borderSize = stream->readInt();
        m_warnInvalid = stream->readBool();
    }

    ref<WorkUnit> createWorkUnit() const {
        return new RectangularWorkUnit();
    }

    ref<WorkResult> createWorkResult() const {
        return new ImageBlock(m_pixelFormat,
            Vector2i(m_blockSize),
            m_sensor->getFilm()->getReconstructionFilter(),
            m_channelCount, m_warnInvalid);
    }

    void prepare() {
        Scene *scene = static_cast<Scene *>(getResource("scene"));
        m_scene = new Scene(scene);
        m_sampler = static_cast<Sampler *>(getResource("sampler"));
        m_sensor = static_cast<Sensor *>(getResource("sensor"));
        m_integrator = static_cast<SamplingIntegrator *>(getResource("integrator"));
        m_scene->removeSensor(scene->getSensor());
        m_scene->addSensor(m_sensor);
        m_scene->setSensor(m_sensor);
        m_scene->setSampler(m_sampler);
        m_scene->setIntegrator(m_integrator);
        m_integrator->wakeup(m_scene, m_resources);
        m_scene->wakeup(m_scene, m_resources);
        m_scene->initializeBidirectional();
    }

    void process(const WorkUnit *workUnit, WorkResult *workResult,
        const bool &stop) {
        const RectangularWorkUnit *rect = static_cast<const RectangularWorkUnit *>(workUnit);
        ImageBlock *block = static_cast<ImageBlock *>(workResult);

#ifdef MTS_DEBUG_FP
        enableFPExceptions();
#endif

        block->setOffset(rect->getOffset());
        block->setSize(rect->getSize());
        m_hilbertCurve.initialize(TVector2<uint8_t>(rect->getSize()));
        m_integrator->renderBlock(m_scene, m_sensor, m_sampler,
            block, stop, m_hilbertCurve.getPoints());

#ifdef MTS_DEBUG_FP
        disableFPExceptions();
#endif
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        stream->writeInt(m_pixelFormat);
        stream->writeInt(m_channelCount);
        stream->writeInt(m_blockSize);
        stream->writeInt(m_borderSize);
        stream->writeBool(m_warnInvalid);
    }

    ref<WorkProcessor> clone() const {
        return new BlockRenderer(m_pixelFormat, m_channelCount,
            m_blockSize, m_borderSize, m_warnInvalid);
    }

    MTS_DECLARE_CLASS()
protected:
    virtual ~BlockRenderer() { }
private:
    ref<Scene> m_scene;
    ref<Sensor> m_sensor;
    ref<Sampler> m_sampler;
    ref<SamplingIntegrator> m_integrator;
    Bitmap::EPixelFormat m_pixelFormat;
    int m_channelCount;
    int m_blockSize;
    int m_borderSize;
    bool m_warnInvalid;
    HilbertCurve2D<uint8_t> m_hilbertCurve;
};

BlockedRenderProcess::BlockedRenderProcess(const RenderJob *parent, RenderQueue *queue,
        int blockSize) : m_queue(queue), m_parent(parent), m_resultCount(0), m_progress(NULL) {
    m_blockSize = blockSize;
    m_resultMutex = new Mutex();
    m_pixelFormat = Bitmap::ESpectrumAlphaWeight;
    m_channelCount = -1;
    m_warnInvalid = true;
}

BlockedRenderProcess::~BlockedRenderProcess() {
    if (m_progress)
        delete m_progress;
}

void BlockedRenderProcess::setPixelFormat(Bitmap::EPixelFormat pixelFormat, int channelCount, bool warnInvalid) {
    m_pixelFormat = pixelFormat;
    m_channelCount = channelCount;
    m_warnInvalid = warnInvalid;
}

ref<WorkProcessor> BlockedRenderProcess::createWorkProcessor() const {
    return new BlockRenderer(m_pixelFormat, m_channelCount,
            m_blockSize, m_borderSize, m_warnInvalid);
}

void BlockedRenderProcess::processResult(const WorkResult *result, bool cancelled) {
    const ImageBlock *block = static_cast<const ImageBlock *>(result);
    UniqueLock lock(m_resultMutex);
    m_film->put(block);
    m_progress->update(++m_resultCount);
    lock.unlock();
    m_queue->signalWorkEnd(m_parent, block, cancelled);
}

ParallelProcess::EStatus BlockedRenderProcess::generateWork(WorkUnit *unit, int worker) {
    EStatus status = BlockedImageProcess::generateWork(unit, worker);
    if (status == ESuccess)
        m_queue->signalWorkBegin(m_parent, static_cast<RectangularWorkUnit *>(unit), worker);
    return status;
}

void BlockedRenderProcess::bindResource(const std::string &name, int id) {
    if (name == "sensor") {
        m_film = static_cast<Sensor *>(Scheduler::getInstance()->getResource(id))->getFilm();
        m_borderSize = m_film->getReconstructionFilter()->getBorderSize();

        Point2i offset = Point2i(0, 0);
        Vector2i size = m_film->getCropSize();

        if (m_film->hasHighQualityEdges()) {
            offset.x -= m_borderSize;
            offset.y -= m_borderSize;
            size.x += 2 * m_borderSize;
            size.y += 2 * m_borderSize;
        }

        if (m_blockSize < m_borderSize)
            Log(EError, "The block size must be larger than the image reconstruction filter radius!");

        BlockedImageProcess::init(offset, size, m_blockSize);
        if (m_progress)
            delete m_progress;
        m_progress = new ProgressReporter("Rendering", m_numBlocksTotal, m_parent);
    }
    BlockedImageProcess::bindResource(name, id);
}

MTS_IMPLEMENT_CLASS(BlockedRenderProcess, false, BlockedImageProcess)
MTS_IMPLEMENT_CLASS_S(BlockRenderer, false, WorkProcessor)
MTS_NAMESPACE_END

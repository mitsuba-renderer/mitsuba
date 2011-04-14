/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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
	BlockRenderer(int blockSize, int borderSize) 
	 : m_blockSize(blockSize), m_borderSize(borderSize) {
	}

	BlockRenderer(Stream *stream, InstanceManager *manager) {
		m_blockSize = stream->readInt();
		m_borderSize = stream->readInt();
		m_collectStatistics = stream->readBool();
	}

	ref<WorkUnit> createWorkUnit() const {
		return new RectangularWorkUnit();
	}

	ref<WorkResult> createWorkResult() const {
		return new ImageBlock(Vector2i(m_blockSize, m_blockSize), m_borderSize, 
			true, true, true, m_collectStatistics);
	}

	void prepare() {
		m_scene = new Scene(static_cast<Scene *>(getResource("scene")));
		/// Variance estimates are required when executing a T-test on the rendered data
		m_collectStatistics = (m_scene->getTestType() == Scene::ETTest);
		m_sampler = static_cast<Sampler *>(getResource("sampler"));
		m_camera = static_cast<Camera *>(getResource("camera"));
		m_integrator = static_cast<SampleIntegrator *>(getResource("integrator"));
		m_scene->setCamera(m_camera);
		m_scene->setSampler(m_sampler);
		m_scene->setIntegrator(m_integrator);
		m_integrator->wakeup(m_resources);
		m_scene->wakeup(m_resources);
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
		m_hilbertCurve.initialize(rect->getSize());
		m_integrator->renderBlock(m_scene, m_camera, m_sampler, 
			block, stop, &m_hilbertCurve.getPoints());

#ifdef MTS_DEBUG_FP
		disableFPExceptions();
#endif
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		stream->writeInt(m_blockSize);
		stream->writeInt(m_borderSize);
		stream->writeBool(m_collectStatistics);
	}

	ref<WorkProcessor> clone() const {
		return new BlockRenderer(m_blockSize, m_borderSize);
	}

	MTS_DECLARE_CLASS()
protected:
	virtual ~BlockRenderer() { }
private:
	ref<Scene> m_scene;
	ref<Camera> m_camera;
	ref<Sampler> m_sampler;
	ref<SampleIntegrator> m_integrator;
	int m_blockSize;
	int m_borderSize;
	int m_collectStatistics;
	HilbertCurve2D<int> m_hilbertCurve;
};


BlockedRenderProcess::BlockedRenderProcess(const RenderJob *parent, RenderQueue *queue,
		int blockSize) : m_queue(queue), m_progress(NULL) {
	m_blockSize = blockSize;
	m_parent = parent;
	m_resultCount = 0;
	m_resultMutex = new Mutex();
}

BlockedRenderProcess::~BlockedRenderProcess() {
	if (m_progress)
		delete m_progress;
}
	
ref<WorkProcessor> BlockedRenderProcess::createWorkProcessor() const {
	return new BlockRenderer(m_blockSize, m_borderSize);
}

void BlockedRenderProcess::processResult(const WorkResult *result, bool cancelled) {
	const ImageBlock *block = static_cast<const ImageBlock *>(result);
	m_resultMutex->lock();
	m_film->putImageBlock(block);
	m_progress->update(++m_resultCount);
	m_resultMutex->unlock();
	m_queue->signalWorkEnd(m_parent, block);
}

ParallelProcess::EStatus BlockedRenderProcess::generateWork(WorkUnit *unit, int worker) {
	EStatus status = BlockedImageProcess::generateWork(unit, worker);
	if (status == ESuccess)
		m_queue->signalWorkBegin(m_parent, static_cast<RectangularWorkUnit *>(unit), worker);
	return status;
}

void BlockedRenderProcess::bindResource(const std::string &name, int id) {
	if (name == "camera") {
		m_film = static_cast<Camera *>(Scheduler::getInstance()->getResource(id))->getFilm();
		const TabulatedFilter *filter = m_film->getTabulatedFilter();
		m_borderSize = (int) std::ceil(std::max(filter->getFilterSize().x,
			filter->getFilterSize().y) - (Float) 0.5);

		Point2i offset = m_film->getCropOffset();
		Vector2i size = m_film->getCropSize();
		if (m_film->hasHighQualityEdges()) {	
			offset.x -= m_borderSize;
			offset.y -= m_borderSize;
			size.x += 2 * m_borderSize;
			size.y += 2 * m_borderSize;
		}
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

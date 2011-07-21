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

#include "ptracer_proc.h"

MTS_NAMESPACE_BEGIN

/**
 * Particle tracer using an adjoint formulation -- meant primarily for
 * verification purposes and test cases. This class follows appendix 
 * 4.A of Eric Veach's PhD thesis and computes the inner product 
 * between emitted radiance and incident importance at the light 
 * source (e.g. I = <Le, Wi>). The importance is recursively 
 * estimated using a Monte Carlo random walk and adjoint BSDFs are
 * used for scattering events. Non-symmetric behavior due to the
 * use of shading normals is handled correctly.
 * For practical reasons, the integral is simultaneously computed
 * for every pixel on the image plane. This is done similarly to 
 * path tracing with next event estimation (PTNEE) by tracing a
 * shadow ray to the camera at every surface/volume interaction.
 * An independent accumulation buffer will be assigned to each
 * processor so that the rendering process can run in parallel.
 * When a perspective camera is used, the importance distribution 
 * on the camera sensor is chosen to be equivalent to a perspective 
 * camera used in a traditional backward ray tracer with uniform 
 * sampling on the image plane. Orthographic cameras are also handled.
 * The approximate number of samples/pixel is specified by the camera's 
 * sampler instance.
 */
class AdjointParticleTracer : public Integrator {
public:
	AdjointParticleTracer(const Properties &props) : Integrator(props) {
		/* Depth to start using russian roulette */
		m_rrDepth = props.getInteger("rrDepth", 10);

		/* Longest visualized path length (<tt>-1</tt>=infinite).
		   A value of <tt>1</tt> will produce a black image, since this integrator
		   does not visualize directly visible light sources, 
		   <tt>2</tt> will lead to single-bounce (direct-only) illumination, and so on. */
		m_maxDepth = props.getInteger("maxDepth", -1);

		/* Granularity of the work units used in parallelizing 
		   the particle tracing task (default: 200K samples).
		   Should be high enough so that sending and accumulating
		   the partially exposed films is not the bottleneck. */
		m_granularity = props.getSize("granularity", 200000);
	}
	
	AdjointParticleTracer(Stream *stream, InstanceManager *manager) 
		: Integrator(stream, manager) {
		m_maxDepth = stream->readInt();
		m_rrDepth = stream->readInt();
		m_granularity = stream->readSize();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Integrator::serialize(stream, manager);
		stream->writeInt(m_maxDepth);
		stream->writeInt(m_rrDepth);
		stream->writeSize(m_granularity);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
			int sceneResID, int cameraResID, int samplerResID) {
		Integrator::preprocess(scene, queue, job, sceneResID, cameraResID, samplerResID);

		Scheduler *sched = Scheduler::getInstance();
		const Camera *camera = static_cast<Camera *>(sched->getResource(cameraResID));
		Vector2i size = camera->getFilm()->getCropSize();

		if (scene->getSubsurfaceIntegrators().size() > 0)
			Log(EError, "Subsurface integrators are not supported by the particle tracer!");
		m_sampleCount = scene->getSampler()->getSampleCount() * size.x * size.y;
		return true;
	}

	void cancel() {
		Scheduler::getInstance()->cancel(m_process);
	}

	bool render(Scene *scene, RenderQueue *queue, 
		const RenderJob *job, int sceneResID, int cameraResID, int samplerResID) {
		ref<Scheduler> scheduler = Scheduler::getInstance();
		ref<Camera> camera = scene->getCamera();
		const Film *film = camera->getFilm();
		size_t sampleCount = scene->getSampler()->getSampleCount();
		size_t nCores = scheduler->getCoreCount();
		Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " samples, " SIZE_T_FMT 
			" %s, " SSE_STR ") ..", film->getCropSize().x, film->getCropSize().y, 
			sampleCount, nCores, nCores == 1 ? "core" : "cores");

		ref<ParallelProcess> process = new CaptureParticleProcess(
			job, queue, m_sampleCount, m_granularity,
			m_maxDepth, m_rrDepth);

		process->bindResource("scene", sceneResID);
		process->bindResource("camera", cameraResID);
		process->bindResource("sampler", samplerResID);
		scheduler->schedule(process);
		m_process = process;
		scheduler->wait(process);
		m_process = NULL;

		return process->getReturnStatus() == ParallelProcess::ESuccess;
	}

	virtual ~AdjointParticleTracer() {
	}

	MTS_DECLARE_CLASS()
protected:
	ref<ParallelProcess> m_process;
	int m_maxDepth, m_rrDepth;
	size_t m_sampleCount, m_granularity;
};

MTS_IMPLEMENT_CLASS_S(AdjointParticleTracer, false, Integrator)
MTS_EXPORT_PLUGIN(AdjointParticleTracer, "Adjoint particle tracer");
MTS_NAMESPACE_END

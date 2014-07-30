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

#include "ptracer_proc.h"

MTS_NAMESPACE_BEGIN

/*! \plugin{ptracer}{Adjoint particle tracer}
 * \order{12}
 * \parameters{
 *     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *         in the generated output image (where \code{-1} corresponds to $\infty$).
 *	       A value of \code{1} will only render directly visible light sources.
 *	       \code{2} will lead to single-bounce (direct-only) illumination,
 *	       and so on. \default{\code{-1}}
 *	   }
 *	   \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after
 *	      which the implementation will start to use the ``russian roulette''
 *	      path termination criterion. \default{\code{5}}
 *	   }
 *	   \parameter{granularity}{\Integer}{
 *        Specifies the work unit granularity used to parallize the particle
 *        tracing task. This should be set high enough so that accumulating
 *        partially exposed images (and potentially sending them over the network)
 *        is not the bottleneck.
 *        \default{200K particles per work unit, i.e. \code{200000}}
 *     }
 *     \parameter{bruteForce}{\Boolean}{
 *        If set to \code{true}, the integrator does not attempt to create
 *        connections to the sensor and purely relies on hitting it via ray
 *        tracing. This is mainly intended for debugging purposes.
 *        \default{\code{false}}
 *     }
 * }
 *
 * This plugin implements a simple adjoint particle tracer. It does
 * essentially the exact opposite of the simple volumetric path tracer
 * (\pluginref[volpathsimple]{volpath\_simple}): instead of tracing rays from
 * the sensor and attempting to connect them to the light source, this
 * integrator shoots particles from the light source and attempts to connect
 * them to the sensor.
 *
 * Usually, this is a relatively useless rendering technique due to
 * its high variance, but there are some cases where it excels.
 * In particular, it does a good job on scenes where most scattering
 * events are directly visible to the camera.
 *
 * When rendering with a finite-aperture sensor (e.g. \pluginref{thinlens})
 * this integrator is able to intersect the actual aperture, which allows
 * it to handle certain caustic paths that would otherwise not be visible.
 *
 * It also supports a specialized ``brute force'' mode, where the integrator
 * does not attempt to create connections to the sensor and purely relies on
 * hitting it via ray tracing. This is one of the worst conceivable rendering
 * and not recommended for any applications. It is mainly included for
 * debugging purposes.
 *
 * The number of traced particles is given by the number of ``samples per
 * pixel'' of the sample generator times the pixel count of the output image.
 * For instance, 16 samples per pixel on a 512$\times$512 image will cause 4M particles
 * to be generated.
 *
 * \remarks{
 *    \item This integrator does not currently work with subsurface scattering
 *    models.
 * }
 */
class AdjointParticleTracer : public Integrator {
public:
	AdjointParticleTracer(const Properties &props) : Integrator(props) {
		/* Depth to start using russian roulette */
		m_rrDepth = props.getInteger("rrDepth", 5);

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

		/* Rely on hitting the sensor via ray tracing? */
		m_bruteForce = props.getBoolean("bruteForce", false);

		if (m_rrDepth <= 0)
			Log(EError, "'rrDepth' must be set to a value than zero!");

		if (m_maxDepth <= 0 && m_maxDepth != -1)
			Log(EError, "'maxDepth' must be set to -1 (infinite) or a value greater than zero!");
	}

	AdjointParticleTracer(Stream *stream, InstanceManager *manager)
		: Integrator(stream, manager) {
		m_maxDepth = stream->readInt();
		m_rrDepth = stream->readInt();
		m_granularity = stream->readSize();
		m_bruteForce = stream->readBool();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Integrator::serialize(stream, manager);
		stream->writeInt(m_maxDepth);
		stream->writeInt(m_rrDepth);
		stream->writeSize(m_granularity);
		stream->writeBool(m_bruteForce);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
			int sceneResID, int sensorResID, int samplerResID) {
		Integrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);

		Scheduler *sched = Scheduler::getInstance();
		const Sensor *sensor = static_cast<Sensor *>(sched->getResource(sensorResID));
		Vector2i size = sensor->getFilm()->getCropSize();

		if (scene->getSubsurfaceIntegrators().size() > 0)
			Log(EError, "Subsurface integrators are not supported by the particle tracer!");
		m_sampleCount = scene->getSampler()->getSampleCount() * size.x * size.y;
		return true;
	}

	void cancel() {
		Scheduler::getInstance()->cancel(m_process);
	}

	bool render(Scene *scene, RenderQueue *queue,
		const RenderJob *job, int sceneResID, int sensorResID, int samplerResID) {
		ref<Scheduler> scheduler = Scheduler::getInstance();
		ref<Sensor> sensor = scene->getSensor();
		const Film *film = sensor->getFilm();
		size_t sampleCount = scene->getSampler()->getSampleCount();
		size_t nCores = scheduler->getCoreCount();
		Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " samples, " SIZE_T_FMT
			" %s, " SSE_STR ") ..", film->getCropSize().x, film->getCropSize().y,
			sampleCount, nCores, nCores == 1 ? "core" : "cores");

		int maxPtracerDepth = m_maxDepth - 1;

		if ((sensor->getType() & (Emitter::EDeltaDirection
			| Emitter::EDeltaPosition)) == 0 && sensor->isOnSurface()) {
			/* The sensor has a finite aperture and a non-degenerate
			   response function -- trace one more bounce, since we
			   can actually try to hit its aperture */
			maxPtracerDepth++;
		}

		ref<ParallelProcess> process = new CaptureParticleProcess(
			job, queue, m_sampleCount, m_granularity,
			maxPtracerDepth, m_maxDepth, m_rrDepth, m_bruteForce);

		process->bindResource("scene", sceneResID);
		process->bindResource("sensor", sensorResID);
		process->bindResource("sampler", samplerResID);
		scheduler->schedule(process);
		m_process = process;
		scheduler->wait(process);
		m_process = NULL;

		return process->getReturnStatus() == ParallelProcess::ESuccess;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "AdjointParticleTracer[" << endl
			<< "  maxDepth = " << m_maxDepth << "," << endl
			<< "  rrDepth = " << m_rrDepth << "," << endl
			<< "  granularity = " << m_granularity << "," << endl
			<< "  bruteForce = " << m_bruteForce << endl
			<< "]";
		return oss.str();
	}


	MTS_DECLARE_CLASS()
protected:
	ref<ParallelProcess> m_process;
	int m_maxDepth, m_rrDepth;
	size_t m_sampleCount, m_granularity;
	bool m_bruteForce;
};

MTS_IMPLEMENT_CLASS_S(AdjointParticleTracer, false, Integrator)
MTS_EXPORT_PLUGIN(AdjointParticleTracer, "Adjoint particle tracer");
MTS_NAMESPACE_END

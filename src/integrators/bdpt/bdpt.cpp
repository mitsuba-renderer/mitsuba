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

#include <mitsuba/bidir/vertex.h>
#include <mitsuba/bidir/edge.h>
#include "bdpt_proc.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{bdpt}{Bidirectional path tracer}
 * \order{5}
 * \parameters{
 *     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *         in the generated output image (where \code{-1} corresponds to $\infty$).
 *	       A value of \code{1} will only render directly visible light sources.
 *	       \code{2} will lead to single-bounce (direct-only) illumination,
 *	       and so on. \default{\code{-1}}
 *	   }
 *	   \parameter{lightImage}{\Boolean}{Include sampling strategies that connect
 *	      paths traced from emitters directly to the camera? (i.e. what \pluginref{ptracer} does)
 *	      This improves the effectiveness of bidirectional path tracing
 *	      but severely increases the local and remote communication
 *	      overhead, since large \emph{light images} must be transferred between threads
 *	      or over the network. See the text below for a more detailed explanation.
 *	      \default{include these strategies, i.e. \code{true}}
 *     }
 *     \parameter{sampleDirect}{\Boolean}{Enable direct sampling strategies? This is a generalization
 *        of direct illumination sampling that works with both emitters and sensors. Usually a good idea.
 *        \default{use direct sampling, i.e. \code{true}}}
 *	   \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after
 *	      which the implementation will start to use the ``russian roulette''
 *	      path termination criterion. \default{\code{5}}
 *	   }
 * }
 *
 ** \renderings{
 *     \rendering{Path tracer, 32 samples/pixel}{integrator_bdpt_path}
 *     \rendering{Bidirectional path tracer, 32 samples/pixel}{integrator_bdpt_bdpt}
 *     \caption{
 *     	The bidirectional path tracer finds light paths by generating partial
 *     	paths starting at the emitters and the sensor and connecting them in every possible way.
 *     	This works particularly well in closed scenes as the one shown above. Here, the unidirectional
 *     	path tracer has severe difficulties finding some of the indirect illumination paths.
 *     	Modeled after a scene by Eric Veach.
 *     }
 * }
 * \renderings{
 *    \includegraphics[width=12cm]{images/integrator_bdpt_sketch.pdf}\hfill\,
 *    \caption{The four different ways in which BDPT can create a direct illumination
 *    path (matching the first row on the next page): \textbf{(a)} Standard path
 *    tracing without direct illumination sampling, \textbf{(b)} path tracing with
 *    direct illumination sampling, \textbf{(c)} Particle tracing with recording of
 *    scattering events observed by the sensor, \textbf{(d)} Particle tracing with
 *    recording of particles that hit the sensor.}\vspace{-3mm}
 * }
 * \renderings{
 *     \unframedbigrendering{The individual sampling strategies that comprise BDPT, but
 *     \emph{without} multiple importance sampling. $s$ denotes the number of steps
 *     taken from the emitters, and $t$ denotes the number of steps from the sensor.
 *     Note how almost every strategy has deficiencies of some kind} {integrator_bdpt_unweighted.pdf}
 *     \unframedbigrendering{The same sampling strategies, but now weighted using multiple importance
 *     sampling---effectively ``turning off'' each strategy where it does not perform well.
 *     The final result is computed by summing all of these images.}{integrator_bdpt_weighted.pdf}
 * }
 *
 * This plugin implements a bidirectional path tracer (short: BDPT) with support for multiple
 * importance sampling, as proposed by Veach and Guibas \cite{Veach1994Bidirectional}.
 *
 * A bidirectional path tracer computes radiance estimates by starting two separate
 * random walks from an emitter and a sensor. The resulting \emph{subpaths} are
 * connected at every possible interaction vertex, creating a large number of complete paths
 * of different lengths. These paths are then used to estimate the amount of
 * radiance that is transferred from the emitter to a pixel on the sensor.
 *
 * Generally, some of the created paths will be undesirable, since they lead to
 * high-variance radiance estimates. To alleviate this situation, BDPT makes
 * use of \emph{multiple importance sampling} which, roughly speaking, weights
 * paths based on their predicted utility.
 *
 * The bidirectional path tracer in Mitsuba is a complete implementation of the
 * technique that handles all sampling strategies, including those that involve
 * direct interactions with the sensor. For this purpose, finite-aperture sensors
 * are explicitly represented by surfaces in the scene so that they can be
 * intersected by random walks started at emitters.
 *
 * Bidirectional path tracing is a relatively ``heavy'' rendering technique---for
 * the same number of samples per pixel, it is easily 3-4 times slower than regular
 * path tracing. However, it usually makes up for this by producing considerably
 * lower-variance radiance estimates (i.e. the output images have less noise).
 *
 * The code parallelizes over multiple cores and machines, but with one caveat:
 * some of the BDPT path sampling strategies are incompatble with the usual
 * approach of rendering an image tile by tile, since they can potentially
 * contribute to \emph{any} pixel on the screen. This means that each
 * rendering work unit must be associated with a full-sized image!
 * When network render nodes are involved or the resolution of this \emph{light image}
 * is very high, a bottleneck can arise where more work is spent  accumulating or
 * transmitting these images than actual rendering.
 *
 * There are two possible resorts should this situation arise: the first one
 * is to reduce the number of work units so that there is approximately one
 * unit per core (and hence one image to transmit per core). This can be done by
 * increasing the block size in the GUI preferences or passing the \code{-b}
 * parameter to the \code{mitsuba} executable. The second option is to simply
 * disable these sampling strategies at the cost of reducing the
 * effectiveness of bidirectional path tracing (particularly, when rendering
 * caustics). For this, set \code{lightImage} to \code{false}.
 * When rendering an image of a reasonable resolution without network nodes,
 * this is not a big concern, hence these strategies are enabled by default.
 *
 * \remarks{
 *    \item This integrator does not work with dipole-style subsurface
 *    scattering models.
 *    \item This integrator does not yet work with certain non-reciprocal BSDFs (i.e.
 *          \pluginref{bumpmap}, but this will be addressed in the future
 * }
 */
class BDPTIntegrator : public Integrator {
public:
	BDPTIntegrator(const Properties &props) : Integrator(props) {
		/* Load the parameters / defaults */
		m_config.maxDepth = props.getInteger("maxDepth", -1);
		m_config.rrDepth = props.getInteger("rrDepth", 5);
		m_config.lightImage = props.getBoolean("lightImage", true);
		m_config.sampleDirect = props.getBoolean("sampleDirect", true);
		m_config.showWeighted = props.getBoolean("showWeighted", false);

		#if BDPT_DEBUG == 1
		if (m_config.maxDepth == -1 || m_config.maxDepth > 6) {
			/* Limit the maximum depth when rendering image
			   matrices in BDPT debug mode (the number of
			   these images grows quadratically with depth) */
			Log(EWarn, "Limiting max. path depth to 6 to avoid an "
				"extremely large number of debug output images");
			m_config.maxDepth = 6;
		}
		#endif

		if (m_config.rrDepth <= 0)
			Log(EError, "'rrDepth' must be set to a value greater than zero!");

		if (m_config.maxDepth <= 0 && m_config.maxDepth != -1)
			Log(EError, "'maxDepth' must be set to -1 (infinite) or a value greater than zero!");
	}

	/// Unserialize from a binary data stream
	BDPTIntegrator(Stream *stream, InstanceManager *manager)
	 : Integrator(stream, manager) {
		m_config = BDPTConfiguration(stream);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Integrator::serialize(stream, manager);
		m_config.serialize(stream);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue,
			const RenderJob *job, int sceneResID, int sensorResID,
			int samplerResID) {
		Integrator::preprocess(scene, queue, job, sceneResID,
				sensorResID, samplerResID);

		if (scene->getSubsurfaceIntegrators().size() > 0)
			Log(EError, "Subsurface integrators are not supported "
				"by the bidirectional path tracer!");

		return true;
	}

	void cancel() {
		Scheduler::getInstance()->cancel(m_process);
	}

	void configureSampler(const Scene *scene, Sampler *sampler) {
		/* Prepare the sampler for tile-based rendering */
		sampler->setFilmResolution(scene->getFilm()->getCropSize(), true);
	}

	bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
			int sceneResID, int sensorResID, int samplerResID) {
		ref<Scheduler> scheduler = Scheduler::getInstance();
		ref<Sensor> sensor = scene->getSensor();
		const Film *film = sensor->getFilm();
		size_t sampleCount = scene->getSampler()->getSampleCount();
		size_t nCores = scheduler->getCoreCount();

		Log(EDebug, "Size of data structures: PathVertex=%i bytes, PathEdge=%i bytes",
			(int) sizeof(PathVertex), (int) sizeof(PathEdge));

		Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " samples, " SIZE_T_FMT
			" %s, " SSE_STR ") ..", film->getCropSize().x, film->getCropSize().y,
			sampleCount, nCores, nCores == 1 ? "core" : "cores");

		m_config.blockSize = scene->getBlockSize();
		m_config.cropSize = film->getCropSize();
		m_config.sampleCount = sampleCount;
		m_config.dump();

		ref<BDPTProcess> process = new BDPTProcess(job, queue, m_config);
		m_process = process;

		process->bindResource("scene", sceneResID);
		process->bindResource("sensor", sensorResID);
		process->bindResource("sampler", samplerResID);
		scheduler->schedule(process);

		scheduler->wait(process);
		m_process = NULL;
		process->develop();

		#if BDPT_DEBUG == 1
			fs::path path = scene->getDestinationFile();
			if (m_config.lightImage)
				process->getResult()->dump(m_config, path.parent_path(), path.stem());
		#endif

		return process->getReturnStatus() == ParallelProcess::ESuccess;
	}

	MTS_DECLARE_CLASS()
private:
	ref<ParallelProcess> m_process;
	BDPTConfiguration m_config;
};

MTS_IMPLEMENT_CLASS_S(BDPTIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(BDPTIntegrator, "Bidirectional path tracer");
MTS_NAMESPACE_END

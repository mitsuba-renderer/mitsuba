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

#include <mitsuba/bidir/util.h>
#include <mitsuba/bidir/pathsampler.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/plugin.h>
#include "erpt_proc.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{erpt}{Energy redistribution path tracing}
 * \order{11}
 * \parameters{
 *     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *         in the generated output image (where \code{-1} corresponds to $\infty$).
 *         A value of \code{1} will only render directly visible light sources.
 *         \code{2} will lead to single-bounce (direct-only) illumination,
 *         and so on. \default{\code{-1}}
 *     }
 *     \parameter{numChains}{\Float}{
 *         On average, how many Markov Chains should be started per
 *         pixel? \default{\code{1}}
 *     }
 *     \parameter{maxChains}{\Float}{
 *         How many Markov Chains should be started \emph{at most} (per
 *         pixel) \default{\code{0}, i.e. this feature is not used}
 *     }
 *     \parameter{chainLength}{\Integer}{
 *         Specifies the number of perturbation steps that are executed per Markov Chain\default{1}.
 *     }
 *     \parameter{directSamples}{\Integer}{
 *         By default, the implementation renders direct illumination component
 *         separately using the \pluginref{direct} plugin, which
 *         uses low-discrepancy number sequences for superior performance
 *         (in other words, it is \emph{not} handled by ERPT). This
 *         parameter specifies the number of samples allocated to that method. To
 *         force MLT to be responsible for the direct illumination
 *         component as well, set this to \code{-1}. \default{\code{16}}
 *     }
 *     \parameter{[lens,multiChain,\!\!\!\!\newline caustic,manifold]\showbreak
 *       \newline Perturbation}{\Boolean}{
 *       These parameters can be used to pick the individual perturbation
 *       strategies that will be used to explore path space. By default, the original set
 *       by Veach and Guibas is enabled (i.e. everything except the manifold
 *       perturbation).
 *     }
 *     \parameter{lambda}{\Float}{
 *         Jump size of the manifold perturbation \default{\code{50}}}
 * }
 * \renderings{
 *  \rendering{A brass chandelier with 24 glass-enclosed bulbs}{integrator_mept_luminaire}
 *  \rendering{Glossy reflective and refractive ableware, lit by the
 *  chandelier on the left}{integrator_mept_tableware}
 *  \vspace{-2mm}
 *  \caption{An interior scene with complex specular and near-specular light paths,
 *  illuminated entirely through caustics. Rendered by this plugin
 *  using the manifold perturbation. This scene was designed by Olesya
 *  Isaenko.\vspace{2mm}}
 * }
 * \renderings{
 * \medrendering{Seed paths generated using bidirectional path tracing.
 *  Note the high variance of paths that involve reflection of sunlight by
 *  the torus.}{integrator_erpt_seeds}
 * \medrendering{Result after running the perturbations of Veach and Guibas
 *  for 800 steps. Some convergence issues remain.}{integrator_erpt_mlt}
 * \medrendering{Result after running the manifold
 * perturbation for the same amount of time}{integrator_erpt_mept}
 * }
 *
 * Energy Redistribution Path Tracing (ERPT) by Cline et al. \cite{Cline2005Energy}
 * combines Path Tracing with the perturbation strategies of Metropolis Light Transport.
 *
 *
 * An initial set of \emph{seed paths} is generated using a standard bidirectional
 * path tracer, and for each one, a MLT-style Markov Chain is subsequently started
 * and executed for some number of steps.
 * This has the effect of redistributing the energy of the individual samples
 * over a larger area, hence the name of this method.
 *
 * \begin{wrapfigure}{R}{0.30\textwidth}
 *  \fbox{\includegraphics[width=.29\textwidth]{images/integrator_mept_egg.jpg}}
 *  \caption{Another view, now with exterior lighting.}
 *  \vspace{-2mm}
 * \end{wrapfigure}
 * This is often a good choice when a (bidirectional) path tracer produces mostly reasonable
 * results except that it finds certain important types of light paths too rarely.
 * ERPT can
 * then explore all of the neighborhing paths as well, to prevent the
 * original sample from showing up as a ``bright pixel'' in the output image.
 *
 * This plugin shares all the perturbation strategies of the \pluginref{mlt} plugin, and
 * the same rules for selecting them apply. In contrast to the original
 * paper by Cline et al., the Mitsuba implementation uses a bidirectional
 * (rather than an unidirectional) bidirectional path tracer to create seed paths.
 * Also, since they add bias to the output, this plugin does not use the image
 * post-processing filters proposed by the authors.
 *
 * The mechanism for selecting Markov Chain seed paths deserves an
 * explanation: when commencing work on a pixel in the output image, the
 * integrator first creates a pool of seed path candidates. The size of this
 * pool is given by the \code{samplesPerPixel} parameter of the sample
 * generator. This should be large enough so that the integrator has a
 * representative set of light paths to work with.
 *
 * Subsequently, one or more of these candidates are chosen (determined by
 * \code{numChains} and \code{maxChains} parameter). For each one, a Markov
 * Chain is created that has an initial configuration matching the seed path.
 * It is simulated for \code{chainLength} iterations, and each intermediate
 * state is recorded in the output image.
 */
class EnergyRedistributionPathTracing : public Integrator {
public:
    EnergyRedistributionPathTracing(const Properties &props) : Integrator(props) {
        m_config.maxDepth = props.getInteger("maxDepth", -1);

        /* Specifies the number of Markov Chains that, on average, are
           started per pixel */
        m_config.numChains = props.getFloat("numChains", 1.0f);
        m_config.maxChains = props.getInteger("maxChains", 0);

        /* Specifies the number of mutations to be performed in each
           Markov Chain */
        m_config.chainLength = props.getInteger("chainLength", 100);

        m_config.directSamples = props.getInteger("directSamples", 16);
        m_config.separateDirect = m_config.directSamples >= 0;

        /* Number of samples used to estimate the average contribution of a
           single sample. In contrast to MLT/PSSMLT, this parameter just
           influences the heuristics that guide how many Markov Chains to start
           for a sample. Usually, there is little reason to change it */
        m_config.luminanceSamples = props.getInteger("luminanceSamples", 15000);

        /* Selectively enable/disable the bidirectional mutation. This is
           probably not desireable for ERPT, since the path tracing stage is
           responsible for generating good seed paths.. */
        m_config.bidirectionalMutation = props.getBoolean("bidirectionalMutation", false);

        /* Selectively enable/disable the lens perturbation */
        m_config.lensPerturbation = props.getBoolean("lensPerturbation", true);

        /* Selectively enable/disable the caustic perturbation */
        m_config.causticPerturbation = props.getBoolean("causticPerturbation", true);

        /* Selectively enable/disable the multi-chain perturbation */
        m_config.multiChainPerturbation = props.getBoolean("multiChainPerturbation", true);

        /* Selectively enable/disable the manifold perturbation */
        m_config.manifoldPerturbation = props.getBoolean("manifoldPerturbation", false);
        m_config.probFactor = props.getFloat("probFactor", props.getFloat("lambda", 50));
        m_config.avgAngleChangeSurface = props.getFloat("avgAngleChangeSurface", 0);
        m_config.avgAngleChangeMedium = props.getFloat("avgAngleChangeMedium", 0);

        if (m_config.maxDepth <= 0 && m_config.maxDepth != -1)
            Log(EError, "'maxDepth' must be set to -1 (infinite) or a value greater than zero!");
    }

    /// Unserialize from a binary data stream
    EnergyRedistributionPathTracing(Stream *stream, InstanceManager *manager)
     : Integrator(stream, manager) {
        m_config = ERPTConfiguration(stream);
    }

    virtual ~EnergyRedistributionPathTracing() { }

    void configureSampler(const Scene *scene, Sampler *sampler) {
        /* Prepare the sampler for tile-based rendering */
        sampler->setFilmResolution(scene->getFilm()->getCropSize(), true);
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
            Log(EError, "Subsurface integrators are not supported by ERPT!");

        return true;
    }

    void cancel() {
        ref<RenderJob> nested = m_nestedJob;
        if (nested)
            nested->cancel();
        Scheduler::getInstance()->cancel(m_process);
    }

    bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID) {
        ref<Scheduler> sched = Scheduler::getInstance();
        ref<Sensor> sensor = scene->getSensor();
        ref<Sampler> sampler = sensor->getSampler();
        const Film *film = sensor->getFilm();
        size_t nCores = sched->getCoreCount();
        size_t sampleCount = sampler->getSampleCount();

        Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT
            " %s, " SSE_STR ", " SIZE_T_FMT " samples/pixel) ..",
            film->getCropSize().x, film->getCropSize().y, nCores,
            nCores == 1 ? "core" : "cores", sampleCount);

        ref<Sampler> indepSampler = static_cast<Sampler *> (PluginManager::getInstance()->
            createObject(MTS_CLASS(Sampler), Properties("independent")));
        indepSampler->configure();

        ref<PathSampler> pathSampler = new PathSampler(PathSampler::EBidirectional, scene,
            indepSampler, indepSampler, indepSampler, m_config.maxDepth, 10,
            m_config.separateDirect, true, true);

        m_config.luminance = pathSampler->computeAverageLuminance(
                m_config.luminanceSamples);
        m_config.blockSize = scene->getBlockSize();

        m_config.dump();

        ref<Bitmap> directImage;
        if (m_config.separateDirect && m_config.directSamples > 0) {
            directImage = BidirectionalUtils::renderDirectComponent(scene,
                sceneResID, sensorResID, queue, job, m_config.directSamples);
            if (directImage == NULL)
                return false;
        }

        ref<ERPTProcess> process = new ERPTProcess(job, queue, m_config, directImage);

        /* Create an independent sampler for use by the MLT chains */
        std::vector<SerializableObject *> indepSamplers(sched->getCoreCount());
        for (size_t i=0; i<sched->getCoreCount(); ++i) {
            ref<Sampler> clonedSampler = indepSampler->clone();
            clonedSampler->incRef();
            indepSamplers[i] = clonedSampler.get();
        }
        int indepSamplerResID = sched->registerMultiResource(indepSamplers);
        for (size_t i=0; i<indepSamplers.size(); ++i)
            indepSamplers[i]->decRef();

        process->bindResource("scene", sceneResID);
        process->bindResource("sensor", sensorResID);
        process->bindResource("sampler", samplerResID);
        process->bindResource("indepSampler", indepSamplerResID);

        m_process = process;
        sched->schedule(process);
        sched->wait(process);
        m_process = NULL;

        sched->unregisterResource(indepSamplerResID);

        return process->getReturnStatus() == ParallelProcess::ESuccess;
    }

    MTS_DECLARE_CLASS()
private:
    ref<ParallelProcess> m_process;
    ref<RenderJob> m_nestedJob;
    ERPTConfiguration m_config;
};

MTS_IMPLEMENT_CLASS_S(EnergyRedistributionPathTracing, false, Integrator)
MTS_EXPORT_PLUGIN(EnergyRedistributionPathTracing, "Energy redistribution path tracing");
MTS_NAMESPACE_END

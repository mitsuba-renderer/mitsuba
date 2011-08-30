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
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/renderproc.h>

MTS_NAMESPACE_BEGIN

Integrator::Integrator(const Properties &props)
 : NetworkedObject(props), m_properties(props) {
}

Integrator::Integrator(Stream *stream, InstanceManager *manager)
 : NetworkedObject(stream, manager) {
}

bool Integrator::preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID) { return true; }
void Integrator::postprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID) { }
void Integrator::configureSampler(Sampler *sampler) { }
void Integrator::serialize(Stream *stream, InstanceManager *manager) const {
	NetworkedObject::serialize(stream, manager);
}
const Integrator *Integrator::getSubIntegrator() const { return NULL; }

SampleIntegrator::SampleIntegrator(const Properties &props)
 : Integrator(props) { }

SampleIntegrator::SampleIntegrator(Stream *stream, InstanceManager *manager)
 : Integrator(stream, manager) { }

void SampleIntegrator::serialize(Stream *stream, InstanceManager *manager) const {
	Integrator::serialize(stream, manager);
}

Spectrum SampleIntegrator::E(const Scene *scene, const Point &p, const Normal &n, Float time,
		const Medium *medium, Sampler *sampler, int nSamples, bool handleIndirect) const {
	Spectrum E(0.0f);
	LuminaireSamplingRecord lRec;
	RadianceQueryRecord rRec(scene, sampler);
	Frame frame(n);

	sampler->generate();
	for (int i=0; i<nSamples; i++) {
		rRec.newQuery(RadianceQueryRecord::ERadianceNoEmission, medium);

		/* Direct */
		if (scene->sampleAttenuatedLuminaire(p, time, medium, lRec, rRec.nextSample2D())) {
			Float dp = dot(lRec.d, n);
			if (dp < 0) 
				E -= lRec.value * dp;
		}

		/* Indirect */
		if (handleIndirect) {
			Vector d = frame.toWorld(squareToHemispherePSA(rRec.nextSample2D()));
			++rRec.depth;
			E += Li(RayDifferential(p, d, time), rRec) * M_PI;
		}
		sampler->advance();
	}
	return E / (Float) nSamples;
}

void SampleIntegrator::cancel() {
	if (m_process)
		Scheduler::getInstance()->cancel(m_process);
}

bool SampleIntegrator::render(Scene *scene,
		RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID) {
	ref<Scheduler> sched = Scheduler::getInstance();
	ref<Camera> camera = static_cast<Camera *>(sched->getResource(cameraResID));
	ref<Film> film = camera->getFilm();

	size_t nCores = sched->getCoreCount();
	const Sampler *sampler = static_cast<const Sampler *>(sched->getResource(samplerResID, 0));
	size_t sampleCount = sampler->getSampleCount();

	Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SIZE_T_FMT 
		" %s, " SSE_STR ") ..", film->getCropSize().x, film->getCropSize().y, 
		sampleCount, sampleCount == 1 ? "sample" : "samples", nCores, 
		nCores == 1 ? "core" : "cores");

	/* This is a sampling-based integrator - parallelize */
	ref<ParallelProcess> proc = new BlockedRenderProcess(job, 
		queue, scene->getBlockSize());
	int integratorResID = sched->registerResource(this);
	proc->bindResource("integrator", integratorResID);
	proc->bindResource("scene", sceneResID);
	proc->bindResource("camera", cameraResID);
	proc->bindResource("sampler", samplerResID);
	scene->bindUsedResources(proc);
	bindUsedResources(proc);
	sched->schedule(proc);

	m_process = proc;
	sched->wait(proc);
	m_process = NULL;
	sched->unregisterResource(integratorResID);

	return proc->getReturnStatus() == ParallelProcess::ESuccess;
}

void SampleIntegrator::bindUsedResources(ParallelProcess *) const {
	/* Do nothing by default */
}

void SampleIntegrator::wakeup(std::map<std::string, SerializableObject *> &) {
	/* Do nothing by default */
}

void SampleIntegrator::renderBlock(const Scene *scene,
	const Camera *camera, Sampler *sampler, ImageBlock *block, 
	const bool &stop, const std::vector<Point2i> *points) const {
	Point2 sample, lensSample;
	RayDifferential eyeRay;
	Float timeSample = 0;
	Spectrum spec;


	block->clear();
	RadianceQueryRecord rRec(scene, sampler);
	bool needsLensSample = camera->needsLensSample();
	bool needsTimeSample = camera->needsTimeSample();
	const TabulatedFilter *filter = camera->getFilm()->getTabulatedFilter();
	Float scaleFactor = 1.0f/std::sqrt((Float) sampler->getSampleCount());

	if (points) {
		/* Use a prescribed traversal order (e.g. using a space-filling curve) */
		if (!block->collectStatistics()) {
			for (size_t i=0; i<points->size(); ++i) {
				Point2i offset = (*points)[i] + Vector2i(block->getOffset());
				if (stop) 
					break;
				sampler->generate();
				for (size_t j = 0; j<sampler->getSampleCount(); j++) {
					rRec.newQuery(RadianceQueryRecord::ECameraRay, camera->getMedium());
					if (needsLensSample)
						lensSample = rRec.nextSample2D();
					if (needsTimeSample)
						timeSample = rRec.nextSample1D();
					sample = rRec.nextSample2D();
					sample.x += offset.x; sample.y += offset.y;
					camera->generateRayDifferential(sample, 
						lensSample, timeSample, eyeRay);
					eyeRay.scaleDifferential(scaleFactor);
					spec = Li(eyeRay, rRec);
					block->putSample(sample, spec, rRec.alpha, filter);
					sampler->advance();
				}
			}
		} else {
			Spectrum mean, meanSqr;
			for (size_t i=0; i<points->size(); ++i) {
				Point2i offset = (*points)[i] + Vector2i(block->getOffset());
				if (stop) 
					break;
				sampler->generate();
				mean = meanSqr = Spectrum(0.0f);
				for (size_t j = 0; j<sampler->getSampleCount(); j++) {
					rRec.newQuery(RadianceQueryRecord::ECameraRay, camera->getMedium());
					if (needsLensSample)
						lensSample = rRec.nextSample2D();
					if (needsTimeSample)
						timeSample = rRec.nextSample1D();
					sample = rRec.nextSample2D();
					sample.x += offset.x; sample.y += offset.y;
					camera->generateRayDifferential(sample, 
						lensSample, timeSample, eyeRay);
					eyeRay.scaleDifferential(scaleFactor);
					spec = Li(eyeRay, rRec);

					/* Numerically robust online variance estimation using an
						algorithm proposed by Donald Knuth (TAOCP vol.2, 3rd ed., p.232) */
					const Spectrum delta = spec - mean;
					mean += delta / ((Float) j+1);
					meanSqr += delta * (spec - mean);
					block->putSample(sample, spec, rRec.alpha, filter);
					block->setVariance(offset.x, offset.y,
						meanSqr / (Float) j, (int) j+1);
					sampler->advance();
				}
			}
		}
	} else {
		/* Simple scanline traversal order */
		const int
			sx = block->getOffset().x,
			sy = block->getOffset().y,
			ex = sx + block->getSize().x,
			ey = sy + block->getSize().y;
		if (!block->collectStatistics()) {
			for (int y = sy; y < ey; y++) {
				for (int x = sx; x < ex; x++) {
					if (stop) 
						break;
					sampler->generate();
					for (size_t j = 0; j<sampler->getSampleCount(); j++) {
						rRec.newQuery(RadianceQueryRecord::ECameraRay, camera->getMedium());
						if (needsLensSample)
							lensSample = rRec.nextSample2D();
						if (needsTimeSample)
							timeSample = rRec.nextSample1D();
						sample = rRec.nextSample2D();
						sample.x += x; sample.y += y;
						camera->generateRayDifferential(sample, 
							lensSample, timeSample, eyeRay);
						eyeRay.scaleDifferential(scaleFactor);
						spec = Li(eyeRay, rRec);
						block->putSample(sample, spec, rRec.alpha, filter);
						sampler->advance();
					}
				}
			}
		} else {
			Spectrum mean, meanSqr;
			for (int y = sy; y < ey; y++) {
				for (int x = sx; x < ex; x++) {
					if (stop) 
						break;
					sampler->generate();
					mean = meanSqr = Spectrum(0.0f);
					for (size_t j = 0; j<sampler->getSampleCount(); j++) {
						rRec.newQuery(RadianceQueryRecord::ECameraRay, camera->getMedium());
						if (needsLensSample)
							lensSample = rRec.nextSample2D();
						if (needsTimeSample)
							timeSample = rRec.nextSample1D();
						sample = rRec.nextSample2D();
						sample.x += x; sample.y += y;
						camera->generateRayDifferential(sample, 
							lensSample, timeSample, eyeRay);
						eyeRay.scaleDifferential(scaleFactor);
						spec = Li(eyeRay, rRec);

						/* Numerically robust online variance estimation using an
						   algorithm proposed by Donald Knuth (TAOCP vol.2, 3rd ed., p.232) */
						const Spectrum delta = spec - mean;
						mean += delta / ((Float) j+1);
						meanSqr += delta * (spec - mean);
						block->putSample(sample, spec, rRec.alpha, filter);
						block->setVariance(x, y, meanSqr / (Float) j, (int) j+1);
						sampler->advance();
					}
				}
			}
		}
	}
}

MonteCarloIntegrator::MonteCarloIntegrator(const Properties &props) : SampleIntegrator(props) {
	/* Depth to begin using russian roulette */
	m_rrDepth = props.getInteger("rrDepth", 10);

	/* Longest visualized path length (\c -1 = infinite). 
	   A value of \c 1 will visualize only directly visible light sources.
	   \c 2 will lead to single-bounce (direct-only) illumination, and so on. */
	m_maxDepth = props.getInteger("maxDepth", -1);

	/**
	 * This parameter specifies the action to be taken when the geometric
	 * and shading normals of a surface inconsistently classify a ray as
	 * being located *both* on the front and back-side.
	 *
	 * When \c strictNormals is set to \c false, the shading normal has 
	 * precedence, and rendering proceeds normally at the risk of
	 * introducing small light leaks (this is the default).
	 *
	 * When \c strictNormals is set to \c true, the random walk is
	 * terminated when encountering such a situation. This may
	 * produce black silhouette edges on badly tesselated meshes.
	 */
	m_strictNormals = props.getBoolean("strictNormals", false);

	if (m_rrDepth <= 0)
		Log(EError, "rrDepth must be greater than zero!");
}

MonteCarloIntegrator::MonteCarloIntegrator(Stream *stream, InstanceManager *manager)
	: SampleIntegrator(stream, manager) {
	m_rrDepth = stream->readInt();
	m_maxDepth = stream->readInt();
	m_strictNormals = stream->readBool();
}

void MonteCarloIntegrator::serialize(Stream *stream, InstanceManager *manager) const {
	SampleIntegrator::serialize(stream, manager);
	stream->writeInt(m_rrDepth);
	stream->writeInt(m_maxDepth);
	stream->writeBool(m_strictNormals);
}

std::string RadianceQueryRecord::toString() const {
	std::ostringstream oss;
	oss << "RadianceQueryRecord[" << std::endl
		<< "  type = { ";
	if (type & EEmittedRadiance) oss << "emitted ";
	if (type & ESubsurfaceRadiance) oss << "subsurface ";
	if (type & EDirectSurfaceRadiance) oss << "direct ";
	if (type & EIndirectSurfaceRadiance) oss << "indirect ";
	if (type & ECausticRadiance) oss << "caustic ";
	if (type & EDirectMediumRadiance) oss << "inscatteredDirect ";
	if (type & EIndirectMediumRadiance) oss << "inscatteredIndirect ";
	if (type & EDistance) oss << "distance ";
	if (type & EOpacity) oss << "opacity ";
	if (type & EIntersection) oss << "intersection ";
	oss << "}," << std::endl
		<< "  depth = " << depth << "," << std::endl
		<< "  its = " << indent(its.toString()) << std::endl
		<< "  alpha = " << alpha << "," << std::endl
		<< "  extra = " << extra << "," << std::endl
		<< "]" << std::endl;
	return oss.str();
}


MTS_IMPLEMENT_CLASS(Integrator, true, NetworkedObject)
MTS_IMPLEMENT_CLASS(SampleIntegrator, true, Integrator)
MTS_IMPLEMENT_CLASS(MonteCarloIntegrator, true, SampleIntegrator)
MTS_NAMESPACE_END

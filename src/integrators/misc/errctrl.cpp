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

#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

static StatsCounter pixelsRendered("Error-controlling integrator", "Pixels rendered");

/**
 * Adaptive integrator - runs a secondary integrator until the 
 * the computed radiance achieves a specifiable relative error 
 * threshold (5% by default) with a certain probability (95% by default). 
 * Internally, it uses a Z-test to decide when to stop collecting samples.
 * While probably repeating a Z-test in this way is likely not entirely 
 * rigorous in the statistical sense, it provides a useful statistically
 * motivated stopping criterion.
 *
 * When used in conjunction with image reconstruction filters, this class
 * ensures that neighboring image regions are not unduly biased by placing
 * many samples at a certain position.
 */
class ErrorControl : public SampleIntegrator {
public:
	ErrorControl(const Properties &props) : SampleIntegrator(props) {
		/* Maximum relative error threshold */
		m_maxError = props.getFloat("maxError", 0.05f);
		/* Minimum numbers of samples (Should be large enough so 
		   that reliable variance estimates can be obtained) */
		m_minSamples = props.getInteger("minSamples", 64);
		/* Absolute maximum number of samples to take. The sample collection 
		   will stop after this many samples even if the variance is still too high.
		   A negative value will be interpreted as infinity. */
		m_maxSamples = props.getInteger("maxSamples", 2048);
		/* Required P-value to accept a sample. */
		m_pval = props.getFloat("pval", 0.05f);
		/* Specifies whether the relative error criterion is applied
		   pixel-by-pixel, or whether the comparison should be made
		   against the overall luminance on the film plane. If this
		   is set to true, the quality in dark regions will be improved
		   at the cost of possibly spending lots of time on them. */
		m_perPixel = props.getBoolean("perPixel", false);
		m_verbose = props.getBoolean("verbose", false);
	}

	ErrorControl(Stream *stream, InstanceManager *manager) 
	 : SampleIntegrator(stream, manager) {
		m_subIntegrator = static_cast<SampleIntegrator *>(manager->getInstance(stream));
		m_minSamples = stream->readInt();
		m_maxSamples = stream->readInt();
		m_maxError = stream->readFloat();
		m_quantile = stream->readFloat();
		m_averageLuminance = stream->readFloat();
		m_pval = stream->readFloat();
		m_perPixel = stream->readBool();
		m_verbose = false;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();

		if (cClass->derivesFrom(MTS_CLASS(Integrator))) {
			if (!cClass->derivesFrom(MTS_CLASS(SampleIntegrator)))
				Log(EError, "The sub-integrator must be derived from the class SampleIntegrator");
			m_subIntegrator = static_cast<SampleIntegrator *>(child);
		} else {
			Integrator::addChild(name, child);
		}
	}

	void configureSampler(Sampler *sampler) {
		m_subIntegrator->configureSampler(sampler);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job, 
			int sceneResID, int cameraResID, int samplerResID) {
		if (!SampleIntegrator::preprocess(scene, queue, job, sceneResID, cameraResID, samplerResID))
			return false;
		if (m_subIntegrator == NULL)
			Log(EError, "No sub-integrator was specified!");
		Sampler *sampler = static_cast<Sampler *>(Scheduler::getInstance()->getResource(samplerResID, 0));
		Camera *camera = static_cast<Camera *>(Scheduler::getInstance()->getResource(cameraResID));
		if (sampler->getClass()->getName() != "IndependentSampler")
			Log(EError, "The error-controlling integrator should only be "
				"used in conjunction with the independent sampler");
		if (!m_subIntegrator->preprocess(scene, queue, job, sceneResID, cameraResID, samplerResID))
			return false;

		Vector2i filmSize = camera->getFilm()->getSize();
		bool needsLensSample = camera->needsLensSample();
		bool needsTimeSample = camera->needsTimeSample();
		const int nSamples = 10000;
		Float luminance = 0, timeSample = 0;
		RadianceQueryRecord rRec(scene, sampler);

		/* Estimate the overall luminance on the image plane */
		for (int i=0; i<nSamples; ++i) {
			Point2 sample, lensSample;
			RayDifferential eyeRay;
			sampler->generate();
			rRec.newQuery(RadianceQueryRecord::ERadiance, camera->getMedium());
			rRec.extra = RadianceQueryRecord::EAdaptiveQuery;
			if (needsLensSample)
				lensSample = rRec.nextSample2D();
			if (needsTimeSample)
				timeSample = rRec.nextSample1D();
			sample = rRec.nextSample2D();
			sample.x *= filmSize.x;
			sample.y *= filmSize.y;
			camera->generateRayDifferential(sample, lensSample, timeSample, eyeRay);

			luminance += m_subIntegrator->Li(eyeRay, rRec).getLuminance();
		}

		m_averageLuminance = luminance / nSamples;
		m_quantile = (Float) normalQuantile(1-m_pval/2);
		Log(EInfo, "Configuring for a %.1f%% confidence interval, quantile=%f, avg. luminance=%f", 
			(1-m_pval)*100, m_quantile, m_averageLuminance);
		return true;
	}

	void renderBlock(const Scene *scene, const Camera *camera, Sampler *sampler, 
			ImageBlock *block, const bool &stop, const std::vector<Point2i> *points) const {
		bool needsLensSample = camera->needsLensSample();
		bool needsTimeSample = camera->needsTimeSample();
		const TabulatedFilter *filter = camera->getFilm()->getTabulatedFilter();

		Float mean, meanSqr;
		Point2 sample, lensSample;
		RayDifferential eyeRay;
		int x, y;
		Float sampleLuminance, timeSample = 0;
		RadianceQueryRecord rRec(scene, sampler);
		int sampleIndex;

		const int sx = block->getOffset().x,
				sy = block->getOffset().y,
				ex = sx + block->getSize().x,
				ey = sy + block->getSize().y;

		block->clear();
		if (points) {
			/* Use a prescribed traversal order (e.g. using a space-filling curve) */
			for (size_t i=0; i<points->size(); ++i) {
				Point2i offset = points->operator[](i) 
					+ Vector2i(block->getOffset());
				sampler->generate();
				mean = meanSqr = 0.0f;
				sampleIndex = 0;

				block->snapshot(offset.x, offset.y);
				while (!stop) {
					rRec.newQuery(RadianceQueryRecord::ECameraRay, camera->getMedium());
					rRec.extra = RadianceQueryRecord::EAdaptiveQuery;

					if (needsLensSample)
						lensSample = rRec.nextSample2D();
					if (needsTimeSample)
						timeSample = rRec.nextSample1D();
					sample = rRec.nextSample2D();
					sample.x += offset.x; sample.y += offset.y;
					camera->generateRayDifferential(sample, 
						lensSample, timeSample, eyeRay);

					Spectrum sampleValue = m_subIntegrator->Li(eyeRay, rRec);

					if (block->putSample(sample, sampleValue, rRec.alpha, filter)) {
						/* Check for problems with the sample */
						sampleLuminance = sampleValue.getLuminance();
					} else {
						sampleLuminance = 0.0f;
					}
					++sampleIndex;
					sampler->advance();

					/* Numerically robust online variance estimation using an
					algorithm proposed by Donald Knuth (TAOCP vol.2, 3rd ed., p.232) */
					const Float delta = sampleLuminance - mean;
					mean += delta / sampleIndex;
					meanSqr += delta * (sampleLuminance - mean);

					if (m_maxSamples >= 0 && sampleIndex >= m_maxSamples) {
						break;
					} else if (sampleIndex >= m_minSamples) {
						/* Variance of the primary estimator */
						const Float variance = meanSqr / (sampleIndex-1);

						Float stdError = std::sqrt(variance/sampleIndex);

						/* Half width of the confidence interval */
						Float ciWidth = stdError * m_quantile;

						if (m_verbose && (sampleIndex % 100) == 0) 
							Log(EDebug, "%i samples, mean=%f, stddev=%f, std error=%f, ci width=%f, max allowed=%f", sampleIndex, mean, 
								std::sqrt(variance), stdError, ciWidth, (m_perPixel ? mean : m_averageLuminance) * m_maxError);

						if (m_perPixel && ciWidth <= m_maxError*mean)
							break;
						else if (!m_perPixel && ciWidth <= m_maxError * m_averageLuminance)
							break;
					}
				}

				/* Ensure that a large amounts of samples in one 
				   pixel do not excessively bias neighboring pixels */
				block->normalize(offset.x, offset.y, 1.0f / sampleIndex);

				if (block->collectStatistics())
					block->setVariance(offset.x, offset.y, 
							Spectrum(meanSqr / (sampleIndex-1)), sampleIndex);
				++pixelsRendered;
			}
		} else {
			/* Use a basic grid traversal */
			for (y = sy; y < ey; y++) {
				for (x = sx; x < ex; x++) {
					sampler->generate();
					mean = meanSqr = 0.0f;
					sampleIndex = 0;

					block->snapshot(x, y);
					while (!stop) {
						rRec.newQuery(RadianceQueryRecord::ECameraRay, camera->getMedium());
						rRec.extra = RadianceQueryRecord::EAdaptiveQuery;

						if (needsLensSample)
							lensSample = rRec.nextSample2D();
						if (needsTimeSample)
							timeSample = rRec.nextSample1D();
						sample = rRec.nextSample2D();
						sample.x += x; sample.y += y;
						camera->generateRayDifferential(sample, 
							lensSample, timeSample, eyeRay);

						Spectrum sampleValue = m_subIntegrator->Li(eyeRay, rRec);

						if (block->putSample(sample, sampleValue, rRec.alpha, filter)) {
							/* Check for problems with the sample */
							sampleLuminance = sampleValue.getLuminance();
						} else {
							sampleLuminance = 0.0f;
						}
						++sampleIndex;
						sampler->advance();

						/* Numerically robust online variance estimation using an
						algorithm proposed by Donald Knuth (TAOCP vol.2, 3rd ed., p.232) */
						const Float delta = sampleLuminance - mean;
						mean += delta / sampleIndex;
						meanSqr += delta * (sampleLuminance - mean);

						if (m_maxSamples >= 0 && sampleIndex >= m_maxSamples) {
							break;
						} else if (sampleIndex >= m_minSamples) {
							/* Variance of the primary estimator */
							const Float variance = meanSqr / (sampleIndex-1);

							Float stdError = std::sqrt(variance/sampleIndex);

							/* Half width of the confidence interval */
							Float ciWidth = stdError * m_quantile;

							if (m_verbose && (sampleIndex % 100) == 0) 
								Log(EDebug, "%i samples, mean=%f, stddev=%f, std error=%f, ci width=%f, max allowed=%f", sampleIndex, mean, 
									std::sqrt(variance), stdError, ciWidth, (m_perPixel ? mean : m_averageLuminance) * m_maxError);

							if (m_perPixel && ciWidth <= m_maxError*mean)
								break;
							else if (!m_perPixel && ciWidth <= m_maxError * m_averageLuminance)
								break;
						}
					}

					/* Ensure that a large amounts of samples in one 
					   pixel do not excessively bias neighboring pixels */
					block->normalize(x, y, 1.0f / sampleIndex);

					if (block->collectStatistics())
						block->setVariance(x, y, Spectrum(meanSqr / (sampleIndex-1)), sampleIndex);
					++pixelsRendered;
				}
			}
		}
	}

	Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
		return m_subIntegrator->Li(ray, rRec);
	}

	Spectrum E(const Scene *scene, const Point &p, const
			Normal &n, Float time, const Medium *medium, Sampler *sampler,
			int nSamples, bool includeIndirect) const { 
		return m_subIntegrator->E(scene, p, n, time, medium,
			sampler, nSamples, includeIndirect);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		SampleIntegrator::serialize(stream, manager);
		manager->serialize(stream, m_subIntegrator.get());

		stream->writeInt(m_minSamples);
		stream->writeInt(m_maxSamples);
		stream->writeFloat(m_maxError);
		stream->writeFloat(m_quantile);
		stream->writeFloat(m_averageLuminance);
		stream->writeFloat(m_pval);
		stream->writeBool(m_perPixel);
	}
	
	void bindUsedResources(ParallelProcess *proc) const {
		m_subIntegrator->bindUsedResources(proc);
	}

	void wakeup(std::map<std::string, SerializableObject *> &params) {
		m_subIntegrator->wakeup(params);
	}

	void cancel() {
		SampleIntegrator::cancel();
		m_subIntegrator->cancel();
	}


	const Integrator *getSubIntegrator() const {
		return m_subIntegrator.get();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "ErrorControl[" << std::endl
			<< "  minSamples = " << m_minSamples << "," << std::endl
			<< "  maxSamples = " << m_maxSamples << "," << std::endl
			<< "  maxError   = " << m_maxError << "," << std::endl
			<< "  quantile   = " << m_quantile << "," << std::endl
			<< "  pValue     = " << m_pval << "," << std::endl
			<< "  subIntegrator = " << indent(m_subIntegrator->toString()) << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<SampleIntegrator> m_subIntegrator;
	int m_minSamples;
	int m_maxSamples;
	Float m_maxError, m_quantile, m_pval, m_averageLuminance;
	bool m_verbose, m_perPixel;
};

MTS_IMPLEMENT_CLASS_S(ErrorControl, false, SampleIntegrator)
MTS_EXPORT_PLUGIN(ErrorControl, "Error-controlling integrator");
MTS_NAMESPACE_END

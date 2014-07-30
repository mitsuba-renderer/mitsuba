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

#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/*! \plugin{direct}{Direct illumination integrator}
 * \order{1}
 * \parameters{
 *     \parameter{shadingSamples}{\Integer}{This convenience parameter can be
 *         used to set both \code{emitterSamples} and \code{bsdfSamples} at
 *         the same time.
 *     }
 *     \parameter{emitterSamples}{\Integer}{Optional more fine-grained
 *        parameter: specifies the number of samples that should be generated
 *        using the direct illumination strategies implemented by the scene's
 *        emitters\default{set to the value of \code{shadingSamples}}
 *     }
 *     \parameter{bsdfSamples}{\Integer}{Optional more fine-grained
 *        parameter: specifies the number of samples that should be generated
 *        using the BSDF sampling strategies implemented by the scene's
 *        surfaces\default{set to the value of \code{shadingSamples}}
 *     }
 *     \parameter{strictNormals}{\Boolean}{Be strict about potential
 *        inconsistencies involving shading normals? See
 *        page~\pageref{sec:strictnormals} for details.
 *        \default{no, i.e. \code{false}}
 *     }
 *     \parameter{hideEmitters}{\Boolean}{Hide directly visible emitters?
 *        See page~\pageref{sec:hideemitters} for details.
 *        \default{no, i.e. \code{false}}
 *     }
 * }
 * \vspace{-1mm}
 * \renderings{
 *     \medrendering{Only BSDF sampling}{integrator_direct_bsdf}
 *     \medrendering{Only emitter sampling}{integrator_direct_lum}
 *     \medrendering{BSDF and emitter sampling}{integrator_direct_both}
 *     \caption{
 *         \label{fig:integrator-direct}
 *         This plugin implements two different strategies for computing the
 *         direct illumination on surfaces. Both of them are dynamically
 *         combined then obtain a robust rendering algorithm.
 *     }
 * }
 *
 * This integrator implements a direct illumination technique that makes use
 * of \emph{multiple importance sampling}: for each pixel sample, the
 * integrator generates a user-specifiable number of BSDF and emitter
 * samples and combines them using the power heuristic. Usually, the BSDF
 * sampling technique works very well on glossy objects but does badly
 * everywhere else (\subfigref{integrator-direct}{a}), while the opposite
 * is true for the emitter sampling technique
 * (\subfigref{integrator-direct}{b}). By combining these approaches, one
 * can obtain a rendering technique that works well in both cases
 * (\subfigref{integrator-direct}{c}).
 *
 * The number of samples spent on either technique is configurable, hence
 * it is also possible to turn this plugin into an emitter sampling-only
 * or BSDF sampling-only integrator.
 *
 * For best results, combine the direct illumination integrator with the
 * low-discrepancy sample generator (\code{ldsampler}). Generally, the number
 * of pixel samples of the sample generator can be kept relatively
 * low (e.g. \code{sampleCount=4}), whereas the \code{shadingSamples}
 * parameter of this integrator should be increased until the variance in
 * the output renderings is acceptable.
 *
 * \remarks{
 *    \item This integrator does not handle participating media or
 *          indirect illumination.
 * }
 */

class MIDirectIntegrator : public SamplingIntegrator {
public:
	MIDirectIntegrator(const Properties &props) : SamplingIntegrator(props) {
		/* Number of shading samples -- this parameter is a shorthand notation
		   to set both 'emitterSamples' and 'bsdfSamples' at the same time*/
		size_t shadingSamples = props.getSize("shadingSamples", 1);

		/* Number of samples to take using the emitter sampling technique */
		m_emitterSamples = props.getSize("emitterSamples", shadingSamples);
		/* Number of samples to take using the BSDF sampling technique */
		m_bsdfSamples = props.getSize("bsdfSamples", shadingSamples);
		/* Be strict about potential inconsistencies involving shading normals? */
		m_strictNormals = props.getBoolean("strictNormals", false);
		/* When this flag is set to true, contributions from directly
		 * visible emitters will not be included in the rendered image */
		m_hideEmitters = props.getBoolean("hideEmitters", false);
		Assert(m_emitterSamples + m_bsdfSamples > 0);
	}

	/// Unserialize from a binary data stream
	MIDirectIntegrator(Stream *stream, InstanceManager *manager)
	 : SamplingIntegrator(stream, manager) {
		m_emitterSamples = stream->readSize();
		m_bsdfSamples = stream->readSize();
		m_strictNormals = stream->readBool();
		m_hideEmitters = stream->readBool();
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		SamplingIntegrator::serialize(stream, manager);
		stream->writeSize(m_emitterSamples);
		stream->writeSize(m_bsdfSamples);
		stream->writeBool(m_strictNormals);
		stream->writeBool(m_hideEmitters);
	}

	void configure() {
		SamplingIntegrator::configure();

		size_t sum = m_emitterSamples + m_bsdfSamples;
		m_weightBSDF = 1 / (Float) m_bsdfSamples;
		m_weightLum = 1 / (Float) m_emitterSamples;
		m_fracBSDF = m_bsdfSamples / (Float) sum;
		m_fracLum = m_emitterSamples / (Float) sum;
	}

	void configureSampler(const Scene *scene, Sampler *sampler) {
		SamplingIntegrator::configureSampler(scene, sampler);
		if (m_emitterSamples > 1)
			sampler->request2DArray(m_emitterSamples);
		if (m_bsdfSamples > 1)
			sampler->request2DArray(m_bsdfSamples);
	}

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		RayDifferential ray(r);
		Spectrum Li(0.0f);
		Point2 sample;

		/* Perform the first ray intersection (or ignore if the
		   intersection has already been provided). */
		if (!rRec.rayIntersect(ray)) {
			/* If no intersection could be found, possibly return
			   radiance from a background emitter */
			if (rRec.type & RadianceQueryRecord::EEmittedRadiance && !m_hideEmitters)
				return scene->evalEnvironment(ray);
			else
				return Spectrum(0.0f);
		}

		/* Possibly include emitted radiance if requested */
		if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance) && !m_hideEmitters)
			Li += its.Le(-ray.d);

		/* Include radiance from a subsurface scattering model if requested */
		if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
			Li += its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

		const BSDF *bsdf = its.getBSDF(ray);

		if (!(rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)
			|| (m_strictNormals && dot(ray.d, its.geoFrame.n)
				* Frame::cosTheta(its.wi) >= 0)) {
			/* Only render the direct illumination component if
			 *
			 * 1. It was requested
			 * 2. The surface has an associated BSDF (i.e. it isn't an index-
			 *    matched medium transition -- this is not supported by 'direct')
			 * 3. If 'strictNormals'=true, when the geometric and shading
			 *    normals classify the incident direction to the same side
			 */
			return Li;
		}

		/* ==================================================================== */
		/*                          Emitter sampling                          */
		/* ==================================================================== */
		bool adaptiveQuery = (rRec.extra & RadianceQueryRecord::EAdaptiveQuery);

		/* Figure out how many BSDF and direct illumination samples to
		   generate, and where the random numbers should come from */
		Point2 *sampleArray;
		size_t numDirectSamples = m_emitterSamples,
			   numBSDFSamples = m_bsdfSamples;
		Float fracLum = m_fracLum, fracBSDF = m_fracBSDF,
		      weightLum = m_weightLum, weightBSDF = m_weightBSDF;

		if (rRec.depth > 1 || adaptiveQuery) {
			/* This integrator is used recursively by another integrator.
			   Be less accurate as this sample will not directly be observed. */
			numBSDFSamples = numDirectSamples = 1;
			fracLum = fracBSDF = .5f;
			weightLum = weightBSDF = 1.0f;
		}

		if (numDirectSamples > 1) {
			sampleArray = rRec.sampler->next2DArray(numDirectSamples);
		} else {
			sample = rRec.nextSample2D(); sampleArray = &sample;
		}

		DirectSamplingRecord dRec(its);
		if (bsdf->getType() & BSDF::ESmooth) {
			/* Only use direct illumination sampling when the surface's
			   BSDF has smooth (i.e. non-Dirac delta) component */
			for (size_t i=0; i<numDirectSamples; ++i) {
				/* Estimate the direct illumination if this is requested */
				Spectrum value = scene->sampleEmitterDirect(dRec, sampleArray[i]);
				if (!value.isZero()) {
					const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

					/* Allocate a record for querying the BSDF */
					BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));

					/* Evaluate BSDF * cos(theta) */
					const Spectrum bsdfVal = bsdf->eval(bRec);

					if (!bsdfVal.isZero() && (!m_strictNormals
							|| dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {
						/* Calculate prob. of sampling that direction using BSDF sampling */
						Float bsdfPdf = emitter->isOnSurface() ? bsdf->pdf(bRec) : 0;

						/* Weight using the power heuristic */
						const Float weight = miWeight(dRec.pdf * fracLum,
								bsdfPdf * fracBSDF) * weightLum;

						Li += value * bsdfVal * weight;
					}
				}
			}
		}

		/* ==================================================================== */
		/*                            BSDF sampling                             */
		/* ==================================================================== */

		if (numBSDFSamples > 1) {
			sampleArray = rRec.sampler->next2DArray(numBSDFSamples);
		} else {
			sample = rRec.nextSample2D(); sampleArray = &sample;
		}

		Intersection bsdfIts;
		for (size_t i=0; i<numBSDFSamples; ++i) {
			/* Sample BSDF * cos(theta) and also request the local density */
			Float bsdfPdf;

			BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
			Spectrum bsdfVal = bsdf->sample(bRec, bsdfPdf, sampleArray[i]);
			if (bsdfVal.isZero())
				continue;

			/* Prevent light leaks due to the use of shading normals */
			const Vector wo = its.toWorld(bRec.wo);
			Float woDotGeoN = dot(its.geoFrame.n, wo);
			if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
				continue;

			/* Trace a ray in this direction */
			Ray bsdfRay(its.p, wo, ray.time);

			Spectrum value;
			if (scene->rayIntersect(bsdfRay, bsdfIts)) {
				/* Intersected something - check if it was an emitter */
				if (!bsdfIts.isEmitter())
					continue;

				value = bsdfIts.Le(-bsdfRay.d);
				dRec.setQuery(bsdfRay, bsdfIts);
			} else {
				/* Intersected nothing -- perhaps there is an environment map? */
				const Emitter *env = scene->getEnvironmentEmitter();

				if (!env || (m_hideEmitters && bRec.sampledType == BSDF::ENull))
					continue;

				value = env->evalEnvironment(RayDifferential(bsdfRay));
				if (!env->fillDirectSamplingRecord(dRec, bsdfRay))
					continue;
			}

			/* Compute the prob. of generating that direction using the
			   implemented direct illumination sampling technique */
			const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
				scene->pdfEmitterDirect(dRec) : 0;

			/* Weight using the power heuristic */
			const Float weight = miWeight(bsdfPdf * fracBSDF,
				lumPdf * fracLum) * weightBSDF;

			Li += value * bsdfVal * weight;
		}

		return Li;
	}

	inline Float miWeight(Float pdfA, Float pdfB) const {
		pdfA *= pdfA; pdfB *= pdfB;
		return pdfA / (pdfA + pdfB);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MIDirectIntegrator[" << endl
			<< "  emitterSamples = " << m_emitterSamples << "," << endl
			<< "  bsdfSamples = " << m_bsdfSamples << "," << endl
			<< "  strictNormals = " << m_strictNormals << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	size_t m_emitterSamples;
	size_t m_bsdfSamples;
	Float m_fracBSDF, m_fracLum;
	Float m_weightBSDF, m_weightLum;
	bool m_strictNormals;
	bool m_hideEmitters;
};

MTS_IMPLEMENT_CLASS_S(MIDirectIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(MIDirectIntegrator, "Direct illumination integrator");
MTS_NAMESPACE_END

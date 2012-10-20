/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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

static StatsCounter avgPathLength("Volumetric path tracer", "Average path length", EAverage);

/*!\plugin{volpath}{Extended volumetric path tracer}
 * \order{4}
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
 *     \parameter{strictNormals}{\Boolean}{Be strict about potential
 *        inconsistencies involving shading normals? See \pluginref{path}
 *        for details.\default{no, i.e. \code{false}}}
 * }
 *
 * This plugin provides a volumetric path tracer that can be used to
 * compute approximate solutions to the radiative transfer equation.
 * Its implementation makes use of multiple importance sampling to
 * combine BSDF and phase function sampling with direct illumination
 * sampling strategies. On surfaces, this integrator behaves exactly
 * like the standard path tracer.
 *
 * \remarks{
 *    \item This integrator will generally perform poorly when rendering
 *      participating media that have a different index of refraction compared
 *      to the surrounding medium.
 *    \item This integrator has poor convergence properties when rendering
 *    caustics and similar effects. In this case, \pluginref{bdpt} or
 *    one of the photon mappers may be preferable.
 * }
 */
class VolumetricPathTracer : public MonteCarloIntegrator {
public:
	VolumetricPathTracer(const Properties &props) : MonteCarloIntegrator(props) { }

	/// Unserialize from a binary data stream
	VolumetricPathTracer(Stream *stream, InstanceManager *manager)
	 : MonteCarloIntegrator(stream, manager) { }

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		MediumSamplingRecord mRec;
		RayDifferential ray(r);
		Spectrum Li(0.0f);
		bool scattered = false;
		Float eta = 1.0f;

		/* Perform the first ray intersection (or ignore if the
		   intersection has already been provided). */
		rRec.rayIntersect(ray);

		Spectrum throughput(1.0f);

		while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
			/* ==================================================================== */
			/*                 Radiative Transfer Equation sampling                 */
			/* ==================================================================== */
			if (rRec.medium && rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler)) {
				/* Sample the integral
				   \int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
				*/
				const PhaseFunction *phase = mRec.getPhaseFunction();

				if (rRec.depth >= m_maxDepth && m_maxDepth != -1) // No more scattering events allowed
					break;

				throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

				/* ==================================================================== */
				/*                          Luminaire sampling                          */
				/* ==================================================================== */

				/* Estimate the single scattering component if this is requested */
				DirectSamplingRecord dRec(mRec.p, mRec.time);

				if (rRec.type & RadianceQueryRecord::EDirectMediumRadiance) {
					int interactions = m_maxDepth - rRec.depth - 1;

					Spectrum value = scene->sampleAttenuatedEmitterDirect(
							dRec, rRec.medium, interactions,
							rRec.nextSample2D(), rRec.sampler);

					if (!value.isZero()) {
						const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

						/* Evaluate the phase function */
						PhaseFunctionSamplingRecord pRec(mRec, -ray.d, dRec.d);
						Float phaseVal = phase->eval(pRec);

						if (phaseVal != 0) {
							/* Calculate prob. of having sampled that direction using
							   phase function sampling */
							Float phasePdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle
									&& interactions == 0)
									? phase->pdf(pRec) : (Float) 0.0f;

							/* Weight using the power heuristic */
							const Float weight = miWeight(dRec.pdf, phasePdf);
							Li += throughput * value * phaseVal * weight;
						}
					}
				}

				/* ==================================================================== */
				/*                         Phase function sampling                      */
				/* ==================================================================== */

				Float phasePdf;
				PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
				Float phaseVal = phase->sample(pRec, phasePdf, rRec.sampler);
				if (phaseVal == 0)
					break;

				bool hitEmitter = false;
				Spectrum value;

				/* Trace a ray in this direction */
				ray = Ray(mRec.p, pRec.wo, ray.time);
				ray.mint = 0;

				if (scene->rayIntersect(ray, its)) {
					/* Intersected something - check if it was a luminaire */
					if (its.isEmitter()) {
						value = its.Le(-ray.d);
						dRec.setQuery(ray, its);
						hitEmitter = true;
					}
				} else {
					/* Intersected nothing -- perhaps there is an environment map? */
					const Emitter *env = scene->getEnvironmentEmitter();

					if (env) {
						value = env->evalEnvironment(ray);
						if (!env->fillDirectSamplingRecord(dRec, ray))
							break;
						hitEmitter = true;
					} else {
						break;
					}
				}

				throughput *= phaseVal;

				/* If a luminaire was hit, estimate the local illumination and
				   weight using the power heuristic */
				if (hitEmitter && (rRec.type & RadianceQueryRecord::EDirectMediumRadiance)) {
					Spectrum transmittance = rRec.medium->evalTransmittance(Ray(ray, 0, its.t));
					const Float emitterPdf = scene->pdfEmitterDirect(dRec);
					Li += throughput * value * transmittance * miWeight(phasePdf, emitterPdf);
				}

				/* ==================================================================== */
				/*                         Multiple scattering                          */
				/* ==================================================================== */

				/* Stop if multiple scattering was not requested */
				if (!(rRec.type & RadianceQueryRecord::EIndirectMediumRadiance))
					break;
				rRec.type = RadianceQueryRecord::ERadianceNoEmission;

				scattered = true;
			} else {
				/* Sample
					tau(x, y) (Surface integral). This happens with probability mRec.pdfFailure
					Account for this and multiply by the proper per-color-channel transmittance.
				*/
				if (rRec.medium)
					throughput *= mRec.transmittance / mRec.pdfFailure;

				if (!its.isValid()) {
					/* If no intersection could be found, possibly return
					   attenuated radiance from a background luminaire */
					if (rRec.type & RadianceQueryRecord::EEmittedRadiance)
						Li += throughput * scene->evalEnvironment(ray);
					break;
				}

				/* Possibly include emitted radiance if requested */
				if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance))
					Li += throughput * its.Le(-ray.d);

				/* Include radiance from a subsurface integrator if requested */
				if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
					Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

				if (rRec.depth >= m_maxDepth && m_maxDepth != -1)
					break;

				/* Prevent light leaks due to the use of shading normals */
				Float wiDotGeoN = -dot(its.geoFrame.n, ray.d),
					  wiDotShN  = Frame::cosTheta(its.wi);
				if (wiDotGeoN * wiDotShN < 0 && m_strictNormals)
					break;

				/* ==================================================================== */
				/*                          Luminaire sampling                          */
				/* ==================================================================== */

				const BSDF *bsdf = its.getBSDF(ray);
				DirectSamplingRecord dRec(its);

				/* Estimate the direct illumination if this is requested */
				if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
						(bsdf->getType() & BSDF::ESmooth)) {
					int interactions = m_maxDepth - rRec.depth - 1;

					Spectrum value = scene->sampleAttenuatedEmitterDirect(
							dRec, its, rRec.medium, interactions,
							rRec.nextSample2D(), rRec.sampler);

					if (!value.isZero()) {
						const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

						/* Evaluate BSDF * cos(theta) */
						BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));
						const Spectrum bsdfVal = bsdf->eval(bRec);

						Float woDotGeoN = dot(its.geoFrame.n, dRec.d);

						/* Prevent light leaks due to the use of shading normals */
						if (!bsdfVal.isZero() && (!m_strictNormals ||
							woDotGeoN * Frame::cosTheta(bRec.wo) > 0)) {
							/* Calculate prob. of having generated that direction
							   using BSDF sampling */
							Float bsdfPdf = (emitter->isOnSurface()
									&& dRec.measure == ESolidAngle
									&& interactions == 0)
									? bsdf->pdf(bRec) : (Float) 0.0f;

							/* Weight using the power heuristic */
							const Float weight = miWeight(dRec.pdf, bsdfPdf);
							Li += throughput * value * bsdfVal * weight;
						}
					}
				}


				/* ==================================================================== */
				/*                            BSDF sampling                             */
				/* ==================================================================== */

				/* Sample BSDF * cos(theta) */
				BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
				Float bsdfPdf;
				Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
				if (bsdfWeight.isZero())
					break;

				/* Prevent light leaks due to the use of shading normals */
				const Vector wo = its.toWorld(bRec.wo);
				Float woDotGeoN = dot(its.geoFrame.n, wo);
				if (woDotGeoN * Frame::cosTheta(bRec.wo) <= 0 && m_strictNormals)
					break;

				/* Trace a ray in this direction */
				ray = Ray(its.p, wo, ray.time);

				/* Keep track of the throughput, medium, and relative
				   refractive index along the path */
				throughput *= bsdfWeight;
				eta *= bRec.eta;
				if (its.isMediumTransition())
					rRec.medium = its.getTargetMedium(ray.d);

				bool hitEmitter = false;
				Spectrum value;

				if (scene->rayIntersect(ray, its)) {
					/* Intersected something - check if it was a luminaire */
					if (its.isEmitter()) {
						value = its.Le(-ray.d);
						dRec.setQuery(ray, its);
						hitEmitter = true;
					}
				} else {
					/* Intersected nothing -- perhaps there is an environment map? */
					const Emitter *env = scene->getEnvironmentEmitter();

					if (env) {
						value = env->evalEnvironment(ray);
						if (!env->fillDirectSamplingRecord(dRec, ray))
							break;
						hitEmitter = true;
					} else {
						break;
					}
				}

				/* If a luminaire was hit, estimate the local illumination and
				   weight using the power heuristic */
				if (hitEmitter && (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)
						&& !((bRec.sampledType & BSDF::ENull) && scattered)) {
					Spectrum transmittance = rRec.medium ?
						rRec.medium->evalTransmittance(Ray(ray, 0, its.t)) : Spectrum(1.0f);
					const Float emitterPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
						scene->pdfEmitterDirect(dRec) : 0;
					Li += throughput * value * transmittance * miWeight(bsdfPdf, emitterPdf);
				}

				/* ==================================================================== */
				/*                         Indirect illumination                        */
				/* ==================================================================== */

				/* Stop if indirect illumination was not requested */
				if (!(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
					break;
				rRec.type = RadianceQueryRecord::ERadianceNoEmission;

				scattered |= bRec.sampledType != BSDF::ENull;
			}

			if (rRec.depth++ >= m_rrDepth) {
				/* Russian roulette: try to keep path weights equal to one,
				   while accounting for the solid angle compression at refractive
				   index boundaries. Stop with at least some probability to avoid
				   getting stuck (e.g. due to total internal reflection) */

				Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
				if (rRec.nextSample1D() >= q)
					break;
				throughput /= q;
			}
		}
		avgPathLength.incrementBase();
		avgPathLength += rRec.depth;
		return Li;
	}

	inline Float miWeight(Float pdfA, Float pdfB) const {
		pdfA *= pdfA; pdfB *= pdfB;
		return pdfA / (pdfA + pdfB);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		MonteCarloIntegrator::serialize(stream, manager);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "VolumetricPathTracer[" << endl
			<< "  maxDepth = " << m_maxDepth << "," << endl
			<< "  rrDepth = " << m_rrDepth << "," << endl
			<< "  strictNormals = " << m_strictNormals << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(VolumetricPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(VolumetricPathTracer, "Volumetric path tracer");
MTS_NAMESPACE_END

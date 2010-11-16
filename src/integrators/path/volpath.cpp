/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/**
 * Volumetric path tracer, which solves the full radiative transfer
 * equation in the presence of participating media. Estimates single
 * scattering using both phase function and luminaire sampling and
 * combines the two with multiple importance sampling and the power
 * heuristic. Afterwards, the phase function sample is reused to
 * recursively estimate the multiple scattering component, which 
 * saves an intersection computation.
 * On surfaces, this integrator behaves exactly like the standard
 * MI path tracer.
 */
class VolumetricPathTracer : public MonteCarloIntegrator {
public:
	VolumetricPathTracer(const Properties &props) : MonteCarloIntegrator(props) {
		/* Beta factor for the power heuristic */
		m_beta = props.getFloat("beta", 2.0f);
	}

	/// Unserialize from a binary data stream
	VolumetricPathTracer(Stream *stream, InstanceManager *manager)
	 : MonteCarloIntegrator(stream, manager) {
		m_beta = stream->readFloat();
	}

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		LuminaireSamplingRecord lRec;
		MediumSamplingRecord mRec;
		RayDifferential ray(r);
		Intersection prevIts;
		Spectrum Li(0.0f);

		/* Perform the first ray intersection (or ignore if the 
		   intersection has already been provided). */
		rRec.rayIntersect(ray);
		Spectrum pathThroughput(1.0f);

		while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
			/* ==================================================================== */
			/*                 Radiative Transfer Equation sampling                 */
			/* ==================================================================== */
			if (scene->sampleDistance(ray, its.t, mRec, rRec.sampler)) {
				const PhaseFunction *phase = mRec.medium->getPhaseFunction();
				Vector wo, wi = -ray.d;

				if (rRec.depth == m_maxDepth && m_maxDepth > 0) // No more scattering events allowed
					break;

				/* Sample the integral
				   \int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
				*/
				pathThroughput *= mRec.sigmaS * mRec.attenuation * mRec.miWeight / mRec.pdf;

				/* ==================================================================== */
				/*                          Luminaire sampling                          */
				/* ==================================================================== */

				/* Estimate the single scattering component if this is requested */
				if (rRec.type & RadianceQueryRecord::EInscatteredDirectRadiance && 
					scene->sampleLuminaireAttenuated(mRec.p, lRec, ray.time, rRec.nextSample2D())) {
					/* Evaluate the phase function */
					Spectrum phaseVal = phase->f(mRec, -ray.d, -lRec.d);

					if (phaseVal.max() > 0) {
						/* Calculate prob. of having sampled that direction using 
						   phase function sampling */
						Float phasePdf = (lRec.luminaire->isIntersectable() 
								|| lRec.luminaire->isBackgroundLuminaire()) ? 
							phase->pdf(mRec, -ray.d, -lRec.d) : 0;

						/* Weight using the power heuristic */
						const Float weight = miWeight(lRec.pdf, phasePdf);
						Li += pathThroughput * lRec.Le * phaseVal * weight;
					}
				}

				/* ==================================================================== */
				/*                         Phase function sampling                      */
				/* ==================================================================== */

				PhaseFunction::ESampledType sampledType;
				Float phasePdf;
				Spectrum phaseVal = phase->sample(mRec, wi, wo, sampledType, phasePdf, rRec.nextSample2D());
				if (phaseVal.max() == 0)
					break;
				phaseVal /= phasePdf;
				prevIts = its;

				/* Trace a ray in this direction */
				ray = Ray(mRec.p, wo, ray.time);
				bool hitLuminaire = false;
				if (scene->rayIntersect(ray, its)) {
					/* Intersected something - check if it was a luminaire */
					if (its.isLuminaire()) {
						lRec = LuminaireSamplingRecord(its, -ray.d);
						hitLuminaire = true;
					}
				} else {
					/* No intersection found. Possibly, there is a background
					   luminaire such as an environment map? */
					if (scene->hasBackgroundLuminaire()) {
						lRec.luminaire = scene->getBackgroundLuminaire();
						lRec.d = -ray.d;
						hitLuminaire = true;
					}
				}
				ray.mint = 0; ray.maxt = its.t; 
				rRec.attenuation = scene->getAttenuation(ray);

				/* If a luminaire was hit, estimate the local illumination and
				   weight using the power heuristic */
				if (hitLuminaire && (rRec.type & RadianceQueryRecord::EInscatteredDirectRadiance)) {
					lRec.Le = lRec.luminaire->Le(lRec);

					/* Prob. of having generated this sample using luminaire sampling */
					const Float lumPdf = (sampledType != PhaseFunction::EDelta) 
						? scene->pdfLuminaire(mRec.p, lRec) : 0;
					Float weight = miWeight(phasePdf, lumPdf);
					Li += pathThroughput * rRec.attenuation * lRec.Le * phaseVal * weight;
				}

				/* ==================================================================== */
				/*                         Multiple scattering                          */
				/* ==================================================================== */

				/* Stop if multiple scattering was not requested */
				if (!(rRec.type & RadianceQueryRecord::EInscatteredIndirectRadiance)) 
					break;
				rRec.type = RadianceQueryRecord::ERadianceNoEmission;

				/* Russian roulette - Possibly stop the recursion */
				if (rRec.depth >= m_rrDepth) {
					if (rRec.nextSample1D() > mRec.albedo)
						break;
					else
						pathThroughput /= mRec.albedo;
				}

				pathThroughput *= phaseVal;
				rRec.depth++;
			} else {
				/* Sample 
					tau(x, y) (Surface integral). This happens with probability mRec.pdf
					Divide this out and multiply with the proper per-color-channel attenuation.
				*/
				pathThroughput *= mRec.attenuation * mRec.miWeight / mRec.pdf;

				if (!its.isValid()) {
					/* If no intersection could be found, possibly return 
					   attenuated radiance from a background luminaire */
					if (rRec.type & RadianceQueryRecord::EEmittedRadiance)
						Li += pathThroughput * scene->LeBackground(ray);
					break;
				}

				const BSDF *bsdf = its.getBSDF(ray);

				/* Possibly include emitted radiance if requested */
				if (its.isLuminaire() && (rRec.type & RadianceQueryRecord::EEmittedRadiance))
					Li += pathThroughput * its.Le(-ray.d);

				/* Include radiance from a subsurface integrator if requested */
				if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
					Li += pathThroughput * its.LoSub(rRec.scene, -ray.d);

				if (rRec.depth == m_maxDepth && m_maxDepth > 0)
					break;

				/* ==================================================================== */
				/*                          Luminaire sampling                          */
				/* ==================================================================== */

				/* Estimate the direct illumination if this is requested */
				if (rRec.type & RadianceQueryRecord::EDirectRadiance && 
					scene->sampleLuminaireAttenuated(its, lRec, rRec.nextSample2D())) {
					/* Allocate a record for querying the BSDF */
					const BSDFQueryRecord bRec(rRec, its, its.toLocal(-lRec.d));

					/* Evaluate BSDF * cos(theta) */
					const Spectrum bsdfVal = bsdf->fCos(bRec);

					if (!bsdfVal.isBlack()) {
						/* Calculate prob. of having sampled that direction
						   using BSDF sampling */
						Float bsdfPdf = (lRec.luminaire->isIntersectable() 
								|| lRec.luminaire->isBackgroundLuminaire()) ? 
							bsdf->pdf(bRec) : 0;

						/* Weight using the power heuristic */
						const Float weight = miWeight(lRec.pdf, bsdfPdf);
						Li += pathThroughput * lRec.Le * bsdfVal * weight;
					}
				}

				/* ==================================================================== */
				/*                            BSDF sampling                             */
				/* ==================================================================== */

				/* Sample BSDF * cos(theta) */
				BSDFQueryRecord bRec(rRec, its, rRec.nextSample2D());
				Float bsdfPdf;
				Spectrum bsdfVal = bsdf->sampleCos(bRec, bsdfPdf);
				if (bsdfVal.isBlack())
					break;
				bsdfVal /= bsdfPdf;
				prevIts = its;

				/* Trace a ray in this direction */
				ray = Ray(its.p, its.toWorld(bRec.wo), ray.time);
				bool hitLuminaire = false;
				if (scene->rayIntersect(ray, its)) {
					/* Intersected something - check if it was a luminaire */
					if (its.isLuminaire()) {
						lRec = LuminaireSamplingRecord(its, -ray.d);
						hitLuminaire = true;
					}
				} else {
					/* No intersection found. Possibly, there is a background
					   luminaire such as an environment map? */
					if (scene->hasBackgroundLuminaire()) {
						lRec.luminaire = scene->getBackgroundLuminaire();
						lRec.d = -ray.d;
						hitLuminaire = true;
					}
				}
				ray.mint = 0; ray.maxt = its.t; 
				rRec.attenuation = scene->getAttenuation(ray);

				/* If a luminaire was hit, estimate the local illumination and
				   weight using the power heuristic */
				if (hitLuminaire && (rRec.type & RadianceQueryRecord::EDirectRadiance)
					&& !(bRec.sampledType & BSDF::EDelta)) {
					lRec.Le = lRec.luminaire->Le(lRec);
					/* Prob. of having generated this sample using luminaire sampling */
					const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
						scene->pdfLuminaire(prevIts, lRec) : 0;
					const Float weight = miWeight(bsdfPdf, lumPdf);
					Li += pathThroughput * rRec.attenuation * lRec.Le * bsdfVal * weight;
				}

				/* ==================================================================== */
				/*                         Indirect illumination                        */
				/* ==================================================================== */
			
				/* Stop if indirect illumination was not requested */
				if (!(rRec.type & RadianceQueryRecord::EIndirectRadiance)) 
					break;
				rRec.type = RadianceQueryRecord::ERadianceNoEmission;

				/* Russian roulette - Possibly stop the recursion */
				if (rRec.depth >= m_rrDepth) {
					/* Assuming that BSDF importance sampling is perfect,
					   'bsdfVal.max()' should equal the maximum albedo
					   over all spectral samples */
					Float approxAlbedo = std::min((Float) 0.9f, bsdfVal.max());
					if (rRec.nextSample1D() > approxAlbedo)
						break;
					else
						pathThroughput /= approxAlbedo;
				}

				pathThroughput *= bsdfVal;
				rRec.depth++;
			}
		}
		return Li;
	}

	inline Float miWeight(Float pdfA, Float pdfB) const {
		pdfA = std::pow(pdfA, m_beta);
		pdfB = std::pow(pdfB, m_beta);
		return pdfA / (pdfA + pdfB);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		MonteCarloIntegrator::serialize(stream, manager);
		stream->writeFloat(m_beta);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "VolumetricPathTracer[" << std::endl
			<< "  beta = " << m_beta << "," << std::endl
			<< "  maxDepth = " << m_maxDepth << "," << std::endl
			<< "  rrDepth = " << m_rrDepth << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Float m_beta;
};

MTS_IMPLEMENT_CLASS_S(VolumetricPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(VolumetricPathTracer, "Volumetric path tracer");
MTS_NAMESPACE_END

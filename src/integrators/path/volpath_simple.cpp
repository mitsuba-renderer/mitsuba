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

static StatsCounter avgPathLength("Volumetric path tracer", "Average path length", EAverage);

/**
 * Volumetric path tracer, which solves the full radiative transfer
 * equation in the presence of participating media. Simplified version
 * without multiple importance sampling - this version can be much
 * faster than the multiple importance sampling version when when 
 * rendering heterogeneous participating media using the 
 * [Coleman et al.] sampling technique, since fewer transmittance
 * evaluations will be required.
 */
class SimpleVolumetricPathTracer : public MonteCarloIntegrator {
public:
	SimpleVolumetricPathTracer(const Properties &props) : MonteCarloIntegrator(props) { }

	/// Unserialize from a binary data stream
	SimpleVolumetricPathTracer(Stream *stream, InstanceManager *manager)
	 : MonteCarloIntegrator(stream, manager) { }

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		LuminaireSamplingRecord lRec;
		MediumSamplingRecord mRec;
		RayDifferential ray(r);
		Spectrum Li(0.0f);

		/* Perform the first ray intersection (or ignore if the 
		   intersection has already been provided). */
		rRec.rayIntersect(ray);
		Spectrum pathThroughput(1.0f);
		bool computeIntersection = false;

		while (rRec.depth < m_maxDepth || m_maxDepth < 0) {
			if (computeIntersection)
				scene->rayIntersect(ray, its);

			/* ==================================================================== */
			/*                 Radiative Transfer Equation sampling                 */
			/* ==================================================================== */
			if (rRec.medium && rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler)) {
				/* Sample the integral
				   \int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
				*/
				const PhaseFunction *phase = rRec.medium->getPhaseFunction();

				pathThroughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

				/* ==================================================================== */
				/*                          Luminaire sampling                          */
				/* ==================================================================== */

				/* Estimate the single scattering component if this is requested */
				if (rRec.type & RadianceQueryRecord::EDirectMediumRadiance && 
					scene->sampleAttenuatedLuminaire(mRec.p, ray.time,
						rRec.medium, lRec, rRec.nextSample2D(), rRec.sampler)) {
					Li += pathThroughput * lRec.value * phase->eval(
							PhaseFunctionQueryRecord(mRec, -ray.d, -lRec.d));
				}

				/* ==================================================================== */
				/*                         Phase function sampling                      */
				/* ==================================================================== */

				PhaseFunctionQueryRecord pRec(mRec, -ray.d);
				Float phaseVal = phase->sample(pRec, rRec.sampler);
				if (phaseVal == 0)
					break;

				/* Trace a ray in this direction */
				ray = Ray(mRec.p, pRec.wo, ray.time);
				ray.mint = 0;
				computeIntersection = true;

				/* ==================================================================== */
				/*                         Multiple scattering                          */
				/* ==================================================================== */

				if (!(rRec.type & RadianceQueryRecord::EIndirectMediumRadiance))
					break; /* Stop if multiple scattering was not requested */

				/* Russian roulette - Possibly stop the recursion */
				if (rRec.depth >= m_rrDepth) {
					if (rRec.nextSample1D() > mRec.albedo)
						break;
					else
						pathThroughput /= mRec.albedo;
				}

				pathThroughput *= phaseVal;
				rRec.depth++;
				rRec.type = RadianceQueryRecord::ERadianceNoEmission;
			} else {
				/* Sample 
					tau(x, y) * (Surface integral). This happens with probability mRec.pdfFailure
					Account for this and multiply by the proper per-color-channel transmittance.
				*/

				if (rRec.medium)
					pathThroughput *= mRec.transmittance / mRec.pdfFailure;

				if (!its.isValid()) {
					/* If no intersection could be found, possibly return 
					   attenuated radiance from a background luminaire */
					if (rRec.type & RadianceQueryRecord::EEmittedRadiance) 
						Li += pathThroughput * scene->LeBackground(ray);
					break;
				}

				computeIntersection = true;

				/* Possibly include emitted radiance if requested */
				if (its.isLuminaire() && (rRec.type & RadianceQueryRecord::EEmittedRadiance))
					Li += pathThroughput * its.Le(-ray.d);

				/* Include radiance from a subsurface integrator if requested */
				if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
					Li += pathThroughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

				const BSDF *bsdf = its.getBSDF(ray);

				if (!bsdf) {
					/* Pass right through the surface (there is no BSDF) */
					if (its.isMediumTransition())
						rRec.medium = its.getTargetMedium(ray.d);
					ray.setOrigin(its.p);
					ray.mint = Epsilon;
					continue;
				}
				
				/* Prevent light leaks due to the use of shading normals */
				Float wiDotGeoN = -dot(its.geoFrame.n, ray.d),
					  wiDotShN  = Frame::cosTheta(its.wi);
				if (wiDotGeoN * wiDotShN < 0 && m_strictNormals) 
					break;

				/* ==================================================================== */
				/*                          Luminaire sampling                          */
				/* ==================================================================== */

				/* Estimate the direct illumination if this is requested */
				if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance && 
					scene->sampleAttenuatedLuminaire(its, rRec.medium, lRec, 
						rRec.nextSample2D(), rRec.sampler)) {
					/* Allocate a record for querying the BSDF */
					const Vector wo = -lRec.d;
					BSDFQueryRecord bRec(its, its.toLocal(wo));
					bRec.sampler = rRec.sampler;
					
					Float woDotGeoN = dot(its.geoFrame.n, wo);
					/* Prevent light leaks due to the use of shading normals */
					if (!m_strictNormals ||
						woDotGeoN * Frame::cosTheta(bRec.wo) > 0)
						Li += pathThroughput * lRec.value * bsdf->eval(bRec);
				}

				/* ==================================================================== */
				/*                            BSDF sampling                             */
				/* ==================================================================== */

				/* Sample BSDF * cos(theta) */
				BSDFQueryRecord bRec(its, rRec.sampler, ERadiance);
				Spectrum bsdfVal = bsdf->sample(bRec, rRec.nextSample2D());
				if (bsdfVal.isZero()) 
					break;
	
				/* Prevent light leaks due to the use of shading normals */
				const Vector wo = its.toWorld(bRec.wo);
				Float woDotGeoN = dot(its.geoFrame.n, wo);
				if (woDotGeoN * Frame::cosTheta(bRec.wo) <= 0 && m_strictNormals)
					break;

				if (its.isMediumTransition())
					rRec.medium = its.getTargetMedium(wo);

				/* In the next iteration, trace a ray in this direction */
				ray = Ray(its.p, wo, ray.time);

				/* ==================================================================== */
				/*                         Indirect illumination                        */
				/* ==================================================================== */
				bool includeEmitted = (bRec.sampledType & BSDF::EDelta) 
					&& (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance);

				if (!(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance)) {
					/* Stop if indirect illumination was not requested (except: just sampled a
					   delta BSDF -- recursively look for emitted radiance to get the
					   direct component of the current interaction) */
					if (includeEmitted) {
						rRec.type = RadianceQueryRecord::EEmittedRadiance;
						if (rRec.depth+1 == m_maxDepth) {
							/* Do one extra recursion to get the emitted radiance component 
							   (e.g. from an envmap). We must reduce rRec.depth 
							   or the loop will terminate */
							--rRec.depth;
						}
					} else {
						break;
					}
				} else {
					if (!includeEmitted) {
						rRec.type = RadianceQueryRecord::ERadianceNoEmission;
					} else {
						if (rRec.depth+1 == m_maxDepth) {
							/* Do one extra recursion to get the emitted radiance component 
							   (e.g. from an envmap). We must reduce rRec.depth 
							   or the loop will terminate */
							--rRec.depth;
							rRec.type = RadianceQueryRecord::EEmittedRadiance;
						} else {
							rRec.type = RadianceQueryRecord::ERadiance;
						}
					}
				}

				/* Russian roulette - Possibly stop the recursion. Don't do this when
				   dealing with a transmission component, since solid angle compression
				   factors cause problems with the heuristic below */
				if (rRec.depth >= m_rrDepth && !(bRec.sampledType & BSDF::ETransmission)) {
					/* Assuming that BSDF importance sampling is perfect,
					   'bsdfVal.max()' should equal the maximum albedo
					   over all spectral samples */
					Float approxAlbedo = std::min((Float) 0.9, bsdfVal.max());
					if (rRec.nextSample1D() > approxAlbedo)
						break;
					else
						pathThroughput /= approxAlbedo;
				}

				pathThroughput *= bsdfVal;
				rRec.depth++;
			}
		}
		avgPathLength.incrementBase();
		avgPathLength += rRec.depth;
		return Li;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		MonteCarloIntegrator::serialize(stream, manager);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SimpleVolumetricPathTracer[" << std::endl
			<< "  maxDepth = " << m_maxDepth << "," << std::endl
			<< "  rrDepth = " << m_rrDepth << "," << std::endl
			<< "  strictNormals = " << m_strictNormals << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(SimpleVolumetricPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(SimpleVolumetricPathTracer, "Simple volumetric path tracer");
MTS_NAMESPACE_END

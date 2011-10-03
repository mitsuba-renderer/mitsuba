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

static StatsCounter avgPathLength("Path tracer", "Average path length", EAverage);

/*! \plugin{path}{Path tracer}
 *
 * \parameters{
 *     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *         in the generated output image (where \code{-1} corresponds to $\infty$).
 *	       A value of \code{1} will only render directly visible light sources.
 *	       \code{2} will lead to single-bounce (direct-only) illumination, 
 *	       and so on. \default{\code{-1}}
 *	   }
 *	   \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after 
 *	      which the implementation will start to use the ``russian roulette'' 
 *	      path termination criterion. \default{\code{10}}
 *	   }
 *     \parameter{strictNormals}{\Boolean}{Be strict about potential
 *        inconsistencies involving shading normals? See \pluginref{path}
 *        for details.\default{no, i.e. \code{false}}}
 * }
 *
 * This integrator implements a path tracer that makes use
 * of \emph{multiple importance sampling}: for each surface interaction, the 
 * integrator generates a single BSDF and emitter sample and combines
 * them using the power heuristic.
 *
 * For best results, combine the path tracer with the
 * low-discrepancy sample generator (\code{ldsampler}).
 *
 * \remarks{
 *    \item This integrator does not handle participating media
 * }
 */

class MIPathTracer : public MonteCarloIntegrator {
public:
	MIPathTracer(const Properties &props)
		: MonteCarloIntegrator(props) { }

	/// Unserialize from a binary data stream
	MIPathTracer(Stream *stream, InstanceManager *manager)
		: MonteCarloIntegrator(stream, manager) { }

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		RayDifferential ray(r);
		Spectrum Li(0.0f);

		/* Perform the first ray intersection (or ignore if the 
		   intersection has already been provided). */
		rRec.rayIntersect(ray);
		ray.mint = Epsilon;

		Spectrum pathThroughput(1.0f);

		while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
			if (!its.isValid()) {
				/* If no intersection could be found, potentially return 
				   radiance from a background luminaire if it exists */
				if (rRec.type & RadianceQueryRecord::EEmittedRadiance)
					Li += pathThroughput * scene->LeBackground(ray);
				break;
			}

			const BSDF *bsdf = its.getBSDF(ray);

			if (EXPECT_NOT_TAKEN(bsdf == NULL)) {
				/* The MI path tracer doesn't support
				   surfaces without a BSDF (e.g. medium transitions)
				   -- give up. */
				break;
			}

			/* Possibly include emitted radiance if requested */
			if (its.isLuminaire() && (rRec.type & RadianceQueryRecord::EEmittedRadiance))
				Li += pathThroughput * its.Le(-ray.d);

			/* Include radiance from a subsurface integrator if requested */
			if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
				Li += pathThroughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

			if (m_maxDepth > 0 && rRec.depth >= m_maxDepth)
				break;

			/* ==================================================================== */
			/*                          Luminaire sampling                          */
			/* ==================================================================== */

			/* Prevent light leaks due to the use of shading normals */
			Float wiDotGeoN = -dot(its.geoFrame.n, ray.d),
				  wiDotShN  = Frame::cosTheta(its.wi);
			if (wiDotGeoN * wiDotShN < 0 && m_strictNormals) 
				break;

			/* Estimate the direct illumination if this is requested */
			LuminaireSamplingRecord lRec;
			if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance && 
				scene->sampleLuminaire(its.p, ray.time, lRec, rRec.nextSample2D())) {
				/* Allocate a record for querying the BSDF */
				const Vector wo = -lRec.d;
				BSDFQueryRecord bRec(its, its.toLocal(wo));
				bRec.sampler = rRec.sampler;
	
				/* Evaluate BSDF * cos(theta) */
				const Spectrum bsdfVal = bsdf->eval(bRec);

				Float woDotGeoN = dot(its.geoFrame.n, wo);

				/* Prevent light leaks due to the use of shading normals */
				if (!bsdfVal.isZero() && (!m_strictNormals
						|| woDotGeoN * Frame::cosTheta(bRec.wo) > 0)) {
					/* Calculate prob. of having sampled that direction
					   using BSDF sampling */
					Float bsdfPdf = (lRec.luminaire->isIntersectable() 
							|| lRec.luminaire->isBackgroundLuminaire()) ? 
						bsdf->pdf(bRec) : 0;

					/* Weight using the power heuristic */
					const Float weight = miWeight(lRec.pdf, bsdfPdf);
					Li += pathThroughput * lRec.value * bsdfVal * weight;
				}
			}

			/* ==================================================================== */
			/*                            BSDF sampling                             */
			/* ==================================================================== */

			/* Sample BSDF * cos(theta) */
			BSDFQueryRecord bRec(its, rRec.sampler, ERadiance);
			Float bsdfPdf;
			Spectrum bsdfVal = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
			if (bsdfVal.isZero()) 
				break;
	
			/* Prevent light leaks due to the use of shading normals */
			const Vector wo = its.toWorld(bRec.wo);
			Float woDotGeoN = dot(its.geoFrame.n, wo);
			if (woDotGeoN * Frame::cosTheta(bRec.wo) <= 0 && m_strictNormals)
				break;

			/* Trace a ray in this direction */
			ray = Ray(its.p, wo, ray.time);

			bool hitLuminaire = false;
			if (scene->rayIntersect(ray, its)) {
				/* Intersected something - check if it was a luminaire */
				if (its.isLuminaire()) {
					lRec = LuminaireSamplingRecord(its, -ray.d);
					lRec.value = its.Le(-ray.d);
					hitLuminaire = true;
				}
			} else {
				/* No intersection found. Possibly, there is a background
					luminaire such as an environment map? */
				if (scene->hasBackgroundLuminaire()) {
					lRec.luminaire = scene->getBackgroundLuminaire();
					lRec.value = lRec.luminaire->Le(ray);
					lRec.d = -ray.d;
					hitLuminaire = true;
				} else {
					rRec.depth++;
					break;
				}
			}

			/* If a luminaire was hit, estimate the local illumination and
			   weight using the power heuristic */
			if (hitLuminaire &&  
				(rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
				/* Prob. of having generated this sample using luminaire sampling */
				const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
					scene->pdfLuminaire(ray.o, lRec) : 0;
				const Float weight = miWeight(bsdfPdf, lumPdf);
				Li += pathThroughput * lRec.value * bsdfVal * weight;
			}

			/* ==================================================================== */
			/*                         Indirect illumination                        */
			/* ==================================================================== */

			/* Set the recursive query type */
			/* Stop if no surface was hit by the BSDF sample or if indirect illumination
			   was not requested */
			if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance)) 
				break;
			rRec.type = RadianceQueryRecord::ERadianceNoEmission;

			/* Russian roulette - Possibly stop the recursion. Don't do this when
			   dealing with a transmission component, since solid angle compression
			   factors cause problems with the heuristic below */
			if (rRec.depth >= m_rrDepth && !(bRec.sampledType & BSDF::ETransmission)) {
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

		/* Store statistics */
		avgPathLength.incrementBase();
		avgPathLength += rRec.depth;

		return Li;
	}

	inline Float miWeight(Float pdfA, Float pdfB) const {
		pdfA *= pdfA;
		pdfB *= pdfB;
		return pdfA / (pdfA + pdfB);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		MonteCarloIntegrator::serialize(stream, manager);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MIPathTracer[" << std::endl
			<< "  maxDepth = " << m_maxDepth << "," << std::endl
			<< "  rrDepth = " << m_rrDepth << "," << std::endl
			<< "  strictNormals = " << m_strictNormals << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(MIPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(MIPathTracer, "MI path tracer");
MTS_NAMESPACE_END

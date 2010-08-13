#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

 /**
  * Extended path tracer -- uses multiple importance sampling to combine 
  * two sampling strategies, namely BSDF and luminaire sampling. 
  * This class does not attempt to solve the full radiative transfer 
  * equation (see <tt>volpath</tt> if this is needed).
  */
class MIPathTracer : public MonteCarloIntegrator {
public:
	MIPathTracer(const Properties &props) : MonteCarloIntegrator(props) {
		/* Beta factor for the power heuristic */
		m_beta = props.getFloat("beta", 2.0f);
	}

	/// Unserialize from a binary data stream
	MIPathTracer(Stream *stream, InstanceManager *manager)
	 : MonteCarloIntegrator(stream, manager) {
		m_beta = stream->readFloat();
	}

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		LuminaireSamplingRecord lRec;
		RayDifferential ray(r);
		Intersection prevIts;
		Spectrum Li(0.0f);

		/* Perform the first ray intersection (or ignore if the 
		   intersection has already been provided). */
		rRec.rayIntersect(ray);
		Spectrum pathThroughput(rRec.attenuation);

		while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
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
				Li += pathThroughput * its.LoSub(-ray.d);

			if (m_maxDepth > 0 && rRec.depth >= m_maxDepth)
				break;

			/* ==================================================================== */
			/*                          Luminaire sampling                          */
			/* ==================================================================== */

			/* Estimate the direct illumination if this is requested */
			if (rRec.type & RadianceQueryRecord::EDirectRadiance && 
				scene->sampleLuminaire(its, lRec, rRec.nextSample2D())) {
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
			ray = Ray(its.p, its.toWorld(bRec.wo));
			bool hitLuminaire = false;
			if (scene->rayIntersect(ray, its)) {
				ray.mint = 0; ray.maxt = its.t; 
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
				} else {
					break;
				}
			}

			/* If a luminaire was hit, estimate the local illumination and
			   weight using the power heuristic */
			if (hitLuminaire &&  
				(rRec.type & RadianceQueryRecord::EDirectRadiance)) {
				lRec.Le = lRec.luminaire->Le(lRec);

				/* Prob. of having generated this sample using luminaire sampling */
				const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ? 
					scene->pdfLuminaire(prevIts, lRec) : 0;
				const Float weight = miWeight(bsdfPdf, lumPdf);
				Li += pathThroughput * lRec.Le * bsdfVal * weight;
			}

			/* ==================================================================== */
			/*                         Indirect illumination                        */
			/* ==================================================================== */

			/* Set the recursive query type */
			/* Stop if no surface was hit by the BSDF sample or if indirect illumination
			   was not requested */
			if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectRadiance)) 
				break;
			rRec.type = RadianceQueryRecord::ERadianceNoEmission;

			/* Russian roulette - Possibly stop the recursion */
			if (rRec.depth >= m_rrDepth) {
				/* Assuming that BSDF importance sampling is perfect,
				   the following should equal the maximum albedo
				   over all spectral samples */
				Float approxAlbedo = std::min((Float) 1, bsdfVal.max());
				if (rRec.nextSample1D() > approxAlbedo)
					break;
				else
					pathThroughput /= approxAlbedo;
			}

			pathThroughput *= bsdfVal;
			rRec.depth++;
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
		oss << "MIPathTracer[" << std::endl
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

MTS_IMPLEMENT_CLASS_S(MIPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(MIPathTracer, "MI path tracer");
MTS_NAMESPACE_END

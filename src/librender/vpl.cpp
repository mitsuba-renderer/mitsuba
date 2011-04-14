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

#include <mitsuba/render/vpl.h>
#include <mitsuba/core/plugin.h>

MTS_NAMESPACE_BEGIN

size_t generateVPLs(const Scene *scene, size_t offset, size_t count, int maxDepth, 
		std::deque<VPL> &vpls) {
	ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
		createObject(MTS_CLASS(Sampler), Properties("halton")));
	EmissionRecord eRec;
	eRec.type = EmissionRecord::EPreview;
	Ray ray;
	Intersection its;
	Spectrum weight, bsdfVal;
	int depth;

	if (maxDepth <= 1)
		return 0;

	const Frame stdFrame(Vector(1,0,0), Vector(0,1,0), Vector(0,0,1));

	while (vpls.size() < count) {
		sampler->setSampleIndex(++offset);
		Point2 areaSample = sampler->next2D(),
		       dirSample  = sampler->next2D();

		/* Sample an emitted particle */
		scene->sampleEmissionArea(eRec, areaSample);
		weight = eRec.value / eRec.pdfArea;
		VPL lumVPL(ELuminaireVPL, weight);
		lumVPL.its.p = eRec.sRec.p;
		lumVPL.its.shFrame = (eRec.luminaire->getType() & Luminaire::EOnSurface)
			? Frame(eRec.sRec.n) : stdFrame;
		lumVPL.luminaire = eRec.luminaire;
		vpls.push_back(lumVPL);

		weight *= eRec.luminaire->sampleEmissionDirection(eRec, dirSample);
		Float cosTheta = (eRec.luminaire->getType() & Luminaire::EOnSurface)
			? absDot(eRec.sRec.n, eRec.d) : 1;
		weight *= cosTheta / eRec.pdfDir;
		ray = Ray(eRec.sRec.p, eRec.d, 0.0f);

		depth = 2;
		while (!weight.isZero() && depth < maxDepth) {
			if (!scene->rayIntersect(ray, its))
				break;

			const BSDF *bsdf = its.shape->getBSDF();
			if (!bsdf) {
				/* Pass right through the surface (there is no BSDF) */
				ray.setOrigin(its.p);
				continue;
			}

			BSDFQueryRecord bRec(its);
			bRec.quantity = EImportance;
			bsdfVal = bsdf->sampleCos(bRec, sampler->next2D());
			if (bsdfVal.isZero())
				break;

			/* Assuming that BSDF importance sampling is perfect,
				the following should equal the maximum albedo
				over all spectral samples */
			Float approxAlbedo = std::min((Float) 0.9f, bsdfVal.max());
			if (sampler->next1D() > approxAlbedo)
				break;
			else
				weight /= approxAlbedo;

			VPL vpl(ESurfaceVPL, weight);
			vpl.its = its;
			vpls.push_back(vpl);
	
			weight *= bsdfVal;
		
			Vector wi = -ray.d, wo = its.toWorld(bRec.wo);
			ray = Ray(its.p, wo, 0.0f);

			/* Prevent light leaks due to the use of shading normals -- [Veach, p. 158] */
			Float wiDotGeoN = dot(its.geoFrame.n, wi),
				  woDotGeoN = dot(its.geoFrame.n, wo);
			if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 || 
				woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
				break;

			/* Adjoint BSDF for shading normals -- [Veach, p. 155] */
			weight *= std::abs(
				(Frame::cosTheta(bRec.wi) * woDotGeoN)/
				(Frame::cosTheta(bRec.wo) * wiDotGeoN));

			++depth;
		}
	}
	return offset;
}

const char *toString(EVPLType type) {
	switch (type) {
		case ELuminaireVPL: return "luminaireVPL";
		case ESurfaceVPL: return "surfaceVPL";
		default:
			SLog(EError, "Unknown VPL type!");
			return NULL;
	}
}

std::string VPL::toString() const {
	std::ostringstream oss;
	oss << "VPL[" << endl
		<< "  type = " << mitsuba::toString(type) << "," << endl
		<< "  P = " << P.toString() << "," << endl;
	if (type == ELuminaireVPL) {
		oss << "  p = " << its.p.toString() << "," << endl;
		oss << "  luminaire = " << indent(luminaire->toString()) << endl;
	} else {
		oss << "  its = " << indent(its.toString()) << endl;
	}
	oss << "]" << endl;
	return oss.str();
}

MTS_NAMESPACE_END

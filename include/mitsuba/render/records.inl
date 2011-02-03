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

#if !defined(__RECORDS_INLINE_H)
#define __RECORDS_INLINE_H

MTS_NAMESPACE_BEGIN
	
inline BSDFQueryRecord::BSDFQueryRecord(RadianceQueryRecord &rRec, 
	const Intersection &its): rRec(&rRec), 
	its(its), wi(its.wi), quantity(ERadiance),
	typeMask(0xFFFFFFFF), sampledType(0), component(-1), sampledComponent(-1) {
}

inline BSDFQueryRecord::BSDFQueryRecord(const Intersection &its)
	: rRec(NULL), its(its), wi(its.wi), quantity(ERadiance),
	typeMask(0xFFFFFFFF), sampledType(0), component(-1), sampledComponent(-1) {
}

inline BSDFQueryRecord::BSDFQueryRecord(RadianceQueryRecord &rRec, 
		const Intersection &its, const Vector &wo)
	: rRec(&rRec), its(its), wi(its.wi), wo(wo), quantity(ERadiance),
	typeMask(0xFFFFFFFF), sampledType(0), component(-1), sampledComponent(-1) {
}
	
inline BSDFQueryRecord::BSDFQueryRecord(const Intersection &its, const Vector &wo)	  
	: rRec(NULL), its(its), wi(its.wi), wo(wo), quantity(ERadiance),
    typeMask(0xFFFFFFFF), sampledType(0), component(-1), sampledComponent(-1) {
}
	
inline BSDFQueryRecord::BSDFQueryRecord(const Intersection &its, const Vector &wi, const Vector &wo) 
  : rRec(NULL), its(its), wi(wi), wo(wo), quantity(ERadiance),
  typeMask(0xFFFFFFFF), sampledType(0), component(-1), sampledComponent(-1) {
}

inline bool Intersection::hasSubsurface() const {
	return shape->hasSubsurface();
}
inline bool Intersection::isLuminaire() const {
	return shape->isLuminaire();
}

inline Spectrum Intersection::Le(const Vector &d) const {
	return shape->getLuminaire()->Le(
		LuminaireSamplingRecord(*this, d));
}

inline Spectrum Intersection::LoSub(const Scene *scene, const Vector &d) const {
	return shape->getSubsurface()->Lo(scene, *this, d);
}
	
inline const BSDF *Intersection::getBSDF(const RayDifferential &ray) {
	const BSDF *bsdf = shape->getBSDF();
	if (bsdf->usesRayDifferentials() && !hasUVPartials)
			computePartials(ray);
	return bsdf;
}

inline LuminaireSamplingRecord::LuminaireSamplingRecord(const Intersection &its, const Vector &dir) {
	sRec.p = its.p;
	sRec.n = its.geoFrame.n;
	d = dir;
	luminaire = its.shape->getLuminaire();
}
	
inline bool RadianceQueryRecord::rayIntersect(const RayDifferential &ray) {
	/* Only search for an intersection if this was explicitly requested */
	if (type & EIntersection) {
		scene->rayIntersect(ray, its);
		if (type & EOpacity)
			alpha = its.isValid() ? 1 : (1 - scene->getAttenuation(
				Ray(ray.o, ray.d, 0, its.t, ray.time)).average());
		if (type & EDistance)
			dist = its.t;
		type ^= EIntersection; // unset the intersection bit
	}
	return its.isValid();
}

inline Point2 RadianceQueryRecord::nextSample2D() {
	return sampler->next2D();
}

inline Float RadianceQueryRecord::nextSample1D() {
	return sampler->next1D();
}

MTS_NAMESPACE_END

#endif /* __RECORDS_INLINE_H */

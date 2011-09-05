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

#if !defined(__RECORDS_INLINE_H)
#define __RECORDS_INLINE_H

MTS_NAMESPACE_BEGIN
	
inline BSDFQueryRecord::BSDFQueryRecord(const Intersection &its, Sampler *sampler, ETransportMode quantity)
	: its(its), sampler(sampler), wi(its.wi), quantity(quantity),
	typeMask(BSDF::EAll), component(-1), sampledType(0), sampledComponent(-1) {
}

inline BSDFQueryRecord::BSDFQueryRecord(const Intersection &its, const Vector &wo, ETransportMode quantity)
	: its(its), sampler(NULL), wi(its.wi), wo(wo), quantity(quantity),
    typeMask(BSDF::EAll), component(-1), sampledType(0), sampledComponent(-1) {
}
	
inline BSDFQueryRecord::BSDFQueryRecord(const Intersection &its, const Vector &wi, const Vector &wo, ETransportMode quantity) 
  : its(its), sampler(NULL), wi(wi), wo(wo), quantity(quantity),
  typeMask(BSDF::EAll), component(-1), sampledType(0), sampledComponent(-1) {
}

void BSDFQueryRecord::reverse() {
	std::swap(wo, wi);
	quantity = (ETransportMode) (1-quantity);
}

inline bool Intersection::hasSubsurface() const {
	return shape->hasSubsurface();
}

inline bool Intersection::isLuminaire() const {
	return shape->isLuminaire();
}

inline Spectrum Intersection::Le(const Vector &d) const {
	return shape->getLuminaire()->Le(ShapeSamplingRecord(*this), d);
}

inline Spectrum Intersection::LoSub(const Scene *scene, 
		Sampler *sampler, const Vector &d, int depth) const {
	return shape->getSubsurface()->Lo(scene, sampler, *this, d, depth);
}

inline const BSDF *Intersection::getBSDF(const RayDifferential &ray) {
	const BSDF *bsdf = shape->getBSDF();

	if (bsdf && bsdf->usesRayDifferentials() && !hasUVPartials)
			computePartials(ray);
	return bsdf;
}

inline bool Intersection::isMediumTransition() const {
	return shape->isMediumTransition();
}

inline const Medium *Intersection::getTargetMedium(const Vector &d) const {
	if (dot(d, geoFrame.n) > 0)
		return shape->getExteriorMedium();
	else
		return shape->getInteriorMedium();
}
	
inline const Medium *Intersection::getTargetMedium(Float cosTheta) const {
	if (cosTheta > 0)
		return shape->getExteriorMedium();
	else
		return shape->getInteriorMedium();
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
		if (type & EOpacity) {
			if (its.isValid())
				alpha = 1.0f;
			else if (medium == NULL)
				alpha = 0.0f;
			else
				alpha = 1-medium->getTransmittance(ray).average();
		}
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

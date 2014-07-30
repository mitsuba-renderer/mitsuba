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

#if !defined(__RECORDS_INLINE_H)
#define __RECORDS_INLINE_H

MTS_NAMESPACE_BEGIN

inline BSDFSamplingRecord::BSDFSamplingRecord(const Intersection &its, Sampler *sampler, ETransportMode mode)
	: its(its), sampler(sampler), wi(its.wi), mode(mode),
	typeMask(BSDF::EAll), component(-1), sampledType(0), sampledComponent(-1) {
}

inline BSDFSamplingRecord::BSDFSamplingRecord(const Intersection &its, const Vector &wo, ETransportMode mode)
	: its(its), sampler(NULL), wi(its.wi), wo(wo), mode(mode),
    typeMask(BSDF::EAll), component(-1), sampledType(0), sampledComponent(-1) {
}

inline BSDFSamplingRecord::BSDFSamplingRecord(const Intersection &its, const Vector &wi, const Vector &wo, ETransportMode mode)
  : its(its), sampler(NULL), wi(wi), wo(wo), mode(mode),
  typeMask(BSDF::EAll), component(-1), sampledType(0), sampledComponent(-1) {
}

void BSDFSamplingRecord::reverse() {
	std::swap(wo, wi);
	mode = (ETransportMode) (1-mode);
}

inline bool Intersection::hasSubsurface() const {
	return shape->hasSubsurface();
}

inline bool Intersection::isEmitter() const {
	return shape->isEmitter();
}

inline bool Intersection::isSensor() const {
	return shape->isSensor();
}

inline Spectrum Intersection::Le(const Vector &d) const {
	return shape->getEmitter()->eval(*this, d);
}

inline Spectrum Intersection::LoSub(const Scene *scene,
		Sampler *sampler, const Vector &d, int depth) const {
	return shape->getSubsurface()->Lo(scene, sampler, *this, d, depth);
}

inline const BSDF *Intersection::getBSDF() const {
	return shape->getBSDF();
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

void Intersection::adjustTime(Float time) {
	if (instance)
		instance->adjustTime(*this, time);
	else if (shape)
		shape->adjustTime(*this, time);
	else
		this->time = time;
}

void Intersection::getNormalDerivative(Vector &dndu, Vector &dndv,
	bool shadingFrame) const {

	if (instance)
		instance->getNormalDerivative(*this, dndu, dndv, shadingFrame);
	else if (shape)
		shape->getNormalDerivative(*this, dndu, dndv, shadingFrame);
}

inline const PhaseFunction *MediumSamplingRecord::getPhaseFunction() const {
	return medium->getPhaseFunction();
}

inline bool RadianceQueryRecord::rayIntersect(const RayDifferential &ray) {
	/* Only search for an intersection if this was explicitly requested */
	if (type & EIntersection) {
		scene->rayIntersect(ray, its);
		if (type & EOpacity) {
			int unused = INT_MAX;

			if (its.isValid()) {
				if (EXPECT_TAKEN(!its.isMediumTransition()))
					alpha = 1.0f;
				else
					alpha = 1-scene->evalTransmittance(its.p, true,
						ray(scene->getBSphere().radius*2), false,
						ray.time, its.getTargetMedium(ray.d), unused).average();
			} else if (medium) {
				alpha = 1-scene->evalTransmittance(ray.o, false,
					ray(scene->getBSphere().radius*2), false,
					ray.time, medium, unused).average();
			} else {
				alpha = 0.0f;
			}
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

inline PositionSamplingRecord::PositionSamplingRecord(const Intersection &its, EMeasure measure)
	: p(its.p), time(its.time), n(its.shFrame.n), measure(measure), uv(its.uv), object(NULL) { }

inline DirectionSamplingRecord::DirectionSamplingRecord(const Intersection &its, EMeasure measure)
	: d(its.toWorld(its.wi)), measure(measure) { }

inline DirectSamplingRecord::DirectSamplingRecord(const Intersection &refIts)
	: PositionSamplingRecord(refIts.time), ref(refIts.p), refN(0.0f) {
	if ((refIts.shape->getBSDF()->getType() & (BSDF::ETransmission | BSDF::EBackSide)) == 0)
		refN = refIts.shFrame.n;
}

inline DirectSamplingRecord::DirectSamplingRecord(const MediumSamplingRecord &refM)
	: PositionSamplingRecord(refM.time), ref(refM.p), refN(0.0f) {
}

void DirectSamplingRecord::setQuery(const Ray &ray, const Intersection &its, EMeasure _measure) {
	p = its.p;
	n = its.shFrame.n;
	measure = _measure;
	uv = its.uv;
	object = its.shape->getEmitter();
	d = ray.d;
	dist = its.t;
}

MTS_NAMESPACE_END

#endif /* __RECORDS_INLINE_H */

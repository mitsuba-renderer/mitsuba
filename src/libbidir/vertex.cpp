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

#include <mitsuba/bidir/path.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

static StatsCounter mediumInconsistencies("Bidirectional layer",
        "Medium inconsistencies in sampleNext()");

void PathVertex::makeEndpoint(const Scene *scene, Float time, ETransportMode mode) {
    memset(this, 0, sizeof(PathVertex));
    type = (mode == EImportance) ? EEmitterSupernode : ESensorSupernode;
    getEndpointRecord() = EndpointRecord(time);
    degenerate = (mode == EImportance)
        ? scene->hasDegenerateEmitters() : scene->hasDegenerateSensor();
}

bool PathVertex::sampleNext(const Scene *scene, Sampler *sampler,
        const PathVertex *pred, const PathEdge *predEdge,
        PathEdge *succEdge, PathVertex *succ,
        ETransportMode mode, bool russianRoulette, Spectrum *throughput) {
    Ray ray;

    memset(succEdge, 0, sizeof(PathEdge));
    memset(succ, 0, sizeof(PathVertex));

    succEdge->medium = (predEdge == NULL) ? NULL : predEdge->medium;
    rrWeight = 1.0f;

    switch (type) {
        case EEmitterSupernode: {
                BDAssert(mode == EImportance && pred == NULL && predEdge == NULL);
                PositionSamplingRecord &pRec = succ->getPositionSamplingRecord();
                const EndpointRecord &eRec = getEndpointRecord();
                pRec = PositionSamplingRecord(eRec.time);
                Spectrum result = scene->sampleEmitterPosition(pRec, sampler->next2D());
                if (result.isZero())
                    return false;

                const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
                weight[EImportance] = result;
                pdf[EImportance] = pRec.pdf;
                measure = pRec.measure;
                succ->type = EEmitterSample;
                succ->degenerate = emitter->getType() & Emitter::EDeltaDirection;

                succEdge->weight[EImportance] = Spectrum(1.0f);
                succEdge->pdf[EImportance] = 1.0f;
                succEdge->medium = emitter->getMedium();

                return true;
            }
            break;

        case ESensorSupernode: {
                BDAssert(mode == ERadiance && pred == NULL && predEdge == NULL);
                PositionSamplingRecord &pRec = succ->getPositionSamplingRecord();
                const EndpointRecord &eRec = getEndpointRecord();
                pRec = PositionSamplingRecord(eRec.time);
                Spectrum result = scene->sampleSensorPosition(pRec, sampler->next2D());
                if (result.isZero())
                    return false;

                const Sensor *sensor = static_cast<const Sensor *>(pRec.object);
                weight[ERadiance] = result;
                pdf[ERadiance] = pRec.pdf;
                measure = pRec.measure;
                succ->type = ESensorSample;
                succ->degenerate = sensor->getType()
                    & Sensor::EDeltaDirection;

                succEdge->weight[ERadiance] = Spectrum(1.0f);
                succEdge->pdf[ERadiance] = 1.0f;
                succEdge->medium = sensor->getMedium();

                return true;
            }
            break;

        case EEmitterSample: {
                BDAssert(mode == EImportance && pred->type == EEmitterSupernode);
                PositionSamplingRecord &pRec = getPositionSamplingRecord();
                const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
                DirectionSamplingRecord dRec;

                Spectrum result = emitter->sampleDirection(dRec, pRec,
                    emitter->needsDirectionSample() ? sampler->next2D() : Point2(0.5f));

                if (result.isZero())
                    return false;

                weight[EImportance] = result;
                weight[ERadiance]   = result * dRec.pdf * (
                    emitter->isOnSurface() ? 1.0f / absDot(dRec.d, pRec.n) : 1.0f);
                pdf[EImportance]    = dRec.pdf;
                pdf[ERadiance]      = 1.0f;

                measure = dRec.measure;
                succEdge->medium = emitter->getMedium();
                ray.time = pRec.time;
                ray.setOrigin(pRec.p);
                ray.setDirection(dRec.d);
            }
            break;

        case ESensorSample: {
                BDAssert(mode == ERadiance && pred->type == ESensorSupernode);
                PositionSamplingRecord &pRec = getPositionSamplingRecord();
                const Sensor *sensor = static_cast<const Sensor *>(pRec.object);
                DirectionSamplingRecord dRec;

                Spectrum result = sensor->sampleDirection(dRec, pRec,
                    sensor->needsDirectionSample() ? sampler->next2D() : Point2(0.5f));

                if (result.isZero())
                    return false;

                weight[EImportance] = result * dRec.pdf * (
                    sensor->isOnSurface() ? 1.0f / absDot(dRec.d, pRec.n) : 1.0f);
                weight[ERadiance]   = result;
                pdf[EImportance]    = 1.0f;
                pdf[ERadiance]      = dRec.pdf;

                measure = dRec.measure;
                succEdge->medium = sensor->getMedium();
                ray.time = pRec.time;
                ray.setOrigin(pRec.p);
                ray.setDirection(dRec.d);
            }
            break;

        case ESurfaceInteraction: {
                const Intersection &its = getIntersection();
                const BSDF *bsdf = its.getBSDF();
                Vector wi = normalize(pred->getPosition() - its.p);
                Vector wo;

                /* Sample the BSDF */
                BSDFSamplingRecord bRec(its, sampler, mode);
                bRec.wi = its.toLocal(wi);
                weight[mode] = bsdf->sample(bRec, pdf[mode], sampler->next2D());
                if (weight[mode].isZero())
                    return false;

                measure = BSDF::getMeasure(bRec.sampledType);
                componentType = (uint16_t) (bRec.sampledType & BSDF::EAll);

                wo = its.toWorld(bRec.wo);

                /* Prevent light leaks due to the use of shading normals */
                Float wiDotGeoN = dot(its.geoFrame.n, wi),
                      woDotGeoN = dot(its.geoFrame.n, wo);
                if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
                    woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                    return false;

                /* Account for medium changes if applicable */
                if (its.isMediumTransition()) {
                    const Medium *expected = its.getTargetMedium(wi);
                    if (expected != predEdge->medium) {
                        #if defined(MTS_BD_TRACE)
                            SLog(EWarn, "Detected an inconsistency: approached "
                                "surface %s within medium %s, but the surface "
                                "states that the ray should have been in medium %s.",
                                its.toString().c_str(), predEdge->medium ?
                                predEdge->medium->toString().c_str() : "null",
                                expected ? expected->toString().c_str() : "null");
                        #endif
                        ++mediumInconsistencies;
                        return false;
                    }
                    succEdge->medium = its.getTargetMedium(wo);
                }

                /* Compute the reverse quantities */
                bRec.reverse();
                pdf[1-mode] = bsdf->pdf(bRec, (EMeasure) measure);
                if (pdf[1-mode] <= RCPOVERFLOW) {
                    /* This can happen rarely due to roundoff errors -- be strict */
                    return false;
                }

                if (!(bsdf->getType() & BSDF::ENonSymmetric)) {
                    /* Make use of symmetry -- no need to re-evaluate
                       everything (only the pdf and cosine factors changed) */
                    weight[1-mode] = weight[mode] * (pdf[mode] / pdf[1-mode]);
                    if (measure == ESolidAngle)
                        weight[1-mode] *=
                            std::abs(Frame::cosTheta(bRec.wo) / Frame::cosTheta(bRec.wi));
                } else {
                    weight[1-mode] = bsdf->eval(bRec, (EMeasure) measure) / pdf[1-mode];
                }
                bRec.reverse();

                /* Adjoint BSDF for shading normals */
                if (mode == EImportance)
                    weight[EImportance] *= std::abs(
                        (Frame::cosTheta(bRec.wi) * woDotGeoN) /
                        (Frame::cosTheta(bRec.wo) * wiDotGeoN));
                else
                    weight[EImportance] *= std::abs(
                        (Frame::cosTheta(bRec.wo) * wiDotGeoN) /
                        (Frame::cosTheta(bRec.wi) * woDotGeoN));

                /// For BDPT & russian roulette, track radiance * eta^2
                if (throughput && mode == ERadiance && bRec.eta != 1)
                    (*throughput) *= bRec.eta * bRec.eta;

                ray.time = its.time;
                ray.setOrigin(its.p);
                ray.setDirection(wo);
            }
            break;

        case EMediumInteraction: {
                const MediumSamplingRecord &mRec = getMediumSamplingRecord();
                const PhaseFunction *phase = succEdge->medium->getPhaseFunction();
                Vector wi = normalize(pred->getPosition() - mRec.p);
                PhaseFunctionSamplingRecord pRec(mRec, wi, mode);

                weight[mode] = mRec.sigmaS * phase->sample(pRec, pdf[mode], sampler);
                if (weight[mode].isZero())
                    return false;

                ray.time = mRec.time;
                ray.mint = 0;
                ray.setOrigin(mRec.p);
                ray.setDirection(pRec.wo);
                measure = ESolidAngle;

                if (!(phase->getType() & PhaseFunction::ENonSymmetric)) {
                    /* Make use of symmetry -- no need to re-evaluate */
                    pdf[1-mode] = pdf[mode];
                    weight[1-mode] = weight[mode];
                } else {
                    pRec.reverse();
                    pdf[1-mode] = phase->pdf(pRec);
                    weight[1-mode] = mRec.sigmaS * (phase->eval(pRec) / pdf[1-mode]);
                }
            }
            break;

        default:
            SLog(EError, "PathVertex::sampleNext(): Encountered an "
                "unsupported vertex type (%i)!", type);
            return false;
    }

    if (throughput) {
        /* For BDPT: keep track of the path throughput to run russian roulette */
        (*throughput) *= weight[mode];

        if (russianRoulette) {
            Float q = std::min(throughput->max(), (Float) 0.95f);

            if (sampler->next1D() > q) {
                measure = EInvalidMeasure;
                return false;
            } else {
                rrWeight = 1.0f / q;
                (*throughput) *= rrWeight;
            }
        }
    }

    if (!succEdge->sampleNext(scene, sampler, this, ray, succ, mode)) {
        /* Sampling a successor edge + vertex failed, hence the vertex
           is not committed to a particular measure yet -- revert. */
        measure = EInvalidMeasure;
        return false;
    } else {
        if (throughput)
            (*throughput) *= succEdge->weight[mode];
    }

    /* Convert from solid angle to area measure */
    if (measure == ESolidAngle) {
        measure = EArea;

        pdf[mode] /= succEdge->length * succEdge->length;
        if (succ->isOnSurface())
            pdf[mode] *= absDot(ray.d, succ->getGeometricNormal());

        if (predEdge->length != 0.0f) {
            pdf[1-mode] /= predEdge->length * predEdge->length;
            if (pred->isOnSurface())
                pdf[1-mode] *= absDot(predEdge->d, pred->getGeometricNormal());
        }
    }

    return true;
}

int PathVertex::sampleSensor(const Scene *scene, Sampler *sampler,
        const Point2i &pixelPosition_, PathEdge *e0, PathVertex *v1,
        PathEdge *e1, PathVertex *v2) {
    BDAssert(type == ESensorSupernode);
    const EndpointRecord &eRec = getEndpointRecord();
    const Sensor *sensor = scene->getSensor();
    Point2 pixelPosition(pixelPosition_);

    memset(e0, 0, sizeof(PathEdge));
    memset(v1, 0, sizeof(PathVertex));

    Point2 pixelSample = sampler->next2D(),
           apertureSample = sensor->needsApertureSample() ? sampler->next2D() : Point2(0.5f);

    PositionSamplingRecord &pRec = v1->getPositionSamplingRecord();
    pRec = PositionSamplingRecord(eRec.time);
    Spectrum result = scene->sampleSensorPosition(pRec,
        (sensor->getType() & Sensor::EPositionSampleMapsToPixels) ? pixelSample
        : apertureSample, &pixelPosition);

    if (result.isZero())
        return 0;

    weight[ERadiance] = result;
    pdf[ERadiance] = pRec.pdf;
    measure = pRec.measure;
    rrWeight = 1.0f;
    v1->type = ESensorSample;
    v1->degenerate = sensor->getType() & Sensor::EDeltaDirection;

    e0->weight[ERadiance] = Spectrum(1.0f);
    e0->pdf[ERadiance] = 1.0f;
    e0->medium = sensor->getMedium();

    DirectionSamplingRecord dRec;
    result = sensor->sampleDirection(dRec, pRec,
        (sensor->getType() & Sensor::EPositionSampleMapsToPixels) ? apertureSample
        : pixelSample, &pixelPosition);
    if (result.isZero())
        return 1;

    memset(e1, 0, sizeof(PathEdge));
    memset(v2, 0, sizeof(PathVertex));

    v1->weight[EImportance] = result * dRec.pdf * (
        sensor->isOnSurface() ? 1.0f / absDot(dRec.d, pRec.n) : 1.0f);
    v1->weight[ERadiance]   = result;
    v1->pdf[EImportance]    = 1.0f;
    v1->pdf[ERadiance]      = dRec.pdf;
    v1->rrWeight = 1.0f;

    v1->measure = dRec.measure;
    e1->medium = sensor->getMedium();

    Ray ray;
    ray.time = pRec.time;
    ray.setOrigin(pRec.p);
    ray.setDirection(dRec.d);

    if (!e1->sampleNext(scene, sampler, v1, ray, v2, ERadiance)) {
        v1->measure = EInvalidMeasure;
        return 1;
    }

    if (v1->measure == ESolidAngle) {
        v1->measure = EArea;
        v1->pdf[ERadiance] /= e1->length * e1->length;
        if (v2->isOnSurface())
            v1->pdf[ERadiance] *= absDot(ray.d, v2->getGeometricNormal());
    }

    return 2;
}

bool PathVertex::perturbPosition(const Scene *scene, Sampler *sampler, Float stddev) {
    Point2 step = warp::squareToStdNormal(sampler->next2D()) * stddev;
    EVertexType type = (EVertexType) this->type;
    Ray ray;

    switch (type) {
        case ESurfaceInteraction: {
                const Intersection &its = getIntersection();

                ray = Ray(its.p + its.geoFrame.s * step.x + its.geoFrame.t * step.y
                        + its.geoFrame.n * Epsilon, -its.geoFrame.n, 0,
                        std::numeric_limits<Float>::infinity(), its.time);
            }
            break;
        case ESensorSample:
        case EEmitterSample: {
                const PositionSamplingRecord &pRec = getPositionSamplingRecord();
                if (pRec.n.isZero())
                    return false;
                Frame frame(pRec.n);
                ray = Ray(pRec.p + frame.s * step.x + frame.t * step.y
                        + frame.n * Epsilon, -frame.n, 0,
                        std::numeric_limits<Float>::infinity(), pRec.time);
            }
            break;
        default:
            SLog(EError, "PathVertex::perturbPosition(): Encountered an "
                "unsupported vertex type (%i)!", type);
            return false;
    }
    Intersection &its = getIntersection();
    Intersection its1, its2;
    scene->rayIntersectAll(ray, its1); ray.d = -ray.d;
    scene->rayIntersectAll(ray, its2);

    bool its1Valid = its1.isValid(), its2Valid = its2.isValid();
    if (its1Valid) {
        if (type == ESurfaceInteraction)
            its1Valid &= its1.getBSDF() == its.getBSDF();
        else if (type == EEmitterSample)
            its1Valid &= its1.isEmitter();
        else if (type == ESensorSample)
            its1Valid &= its1.isSensor();
    }

    if (its2Valid) {
        if (type == ESurfaceInteraction)
            its2Valid &= its2.getBSDF() == its.getBSDF();
        else if (type == EEmitterSample)
            its2Valid &= its2.isEmitter();
        else if (type == ESensorSample)
            its2Valid &= its2.isSensor();
    }

    if (its1Valid && its2Valid) {
        if (its1.t < its2.t)
            its = its1;
        else
            its = its2;
    } else if (its1Valid) {
        its = its1;
    } else if (its2Valid) {
        its = its2;
    } else {
        return false;
    }

    this->type = ESurfaceInteraction;
    return cast(scene, type);
}

Float PathVertex::perturbPositionPdf(const PathVertex *target, Float stddev) const {
    BDAssert(type == target->type);
    switch (type) {
        case ESurfaceInteraction: {
                const Intersection &itsOld = getIntersection();
                const Intersection &itsNew = target->getIntersection();
                Vector rel = itsOld.geoFrame.toLocal(itsOld.p - itsNew.p);
                Point2 rel2 = Point2(rel.x, rel.y) / stddev;

                return warp::squareToStdNormalPdf(rel2) * absDot(itsOld.geoFrame.n, itsNew.geoFrame.n) / (stddev*stddev);
            }
            break;

        case ESensorSample:
        case EEmitterSample: {
                const PositionSamplingRecord &prOld = getPositionSamplingRecord();
                const PositionSamplingRecord &prNew = target->getPositionSamplingRecord();
                Vector rel = Frame(prOld.n).toLocal(prOld.p - prNew.p);
                Point2 rel2 = Point2(rel.x, rel.y) / stddev;

                return warp::squareToStdNormalPdf(rel2) * absDot(prOld.n, prNew.n) / (stddev*stddev);
            }
            break;

        default:
            SLog(EError, "PathVertex::perturbPositionPdf(): Encountered an "
                "unsupported vertex type (%i)!", type);
            return 0.0f;
    }
}

bool PathVertex::perturbDirection(const Scene *scene, const PathVertex *pred,
    const PathEdge *predEdge, PathEdge *succEdge, PathVertex *succ,
    const Vector &d, Float dist, EVertexType desiredType, ETransportMode mode) {
    Ray ray(getPosition(), d, pred->getTime());

    memset(succEdge, 0, sizeof(PathEdge));
    memset(succ, 0, sizeof(PathVertex));

    succ->measure = EInvalidMeasure;
    succEdge->medium = (predEdge == NULL) ? NULL : predEdge->medium;
    BDAssert(!isSupernode());
    if (isDegenerate())
        return false;

    switch (type) {
        case EEmitterSample: {
                BDAssert(mode == EImportance && pred->type == EEmitterSupernode);
                PositionSamplingRecord &pRec = getPositionSamplingRecord();
                const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
                DirectionSamplingRecord dRec(d);

                Spectrum value = emitter->evalDirection(dRec, pRec);
                Float prob = emitter->pdfDirection(dRec, pRec);

                if (value.isZero() || prob <= RCPOVERFLOW)
                    return false;

                weight[EImportance] = value/prob;
                weight[ERadiance]   = value * (
                    emitter->isOnSurface() ? 1.0f / absDot(dRec.d, pRec.n) : 1.0f);
                pdf[EImportance]    = prob;
                pdf[ERadiance]      = 1.0f;

                measure = dRec.measure;
                succEdge->medium = emitter->getMedium();
            }
            break;

        case ESensorSample: {
                BDAssert(mode == ERadiance && pred->type == ESensorSupernode);
                PositionSamplingRecord &pRec = getPositionSamplingRecord();
                const Sensor *sensor = static_cast<const Sensor *>(pRec.object);
                DirectionSamplingRecord dRec(d);

                Spectrum value = sensor->evalDirection(dRec, pRec);
                Float prob = sensor->pdfDirection(dRec, pRec);

                if (value.isZero() || prob <= RCPOVERFLOW)
                    return false;

                weight[EImportance] = value * (
                    sensor->isOnSurface() ? 1.0f / absDot(dRec.d, pRec.n) : 1.0f);
                weight[ERadiance]   = value / prob;
                pdf[EImportance]    = 1.0f;
                pdf[ERadiance]      = prob;

                measure = dRec.measure;
                succEdge->medium = sensor->getMedium();
            }
            break;

        case ESurfaceInteraction: {
                const Intersection &its = getIntersection();
                const BSDF *bsdf = its.getBSDF();
                Vector wi = normalize(pred->getPosition() - its.p);
                Vector wo(d);

                BSDFSamplingRecord bRec(its, its.toLocal(wi), its.toLocal(wo), mode);

                Spectrum value = bsdf->eval(bRec);
                Float prob = bsdf->pdf(bRec);

                if (value.isZero() || prob <= RCPOVERFLOW)
                    return false;

                weight[mode] = value/prob;
                pdf[mode] = prob;

                /* Prevent light leaks due to the use of shading normals */
                Float wiDotGeoN = dot(its.geoFrame.n, wi),
                      woDotGeoN = dot(its.geoFrame.n, wo);
                if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
                    woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                    return false;

                /* Account for medium changes if applicable */
                if (its.isMediumTransition()) {
                    const Medium *expected = its.getTargetMedium(wi);
                    if (expected != predEdge->medium) {
                        #if defined(MTS_BD_TRACE)
                            SLog(EWarn, "Detected an inconsistency: approached "
                                "surface %s within medium %s, but the surface "
                                "states that the ray should be in medium %s.",
                                its.toString().c_str(), predEdge->medium ?
                                predEdge->medium->toString().c_str() : "null",
                                expected ? expected->toString().c_str() : "null");
                        #endif
                        ++mediumInconsistencies;
                        return false;
                    }
                    succEdge->medium = its.getTargetMedium(wo);
                }
                measure = ESolidAngle;
                componentType = BSDF::ESmooth;

                /* Compute the reverse quantities */
                bRec.reverse();
                pdf[1-mode] = bsdf->pdf(bRec, ESolidAngle);
                if (pdf[1-mode] <= RCPOVERFLOW) {
                    /* This can happen rarely due to roundoff errors -- be strict */
                    return false;
                }
                if (!(bsdf->getType() & BSDF::ENonSymmetric)) {
                    /* Make use of symmetry -- no need to re-evaluate
                       everything (only the pdf and cosine factors changed) */
                    weight[1-mode] = weight[mode] * std::abs(
                        (pdf[mode] * Frame::cosTheta(bRec.wo)) /
                        (pdf[1-mode] * Frame::cosTheta(bRec.wi)));
                } else {
                    weight[1-mode] = bsdf->eval(bRec, ESolidAngle) / pdf[1-mode];
                }
                bRec.reverse();

                /* Adjoint BSDF for shading normals */
                if (mode == EImportance)
                    weight[EImportance] *= std::abs(
                        (Frame::cosTheta(bRec.wi) * woDotGeoN) /
                        (Frame::cosTheta(bRec.wo) * wiDotGeoN));
                else
                    weight[EImportance] *= std::abs(
                        (Frame::cosTheta(bRec.wo) * wiDotGeoN) /
                        (Frame::cosTheta(bRec.wi) * woDotGeoN));
            }
            break;

        case EMediumInteraction: {
                const MediumSamplingRecord &mRec = getMediumSamplingRecord();
                const PhaseFunction *phase = succEdge->medium->getPhaseFunction();
                Vector wi = normalize(pred->getPosition() - mRec.p);
                PhaseFunctionSamplingRecord pRec(mRec, wi, d, mode);

                Float value = phase->eval(pRec);
                Float prob = phase->pdf(pRec);

                if (value == 0 || prob <= RCPOVERFLOW)
                    return false;

                pdf[mode] = prob;
                weight[mode] = mRec.sigmaS * (value/prob);

                measure = ESolidAngle;

                if (!(phase->getType() & PhaseFunction::ENonSymmetric)) {
                    /* Make use of symmetry -- no need to re-evaluate */
                    pdf[1-mode] = pdf[mode];
                    weight[1-mode] = weight[mode];
                } else {
                    pRec.reverse();
                    pdf[1-mode] = phase->pdf(pRec);
                    weight[1-mode] = mRec.sigmaS * (phase->eval(pRec) / pdf[1-mode]);
                }
            }
            break;

        default:
            SLog(EError, "PathVertex::perturbDirection(): Encountered an "
                "unsupported vertex type (%i)!", type);
            return false;
    }

    if (!succEdge->perturbDirection(scene, this, ray, dist, desiredType, succ, mode)) {
        measure = EInvalidMeasure;
        return false;
    }

    /* Convert from solid angle to area measure */
    if (measure == ESolidAngle) {
        measure = EArea;

        pdf[mode] /= succEdge->length * succEdge->length;
        if (succ->isOnSurface())
            pdf[mode] *= absDot(ray.d, succ->getGeometricNormal());

        if (predEdge->length != 0.0f) {
            pdf[1-mode] /= predEdge->length * predEdge->length;
            if (pred->isOnSurface())
                pdf[1-mode] *= absDot(predEdge->d, pred->getGeometricNormal());
        }
    }

    return true;
}

bool PathVertex::propagatePerturbation(const Scene *scene, const PathVertex *pred,
        const PathEdge *predEdge, PathEdge *succEdge, PathVertex *succ,
        unsigned int componentType_, Float dist, EVertexType desiredType, ETransportMode mode) {
    BDAssert(isSurfaceInteraction());

    const Intersection &its = getIntersection();
    const BSDF *bsdf = its.getBSDF();
    if (!(bsdf->getType() & BSDF::EDelta))
        return false;

    memset(succEdge, 0, sizeof(PathEdge));
    memset(succ, 0, sizeof(PathVertex));

    Vector wi = normalize(pred->getPosition() - its.p);

    BSDFSamplingRecord bRec(its, NULL, mode);
    bRec.typeMask = componentType_;

    bRec.wi = its.toLocal(wi);
    if (bsdf->sample(bRec, Point2(0.5f)).isZero())
        return false;

    Vector wo = its.toWorld(bRec.wo);

    /* Prevent light leaks due to the use of shading normals */
    Float wiDotGeoN = dot(its.geoFrame.n, wi),
          woDotGeoN = dot(its.geoFrame.n, wo);
    if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
        woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
        return false;

    bRec.typeMask = BSDF::EAll;
    Float prob = bsdf->pdf(bRec, EDiscrete);
    if (prob <= RCPOVERFLOW) {
        SLog(EWarn, "Unable to recreate specular vertex in perturbation (bsdf=%s)",
            bsdf->toString().c_str());
        return false;
    }

    weight[mode] = bsdf->eval(bRec, EDiscrete) / prob;
    pdf[mode] = prob;
    measure = EDiscrete;
    componentType = componentType_;

    if (weight[mode].isZero() || prob <= RCPOVERFLOW)
        return false;

    /* Account for medium changes if applicable */
    if (its.isMediumTransition()) {
        const Medium *expected = its.getTargetMedium(wi);
        if (expected != predEdge->medium) {
            #if defined(MTS_BD_TRACE)
                SLog(EWarn, "Detected an inconsistency: approached "
                    "surface %s within medium %s, but the surface "
                    "states that the ray should be in medium %s.",
                    its.toString().c_str(), predEdge->medium ?
                    predEdge->medium->toString().c_str() : "null",
                    expected ? expected->toString().c_str() : "null");
            #endif
            ++mediumInconsistencies;
            return false;
        }
        succEdge->medium = its.getTargetMedium(wo);
    } else {
        succEdge->medium = predEdge->medium;
    }

    /* Compute the reverse quantities */
    bRec.reverse();
    pdf[1-mode] = bsdf->pdf(bRec, EDiscrete);
    if (pdf[1-mode] <= RCPOVERFLOW) {
        /* This can happen rarely due to roundoff errors -- be strict */
        return false;
    }

    if (!(bsdf->getType() & BSDF::ENonSymmetric))
        weight[1-mode] = weight[mode];
    else
        weight[1-mode] = bsdf->eval(bRec, EDiscrete) / pdf[1-mode];
    bRec.reverse();

    /* Adjoint BSDF for shading normals */
    if (mode == EImportance)
        weight[EImportance] *= std::abs(
            (Frame::cosTheta(bRec.wi) * woDotGeoN) /
            (Frame::cosTheta(bRec.wo) * wiDotGeoN));
    else
        weight[EImportance] *= std::abs(
            (Frame::cosTheta(bRec.wo) * wiDotGeoN) /
            (Frame::cosTheta(bRec.wi) * woDotGeoN));

    Ray ray(its.p, wo, its.time);
    if (!succEdge->perturbDirection(scene, this, ray, dist, desiredType, succ, mode)) {
        measure = EInvalidMeasure;
        return false;
    }

    return true;
}

Spectrum PathVertex::eval(const Scene *scene, const PathVertex *pred,
        const PathVertex *succ, ETransportMode mode, EMeasure measure) const {
    Spectrum result(0.0f);
    Vector wo(0.0f);

    switch (type) {
        case EEmitterSupernode: {
                if (mode != EImportance || pred != NULL || succ->type != EEmitterSample)
                    return Spectrum(0.0f);
                PositionSamplingRecord pRec = succ->getPositionSamplingRecord();
                const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
                pRec.measure = measure;
                return emitter->evalPosition(pRec);
            }
            break;

        case ESensorSupernode: {
                if (mode != ERadiance || pred != NULL || succ->type != ESensorSample)
                    return Spectrum(0.0f);
                PositionSamplingRecord pRec = succ->getPositionSamplingRecord();
                const Sensor *sensor = static_cast<const Sensor *>(pRec.object);
                pRec.measure = measure;
                return sensor->evalPosition(pRec);
            }
            break;

        case EEmitterSample: {
                Point target;
                if (mode == EImportance && pred->type == EEmitterSupernode)
                    target = succ->getPosition();
                else if (mode == ERadiance && succ->type == EEmitterSupernode)
                    target = pred->getPosition();
                else
                    return Spectrum(0.0f);

                const PositionSamplingRecord &pRec = getPositionSamplingRecord();
                const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
                wo = normalize(target - pRec.p);
                DirectionSamplingRecord dRec(wo, measure == EArea ? ESolidAngle : measure);
                result = emitter->evalDirection(dRec, pRec);
                Float dp = absDot(pRec.n, wo);
                if (measure != EDiscrete && dp != 0)
                    result /= dp;
            }
            break;

        case ESensorSample: {
                Point target;
                if (mode == ERadiance && pred->type == ESensorSupernode)
                    target = succ->getPosition();
                else if (mode == EImportance && succ->type == ESensorSupernode)
                    target = pred->getPosition();
                else
                    return Spectrum(0.0f);

                const PositionSamplingRecord &pRec = getPositionSamplingRecord();
                const Sensor *sensor = static_cast<const Sensor *>(pRec.object);
                wo = normalize(target - pRec.p);
                DirectionSamplingRecord dRec(wo, measure == EArea ? ESolidAngle : measure);
                result = sensor->evalDirection(dRec, pRec);
                Float dp = absDot(pRec.n, wo);
                if (measure != EDiscrete && dp != 0)
                    result /= dp;

                return result;
            }
            break;

        case ESurfaceInteraction: {
                const Intersection &its = getIntersection();
                const BSDF *bsdf = its.getBSDF();

                Point predP = pred->getPosition(),
                      succP = succ->getPosition();

                Vector wi = normalize(predP - its.p);
                wo = normalize(succP - its.p);

                BSDFSamplingRecord bRec(its, its.toLocal(wi),
                        its.toLocal(wo), mode);

                if (measure == EArea)
                    measure = ESolidAngle;

                result = bsdf->eval(bRec, measure);

                /* Prevent light leaks due to the use of shading normals */
                Float wiDotGeoN = dot(its.geoFrame.n, wi),
                      woDotGeoN = dot(its.geoFrame.n, wo);

                if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
                    woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                    return Spectrum(0.0f);

                if (mode == EImportance) {
                    /* Adjoint BSDF for shading normals */
                    result *= std::abs(
                        (Frame::cosTheta(bRec.wi) * woDotGeoN) /
                        (Frame::cosTheta(bRec.wo) * wiDotGeoN));
                }

                if (measure != EDiscrete && Frame::cosTheta(bRec.wo) != 0)
                    result /= std::abs(Frame::cosTheta(bRec.wo));
            }
            break;

        case EMediumInteraction: {
                if (measure != ESolidAngle && measure != EArea)
                    return Spectrum(0.0f);

                const MediumSamplingRecord &mRec = getMediumSamplingRecord();

                Point predP = pred->getPosition(),
                      succP = succ->getPosition();

                Vector wi = normalize(predP - mRec.p);
                wo = normalize(succP - mRec.p);

                const PhaseFunction *phase = mRec.medium->getPhaseFunction();
                PhaseFunctionSamplingRecord pRec(mRec, wi, wo, mode);

                result = mRec.sigmaS * phase->eval(pRec);
            }
            break;

        default:
            SLog(EError, "PathVertex::eval(): Encountered an "
                "unsupported vertex type (%i)!", type);
            return Spectrum(0.0f);
    }

    return result;
}

Float PathVertex::evalPdf(const Scene *scene, const PathVertex *pred,
        const PathVertex *succ, ETransportMode mode, EMeasure measure) const {
    Vector wo(0.0f);
    Float dist = 0.0f, result = 0.0f;

    switch (type) {
        case EEmitterSupernode: {
                if (mode != EImportance || pred != NULL || succ->type != EEmitterSample)
                    return 0.0f;
                PositionSamplingRecord pRec = succ->getPositionSamplingRecord();
                pRec.measure = measure;
                return scene->pdfEmitterPosition(pRec);
            }
            break;

        case ESensorSupernode: {
                if (mode != ERadiance || pred != NULL || succ->type != ESensorSample)
                    return 0.0f;
                PositionSamplingRecord pRec = succ->getPositionSamplingRecord();
                pRec.measure = measure;
                return scene->pdfSensorPosition(pRec);
            }
            break;

        case EEmitterSample: {
                if (mode == ERadiance && succ->type == EEmitterSupernode)
                    return 1.0f;
                else if (mode != EImportance || pred->type != EEmitterSupernode)
                    return 0.0f;

                const PositionSamplingRecord &pRec = getPositionSamplingRecord();
                const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
                wo = succ->getPosition() - pRec.p;
                dist = wo.length(); wo /= dist;
                DirectionSamplingRecord dRec(wo, measure == EArea ? ESolidAngle : measure);
                result = emitter->pdfDirection(dRec, pRec);
            }
            break;

        case ESensorSample: {
                if (mode == EImportance && succ->type == ESensorSupernode)
                    return 1.0f;
                else if (mode != ERadiance || pred->type != ESensorSupernode)
                    return 0.0f;

                const PositionSamplingRecord &pRec = getPositionSamplingRecord();
                const Sensor *sensor = static_cast<const Sensor *>(pRec.object);
                wo = succ->getPosition() - pRec.p;
                dist = wo.length(); wo /= dist;
                DirectionSamplingRecord dRec(wo, measure == EArea ? ESolidAngle : measure);
                result = sensor->pdfDirection(dRec, pRec);
            }
            break;

        case ESurfaceInteraction: {
                const Intersection &its = getIntersection();
                const BSDF *bsdf = its.getBSDF();
                wo = succ->getPosition() - its.p;
                dist = wo.length(); wo /= dist;

                Point predP = pred->getPosition();
                Vector wi = normalize(predP - its.p);

                BSDFSamplingRecord bRec(its, its.toLocal(wi), its.toLocal(wo), mode);
                result = bsdf->pdf(bRec, measure == EArea ? ESolidAngle : measure);

                /* Prevent light leaks due to the use of shading normals */
                Float wiDotGeoN = dot(its.geoFrame.n, wi),
                      woDotGeoN = dot(its.geoFrame.n, wo);

                if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
                    woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                    return 0.0f;
            }
            break;

        case EMediumInteraction: {
                if (measure != ESolidAngle && measure != EArea)
                    return 0.0f;

                const MediumSamplingRecord &mRec = getMediumSamplingRecord();
                wo = succ->getPosition() - mRec.p;
                dist = wo.length(); wo /= dist;

                Point predP = pred->getPosition();
                Vector wi = normalize(predP - mRec.p);

                const PhaseFunction *phase = mRec.medium->getPhaseFunction();
                PhaseFunctionSamplingRecord pRec(mRec, wi, wo, mode);

                result = phase->pdf(pRec);
            }
            break;

        default:
            SLog(EError, "PathVertex::evalPdf(): Encountered an "
                "unsupported vertex type (%i)!", type);
            return 0.0f;
    }

    if (measure == EArea) {
        result /= dist * dist;
        if (succ->isOnSurface())
            result *= absDot(wo, succ->getGeometricNormal());
    }

    return result;
}

Spectrum PathVertex::sampleDirect(const Scene *scene, Sampler *sampler,
        PathVertex *endpoint, PathEdge *edge, PathVertex *sample, ETransportMode mode) const {
    if (isDegenerate() || isAbsorbing())
        return Spectrum(0.0f);

    memset(edge, 0, sizeof(PathEdge));
    memset(endpoint, 0, sizeof(PathVertex));
    memset(sample, 0, sizeof(PathVertex));

    bool emitter = (mode == EImportance);
    DirectSamplingRecord dRec;
    if (isSurfaceInteraction())
        dRec = DirectSamplingRecord(getIntersection());
    else
        dRec = DirectSamplingRecord(getPosition(), getTime());

    Spectrum value;
    Point2 rsamp(0.5f);
    if (sampler)
        rsamp = sampler->next2D();

    if (emitter)
        value = scene->sampleEmitterDirect(dRec, rsamp, false);
    else
        value = scene->sampleSensorDirect(dRec, rsamp, false);

    if (value.isZero())
        return Spectrum(0.0f);

    const AbstractEmitter *ae = static_cast<const AbstractEmitter *>(dRec.object);
    bool degenPos = ae->getType() & AbstractEmitter::EDeltaPosition;
    bool degenDir = ae->getType() & AbstractEmitter::EDeltaDirection;

    endpoint->type = emitter ? EEmitterSupernode : ESensorSupernode;
    endpoint->measure = degenPos ? EDiscrete : EArea;
    endpoint->degenerate = degenPos;

    endpoint->getEndpointRecord() = EndpointRecord(dRec.time);

    /* Be resilient to FP issues */
    if ((dRec.ref - dRec.p).lengthSquared() <= 0)
        return Spectrum(0.0f);

    edge->medium = ae->getMedium();
    edge->pdf[mode] = 1.0f;
    edge->weight[mode] = Spectrum(1.0f);

    sample->type = emitter ? EEmitterSample : ESensorSample;
    sample->degenerate = degenDir;
    sample->measure = degenDir ? EDiscrete : EArea;
    sample->getPositionSamplingRecord() = PositionSamplingRecord(dRec);

    dRec.measure = (EMeasure) endpoint->measure;
    endpoint->pdf[mode] = emitter ?
        scene->pdfEmitterPosition(dRec) :
        scene->pdfSensorPosition(dRec);

    endpoint->weight[mode] = ae->evalPosition(dRec) / endpoint->pdf[mode];

    return value;
}

Float PathVertex::evalPdfDirect(const Scene *scene,
        const PathVertex *sample, ETransportMode mode, EMeasure measure) const {
    BDAssert((mode == EImportance && sample->type == EEmitterSample) ||
        (mode == ERadiance && sample->type == ESensorSample));
    bool emitter = (mode == EImportance);

    DirectSamplingRecord dRec;
    if (isSurfaceInteraction())
        dRec = DirectSamplingRecord(getIntersection());
    else
        dRec = DirectSamplingRecord(getPosition(), getTime());

    const PositionSamplingRecord &pRec = sample->getPositionSamplingRecord();

    dRec.p = pRec.p;
    dRec.n = pRec.n;
    dRec.uv = pRec.uv;
    dRec.measure = measure;
    dRec.object = pRec.object;
    dRec.d = sample->getPosition() - getPosition();
    dRec.dist = dRec.d.length();
    dRec.d /= dRec.dist;

    if (emitter)
        return scene->pdfEmitterDirect(dRec);
    else
        return scene->pdfSensorDirect(dRec);
}

bool PathVertex::cast(const Scene *scene, EVertexType desired) {
    if (desired == type) {
        return true;
    } else if (desired == EEmitterSample) {
        if (type != ESurfaceInteraction)
            return false;

        const Intersection &its = getIntersection();
        const Emitter *emitter = its.shape->getEmitter();
        if (!emitter)
            return false;

        type = desired;
        PositionSamplingRecord pRec(its);
        pRec.object = emitter;
        pRec.pdf = 0.0f;
        getPositionSamplingRecord() = pRec;
        measure = pRec.measure;
        degenerate = emitter->getType() & Emitter::EDeltaDirection;

        return true;
    } else if (desired == ESensorSample) {
        if (type != ESurfaceInteraction)
            return false;

        const Intersection &its = getIntersection();
        const Sensor *sensor = its.shape->getSensor();

        if (sensor != scene->getSensor())
            return false;

        type = desired;
        PositionSamplingRecord pRec(its);
        pRec.object = sensor;
        pRec.pdf = 0.0f;

        Vector2i size = sensor->getFilm()->getSize();
        pRec.uv.x *= size.x; pRec.uv.y *= size.y;
        getPositionSamplingRecord() = pRec;
        measure = pRec.measure;
        degenerate = sensor->getType() & Sensor::EDeltaDirection;

        return true;
    } else {
        SLog(EError, "Unsupported conversion request from type %i->%i!",
                type, desired);
        return false;
    }
}

bool PathVertex::update(const Scene *scene, const PathVertex *pred,
        const PathVertex *succ, ETransportMode mode, EMeasure measure) {

    pdf[mode]       = evalPdf(scene, pred, succ, mode, measure);
    pdf[1-mode]     = evalPdf(scene, succ, pred, (ETransportMode) (1-mode), measure);
    weight[mode]    = eval(scene, pred, succ, mode, measure);
    weight[1-mode]  = eval(scene, succ, pred, (ETransportMode) (1-mode), measure);

    if (weight[mode].isZero() || pdf[mode] <= RCPOVERFLOW)
        return false;

    Float weightFwd = pdf[mode]   <= RCPOVERFLOW ? 0 : 1 / pdf[mode],
          weightBkw = pdf[1-mode] <= RCPOVERFLOW ? 0 : 1 / pdf[1-mode];

    this->measure = measure;

    if (!isSupernode() && measure == EArea) {
        if (!pred->isSupernode()) {
            Vector d = pred->getPosition() - getPosition();
            Float invDistSqr = 1.0f / d.lengthSquared();
            weightBkw *= invDistSqr;
            d *= std::sqrt(invDistSqr);
            if (isOnSurface() && isConnectable())
                weightBkw *= absDot(getShadingNormal(), d);
            if (pred->isOnSurface())
                weightBkw *= absDot(pred->getGeometricNormal(), d);
        }

        if (!succ->isSupernode()) {
            Vector d = succ->getPosition() - getPosition();
            Float invDistSqr = 1.0f / d.lengthSquared();
            weightFwd *= invDistSqr;
            d *= std::sqrt(invDistSqr);
            if (isOnSurface() && isConnectable())
                weightFwd *= absDot(getShadingNormal(), d);
            if (succ->isOnSurface())
                weightFwd *= absDot(succ->getGeometricNormal(), d);
        }
        if (isSurfaceInteraction())
            componentType = BSDF::ESmooth;
    }

    weight[mode]   *= weightFwd;
    weight[1-mode] *= weightBkw;

    return true;
}

Point PathVertex::getPosition() const {
    switch (type) {
        case ESurfaceInteraction:
            return getIntersection().p;
        case EMediumInteraction:
            return getMediumSamplingRecord().p;
        case EEmitterSample:
        case ESensorSample:
            return getPositionSamplingRecord().p;
        default:
            SLog(EError, "PathVertex::getPosition(): Encountered an "
                "unsupported vertex type (%i)!", type);
            return Point(0.0f);
    }
}

Normal PathVertex::getShadingNormal() const {
    switch (type) {
        case ESurfaceInteraction:
            return getIntersection().shFrame.n;
        case EEmitterSample:
        case ESensorSample:
            return getPositionSamplingRecord().n;
        default:
            SLog(EError, "PathVertex::getShadingNormal(): Encountered an "
                "unsupported vertex type (%i)!", type);
            return Normal(0.0f);
    }
}

Normal PathVertex::getGeometricNormal() const {
    switch (type) {
        case ESurfaceInteraction:
            return getIntersection().geoFrame.n;
        case EEmitterSample:
        case ESensorSample:
            return getPositionSamplingRecord().n;
        default:
            SLog(EError, "PathVertex::getGeometricNormal(): Encountered an "
                "unsupported vertex type (%i)!", type);
            return Normal(0.0f);
    }
}

Float PathVertex::getTime() const {
    switch (type) {
        case ESurfaceInteraction:
            return getIntersection().time;
        case EMediumInteraction:
            return getMediumSamplingRecord().time;
        case EEmitterSample:
        case ESensorSample:
            return getPositionSamplingRecord().time;
        case EEmitterSupernode:
        case ESensorSupernode:
            return getEndpointRecord().time;
        default:
            SLog(EError, "PathVertex::getTime(): Encountered an "
                "unsupported vertex type (%i)!", type);
            return 0.0f;
    }
}

const Medium *PathVertex::getTargetMedium(const PathEdge *predEdge,
    const PathVertex *succ) const {
    if (isSurfaceInteraction()) {
        const Intersection &its = getIntersection();
        if (its.isMediumTransition())
            return its.getTargetMedium(succ->getPosition() - its.p);
    }
    return predEdge->medium;
}

const Medium *PathVertex::getTargetMedium(const PathEdge *predEdge,
    const Vector &d) const {
    if (isSurfaceInteraction()) {
        const Intersection &its = getIntersection();
        if (its.isMediumTransition())
            return its.getTargetMedium(d);
    }
    return predEdge->medium;
}

bool PathVertex::updateSamplePosition(const PathVertex *v) {
    BDAssert(isSensorSample());

    PositionSamplingRecord &pRec = getPositionSamplingRecord();
    const Sensor *sensor = static_cast<const Sensor *>(pRec.object);
    DirectionSamplingRecord dRec(v->getPosition() - getPosition());

    return sensor->getSamplePosition(pRec, dRec, pRec.uv);
}

bool PathVertex::getSamplePosition(const PathVertex *v, Point2 &result) const {
    BDAssert(isSensorSample());

    const PositionSamplingRecord &pRec = getPositionSamplingRecord();
    const Sensor *sensor = static_cast<const Sensor *>(pRec.object);
    const DirectionSamplingRecord dRec(v->getPosition() - getPosition());

    return sensor->getSamplePosition(pRec, dRec, result);
}

bool PathVertex::connect(const Scene *scene,
        const PathVertex *pred, const PathEdge *predEdge,
        PathVertex *vs, PathEdge *edge, PathVertex *vt,
        const PathEdge *succEdge, const PathVertex *succ) {

    if (vs->isEmitterSupernode()) {
        if (!vt->cast(scene, PathVertex::EEmitterSample))
            return false;
    } else if (vt->isSensorSupernode()) {
        /* If possible, convert 'vs' into an sensor sample */
        if (!vs->cast(scene, PathVertex::ESensorSample))
            return false;
    } else if (vs->getPosition() == vt->getPosition()) {
        /* Check for this here to avoid dividing by zero
           when computing the direction vs->vt. */
        return false;
    }

    if (vs->isDegenerate() || vt->isDegenerate())
        return false;

    if (!vs->update(scene, pred, vt, EImportance))
        return false;

    if (!vt->update(scene, succ, vs, ERadiance))
        return false;

    return edge->connect(scene, predEdge, vs, vt, succEdge);
}

bool PathVertex::connect(const Scene *scene,
        const PathVertex *pred, const PathEdge *predEdge,
        PathVertex *vs, PathEdge *edge, PathVertex *vt,
        const PathEdge *succEdge, const PathVertex *succ,
        EMeasure vsMeasure, EMeasure vtMeasure) {

    if (vs->isEmitterSupernode()) {
        if (!vt->cast(scene, PathVertex::EEmitterSample))
            return false;
    } else if (vt->isSensorSupernode()) {
        /* If possible, convert 'vs' into an sensor sample */
        if (!vs->cast(scene, PathVertex::ESensorSample))
            return false;
    }

    if (!vs->update(scene, pred, vt, EImportance, vsMeasure))
        return false;

    if (!vt->update(scene, succ, vs, ERadiance, vtMeasure))
        return false;

    return edge->connect(scene, predEdge, vs, vt, succEdge);
}


PathVertex *PathVertex::clone(MemoryPool &pool) const {
    PathVertex *result = pool.allocVertex();
    *result = *this;
    return result;
}

std::string PathVertex::toString() const {
    std::ostringstream oss;
    oss << "PathVertex[" << endl
        << "  type = " << (EVertexType) type << "," << endl;

    switch (type) {
        case ESensorSupernode:
            oss << "  data = " << indent(getEndpointRecord().toString()) << "," << endl;
            break;
        case EEmitterSupernode:
            oss << "  data = " << indent(getEndpointRecord().toString()) << "," << endl;
            break;
        case ESensorSample:
            oss << "  data = " << indent(getPositionSamplingRecord().toString()) << "," << endl;
            break;
        case EEmitterSample:
            oss << "  data = " << indent(getPositionSamplingRecord().toString()) << "," << endl;
            break;
        case ESurfaceInteraction:
            oss << "  data = " << indent(getIntersection().toString()) << "," << endl
                << "  componentType = " << (int) componentType << "," << endl;
            break;
        case EMediumInteraction:
            oss << "  data = " << indent(getMediumSamplingRecord().toString()) << "," << endl;
            break;
        default:
            break;
    }

    oss << "  degenerate = " << (degenerate ? "true" : "false") << "," << endl
        << "  measure = " << (EMeasure) measure << "," << endl
        << "  weight[importance] = " << weight[EImportance].toString() << "," << endl
        << "  weight[radiance] = " << weight[ERadiance].toString() << "," << endl
        << "  pdf[importance] = " << pdf[EImportance] << "," << endl
        << "  pdf[radiance] = " << pdf[ERadiance] << endl
        << "]";
    return oss.str();
}

bool PathVertex::operator==(const PathVertex &other) const {
    return other.type == type &&
        other.degenerate == degenerate &&
        other.measure == measure &&
        other.componentType == componentType &&
        other.weight[EImportance] == weight[EImportance] &&
        other.weight[ERadiance] == weight[ERadiance] &&
        other.pdf[EImportance] == pdf[EImportance] &&
        other.pdf[ERadiance] == pdf[ERadiance] &&
        memcmp(other.data, data, EDataSize) == 0;
}

std::ostream &operator<<(std::ostream &os, PathVertex::EVertexType type) {
    switch (type) {
        case PathVertex::ESensorSupernode: os << "sensorSupernode"; break;
        case PathVertex::EEmitterSupernode: os << "emitterSupernode"; break;
        case PathVertex::ESensorSample: os << "sensorSample"; break;
        case PathVertex::EEmitterSample: os << "emitterSample"; break;
        case PathVertex::ESurfaceInteraction: os << "surfaceInteraction"; break;
        case PathVertex::EMediumInteraction: os << "mediumInteraction"; break;
        case PathVertex::EInvalid: default: os << "invalidType"; break;
    }
    return os;
}
MTS_NAMESPACE_END


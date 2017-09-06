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
        "Medium inconsistencies in connect()");

bool PathEdge::sampleNext(const Scene *scene, Sampler *sampler,
        const PathVertex *pred, const Ray &ray, PathVertex *succ,
        ETransportMode mode) {
    /* First, check if there is a surface in the sampled direction */
    Intersection &its = succ->getIntersection();
    bool surface = scene->rayIntersectAll(ray, its);

    /* Sample the RTE in-scattering integral -- this determines whether the
       next vertex is invalid or a surface or medium scattering event */
    MediumSamplingRecord mRec;
    if (medium && medium->sampleDistance(Ray(ray, 0, its.t), mRec, sampler)) {
        succ->type = PathVertex::EMediumInteraction;
        succ->degenerate = false;
        succ->getMediumSamplingRecord() = mRec;
        length = mRec.t;
    } else if (surface) {
        succ->type = PathVertex::ESurfaceInteraction;
        succ->degenerate = !(its.getBSDF()->hasComponent(BSDF::ESmooth) ||
                its.shape->isEmitter() || its.shape->isSensor());
        length = its.t;
    } else {
        return false;
    }

    if (length == 0)
        return false;

    if (!medium) {
        weight[ERadiance] = weight[EImportance] = Spectrum(1.0f);
        pdf[ERadiance] = pdf[EImportance] = 1.0f;
    } else {
        if (mRec.transmittance.isZero())
            return false;
        pdf[mode]   = succ->isMediumInteraction() ? mRec.pdfSuccess    : mRec.pdfFailure;
        pdf[1-mode] = pred->isMediumInteraction() ? mRec.pdfSuccessRev : mRec.pdfFailure;
        weight[mode]   = mRec.transmittance / pdf[mode];
        weight[1-mode] = mRec.transmittance / pdf[1-mode];
    }
    d = ray.d;
    /* Direction always points along the light path (from the light source along the path) */
    if (mode == ERadiance)
        d = -d;

    return true;
}

bool PathEdge::perturbDirection(const Scene *scene,
        const PathVertex *pred, const Ray &ray, Float dist,
        PathVertex::EVertexType desiredType, PathVertex *succ,
        ETransportMode mode) {
    /* First, check if there is a surface in the sampled direction */
    Intersection &its = succ->getIntersection();
    bool surface = scene->rayIntersectAll(ray, its);

    /* Sample the RTE in-scattering integral -- this determines whether the
       next vertex is invalid or a surface or medium scattering event */
    MediumSamplingRecord mRec;

    bool wantMedium =
        desiredType == PathVertex::EMediumInteraction;

    if ((wantMedium && dist > its.t) || dist <= 0)
        return false;

    if (medium)
        medium->eval(Ray(ray, 0,
            wantMedium ? std::min(dist, its.t) : its.t), mRec);

    if (medium && wantMedium) {
        succ->type = PathVertex::EMediumInteraction;
        succ->degenerate = false;
        length = dist;
        mRec.p = ray(dist);
        mRec.t = dist;
        succ->getMediumSamplingRecord() = mRec;
    } else if (surface && !wantMedium) {
        succ->type = PathVertex::ESurfaceInteraction;
        succ->degenerate = !(its.getBSDF()->hasComponent(BSDF::ESmooth) ||
                its.shape->isEmitter() || its.shape->isSensor());
        length = its.t;
    } else {
        return false;
    }
    d = ray.d;
    /* Direction always points along the light path (from the light source along the path) */
    if (mode == ERadiance)
        d = -d;

    if (length == 0)
        return false;

    if (!medium) {
        weight[ERadiance] = weight[EImportance] = Spectrum(1.0f);
        pdf[ERadiance] = pdf[EImportance] = 1.0f;
    } else {
        pdf[mode]   = succ->isMediumInteraction() ? mRec.pdfSuccess    : mRec.pdfFailure;
        pdf[1-mode] = pred->isMediumInteraction() ? mRec.pdfSuccessRev : mRec.pdfFailure;
        if (pdf[mode] == 0 || pdf[1-mode] == 0)
            return false;
        weight[mode]   = mRec.transmittance / pdf[mode];
        weight[1-mode] = mRec.transmittance / pdf[1-mode];
    }

    return true;
}

Spectrum PathEdge::evalTransmittance(const PathVertex *pred, const PathVertex *succ) const {
    if (succ->isSupernode())
        return Spectrum(0.0f);
    else if (!medium || pred->isSupernode())
        return Spectrum(1.0f);

    Point a = pred->getPosition(),
          b = succ->getPosition();
    Vector d(b-a);

    Float length = d.length();
    return medium->evalTransmittance(
        Ray(a, d/length, 0, length, pred->getTime()));
}

Float PathEdge::evalPdf(const PathVertex *pred,
        const PathVertex *succ) const {
    if (succ->isSupernode())
        return 0.0f;
    else if (!medium || pred->isSupernode())
        return 1.0f;

    Point a = pred->getPosition(),
          b = succ->getPosition();
    Vector d(b-a);

    Float length = d.length();
    Ray ray(a, d/length, 0, length, pred->getTime());

    MediumSamplingRecord mRec;
    medium->eval(ray, mRec);

    return succ->isMediumInteraction() ?
        mRec.pdfSuccess : mRec.pdfFailure;
}

Spectrum PathEdge::evalCached(const PathVertex *pred, const PathVertex *succ,
        unsigned int what) const {
    /* Extract the requested information based on what is currently cached in the
       vertex. The actual computation that has to happen here is pretty awful, but
       it works. It might be worth to change the caching scheme to make this function
       simpler in a future revision */
    Spectrum result(1.0f);

    if (length == 0) {
        if (what & EValueImp)
            result *= pred->weight[EImportance] * pred->pdf[EImportance];
        if (what & EValueRad)
            result *= succ->weight[ERadiance] * succ->pdf[ERadiance];
    } else {
        if (what & EValueImp) {
            Float tmp = pred->pdf[EImportance];
            if (pred->isConnectable()) {
                tmp *= length * length;
                if (succ->isOnSurface())
                    tmp /= dot(succ->getGeometricNormal(), d);
                if (pred->isOnSurface() && !(what & ECosineImp))
                    tmp /= dot(pred->getShadingNormal(), d);
            }
            result *= pred->weight[EImportance] * std::abs(tmp);
        } else if ((what & ECosineImp) && pred->isOnSurface() && pred->isConnectable()) {
            result *= absDot(pred->getShadingNormal(), d);
        }

        if (what & EValueRad) {
            Float tmp = succ->pdf[ERadiance];
            if (succ->isConnectable()) {
                tmp *= length * length;
                if (pred->isOnSurface())
                    tmp /= dot(pred->getGeometricNormal(), d);
                if (succ->isOnSurface() && !(what & ECosineRad))
                    tmp /= dot(succ->getShadingNormal(), d);
            }
            result *= succ->weight[ERadiance] * std::abs(tmp);
        } else if ((what & ECosineRad) && succ->isOnSurface() && succ->isConnectable()) {
            result *= absDot(succ->getShadingNormal(), d);
        }

        if (what & EInverseSquareFalloff)
            result /= length * length;

        if (what & ETransmittance)
            result *= weight[EImportance] * pdf[EImportance];
    }

    return result;
}

bool PathEdge::connect(const Scene *scene,
            const PathEdge *predEdge, const PathVertex *vs,
            const PathVertex *vt, const PathEdge *succEdge) {

    if (vs->isEmitterSupernode() || vt->isSensorSupernode()) {
        Float radianceTransport   = vt->isSensorSupernode() ? 1.0f : 0.0f,
              importanceTransport = 1-radianceTransport;

        medium = NULL;
        d = Vector(0.0f);
        length = 0.0f;
        pdf[ERadiance]   = radianceTransport;
        pdf[EImportance] = importanceTransport;
        weight[ERadiance] = Spectrum(radianceTransport);
        weight[EImportance] = Spectrum(importanceTransport);
    } else {
        Point vsp = vs->getPosition(), vtp = vt->getPosition();
        d = vsp-vtp;
        length = d.length();
        d /= length;

        Ray ray(vtp, d, vt->isOnSurface() ? Epsilon : 0, length *
            (vs->isOnSurface() ? (1-ShadowEpsilon) : 1), vs->getTime());

        /* Check for occlusion */
        if (scene->rayIntersectAll(ray))
            return false;

        const Medium *vtMedium = vt->getTargetMedium(succEdge, d);
        const Medium *vsMedium = vs->getTargetMedium(predEdge, -d);

        if (vsMedium != vtMedium) {
            #if defined(MTS_BD_TRACE)
                SLog(EWarn, "PathEdge::connect(): attempted two connect "
                    "two vertices that disagree about the medium in between! "
                    "Please check your scene for leaks.");
            #endif
            ++mediumInconsistencies;
            return false;
        }

        medium = vtMedium;

        if (medium) {
            MediumSamplingRecord mRec;
            medium->eval(ray, mRec);

            pdf[EImportance] = vt->isMediumInteraction() ? mRec.pdfSuccessRev : mRec.pdfFailure;
            pdf[ERadiance]   = vs->isMediumInteraction() ? mRec.pdfSuccess    : mRec.pdfFailure;

            /* Fail if there is no throughput */
            if (mRec.transmittance.isZero() || pdf[EImportance] == 0 || pdf[ERadiance] == 0)
                return false;

            weight[EImportance] = mRec.transmittance / pdf[EImportance];
            weight[ERadiance]   = mRec.transmittance / pdf[ERadiance];
        } else {
            weight[ERadiance] = weight[EImportance] = Spectrum(1.0f);
            pdf[ERadiance] = pdf[EImportance] = 1.0f;
        }
    }

    /* Direction always points along the light path (from the light source along the path) */
    d = -d;

    return true;
}

bool PathEdge::pathConnect(const Scene *scene, const PathEdge *predEdge,
        const PathVertex *vs, Path &result, const PathVertex *vt,
        const PathEdge *succEdge, int maxInteractions, MemoryPool &pool) {
    BDAssert(result.edgeCount() == 0 && result.vertexCount() == 0);

    if (vs->isEmitterSupernode() || vt->isSensorSupernode()) {
        Float radianceTransport   = vt->isSensorSupernode() ? 1.0f : 0.0f,
              importanceTransport = 1-radianceTransport;
        PathEdge *edge = pool.allocEdge();
        edge->medium = NULL;
        edge->length = 0.0f;
        edge->d = Vector(0.0f);
        edge->pdf[ERadiance]   = radianceTransport;
        edge->pdf[EImportance] = importanceTransport;
        edge->weight[ERadiance] = Spectrum(radianceTransport);
        edge->weight[EImportance] = Spectrum(importanceTransport);
        result.append(edge);
    } else {
        Point vsp = vs->getPosition(), vtp = vt->getPosition();
        Vector d(vsp-vtp);
        Float remaining = d.length();
        d /= remaining;
        if (remaining == 0) {
            #if defined(MTS_BD_DEBUG)
                SLog(EWarn, "Tried to connect %s and %s, which are located at exactly the same position!",
                    vs->toString().c_str(), vt->toString().c_str());
            #endif
            return false;
        }

        Float lengthFactor = vs->isOnSurface() ? (1-ShadowEpsilon) : 1;
        Ray ray(vtp, d, vt->isOnSurface() ? Epsilon : 0,
                remaining * lengthFactor, vs->getTime());
        const Medium *medium = vt->getTargetMedium(succEdge,  d);

        int interactions = 0;

        Intersection its;
        while (true) {
            bool surface = scene->rayIntersectAll(ray, its);

            if (surface && (interactions == maxInteractions ||
                !(its.getBSDF()->getType() & BSDF::ENull))) {
                /* Encountered an occluder -- zero transmittance. */
                result.release(pool);
                return false;
            }

            /* Construct an edge */
            PathEdge *edge = pool.allocEdge();
            result.append(edge);
            edge->length = std::min(its.t, remaining);
            edge->medium = medium;
            /* Direction always points along the light path (from the light source along the path) */
            edge->d = -d;

            if (medium) {
                MediumSamplingRecord mRec;
                medium->eval(Ray(ray, 0, edge->length), mRec);
                edge->pdf[ERadiance] = (surface || !vs->isMediumInteraction())
                    ? mRec.pdfFailure : mRec.pdfSuccess;
                edge->pdf[EImportance] = (interactions > 0 || !vt->isMediumInteraction())
                    ? mRec.pdfFailure : mRec.pdfSuccessRev;

                if (edge->pdf[ERadiance] == 0 || edge->pdf[EImportance] == 0
                        || mRec.transmittance.isZero()) {
                    /* Zero transmittance */
                    result.release(pool);
                    return false;
                }
                edge->weight[EImportance] = mRec.transmittance / edge->pdf[EImportance];
                edge->weight[ERadiance]   = mRec.transmittance / edge->pdf[ERadiance];
            } else {
                edge->weight[ERadiance] = edge->weight[EImportance] = Spectrum(1.0f);
                edge->pdf[ERadiance] = edge->pdf[EImportance] = 1.0f;
            }

            if (!surface || remaining - its.t < 0)
                break;

            /* Advance the ray */
            ray.o = ray(its.t);
            remaining -= its.t;
            ray.mint = Epsilon;
            ray.maxt = remaining * lengthFactor;

            const BSDF *bsdf = its.getBSDF();

            /* Account for the ENull interaction */
            Vector wo = its.toLocal(ray.d);
            BSDFSamplingRecord bRec(its, -wo, wo, ERadiance);
            bRec.component = BSDF::ENull;
            Float nullPdf = bsdf->pdf(bRec, EDiscrete);
            if (nullPdf == 0) {
                result.release(pool);
                return false;
            }

            PathVertex *vertex = pool.allocVertex();
            vertex->type = PathVertex::ESurfaceInteraction;
            vertex->degenerate = !(bsdf->hasComponent(BSDF::ESmooth)
                || its.shape->isEmitter() || its.shape->isSensor());
            vertex->measure = EDiscrete;
            vertex->componentType = BSDF::ENull;
            vertex->pdf[EImportance] = vertex->pdf[ERadiance] = nullPdf;
            vertex->weight[EImportance] = vertex->weight[ERadiance]
                = bsdf->eval(bRec, EDiscrete) / nullPdf;
            vertex->rrWeight = 1.0f;
            vertex->getIntersection() = its;
            result.append(vertex);

            if (its.isMediumTransition()) {
                const Medium *expected = its.getTargetMedium(-ray.d);
                if (medium != expected) {
                    #if defined(MTS_BD_TRACE)
                        SLog(EWarn, "PathEdge::pathConnect(): attempted two connect "
                            "two vertices that disagree about the medium in between! "
                            "Please check your scene for leaks.");
                    #endif
                    ++mediumInconsistencies;
                    result.release(pool);
                    return false;
                }
                medium = its.getTargetMedium(ray.d);
            }

            if (++interactions > 100) { /// Just a precaution..
                SLog(EWarn, "pathConnect(): round-off error issues?");
                result.release(pool);
                return false;
            }
        }

        if (medium != vs->getTargetMedium(predEdge, -d)) {
            #if defined(MTS_BD_TRACE)
                SLog(EWarn, "PathEdge::pathConnect(): attempted two connect "
                    "two vertices that disagree about the medium in between! "
                    "Please check your scene for leaks.");
            #endif
            ++mediumInconsistencies;
            result.release(pool);
            return false;
        }
    }

    result.reverse();

    BDAssert(result.edgeCount() == result.vertexCount() + 1);
    BDAssert((int) result.vertexCount() <= maxInteractions || maxInteractions < 0);

    return true;
}

bool PathEdge::pathConnectAndCollapse(const Scene *scene, const PathEdge *predEdge,
        const PathVertex *vs, const PathVertex *vt,
        const PathEdge *succEdge, int &interactions) {
    if (vs->isEmitterSupernode() || vt->isSensorSupernode()) {
        Float radianceTransport   = vt->isSensorSupernode() ? 1.0f : 0.0f,
              importanceTransport = 1-radianceTransport;
        medium = NULL;
        length = 0.0f;
        d = Vector(0.0f);
        pdf[ERadiance]   = radianceTransport;
        pdf[EImportance] = importanceTransport;
        weight[ERadiance] = Spectrum(radianceTransport);
        weight[EImportance] = Spectrum(importanceTransport);
        interactions = 0;
    } else {
        Point vsp = vs->getPosition(), vtp = vt->getPosition();
        d = vsp-vtp;
        length = d.length();
        int maxInteractions = interactions;
        interactions = 0;

        if (length == 0) {
            #if defined(MTS_BD_DEBUG)
                SLog(EWarn, "Tried to connect %s and %s, which are located at exactly the same position!",
                    vs->toString().c_str(), vt->toString().c_str());
            #endif
            return false;
        }

        d /= length;
        Float lengthFactor = vs->isOnSurface() ? (1-ShadowEpsilon) : 1;
        Ray ray(vtp, d, vt->isOnSurface() ? Epsilon : 0, length * lengthFactor, vs->getTime());

        weight[ERadiance] = Spectrum(1.0f);
        weight[EImportance] = Spectrum(1.0f);
        pdf[ERadiance] = 1.0f;
        pdf[EImportance] = 1.0f;

        Intersection its;
        Float remaining = length;
        medium = vt->getTargetMedium(succEdge, d);

        while (true) {
            bool surface = scene->rayIntersectAll(ray, its);

            if (surface && (interactions == maxInteractions ||
                !(its.getBSDF()->getType() & BSDF::ENull))) {
                /* Encountered an occluder -- zero transmittance. */
                return false;
            }

            if (medium) {
                Float segmentLength = std::min(its.t, remaining);
                MediumSamplingRecord mRec;
                medium->eval(Ray(ray, 0, segmentLength), mRec);

                Float pdfRadiance = (surface || !vs->isMediumInteraction())
                    ? mRec.pdfFailure : mRec.pdfSuccess;
                Float pdfImportance = (interactions > 0 || !vt->isMediumInteraction())
                    ? mRec.pdfFailure : mRec.pdfSuccessRev;

                if (pdfRadiance == 0 || pdfImportance == 0 || mRec.transmittance.isZero()) {
                    /* Zero transmittance */
                    return false;
                }

                weight[EImportance] *= mRec.transmittance / pdfImportance;
                weight[ERadiance] *= mRec.transmittance / pdfRadiance;
                pdf[EImportance] *= pdfImportance;
                pdf[ERadiance] *= pdfRadiance;
            }

            if (!surface || remaining - its.t < 0)
                break;

            /* Advance the ray */
            ray.o = ray(its.t);
            remaining -= its.t;
            ray.mint = Epsilon;
            ray.maxt = remaining * lengthFactor;

            /* Account for the ENull interaction */
            const BSDF *bsdf = its.getBSDF();
            Vector wo = its.toLocal(ray.d);
            BSDFSamplingRecord bRec(its, -wo, wo, ERadiance);
            bRec.component = BSDF::ENull;
            Float nullPdf = bsdf->pdf(bRec, EDiscrete);
            if (nullPdf == 0)
                return false;

            Spectrum nullWeight = bsdf->eval(bRec, EDiscrete) / nullPdf;

            weight[EImportance] *= nullWeight;
            weight[ERadiance] *= nullWeight;
            pdf[EImportance] *= nullPdf;
            pdf[ERadiance] *= nullPdf;

            if (its.isMediumTransition()) {
                const Medium *expected = its.getTargetMedium(-ray.d);
                if (medium != expected) {
                    #if defined(MTS_BD_TRACE)
                        SLog(EWarn, "PathEdge::pathConnectAndCollapse(): attempted two connect "
                            "two vertices that disagree about the medium in between! "
                            "Please check your scene for leaks.");
                    #endif
                    ++mediumInconsistencies;
                    return false;
                }
                medium = its.getTargetMedium(ray.d);
            }

            if (++interactions > 100) { /// Just a precaution..
                SLog(EWarn, "pathConnectAndCollapse(): round-off error issues?");
                return false;
            }
        }

        if (medium != vs->getTargetMedium(predEdge, -d)) {
            #if defined(MTS_BD_TRACE)
                SLog(EWarn, "PathEdge::pathConnectAndCollapse(): attempted two connect "
                    "two vertices that disagree about the medium in between! "
                    "Please check your scene for leaks.");
            #endif
            ++mediumInconsistencies;
            return false;
        }
    }

    /* Direction always points along the light path (from the light source along the path) */
    d = -d;

    return true;
}

PathEdge *PathEdge::clone(MemoryPool &pool) const {
    PathEdge *result = pool.allocEdge();
    *result = *this;
    return result;
}

bool PathEdge::operator==(const PathEdge &other) const {
    return other.medium == medium &&
        other.d == d &&
        other.length == length &&
        other.weight[EImportance] == weight[EImportance] &&
        other.weight[ERadiance] == weight[ERadiance] &&
        other.pdf[EImportance] == pdf[EImportance] &&
        other.pdf[ERadiance] == pdf[ERadiance];
}

std::string PathEdge::toString() const {
    std::ostringstream oss;
    oss << "PathEdge[" << endl
        << "  medium = " << indent(medium ? medium->toString().c_str() : "null") << "," << endl
        << "  d = " << d.toString() << "," << endl
        << "  length = " << length << "," << endl
        << "  weight[importance] = " << weight[EImportance].toString() << "," << endl
        << "  weight[radiance] = " << weight[ERadiance].toString() << "," << endl
        << "  pdf[importance] = " << pdf[EImportance] << "," << endl
        << "  pdf[radiance] = " << pdf[ERadiance] << endl
        << "]";
    return oss.str();
}

MTS_NAMESPACE_END

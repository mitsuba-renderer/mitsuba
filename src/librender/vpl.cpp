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

#include <mitsuba/render/vpl.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

static StatsCounter prunedVPLs("VPL renderer", "Pruned VPLs", EPercentage);

static void appendVPL(const Scene *scene, Random *random,
    VPL &vpl, bool prune, std::deque<VPL> &vpls) {
    prunedVPLs.incrementBase();

    const Sensor *sensor = scene->getSensor();
    Float time = random->nextFloat();

    if (prune) {
        /* Possibly reject VPLs if they are unlikely to be
           visible from the camera */
        int nSuccesses = 0, nSamples = 50;
        const Shape *shape;
        Normal n;
        Point2 uv;
        Ray ray;
        Vector2i size = sensor->getFilm()->getCropSize();
        for (int i=0; i<nSamples; ++i) {
            if (sensor->needsTimeSample())
                time = random->nextFloat();

            sensor->sampleRay(ray, Point2(random->nextFloat() * size.x,
                    random->nextFloat() * size.y), Point2(0.5f), time);

            Float t;
            if (scene->rayIntersect(ray, t, shape, n, uv)) {
                Point p = ray(t);
                Vector d = vpl.its.p - p;
                Float length = d.length();
                Ray shadowRay(p, d/length, Epsilon, length*(1-ShadowEpsilon), time);

                if (!scene->rayIntersect(shadowRay))
                    ++nSuccesses;
            } else {
                ++nSuccesses; // be conservative
            }
        }
        /// Have a small chance of acceptance in any case
        Float acceptanceProb = (nSuccesses+1) / (Float) (nSamples+1);
        if (random->nextFloat() < acceptanceProb) {
            vpl.P /= acceptanceProb;
            vpls.push_back(vpl);
        } else {
            ++prunedVPLs;
        }
    } else {
        vpls.push_back(vpl);
    }
}

size_t generateVPLs(const Scene *scene, Random *random,
        size_t offset, size_t count, int maxDepth, bool prune, std::deque<VPL> &vpls) {
    if (maxDepth <= 1)
        return 0;

    static Sampler *sampler = NULL;
    if (!sampler) {
        Properties props("halton");
        props.setInteger("scramble", 0);
        sampler = static_cast<Sampler *> (PluginManager::getInstance()->
            createObject(MTS_CLASS(Sampler), props));
        sampler->configure();
        sampler->generate(Point2i(0));
    }

    const Sensor *sensor = scene->getSensor();
    Float time = sensor->getShutterOpen()
        + sensor->getShutterOpenTime() * sampler->next1D();

    const Frame stdFrame(Vector(1,0,0), Vector(0,1,0), Vector(0,0,1));
    int retries = 0;

    while (vpls.size() < count) {
        sampler->setSampleIndex(++offset);

        if (vpls.empty() && ++retries > 10000) {
            /* Unable to generate VPLs in this scene -- give up. */
            return 0;
        }

        PositionSamplingRecord pRec(time);
        DirectionSamplingRecord dRec;
        Spectrum weight = scene->sampleEmitterPosition(pRec,
            sampler->next2D());

        size_t start = vpls.size();

        /* Sample an emitted particle */
        const Emitter *emitter = static_cast<const Emitter *>(pRec.object);

        if (!emitter->isEnvironmentEmitter() && emitter->needsDirectionSample()) {
            VPL lumVPL(EPointEmitterVPL, weight);
            lumVPL.its.p = pRec.p;
            lumVPL.its.time = time;
            lumVPL.its.shFrame = pRec.n.isZero() ? stdFrame : Frame(pRec.n);
            lumVPL.emitter = emitter;
            appendVPL(scene, random, lumVPL, prune, vpls);

            weight *= emitter->sampleDirection(dRec, pRec, sampler->next2D());
        } else {
            /* Hack to get the proper information for directional VPLs */
            DirectSamplingRecord diRec(
                scene->getKDTree()->getAABB().getCenter(), pRec.time);

            Spectrum weight2 = emitter->sampleDirect(diRec, sampler->next2D())
                / scene->pdfEmitterDiscrete(emitter);

            if (weight2.isZero())
                continue;

            VPL lumVPL(EDirectionalEmitterVPL, weight2);
            lumVPL.its.p = Point(0.0);
            lumVPL.its.time = time;
            lumVPL.its.shFrame = Frame(-diRec.d);
            lumVPL.emitter = emitter;
            appendVPL(scene, random, lumVPL, false, vpls);
            dRec.d = -diRec.d;

            Point2 offset = warp::squareToUniformDiskConcentric(sampler->next2D());
            Vector perpOffset = Frame(diRec.d).toWorld(Vector(offset.x, offset.y, 0));
            BSphere geoBSphere = scene->getKDTree()->getAABB().getBSphere();
            pRec.p = geoBSphere.center + (perpOffset - dRec.d) * geoBSphere.radius;
            weight = weight2 * M_PI * geoBSphere.radius * geoBSphere.radius;
        }

        int depth = 2;
        Ray ray(pRec.p, dRec.d, time);
        Intersection its;

        while (!weight.isZero() && (depth < maxDepth || maxDepth == -1)) {
            if (!scene->rayIntersect(ray, its))
                break;

            const BSDF *bsdf = its.getBSDF();
            BSDFSamplingRecord bRec(its, sampler, EImportance);
            Spectrum bsdfVal = bsdf->sample(bRec, sampler->next2D());
            if (bsdfVal.isZero())
                break;

            /* Assuming that BSDF importance sampling is perfect,
                the following should equal the maximum albedo
                over all spectral samples */
            Float approxAlbedo = std::min((Float) 0.95f, bsdfVal.max());
            if (sampler->next1D() > approxAlbedo)
                break;
            else
                weight /= approxAlbedo;

            VPL vpl(ESurfaceVPL, weight);
            vpl.its = its;

            if (BSDF::getMeasure(bRec.sampledType) == ESolidAngle)
                appendVPL(scene, random, vpl, prune, vpls);

            weight *= bsdfVal;

            Vector wi = -ray.d, wo = its.toWorld(bRec.wo);
            ray = Ray(its.p, wo, 0.0f);

            /* Prevent light leaks due to the use of shading normals -- [Veach, p. 158] */
            Float wiDotGeoN = dot(its.geoFrame.n, wi),
                woDotGeoN = dot(its.geoFrame.n, wo);
            if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
                woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                break;

            /* Disabled for now -- this increases VPL weights
               and accuracy is not really a big requirement */
            #if 0
                /* Adjoint BSDF for shading normals -- [Veach, p. 155] */
                weight *= std::abs(
                    (Frame::cosTheta(bRec.wi) * woDotGeoN)/
                    (Frame::cosTheta(bRec.wo) * wiDotGeoN));
            #endif

            ++depth;
        }

        size_t end = vpls.size();
        for (size_t i=start; i<end; ++i)
            vpls[i].emitterScale = 1.0f / (end - start);
    }

    return offset;
}

const char *toString(EVPLType type) {
    switch (type) {
        case EPointEmitterVPL: return "emitterVPL";
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
    if (type == EPointEmitterVPL) {
        oss << "  p = " << its.p.toString() << "," << endl;
        oss << "  emitter = " << indent(emitter->toString()) << endl;
    } else {
        oss << "  its = " << indent(its.toString()) << endl;
    }
    oss << "]" << endl;
    return oss.str();
}

MTS_NAMESPACE_END

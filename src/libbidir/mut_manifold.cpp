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

#include <mitsuba/bidir/mut_manifold.h>
#include <mitsuba/bidir/manifold.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/vmf.h>

MTS_NAMESPACE_BEGIN

#define DIFF_SAMPLES 50000

static StatsCounter statsAcceptedRad("Manifold perturbation",
        "Acceptance rate (rad. transport)", EPercentage);
static StatsCounter statsGeneratedRad("Manifold perturbation",
        "Successful generation rate (rad. transport)", EPercentage);

static StatsCounter statsAcceptedImp("Manifold perturbation",
        "Acceptance rate (imp. transport)", EPercentage);
static StatsCounter statsGeneratedImp("Manifold perturbation",
        "Successful generation rate (imp. transport)", EPercentage);
static StatsCounter statsUsedManifold("Manifold perturbation",
        "Perturbations involving manifold walks", EPercentage);
static StatsCounter statsNonReversible("Manifold perturbation",
        "Non-reversible walks", EPercentage);
static StatsCounter statsRoughMediumSpecular("Manifold perturbation",
        "Medium treated as specular", EPercentage);
static StatsCounter statsRoughSurfaceSpecular("Manifold perturbation",
        "Rough material treated as specular", EPercentage);

Float ManifoldPerturbation::m_thetaDiffSurface;
Float ManifoldPerturbation::m_thetaDiffMedium;
int ManifoldPerturbation::m_thetaDiffSurfaceSamples;
int ManifoldPerturbation::m_thetaDiffMediumSamples;
Mutex *ManifoldPerturbation::m_thetaDiffMutex = new Mutex();

ManifoldPerturbation::ManifoldPerturbation(const Scene *scene, Sampler *sampler,
          MemoryPool &pool, Float probFactor, bool enableOffsetManifolds,
          bool enableSpecularMedia, Float avgAngleChangeSurface,
          Float avgAngleChangeMedium) : m_scene(scene),
      m_sampler(sampler), m_pool(pool),
      m_probFactor(probFactor),
      m_enableOffsetManifolds(enableOffsetManifolds),
      m_enableSpecularMedia(enableSpecularMedia) {
    m_manifold = new SpecularManifold(scene);

    if (avgAngleChangeSurface != 0) {
        Log(EInfo, "Using avg. angle change (surface) from configuration: %f", avgAngleChangeSurface);
        m_thetaDiffSurface = degToRad(avgAngleChangeSurface)*DIFF_SAMPLES;
        m_thetaDiffSurfaceSamples = DIFF_SAMPLES;
    } else {
        /* Initial guess: 1 degree change */
        m_thetaDiffSurface = degToRad(1.0f);
        m_thetaDiffSurfaceSamples = 50;
    }

    if (avgAngleChangeMedium != 0) {
        Log(EInfo, "Using avg. angle change (medium) from configuration: %f", avgAngleChangeMedium);
        m_thetaDiffMedium = degToRad(avgAngleChangeMedium)*DIFF_SAMPLES;
        m_thetaDiffMediumSamples = DIFF_SAMPLES;
    } else {
        /* Initial guess: 1 degree change */
        m_thetaDiffMedium = degToRad(1.0f);
        m_thetaDiffMediumSamples = 50;
    }

    m_mediumDensityMultiplier = m_probFactor;
}

ManifoldPerturbation::~ManifoldPerturbation() {
}

Mutator::EMutationType ManifoldPerturbation::getType() const {
    return EManifoldPerturbation;
}

Float ManifoldPerturbation::suitability(const Path &path) const {
    return path.length() >= 4 ? 1.0f : 0.0f;
}

Float ManifoldPerturbation::nonspecularProbSurface(Float alpha) const {
    if (alpha == std::numeric_limits<Float>::infinity())
        return 1.0f;
    else if (!m_enableOffsetManifolds)
        return alpha == 0 ? 0.0f : 1.0f;

    Float q = MTS_MANIFOLD_QUANTILE_SURFACE;
    Float theta_domain = std::atan(-math::fastlog(1-q) * alpha*alpha);
    Float theta_diff = m_thetaDiffSurfaceSamples > 0 ? (m_thetaDiffSurface
            / (Float) m_thetaDiffSurfaceSamples) : (Float) 0.0f;

    return (1-std::cos(theta_domain))
         / (1-std::cos(theta_domain + theta_diff));
}

Float ManifoldPerturbation::nonspecularProbMedium(Float g_) const {
    if (g_ == 0 || !m_enableOffsetManifolds || !m_enableSpecularMedia)
        return 1.0f;

    Float
        g = std::abs(g_),
        q = MTS_MANIFOLD_QUANTILE_MEDIUM,
        t0 = 1+g,
        t1 = 1+g*g,
        t2 = t0 - 2*g*q;

    Float theta_domain = math::safe_acos(
        (t0*t0 - 2*t0*t1*q + 2*g*t1*q*q) / (t2*t2));

    Float theta_diff = m_thetaDiffMediumSamples > 0 ? (m_thetaDiffMedium
            / (Float) m_thetaDiffMediumSamples) : (Float) 0.0f;
    Float theta_newdomain = std::min(M_PI, theta_domain + theta_diff);

    return (1-std::cos(theta_domain))
         / (1-std::cos(theta_newdomain));
}

Float ManifoldPerturbation::nonspecularProb(const PathVertex *vertex) const {
    if (!vertex->isConnectable())
        return 0.0f;

    if (vertex->isSurfaceInteraction()) {
        const Intersection &its = vertex->getIntersection();
        const BSDF *bsdf = its.getBSDF();
        Float nonspecProb = 0;
        int nonspecProbSamples = 0;
        for (int i=0; i<bsdf->getComponentCount(); ++i) {
            if (bsdf->getType(i) & BSDF::ESmooth) {
                nonspecProb += nonspecularProbSurface(bsdf->getRoughness(its, i));
                nonspecProbSamples++;
            }
        }
        BDAssert(nonspecProbSamples > 0);
        if (nonspecProbSamples > 1)
            nonspecProb /= nonspecProbSamples;
        return nonspecProb;
    } else if (vertex->isMediumInteraction()) {
        const MediumSamplingRecord &mRec = vertex->getMediumSamplingRecord();
        return nonspecularProbMedium(mRec.getPhaseFunction()->getMeanCosine());
    } else {
        return 1.0f;
    }
}

int ManifoldPerturbation::getSpecularChainEnd(const Path &path, int pos, int step) {
    while (true) {
        if (pos < 0 || pos > path.length())
            return -1;

        const PathVertex *vertex = path.vertex(pos);
        Float prob = nonspecularProb(vertex);

        if (vertex->isSurfaceInteraction() && vertex->isConnectable())
            statsRoughSurfaceSpecular.incrementBase();
        else if (vertex->isMediumInteraction())
            statsRoughMediumSpecular.incrementBase();

        if (prob == 1 || (prob > 0 && m_sampler->next1D() <= prob)) {
            break;
        } else {
            if (vertex->isSurfaceInteraction() && vertex->isConnectable())
                ++statsRoughSurfaceSpecular;
            else if (vertex->isMediumInteraction())
                ++statsRoughMediumSpecular;
        }

        pos += step;
    }

    return pos;
}

bool ManifoldPerturbation::sampleMutationRecord(
        const Path &source, int &a, int &b, int &c, int &step) {
    int k = source.length();
    Float sample = m_sampler->next1D();
    a = -1;

    if (source.vertex(k-1)->isConnectable()) {
        /* Extra optimization: slightly prefer perturbations from the sensor */
        #define SENSOR_PROB (Float) 0.25f

        if (sample < SENSOR_PROB) {
            a = k-1;
            step = -1;
        } else {
            sample = (sample - SENSOR_PROB) * (1 / (1-SENSOR_PROB));
        }
    }

    if (a < 0) {
        step = sample < 0.5f ? 1 : -1;

        /* Sample the starting vertex of a subpath perturbation */
        a = std::min((int) ((k+1) * m_sampler->next1D()), k);

        /* Probabilistically treat as non-specular */
        Float nonspecProb = nonspecularProb(source.vertex(a));
        if (nonspecProb == 0 || m_sampler->next1D() > nonspecProb) {
            /* Don't start perturbations at specular vertices */
            return false;
        }
    }

    if ((b = getSpecularChainEnd(source, a + step, step)) == -1)
        return false;

    if ((c = getSpecularChainEnd(source, b + step, step)) == -1)
        return false;

    return true;
}

bool ManifoldPerturbation::sampleMutation(
        Path &source, Path &proposal, MutationRecord &muRec, const MutationRecord& sourceMuRec) {
    int k = source.length();

    int a, b, c, step, tries = 0;
    while (!sampleMutationRecord(source, a, b, c, step)) {
        if (tries++ > 1000) {
            SLog(EWarn, "Internal error -- can't decide on a mutation strategy!");
            return false;
        }
    }

    ETransportMode mode = (step == 1) ? EImportance : ERadiance;
    int l = std::min(a, c), m = std::max(a, c);
    int q = std::min(b, b+step);

    if (mode == EImportance) {
        statsAcceptedImp.incrementBase();
        statsGeneratedImp.incrementBase();
    } else {
        statsAcceptedRad.incrementBase();
        statsGeneratedRad.incrementBase();
    }

    muRec = MutationRecord(EManifoldPerturbation, l, m, m-l,
        source.getPrefixSuffixWeight(l, m));

    #if MTS_MANIFOLD_DEBUG == 1
        cout << "Sampled manifold perturbation " << a << " -> " << b << " -> " << c
            << " (k=" << k << ") for path " << source.summarize() << endl;
    #endif

    muRec.extra[0] = a;
    muRec.extra[1] = b;
    muRec.extra[2] = c;
    muRec.extra[3] = step;
    muRec.extra[4] = mode;

    /* Allocate memory for the proposed path */
    proposal.clear();
    proposal.append(source, 0, l+1);
    proposal.append(m_pool.allocEdge());
    for (int i=l+1; i<m; ++i) {
        proposal.append(m_pool.allocVertex());
        proposal.append(m_pool.allocEdge());
    }
    proposal.append(source, m, k+1);

    proposal.vertex(a) = proposal.vertex(a)->clone(m_pool);
    proposal.vertex(c) = proposal.vertex(c)->clone(m_pool);

    const PathVertex
        *vb_old = source.vertex(b),
        *vb_new = proposal.vertex(b);

    if (a != 0 && a != k) {
        /* Sample the first vertex */
        const PathVertex
            *pred_old     = source.vertex(a-step),
            *vertex_old   = source.vertex(a),
            *succ_old     = source.vertex(a+step);
        const PathEdge
            *succEdge_old = source.edge(mode == EImportance ? a : a-1);
        PathVertex
            *pred         = proposal.vertex(a-step),
            *vertex       = proposal.vertex(a),
            *succ         = proposal.vertex(a+step);
        PathEdge
            *predEdge     = proposal.edge(mode == EImportance ? a-step : a-1-step),
            *succEdge     = proposal.edge(mode == EImportance ? a : a-1);

        Float prob_old = std::max(INV_FOURPI, vertex_old->evalPdf(
                m_scene, pred_old, succ_old, mode, ESolidAngle));

        VonMisesFisherDistr vMF(
            VonMisesFisherDistr::forPeakValue(prob_old * m_probFactor * m_probFactor));

        Vector sampled = vMF.sample(m_sampler->next2D());
        Vector wo_old = normalize(succ_old->getPosition()
                - vertex_old->getPosition());
        Vector wo_new = Frame(wo_old).toWorld(sampled);

        if (!vertex->perturbDirection(m_scene,
                pred, predEdge, succEdge, succ, wo_new,
                succEdge_old->length, succ_old->getType(), mode)) {
            goto fail;
        }
    } else {
        const PathVertex *vertex_old = source.vertex(a);
        const PathVertex *succ_old = source.vertex(a+step);
        if (!succ_old->isConnectable())
            goto fail;

        PathVertex *vertex = proposal.vertex(a);
        PathVertex *succ = proposal.vertex(a+step);
        const PathEdge *succEdge_old = source.edge(mode == EImportance ? a : a-1);
        PathEdge *succEdge = proposal.edge(mode == EImportance ? a : a-1);

        *succ = *succ_old;
        *succEdge = *succEdge_old;

        Float pdf = vertex_old->pdf[mode] * m_probFactor * m_probFactor;
        Float stddev = 1.0f / std::sqrt(2*M_PI * pdf);
        if (!succ->perturbPosition(m_scene, m_sampler, stddev))
            goto fail;

        vertex->update(m_scene, NULL, succ, mode, (EMeasure) succ_old->measure);
    }

    /* Generate subsequent vertices between a .. b deterministically */
    for (int i = a + step; i != b; i += step) {
        const PathVertex
            *pred_old     = source.vertex(i-step),
            *vertex_old   = source.vertex(i),
            *succ_old     = source.vertex(i+step);
        const PathEdge
            *succEdge_old = source.edge(mode == EImportance ? i : i-1);
        PathVertex
            *pred         = proposal.vertex(i-step),
            *vertex       = proposal.vertex(i),
            *succ         = proposal.vertex(i+step);
        PathEdge
            *predEdge     = proposal.edge(mode == EImportance ? i-step : i-1-step),
            *succEdge     = proposal.edge(mode == EImportance ? i : i-1);

        if (vertex_old->isSurfaceInteraction()) {
            const Intersection
                &its_old = vertex_old->getIntersection(),
                &its_new = vertex->getIntersection();

            Vector
                wi_old = its_old.toLocal(normalize(pred_old->getPosition() - its_old.p)),
                wo_old = its_old.toLocal(normalize(succ_old->getPosition() - its_old.p));

            bool reflection = Frame::cosTheta(wi_old) * Frame::cosTheta(wo_old) > 0;
            Float eta = vertex_old->getIntersection().getBSDF()->getEta();
            Vector wi_world = normalize(pred->getPosition() - vertex->getPosition()),
                wo_world(0.0f);

            /// todo: this is perhaps a bit drastic
            if (its_old.getBSDF() != its_new.getBSDF())
                goto fail;

            if (vertex_old->isConnectable()) {
                Vector m(0.0f);
                if (reflection)
                    m = normalize(wi_old + wo_old);
                else if (eta != 1)
                    m = normalize(wi_old.z < 0 ? (wi_old*eta + wo_old)
                        : (wi_old + wo_old*eta));
                m = its_new.toWorld(m.z > 0 ? m : -m);

                if (reflection) {
                    wo_world = reflect(wi_world, m);
                } else {
                    if (eta != 1) {
                        wo_world = refract(wi_world, m, eta);
                        if (wo_world.isZero())
                            goto fail;
                    } else {
                        wo_world = -wi_world;
                    }
                }

                Float dist = succEdge_old->length;
                if (i+step == b && succ_old->isMediumInteraction())
                    dist += perturbMediumDistance(m_sampler, succ_old);

                if (!vertex->perturbDirection(m_scene,
                        pred, predEdge, succEdge, succ, wo_world,
                        dist, succ_old->getType(), mode)) {
                    goto fail;
                }
            } else {
                int component = reflection ? BSDF::EDeltaReflection :
                    (BSDF::EDeltaTransmission | BSDF::ENull);

                Float dist = succEdge_old->length;
                if (i+step == b && succ_old->isMediumInteraction())
                    dist += perturbMediumDistance(m_sampler, succ_old);

                if (!vertex->propagatePerturbation(m_scene,
                        pred, predEdge, succEdge, succ, component, dist,
                        succ_old->getType(), mode)) {
                    goto fail;
                }
            }
        } else if (vertex_old->isMediumInteraction()) {
            Point p_old = vertex_old->getPosition(),
                  p_new = vertex->getPosition();

            Normal
                n_old(normalize(p_old - pred_old->getPosition())),
                n_new(normalize(p_new - pred->getPosition()));

            Vector
                dpdu_old = Vector(p_old) - dot(Vector(p_old), n_old) * n_old,
                dpdu_new = Vector(p_new) - dot(Vector(p_new), n_new) * n_new;

            Vector dpdv_old, dpdv_new, wo_old, wo_new;
            Float cosTheta, cosPhi, sinPhi;

            if (dpdu_old.isZero() || dpdu_new.isZero())
                goto fail;

            dpdu_old = normalize(dpdu_old);
            dpdu_new = normalize(dpdu_new);
            dpdv_old = cross(n_old, dpdu_old);
            dpdv_new = cross(n_new, dpdu_new);

            wo_old = normalize(succ_old->getPosition() - p_old);

            cosTheta = dot(wo_old, n_old);

            Float dTheta = warp::squareToStdNormal(m_sampler->next2D()).x
                * 0.5f * M_PI / m_probFactor;
            math::sincos(dTheta, &sinPhi, &cosPhi);

            Float x = dot(wo_old, dpdu_old), y = dot(wo_old, dpdv_old);
            Float x_new = x * cosPhi - y*sinPhi,
                  y_new = x * sinPhi + y*cosPhi;

            wo_new = dpdu_new * x_new + dpdv_new * y_new + n_new * cosTheta;

            Float dist = succEdge_old->length;
            if (i+step == b && succ_old->isMediumInteraction())
                dist += perturbMediumDistance(m_sampler, succ_old);

            if (!vertex->perturbDirection(m_scene,
                    pred, predEdge, succEdge, succ, wo_new,
                    dist, succ_old->getType(), mode))
                goto fail;
        } else {
            Log(EError, "Unsupported vertex type!");
        }
    }

    if (!vb_new->isConnectable())
        goto fail;

    statsUsedManifold.incrementBase();
    statsNonReversible.incrementBase();
    if (std::abs(b-c) > 1) {
        /* Choose a local parameterization of the specular manifold using
           a plane with the following normal (which is computed in a
           reversible manner) */

        Point p0;
        Normal n, n1, n2;
        ++statsUsedManifold;

        if (!vb_old->isMediumInteraction()) {
            n1 = vb_old->getGeometricNormal();
            n2 = vb_new->getGeometricNormal();
        } else {
            n1 = normalize(vb_old->getPosition() - source.vertex(b-step)->getPosition());
            n2 = normalize(vb_new->getPosition() - proposal.vertex(b-step)->getPosition());
        }

        Vector rel = vb_new->getPosition() - vb_old->getPosition();
        Float len = rel.length();
        if (len == 0)
            goto fail;
        rel /= len;

        if (dot(n1, n2) < 0)
            n1 = -n1;
        n = n1 + n2;
        n = n - dot(rel, n)*rel;
        len = n.length();
        if (len == 0)
            goto fail;
        n /= len;

        if (!m_manifold->init(source, c, b))
            goto fail;
        p0 = m_manifold->getPosition(1);
        if (!m_manifold->move(vb_new->getPosition(), n))
            goto fail;
        if (!m_manifold->update(proposal, c, b))
            goto fail;
        if (!m_manifold->move(vb_old->getPosition(), n)) {
            ++statsNonReversible;
            goto fail;
        }

        Point p1 = m_manifold->getPosition(1);
        Float relerr = (p0-p1).length() / std::max(std::max(std::abs(p0.x),
            std::abs(p0.y)), std::abs(p0.z));
        if (relerr > ShadowEpsilon) {
            ++statsNonReversible;
            goto fail;
        }
    }

    if (((vb_old->isSurfaceInteraction() && m_thetaDiffSurfaceSamples < DIFF_SAMPLES) ||
        (vb_old->isMediumInteraction() && m_thetaDiffMediumSamples < DIFF_SAMPLES)) && b+1 != k && b-1 != 0) {
        LockGuard guard(m_thetaDiffMutex);

        if ((vb_old->isSurfaceInteraction() && m_thetaDiffSurfaceSamples < DIFF_SAMPLES) ||
            (vb_old->isMediumInteraction() && m_thetaDiffMediumSamples < DIFF_SAMPLES)) {
            /* Compute the half direction-vector change */
            const PathVertex
                *pred_old     = source.vertex(b-step),
                *vertex_old   = source.vertex(b),
                *succ_old     = source.vertex(b+step);

            const PathVertex
                *pred_new     = proposal.vertex(b-step),
                *vertex_new   = proposal.vertex(b),
                *succ_new     = proposal.vertex(b+step);

            Vector wi_old = normalize(pred_old->getPosition() - vertex_old->getPosition());
            Vector wo_old = normalize(succ_old->getPosition() - vertex_old->getPosition());
            Vector wi_new = normalize(pred_new->getPosition() - vertex_new->getPosition());
            Vector wo_new = normalize(succ_new->getPosition() - vertex_new->getPosition());

            if (vb_old->isSurfaceInteraction()) {
                const BSDF *bsdf_old = vertex_old->getIntersection().getBSDF();
                Vector n_old = vertex_old->getShadingNormal();
                bool reflection = dot(wi_old, n_old) * dot(wo_old, n_old) > 0;
                Vector m_old(0.0f), m_new(0.0f);

                if (reflection) {
                    m_old = wi_old + wo_old;
                    m_new = wi_new + wo_new;
                } else {
                    Float eta = bsdf_old->getEta();
                    if (eta != 1) {
                        if (dot(wi_old, n_old) < 0)
                            eta = 1/eta;

                        m_old = wi_old + wo_old * eta;
                        m_new = wi_new + wo_new * eta;
                    }
                }

                if (!m_old.isZero()) {
                    m_thetaDiffSurface += unitAngle(normalize(m_old), normalize(m_new));
                    m_thetaDiffSurfaceSamples++;
                }

                if (m_thetaDiffSurfaceSamples == DIFF_SAMPLES)
                    Log(EDebug, "Average angle change (surface): %f, p(.2)=%f, p(.1)=%f, p(.01)=%f",
                        radToDeg(m_thetaDiffSurface/m_thetaDiffSurfaceSamples),
                        nonspecularProbSurface(0.2f),
                        nonspecularProbSurface(0.1f),
                        nonspecularProbSurface(0.01f));
            } else {
                m_thetaDiffMedium += std::abs(unitAngle(wi_old, wo_old) - unitAngle(wi_new, wo_new));
                m_thetaDiffMediumSamples++;

                if (m_thetaDiffMediumSamples == DIFF_SAMPLES)
                    Log(EDebug, "Average angle change (medium): %f, p(.9)=%f, p(.99)=%f",
                        radToDeg(m_thetaDiffMedium/m_thetaDiffMediumSamples),
                        nonspecularProbMedium(0.9f),
                        nonspecularProbMedium(0.99f));
            }
        }
    }

    if (!PathVertex::connect(m_scene,
            proposal.vertexOrNull(q-1),
            proposal.edgeOrNull(q-1),
            proposal.vertex(q),
            proposal.edge(q),
            proposal.vertex(q+1),
            proposal.edgeOrNull(q+1),
            proposal.vertexOrNull(q+2),
            source.vertex(q)->isConnectable() ? EArea : EDiscrete,
            source.vertex(q+1)->isConnectable() ? EArea : EDiscrete))
        goto fail;

    if (m >= k-1)
        proposal.vertex(k-1)->updateSamplePosition(
            proposal.vertex(k-2));
    BDAssert(source.matchesConfiguration(proposal));

    if (mode == EImportance)
        ++statsGeneratedImp;
    else
        ++statsGeneratedRad;
    return true;

fail:
    proposal.release(l, m+1, m_pool);
    return false;
}

Float ManifoldPerturbation::Q(const Path &source, const Path &proposal,
        const MutationRecord &muRec) const {
    int a = muRec.extra[0],
        b = muRec.extra[1],
        c = muRec.extra[2],
        k = source.length(),
        step = muRec.extra[3];

    const PathVertex
        *vb_old = source.vertex(b),
        *vb_new = proposal.vertex(b);

    ETransportMode mode = (ETransportMode) muRec.extra[4];

    Spectrum weight = muRec.weight;

    if (a != 0 && a != k) {
        /* Compute the density of the first vertex */
        const PathVertex
            *pred         = source.vertex(a-step),
            *vertex       = source.vertex(a),
            *succ_old     = source.vertex(a+step),
            *succ_new     = proposal.vertex(a+step);

        Float prob_old = std::max(INV_FOURPI, vertex->evalPdf(
                m_scene.get(), pred, succ_old, mode, ESolidAngle));

        Vector wo_old = normalize(succ_old->getPosition()
                - vertex->getPosition());
        Vector wo_new = normalize(succ_new->getPosition()
                - vertex->getPosition());

        Float dp = dot(wo_old, wo_new);

        VonMisesFisherDistr vMF(
            VonMisesFisherDistr::forPeakValue(prob_old * m_probFactor * m_probFactor));

        /* Compute outgoing density wrt. proj. SA measure */
        Float prob = vMF.eval(dp);
        if (vertex->isOnSurface())
            prob /= absDot(wo_new, vertex->getShadingNormal());

        /* Convert to area density at x_b */
        prob *= m_manifold->G(proposal, a, b);

        if (prob <= RCPOVERFLOW)
            return 0.0f;

#if defined(MTS_DEBUG_FP)
        disableFPExceptions();
#endif

        weight /= prob;

#if defined(MTS_DEBUG_FP)
        enableFPExceptions();
#endif

        /* Catch very low probabilities which round to +inf in the above division operation */
        if (!std::isfinite(weight.average()))
            return 0.0f;
    } else {
        Frame frame(source.vertex(a+step)->getGeometricNormal());

        Float stddev = 1.0f / std::sqrt(2*M_PI *
            source.vertex(a)->pdf[mode] * m_probFactor * m_probFactor);

        Float pdf = source.vertex(a+step)->perturbPositionPdf(proposal.vertex(a+step), stddev);
        if (pdf <= RCPOVERFLOW)
            return 0.0f;

        weight /= pdf;
    }

    weight *=
          m_manifold->multiG(proposal, a, c) *
          m_manifold->det(proposal, a, b, c);

    for (int i = a; i != b; i += step) {
        int l = std::min(i, i+step),
            r = std::max(i, i+step);

        if (!proposal.vertex(i)->isConnectable())
            weight *= proposal.edge(l)->evalCached(proposal.vertex(l),
                proposal.vertex(r), PathEdge::ETransmittance | ((mode == EImportance)
                    ? PathEdge::EValueCosineImp : PathEdge::EValueCosineRad));
        else
            weight *= proposal.edge(l)->evalCached(proposal.vertex(l),
                proposal.vertex(r), PathEdge::ETransmittance | ((mode == EImportance)
                    ? PathEdge::EValueImp : PathEdge::EValueRad));

        if (i != a)
            weight /= specularProb(source.vertex(i));
    }

    Float nonspec = nonspecularProb(vb_old);
    if (nonspec == 0)
        return 0.0f;
    weight /= nonspec;

    for (int i = c; i != b; i -= step) {
        int l = std::min(i, i-step),
            r = std::max(i, i-step);

        if (!proposal.vertex(i)->isConnectable())
            weight *= proposal.edge(l)->evalCached(proposal.vertex(l),
                proposal.vertex(r), PathEdge::ETransmittance | ((mode == EImportance)
                    ? PathEdge::EValueCosineRad : PathEdge::EValueCosineImp));
        else
            weight *= proposal.edge(l)->evalCached(proposal.vertex(l),
                proposal.vertex(r), PathEdge::ETransmittance | ((mode == EImportance)
                    ? PathEdge::EValueRad : PathEdge::EValueImp));

        if (i == c)
            weight /= nonspecularProb(source.vertex(i));
        else
            weight /= specularProb(source.vertex(i));
    }

    if (mode == EImportance)
        weight *= proposal.edge(b)->evalCached(
                vb_new, proposal.vertex(b+1), PathEdge::EValueImp);
    else
        weight *= proposal.edge(b-1)->evalCached(
                proposal.vertex(b-1), vb_new, PathEdge::EValueRad);

    if (vb_old->isMediumInteraction()) {
        int vbEdge = std::min(b, b-step);
        weight /= pdfMediumPerturbation(vb_old,
            source.edge(vbEdge), proposal.edge(vbEdge));
    }

    Float lum = weight.getLuminance();

    if (lum <= RCPOVERFLOW || !std::isfinite(lum)) {
        Log(EWarn, "Internal error in manifold perturbation: luminance = %f!", lum);
        return 0.f;
    }

    return 1.0f / lum;
}

void ManifoldPerturbation::accept(const MutationRecord &muRec) {
    ETransportMode mode = (ETransportMode) muRec.extra[4];
    if (mode == EImportance)
        ++statsAcceptedImp;
    else
        ++statsAcceptedRad;

}

MTS_IMPLEMENT_CLASS(ManifoldPerturbation, false, Mutator)
MTS_NAMESPACE_END

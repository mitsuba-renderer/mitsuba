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

#include <mitsuba/bidir/manifold.h>
#include <mitsuba/bidir/path.h>
#include <mitsuba/core/statistics.h>
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_NO_DEBUG
#include <Eigen/LU>
#include <Eigen/Geometry>

MTS_NAMESPACE_BEGIN

/* Some statistics counters */
static StatsCounter statsStepFailed(
        "Specular manifold", "Retries (step failed)");
static StatsCounter statsStepTooFar(
        "Specular manifold", "Retries (step increased distance)");
static StatsCounter statsStepSuccess(
        "Specular manifold", "Successful steps");
static StatsCounter statsAvgIterations(
        "Specular manifold", "Avg. iterations per walk", EAverage);
static StatsCounter statsAvgIterationsSuccess(
        "Specular manifold", "Avg. iterations per successful walk", EAverage);
static StatsCounter statsAvgManifoldSize(
        "Specular manifold", "Avg. manifold size", EAverage);
static StatsCounter statsSuccessfulWalks(
        "Specular manifold", "Successful walks", EPercentage);
static StatsCounter statsMediumSuccess(
        "Specular manifold", "Successful walks w/ media", EPercentage);
static StatsCounter statsNonManifold(
        "Specular manifold", "Non-manifold", EPercentage);
static StatsCounter statsUpdateFailed(
        "Specular manifold", "Update failed");
static StatsCounter statsMaxManifold(
        "Specular manifold", "Max. manifold size", EMaximumValue);

SpecularManifold::SpecularManifold(const Scene *scene, int maxIterations)
  : m_scene(scene) {
    m_maxIterations = maxIterations > 0 ? maxIterations :
        MTS_MANIFOLD_MAX_ITERATIONS;
}

bool SpecularManifold::init(const Path &path, int start, int end) {
    int step = start < end ? 1 : -1;
    if (path.vertex(start)->isSupernode())
        start += step;
    if (path.vertex(end)->isSupernode())
        end -= step;

    const PathVertex
        *vs = path.vertex(start),
        *ve = path.vertex(end);

    /* Create the initial vertex that is pinned in position by default */
    SimpleVertex v(EPinnedPosition, vs->getPosition());

    /* When the endpoint is on an orthographic camera or directional light
       source, switch to a directionally pinned vertex instead */
    if (vs->getType() & (PathVertex::ESensorSample | PathVertex::EEmitterSample)) {
        const PositionSamplingRecord &pRec
            = vs->getPositionSamplingRecord();
        uint32_t type = static_cast<const AbstractEmitter *>(pRec.object)->getType()
            & (AbstractEmitter::EDeltaDirection | AbstractEmitter::EDeltaPosition);
        if (type == AbstractEmitter::EDeltaDirection) {
            v.type = EPinnedDirection;
            v.gn = v.n = pRec.n;
            coordinateSystem(pRec.n, v.dpdu, v.dpdv);
        }
    }

    m_time = vs->getTime();
    m_vertices.clear();
    m_vertices.push_back(v);

    for (int i=start + step; i != end; i += step) {
        const PathVertex
            *pred = path.vertex(i-step),
            *vertex = path.vertex(i),
            *succ = path.vertex(i+step);

        if (vertex->isSurfaceInteraction()) {
            const Intersection &its = vertex->getIntersection();
            const BSDF *bsdf = its.getBSDF();

            v.p = its.p;
            v.gn = its.geoFrame.n;
            v.n = its.shFrame.n;
            v.dpdu = its.dpdu;
            v.dpdv = its.dpdv;
            v.object = bsdf;
            v.degenerate = !vertex->isConnectable();
            const Shape *shape = its.instance != NULL ? its.instance : its.shape;
            shape->getNormalDerivative(its, v.dndu, v.dndv);

            /* Turn into an orthonormal parameterization at 'p' */
            Float invLen = 1 / v.dpdu.length();
            v.dpdu *= invLen;
            v.dndu *= invLen;
            Float dp = dot(v.dpdu, v.dpdv);
            Vector dpdv = v.dpdv - dp * v.dpdu;
            Vector dndv = v.dndv - dp * v.dndu;
            invLen = 1 / dpdv.length();
            v.dpdv = dpdv * invLen;
            v.dndv = dndv * invLen;

            Vector wPred = pred->getPosition() - v.p;
            Vector wSucc = succ->getPosition() - v.p;

            if (dot(v.gn, wPred) * dot(v.gn, wSucc) < 0) {
                v.type = ERefraction;
                v.eta = bsdf->getEta();
            } else {
                v.type = EReflection;
                v.eta = 1.0f;
            }
        } else if (vertex->isMediumInteraction()) {
            const MediumSamplingRecord &mRec = vertex->getMediumSamplingRecord();

            Vector wi = pred->getPosition() - mRec.p;
            Float invLength = 1.0f / wi.length();
            wi *= invLength;

            v.p = mRec.p;
            v.gn = v.n = Normal(0.0f);

            Vector s, t;
            coordinateSystem(wi, s, t);

            v.dpdu = s;
            v.dpdv = t;
            v.dndu = s * invLength;
            v.dndv = t * invLength;

            v.object = mRec.getPhaseFunction();
            v.eta = 1.0f;
            v.degenerate = false;
            v.type = EMedium;
        } else {
            Log(EError, "Unknown vertex type! : %s", vertex->toString().c_str());
        }

        m_vertices.push_back(v);
    }

    v = SimpleVertex(EMovable, ve->getPosition());
    m_vertices.push_back(v);

    #if MTS_MANIFOLD_DEBUG == 1
        cout << "==========================================" << endl;
        cout << "Initialized specular manifold: " << toString() << endl;
    #endif

    return true;
}

bool SpecularManifold::computeTangents() {
    const int n = static_cast<int>(m_vertices.size() - 1);

    m_vertices[0].Tp.setZero();
    m_vertices[m_vertices.size()-1].Tp.setIdentity();

    if (m_vertices.size() == 2) /* Nothing to do */
        return true;

    /* Matrix assembly stage */
    for (int i=0; i<n; ++i) {
        SimpleVertex *v = &m_vertices[i];

        Vector wo = v[1].p - v[0].p;
        Float ilo = wo.length();

        if (ilo == 0)
            return false;
        ilo = 1/ilo; wo *= ilo;

        if (v[0].type == EPinnedPosition) {
            v[0].a.setZero();
            v[0].b.setIdentity();
            v[0].c.setZero();
            continue;
        } else if (v[0].type == EPinnedDirection) {

            Vector dC_dnext_u = (v[1].dpdu - wo * dot(wo, v[1].dpdu)) * ilo;
            Vector dC_dnext_v = (v[1].dpdv - wo * dot(wo, v[1].dpdv)) * ilo;
            Vector dC_dcur_u = (wo * dot(wo, v[0].dpdu) - v[0].dpdu) * ilo;
            Vector dC_dcur_v = (wo * dot(wo, v[0].dpdv) - v[0].dpdv) * ilo;

            v[0].a.setZero();
            v[0].b = Matrix2x2(
                Vector2(dot(dC_dcur_u, v[0].dpdu), dot(dC_dcur_u, v[0].dpdv)),
                Vector2(dot(dC_dcur_v, v[0].dpdu), dot(dC_dcur_v, v[0].dpdv))
            );
            v[0].c = Matrix2x2(
                Vector2(dot(dC_dnext_u, v[0].dpdu), dot(dC_dnext_u, v[0].dpdv)),
                Vector2(dot(dC_dnext_v, v[0].dpdu), dot(dC_dnext_v, v[0].dpdv))
            );
            continue;
        }

        Vector wi = v[-1].p - v[0].p;
        Float ili = wi.length();

        if (ili == 0)
            return false;

        ili = 1/ili; wi *= ili;

        if (v[0].type == EReflection || v[0].type == ERefraction) {
            Float eta = v[0].eta;
            bool normalizeH = !(v[0].type == ERefraction && eta == 1);

            /* Compute the half vector and a few useful projections */
            Vector H;
            Float ilh;
            if (normalizeH) {
                /* Generally compute derivatives with respect to the normalized
                   half-vector. When given an index-matched refraction event,
                   don't perform this normalization, since the desired vertex
                   configuration is actually where H = 0. */

                if (dot(wi, v[0].gn) < 0)
                    eta = 1 / eta;

                H = wi + eta * wo;
                ilh = 1 / H.length();
                H *= ilh;
            } else {
                H = wi + wo;
                ilh = 1.0f;
            }

            /* Orient the half-vector so that it points in the same
               hemisphere as the geometric surface normal */

            Float dot_H_n    = dot(v[0].n, H),
                  dot_H_dndu = dot(v[0].dndu, H),
                  dot_H_dndv = dot(v[0].dndv, H),
                  dot_u_n    = dot(v[0].dpdu, v[0].n),
                  dot_v_n    = dot(v[0].dpdv, v[0].n);

            /* Local shading tangent frame */
            Vector s = v[0].dpdu - dot_u_n * v[0].n;
            Vector t = v[0].dpdv - dot_v_n * v[0].n;

            ilo *= eta * ilh; ili *= ilh;

            /* Derivatives of C with respect to x_{i-1} */
            Vector
                dH_du = (v[-1].dpdu - wi * dot(wi, v[-1].dpdu)) * ili,
                dH_dv = (v[-1].dpdv - wi * dot(wi, v[-1].dpdv)) * ili;

            if (normalizeH) {
                dH_du -= H * dot(dH_du, H);
                dH_dv -= H * dot(dH_dv, H);
            }

            v[0].a = Matrix2x2(
                dot(dH_du, s), dot(dH_dv, s),
                dot(dH_du, t), dot(dH_dv, t));

            /* Derivatives of C with respect to x_i */
            dH_du = -v[0].dpdu * (ili + ilo) + wi * (dot(wi, v[0].dpdu) * ili)
                                             + wo * (dot(wo, v[0].dpdu) * ilo);
            dH_dv = -v[0].dpdv * (ili + ilo) + wi * (dot(wi, v[0].dpdv) * ili)
                                             + wo * (dot(wo, v[0].dpdv) * ilo);

            if (normalizeH) {
                dH_du -= H * dot(dH_du, H);
                dH_dv -= H * dot(dH_dv, H);
            }

            v[0].b = Matrix2x2(
                dot(dH_du, s) - dot(v[0].dpdu, v[0].dndu) * dot_H_n - dot_u_n * dot_H_dndu,
                dot(dH_dv, s) - dot(v[0].dpdu, v[0].dndv) * dot_H_n - dot_u_n * dot_H_dndv,
                dot(dH_du, t) - dot(v[0].dpdv, v[0].dndu) * dot_H_n - dot_v_n * dot_H_dndu,
                dot(dH_dv, t) - dot(v[0].dpdv, v[0].dndv) * dot_H_n - dot_v_n * dot_H_dndv);

            /* Derivatives of C with respect to x_{i+1} */
            dH_du = (v[1].dpdu - wo * dot(wo, v[1].dpdu)) * ilo;
            dH_dv = (v[1].dpdv - wo * dot(wo, v[1].dpdv)) * ilo;

            if (normalizeH) {
                dH_du -= H * dot(dH_du, H);
                dH_dv -= H * dot(dH_dv, H);
            }

            v[0].c = Matrix2x2(
                dot(dH_du, s), dot(dH_dv, s),
                dot(dH_du, t), dot(dH_dv, t));

            /* Store the microfacet normal wrt. the local (orthonormal) shading frame */
            s = normalize(s);
            t = cross(v[0].n, s);
            v[0].m = Vector(dot(s, H), dot(t, H), dot(v[0].n, H));
            if (dot(H, v[0].gn) < 0)
                v[0].m = -v[0].m;
        } else if (v[0].type == EMedium) {
            Vector dwi_dpred_u = (v[-1].dpdu - wi * dot(wi, v[-1].dpdu)) * ili;
            Vector dwi_dpred_v = (v[-1].dpdv - wi * dot(wi, v[-1].dpdv)) * ili;
            Vector dwi_dcur_u  = (-v[0].dpdu + wi * dot(wi, v[ 0].dpdu)) * ili;
            Vector dwi_dcur_v  = (-v[0].dpdv + wi * dot(wi, v[ 0].dpdv)) * ili;

            Vector t, dt_dpred_u, dt_dpred_v, dt_dcur_u, dt_dcur_v;

            /* Compute the local frame and derivatives thereof */
            if (std::abs(wi.x) > std::abs(wi.y)) {
                Float tl = 1.0f / std::sqrt(wi.x * wi.x + wi.z * wi.z);
                t = Vector(wi.z * tl, 0.0f, -wi.x * tl);

                dt_dpred_u = Vector(dwi_dpred_u.z*tl, 0.0f, -dwi_dpred_u.x*tl);
                dt_dpred_v = Vector(dwi_dpred_v.z*tl, 0.0f, -dwi_dpred_v.x*tl);
                dt_dcur_u  = Vector(dwi_dcur_u.z*tl,  0.0f, -dwi_dcur_u.x*tl);
                dt_dcur_v  = Vector(dwi_dcur_v.z*tl,  0.0f, -dwi_dcur_v.x*tl);
            } else {
                Float tl = 1.0f / std::sqrt(wi.y * wi.y + wi.z * wi.z);
                t = Vector(0.0f, wi.z * tl, -wi.y * tl);

                dt_dpred_u = Vector(0.0f, dwi_dpred_u.z*tl, -dwi_dpred_u.y*tl);
                dt_dpred_v = Vector(0.0f, dwi_dpred_v.z*tl, -dwi_dpred_v.y*tl);
                dt_dcur_u  = Vector(0.0f, dwi_dcur_u.z*tl,  -dwi_dcur_u.y*tl);
                dt_dcur_v  = Vector(0.0f, dwi_dcur_v.z*tl,  -dwi_dcur_v.y*tl);
            }

            dt_dpred_u -= t * dot(t, dt_dpred_u);
            dt_dpred_v -= t * dot(t, dt_dpred_v);
            dt_dcur_u  -= t * dot(t, dt_dcur_u);
            dt_dcur_v  -= t * dot(t, dt_dcur_v);

            Vector s = cross(t, wi);
            Vector ds_dpred_u = cross(dt_dpred_u, wi) + cross(t, dwi_dpred_u);
            Vector ds_dpred_v = cross(dt_dpred_v, wi) + cross(t, dwi_dpred_v);
            Vector ds_dcur_u  = cross(dt_dcur_u, wi)  + cross(t, dwi_dcur_u);
            Vector ds_dcur_v  = cross(dt_dcur_v, wi)  + cross(t, dwi_dcur_v);

            /* Some tangential projections */
            Vector2
                t_cur_dpdu (dot(v[ 0].dpdu, s), dot(v[ 0].dpdu, t)),
                t_cur_dpdv (dot(v[ 0].dpdv, s), dot(v[ 0].dpdv, t)),
                t_next_dpdu(dot(v[ 1].dpdu, s), dot(v[ 1].dpdu, t)),
                t_next_dpdv(dot(v[ 1].dpdv, s), dot(v[ 1].dpdv, t)),
                t_wo = Vector2(dot(wo, s), dot(wo, t));

            v[0].a = Matrix2x2(
                Vector2(dot(ds_dpred_u, wo), dot(dt_dpred_u, wo)),
                Vector2(dot(ds_dpred_v, wo), dot(dt_dpred_v, wo))
            );

            v[0].b = Matrix2x2(
                (t_wo * dot(wo, v[0].dpdu) - t_cur_dpdu) * ilo +
                Vector2(dot(ds_dcur_u, wo), dot(dt_dcur_u, wo)),
                (t_wo * dot(wo, v[0].dpdv) - t_cur_dpdv) * ilo +
                Vector2(dot(ds_dcur_v, wo), dot(dt_dcur_v, wo)));

            v[0].c = Matrix2x2(
                (t_next_dpdu - t_wo * dot(wo, v[1].dpdu)) * ilo,
                (t_next_dpdv - t_wo * dot(wo, v[1].dpdv)) * ilo);

            v[0].m = Vector(dot(s, wo), dot(t, wo), dot(wi, wo));
        } else {
            Log(EError, "Unknown vertex type!");
        }
    }

    /* Find the tangent space with respect to translation of the last
       vertex. For this, we must solve a tridiagonal system. The following is
       simplified version of the block tridiagonal LU factorization algorithm
       for this specific problem */
    Matrix2x2 Li;
    if (!m_vertices[0].b.invert(Li))
        return false;

    for (int i=0; i < n - 1; ++i) {
        m_vertices[i].u = Li * m_vertices[i].c;
        Matrix2x2 temp = m_vertices[i+1].b - m_vertices[i+1].a * m_vertices[i].u;
        if (!temp.invert(Li))
            return false;
    }

    m_vertices[n-1].Tp = -Li * m_vertices[n-1].c;

    for (int i=n-2; i>=0; --i)
        m_vertices[i].Tp = -m_vertices[i].u * m_vertices[i+1].Tp;
    return true;
}

bool SpecularManifold::project(const Vector &d) {
    const SimpleVertex &last = m_vertices[m_vertices.size()-1];
    Float du = dot(d, last.dpdu), dv = dot(d, last.dpdv);

    Ray ray(Point(0.0f), Vector(1.0f), 0); // make gcc happy
    Intersection its;

    m_proposal.clear();
    for (size_t i=0; i<m_vertices.size(); ++i) {
        m_proposal.push_back(m_vertices[i]);
        SimpleVertex &vertex = m_proposal[i];

        if (i == 0) {
            Point p0 = m_vertices[0].p + m_vertices[0].map(du, dv);
            Point p1 = m_vertices[1].p + m_vertices[1].map(du, dv);

            ray = Ray(p0, normalize(p1 - p0), m_time);
            vertex.p = ray.o;
            continue;
        } else if (vertex.type == EMovable) {
            Float dp = dot(ray.d, vertex.n);
            if (std::abs(dp) < Epsilon)
                return false;

            Float t = dot(vertex.p - ray.o, vertex.n) / dp;
            vertex.p = ray(t);
            break;
        } else if (vertex.type == EReflection) {
            if (!m_scene->rayIntersect(ray, its))
                return false;

            Vector n = its.shFrame.n,
                   s = its.dpdu, t;
            s = normalize(s - n * dot(n, s));
            t = cross(n, s);

            Vector m = s * vertex.m[0] + t * vertex.m[1] + n * vertex.m[2];

            ray.setOrigin(its.p);
            ray.setDirection(reflect(-ray.d, m));
        } else if (vertex.type == ERefraction) {
            if (!m_scene->rayIntersect(ray, its))
                return false;

            Vector n = its.shFrame.n,
                   s = its.dpdu, t;
            s = normalize(s - n * dot(n, s));
            t = cross(n, s);

            Vector m = s * vertex.m[0] + t * vertex.m[1] + n * vertex.m[2];
            Vector refracted = refract(-ray.d, m, its.shape->getBSDF()->getEta());

            if (refracted.isZero())
                return false;

            ray.setOrigin(its.p);
            ray.setDirection(refracted);
        } else if (vertex.type == EMedium) {
            Float length = (m_vertices[i].p - m_vertices[i-1].p).length(),
                  invLength = 1.0f / length;

            /* Check for occlusion */
            if (m_scene->rayIntersect(Ray(ray, Epsilon, length)))
                return false;

            vertex.p = ray(length);
            vertex.n = Vector(0.0f);

            Vector wi = -ray.d, s, t;
            coordinateSystem(wi, s, t);

            vertex.dpdu = s;
            vertex.dpdv = t;
            vertex.dndu = s * invLength;
            vertex.dndv = t * invLength;

            ray.setOrigin(vertex.p);
            ray.setDirection(s * vertex.m[0] + t * vertex.m[1] + wi * vertex.m[2]);
        } else {
            Log(EError, "Unsupported vertex type!");
        }

        if (vertex.type != EMedium) {
            if (vertex.object != its.shape->getBSDF())
                return false;

            vertex.p = its.p;
            vertex.n = its.shFrame.n;
            vertex.gn = its.geoFrame.n;
            vertex.dpdu = its.dpdu;
            vertex.dpdv = its.dpdv;

            const Shape *shape = its.instance != NULL ? its.instance : its.shape;
            shape->getNormalDerivative(its,
                vertex.dndu, vertex.dndv);

            /* Turn into an orthonormal parameterization at 'p' */
            Float invLen = 1 / vertex.dpdu.length();
            vertex.dpdu *= invLen; vertex.dndu *= invLen;
            Float dp = dot(vertex.dpdu, vertex.dpdv);
            Vector dpdv = vertex.dpdv - dp * vertex.dpdu;
            Vector dndv = vertex.dndv - dp * vertex.dndu;
            invLen = 1 / dpdv.length();
            vertex.dpdv = dpdv * invLen;
            vertex.dndv = dndv * invLen;
        }
    }
    return true;
}

bool SpecularManifold::move(const Point &target, const Normal &n) {
    SimpleVertex &last = m_vertices[m_vertices.size()-1];

    #if MTS_MANIFOLD_DEBUG == 1
        cout << "moveTo(" << last.p.toString() << " => " << target.toString() << ", n=" << n.toString() << ")" << endl;
    #endif

    if (m_vertices.size() == 2 && m_vertices[0].type == EPinnedPosition) {
        /* Nothing to do */
        return true;
    }

    bool medium = false;
    for (size_t i=0; i<m_vertices.size(); ++i) {
        if (m_vertices[i].type == EMedium)
            medium = true;
    }

    if (medium)
        statsMediumSuccess.incrementBase();

    statsAvgManifoldSize.incrementBase();
    statsAvgManifoldSize += m_vertices.size();
    statsMaxManifold.recordMaximum(m_vertices.size());

    statsSuccessfulWalks.incrementBase();

    Float invScale = 1.0f / std::max(std::max(std::abs(target.x),
            std::abs(target.y)), std::abs(target.z));
    Float stepSize = 1;

    BDAssert(last.type == EMovable);
    coordinateSystem(n, last.dpdu, last.dpdv);
    last.n = n;

    m_proposal.reserve(m_vertices.size());
    m_iterations = 0;
    statsAvgIterations.incrementBase();
    while (m_iterations < m_maxIterations) {
        Vector rel = target - m_vertices[m_vertices.size()-1].p;
        Float dist = rel.length(), newDist;
        if (dist * invScale < MTS_MANIFOLD_EPSILON) {
            /* Check for an annoying corner-case where the last
               two vertices converge to the same point (this can
               happen e.g. on rough planar reflectors) */
            dist = (m_vertices[m_vertices.size()-1].p
                  - m_vertices[m_vertices.size()-2].p).length();
            if (dist * invScale < Epsilon) {
                return false;
            }

            /* The manifold walk converged. */
            ++statsSuccessfulWalks;
            statsAvgIterationsSuccess.incrementBase();
            statsAvgIterationsSuccess += m_iterations;
            if (medium)
                ++statsMediumSuccess;
            #if MTS_MANIFOLD_DEBUG == 1
                cout << "move(): converged after " << m_iterations << " iterations" << endl;
                cout << "Final configuration:" << toString() << endl;
            #endif
            return true;
        }
        m_iterations++;
        ++statsAvgIterations;

        /* Compute the tangent vectors for the current path */
        statsNonManifold.incrementBase();
        if (!computeTangents()) {
            ++statsNonManifold;
            #if MTS_MANIFOLD_DEBUG == 1
                cout << "move(): unable to compute tangents!" << endl;
            #endif
            return false;
        }

        /* Take a step using the computed tangents and project
           back on the manifold */
        #if MTS_MANIFOLD_DEBUG == 1
            const SimpleVertex &last = m_vertices[m_vertices.size()-1];
            Float du = dot(rel, last.dpdu), dv = dot(rel, last.dpdv);
            cout << "project(du=" << du << ", dv=" << dv << ", stepSize=" << stepSize << ")" << endl;
        #endif

        if (!project(rel * stepSize)) {
            #if MTS_MANIFOLD_DEBUG == 1
                cout << "project failed!" << endl;
            #endif
            ++statsStepFailed;
            goto failure;
        }

        /* Reject if the step increased the distance */
        newDist = (target - m_proposal[m_proposal.size()-1].p).length();
        #if MTS_MANIFOLD_DEBUG == 1
            cout << "Distance: " << dist << " -> " << newDist << endl;
        #endif
        if (newDist > dist) {
            ++statsStepTooFar;
            #if MTS_MANIFOLD_DEBUG == 1
                cout << "-> Rejecting!" << endl;
            #endif
            goto failure;
        }
        #if MTS_MANIFOLD_DEBUG == 1
            cout << "-> Accepting!" << endl;
        #endif
        ++statsStepSuccess;

        m_proposal.swap(m_vertices);

        /* Increase the step size */
        stepSize = std::min((Float) 1.0f, stepSize * 2.0f);
        continue;
    failure:
        /* Reduce the step size */
        stepSize /= 2.0f;
    }
    #if MTS_MANIFOLD_DEBUG == 1
        cout << "Exceeded the max. iteration count!" << endl;
    #endif

    return false;
}

bool SpecularManifold::update(Path &path, int start, int end) {
    int step;
    ETransportMode mode;

    if (start < end) {
        step = 1; mode = EImportance;
    } else {
        step = -1; mode = ERadiance;
    }

    int last = (int) m_vertices.size() - 2;
    if (m_vertices[0].type == EPinnedDirection)
        last = std::max(last, 1);

    for (int j=0, i=start; j < last; ++j, i += step) {
        const SimpleVertex
            &v = m_vertices[j],
            &vn = m_vertices[j+1];

        PathVertex
            *pred   = path.vertexOrNull(i-step),
            *vertex = path.vertex(i),
            *succ   = path.vertex(i+step);

        int predEdgeIdx = (mode == EImportance) ? i-step : i-step-1;
        PathEdge *predEdge = path.edgeOrNull(predEdgeIdx),
                 *succEdge = path.edge(predEdgeIdx + step);

        Vector d = vn.p - v.p;
        Float length = d.length();
        d /= length;
        PathVertex::EVertexType desiredType = vn.type == EMedium ?
            PathVertex::EMediumInteraction : PathVertex::ESurfaceInteraction;

        if (v.type == EPinnedDirection) {
            /* Create a fake vertex and use it to call sampleDirect(). This is
               kind of terrible -- a nicer API is needed to cleanly support this */
            PathVertex temp;
            temp.type = PathVertex::EMediumInteraction;
            temp.degenerate = false;
            temp.measure = EArea;
            MediumSamplingRecord &mRec = temp.getMediumSamplingRecord();
            mRec.time = m_time;
            mRec.p = vn.p;

            if (temp.sampleDirect(m_scene, NULL, vertex, succEdge, succ, mode).isZero()) {
                #if MTS_MANIFOLD_DEBUG == 1
                    cout << "update(): failed in sampleDirect()!" << endl;
                #endif
                ++statsUpdateFailed;
                return false;
            }

            if (m_vertices.size() >= 3) {
                PathVertex *succ2 = path.vertex(i+2*step);
                PathEdge *succ2Edge = path.edge(predEdgeIdx + 2*step);
                if (!succ->sampleNext(m_scene, NULL, vertex, succEdge, succ2Edge, succ2, mode)) {
                    #if MTS_MANIFOLD_DEBUG == 1
                        cout << "update(): failed in sampleNext() / pinned direction!" << endl;
                    #endif
                    ++statsUpdateFailed;
                    return false;
                }
            }
            i += step;
        } else if (!v.degenerate) {
            if (!vertex->perturbDirection(m_scene,
                    pred, predEdge, succEdge, succ, d,
                    length, desiredType, mode)) {
                #if MTS_MANIFOLD_DEBUG == 1
                    cout << "update(): failed in perturbDirection()" << endl;
                #endif
                ++statsUpdateFailed;
                return false;
            }

            Float relerr = (vn.p - succ->getPosition()).length() /
                std::max(std::max(std::abs(vn.p.x),
                    std::abs(vn.p.y)), std::abs(vn.p.z));

            if (relerr > 1e-3f) {
                // be extra-cautious
                #if MTS_MANIFOLD_DEBUG == 1
                    cout << "update(): failed, relative error of perturbDirection() too high:" << relerr << endl;
                #endif
                ++statsUpdateFailed;
                return false;
            }
        } else {
            unsigned int compType;
            if (v.type == ERefraction)
                compType = v.eta != 1 ? BSDF::EDeltaTransmission : (BSDF::ENull | BSDF::EDeltaTransmission);
            else
                compType = BSDF::EDeltaReflection;

            if (!vertex->propagatePerturbation(m_scene,
                    pred, predEdge, succEdge, succ, compType,
                    length, desiredType, mode)) {
                #if MTS_MANIFOLD_DEBUG == 1
                    cout << "update(): failed in propagatePerturbation()" << endl;
                #endif
                ++statsUpdateFailed;
                return false;
            }

            Float relerr = (vn.p - succ->getPosition()).length() /
                std::max(std::max(std::abs(vn.p.x),
                    std::abs(vn.p.y)), std::abs(vn.p.z));
            if (relerr > 1e-3f) {
                // be extra-cautious
                #if MTS_MANIFOLD_DEBUG == 1
                    cout << "update(): failed, relative error of propagatePerturbation() too high:" << relerr << endl;
                #endif
                ++statsUpdateFailed;
                return false;
            }
        }
    }

    return true;
}

Float SpecularManifold::det(const Path &path, int a, int b, int c) {
    int k = path.length();

    if (a == 0 || a == k)
        std::swap(a, c);

    int step = b > a ? 1 : -1, nGlossy = 0, nSpecular = 0;

    for (int i=a + step; i != c; i += step) {
        if (path.vertex(i)->isConnectable())
            ++nGlossy;
        else
            ++nSpecular;
    }

    if (nGlossy <= 1) /* No glossy materials -- we don't need this derivative */
        return 1.0f;

    bool success = init(path, a, c);
    BDAssert(success);

    int b_idx = std::abs(b-a);
    SimpleVertex &vb = m_vertices[b_idx];
    const PathVertex *pb = path.vertex(b);

    if (pb->isMediumInteraction()) {
        vb.n = Vector(path.edge(a < b ? (b-1) : b)->d);
    } else {
        vb.n = pb->getShadingNormal();
    }
    coordinateSystem(vb.n, vb.dpdu, vb.dpdv);

    if (!computeTangents()) {
        Log(EWarn, "Could not compute tangents!");
        return 0.0f;
    }

    m_vertices[b_idx].a.setZero();
    m_vertices[b_idx].b.setIdentity();
    m_vertices[b_idx].c.setZero();

    if (nSpecular == 0) {
        /* The chain only consists of glossy vertices -- simply compute the
           determinant of the block tridiagonal matrix A.

           See D.K. Salkuyeh, Comments on "A note on a three-term recurrence for a
           tridiagonal matrix", Appl. Math. Comput. 176 (2006) 442-444. */

        Matrix2x2 Di(0.0f), D = m_vertices[1].b;

        Float det = D.det();
        for (size_t i=2; i<m_vertices.size()-1; ++i) {
            if (!D.invert(Di)) {
                Log(EWarn, "Could not invert matrix!");
                return 0.0f;
            }

            D = m_vertices[i].b - m_vertices[i].a * Di * m_vertices[i-1].c;
            det *= D.det();
        }

        return std::abs(1 / det);
    } else {
        /* The chain contains both glossy and specular materials. Compute the
           determinant of A^-1, where rows corresponding to specular vertices
           have been crossed out. The performance of the following is probably
           terrible (lots of dynamic memory allocation), but it works and
           this case happens rarely enough .. */

        Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> A(2*(nGlossy + nSpecular), 2*(nGlossy + nSpecular));
        A.setZero();

        for (int j=0, i=0; j<nGlossy+nSpecular; ++j) {
            if (j-1 >= 0) {
                A(2*i,   2*(j-1))   = m_vertices[j+1].a(0,0);
                A(2*i,   2*(j-1)+1) = m_vertices[j+1].a(0,1);
                A(2*i+1, 2*(j-1))   = m_vertices[j+1].a(1,0);
                A(2*i+1, 2*(j-1)+1) = m_vertices[j+1].a(1,1);
            }
            A(2*i,   2*j)   = m_vertices[j+1].b(0,0);
            A(2*i,   2*j+1) = m_vertices[j+1].b(0,1);
            A(2*i+1, 2*j)   = m_vertices[j+1].b(1,0);
            A(2*i+1, 2*j+1) = m_vertices[j+1].b(1,1);

            if (j+1 < nGlossy + nSpecular) {
                A(2*i,   2*(j+1))   = m_vertices[j+1].c(0,0);
                A(2*i,   2*(j+1)+1) = m_vertices[j+1].c(0,1);
                A(2*i+1, 2*(j+1))   = m_vertices[j+1].c(1,0);
                A(2*i+1, 2*(j+1)+1) = m_vertices[j+1].c(1,1);
            }
            ++i;
        }

        /* Compute the inverse and "cross out" irrelevant columns and rows */
        Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> Ai = A.inverse();

        for (int i=0; i<nGlossy+nSpecular; ++i) {
            if (!m_vertices[i+1].degenerate)
                continue;

            Ai.row(2*i).setZero();
            Ai.col(2*i).setZero();
            Ai.row(2*i+1).setZero();
            Ai.col(2*i+1).setZero();

            Ai.block<2,2>(2*i, 2*i).setIdentity();
        }

        return std::abs(Ai.determinant());
    }
}

Float SpecularManifold::multiG(const Path &path, int a, int b) {
    if (a == 0)
        ++a;
    else if (a == path.length())
        --a;
    if (b == 0)
        ++b;
    else if (b == path.length())
        --b;

    int step = b > a ? 1 : -1;
    while (!path.vertex(b)->isConnectable())
        b -= step;
    while (!path.vertex(a)->isConnectable())
        a += step;

    Float result = 1;

    BDAssert(path.vertex(a)->isConnectable() && path.vertex(b)->isConnectable());
    for (int i = a + step, start = a; i != b + step; i += step) {
        if (path.vertex(i)->isConnectable()) {
            result *= G(path, start, i);
            start = i;
        }
    }

    return result;
}

Float SpecularManifold::G(const Path &path, int a, int b) {
    if (std::abs(a-b) == 1) {
        if (a > b)
            std::swap(a, b);
        return path.edge(a)->evalCached(path.vertex(a),
            path.vertex(b), PathEdge::EGeometricTerm)[0];
    }

    Assert(path.vertex(a)->isConnectable());
    Assert(path.vertex(b)->isConnectable());
    int step = b > a ? 1 : -1;

    bool success = init(path, a, b);
    BDAssert(success);

    SimpleVertex &last = m_vertices[m_vertices.size()-1];
    const PathVertex *vb = path.vertex(b);
    if (!vb->isOnSurface()) {
        last.n = Vector(path.edge(a < b ? (b-1) : b)->d);
    } else {
        last.n = vb->getShadingNormal();
    }
    coordinateSystem(last.n, last.dpdu, last.dpdv);

    statsNonManifold.incrementBase();
    if (!computeTangents()) {
        ++statsNonManifold;
        Log(EWarn, "SpecularManifold::evalG(): non-manifold configuration!");
        return 0;
    }

    Float result;
    if (m_vertices[0].type == EPinnedDirection) {
        result = cross(m_vertices[0].map(1, 0), m_vertices[0].map(0, 1)).length();
    } else if (m_vertices[0].type == EPinnedPosition) {
        Vector d = m_vertices[1].p - m_vertices[0].p;
        Float lengthSqr = d.lengthSquared(), invLength = 1/std::sqrt(lengthSqr);

        result = cross(m_vertices[1].map(1, 0), m_vertices[1].map(0, 1)).length() / lengthSqr;

        if (path.vertex(a)->isOnSurface())
            result *= absDot(d, path.vertex(a)->getShadingNormal()) * invLength;

        if (path.vertex(a+step)->isOnSurface())
            result *= absDot(d, path.vertex(a+step)->getShadingNormal()) * invLength;
    } else {
        Log(EError, "Invalid vertex type!");
        return 0;
    }

    return result;
}

std::string SpecularManifold::SimpleVertex::toString() const {
    std::ostringstream oss;

    oss << "SimpleVertex[" << endl
        << "  type = ";

    switch (type) {
        case EPinnedPosition: oss << "pinnedPosition"; break;
        case EPinnedDirection: oss << "pinnedDirection"; break;
        case EReflection: oss << "reflection"; break;
        case ERefraction: oss << "refraction"; break;
        case EMedium: oss << "medium"; break;
        case EMovable: oss << "movable"; break;
        default: SLog(EError, "Unknown vertex type!");
    }

    oss << "," << endl
        << "  p = " << p.toString() << "," << endl
        << "  n = " << n.toString() << "," << endl
        << "  m = " << m.toString() << "," << endl
        << "  dpdu = " << dpdu.toString() << "," << endl
        << "  dpdv = " << dpdv.toString() << "," << endl
        << "  dndu = " << dndu.toString() << "," << endl
        << "  dndv = " << dndv.toString() << "," << endl
        << "  eta = " << eta << "," << endl
        << "  object = " << (object ? indent(object->toString()).c_str() : "null") << endl
        << "]";

    return oss.str();
}

std::string SpecularManifold::toString() const {
    std::ostringstream oss;

    oss << "SpecularManifold[" << endl;
    for (size_t i=0; i<m_vertices.size(); ++i) {
        oss << "  " << i << " => " << indent(m_vertices[i].toString());
        if (i+1 < m_vertices.size())
            oss << ",";
        oss << endl;
    }
    oss << "]";

    return oss.str();
}

MTS_IMPLEMENT_CLASS(SpecularManifold, false, Object)
MTS_NAMESPACE_END

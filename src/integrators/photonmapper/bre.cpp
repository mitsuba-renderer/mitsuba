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

#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/core/timer.h>
#if defined(MTS_OPENMP)
# include <omp.h>
#endif
#include "bre.h"

MTS_NAMESPACE_BEGIN

BeamRadianceEstimator::BeamRadianceEstimator(const PhotonMap *pmap, size_t lookupSize) {
    /* Use an optimization proposed by Jarosz et al, which accelerates
       the radius computation by extrapolating radius information obtained
       from a kd-tree lookup of a smaller size */
    size_t reducedLookupSize = (size_t) std::sqrt((Float) lookupSize);
    Float sizeFactor = (Float) lookupSize / (Float) reducedLookupSize;

    m_photonCount = pmap->size();
    m_scaleFactor = pmap->getScaleFactor();
    m_depth = pmap->getDepth();

    Log(EInfo, "Allocating %s of memory for the BRE acceleration data structure",
        memString(sizeof(BRENode) * m_photonCount).c_str());
    m_nodes = new BRENode[m_photonCount];

    Log(EInfo, "Computing photon radii ..");
    #if defined(MTS_OPENMP)
        int tcount = mts_omp_get_max_threads();
    #else
        int tcount = 1;
    #endif
    PhotonMap::SearchResult **resultsPerThread = new PhotonMap::SearchResult*[tcount];
    for (int i=0; i<tcount; ++i)
        resultsPerThread[i] = new PhotonMap::SearchResult[reducedLookupSize+1];

    ref<Timer> timer = new Timer();
    #if defined(MTS_OPENMP)
        #pragma omp parallel for
    #endif
    for (int i=0; i<(int) m_photonCount; ++i) {
        #if defined(MTS_OPENMP)
            int tid = mts_omp_get_thread_num();
        #else
            int tid = 0;
        #endif

        PhotonMap::SearchResult *results = resultsPerThread[tid];
        const Photon &photon = pmap->operator[](i);
        BRENode &node = m_nodes[i];
        node.photon = photon;

        Float searchRadiusSqr = std::numeric_limits<Float>::infinity();
        pmap->nnSearch(photon.getPosition(), searchRadiusSqr, reducedLookupSize, results);

        /* Compute photon radius based on a locally uniform density assumption */
        node.radius = std::sqrt(searchRadiusSqr * sizeFactor);
    }
    Log(EInfo, "Done (took %i ms)", timer->getMilliseconds());

    Log(EInfo, "Generating a hierarchy for the beam radiance estimate");
    timer->reset();

    buildHierarchy(0);
    Log(EInfo, "Done (took %i ms)", timer->getMilliseconds());

    for (int i=0; i<tcount; ++i)
        delete[] resultsPerThread[i];
    delete[] resultsPerThread;
}

BeamRadianceEstimator::BeamRadianceEstimator(Stream *stream, InstanceManager *manager) {
    m_photonCount = stream->readSize();
    m_depth = stream->readSize();
    m_scaleFactor = stream->readFloat();
    m_nodes = new BRENode[m_photonCount];
    for (size_t i=0; i<m_photonCount; ++i) {
        BRENode &node = m_nodes[i];
        node.aabb = AABB(stream);
        node.photon = Photon(stream);
        node.radius = stream->readFloat();
    }
}

void BeamRadianceEstimator::serialize(Stream *stream, InstanceManager *manager) const {
    Log(EDebug, "Serializing a BRE data structure (%s)",
            memString(m_photonCount * sizeof(BRENode)).c_str());
    stream->writeSize(m_photonCount);
    stream->writeSize(m_depth);
    stream->writeFloat(m_scaleFactor);
    for (size_t i=0; i<m_photonCount; ++i) {
        BRENode &node = m_nodes[i];
        node.aabb.serialize(stream);
        node.photon.serialize(stream);
        stream->writeFloat(node.radius);
    }
}

AABB BeamRadianceEstimator::buildHierarchy(IndexType index) {
    BRENode &node = m_nodes[index];

    Point center = node.photon.getPosition();
    Float radius = node.radius;
    node.aabb = AABB(
        center - Vector(radius, radius, radius),
        center + Vector(radius, radius, radius)
    );

    if (!node.photon.isLeaf()) {
        IndexType left = node.photon.getLeftIndex(index);
        IndexType right = node.photon.getRightIndex(index);
        if (left)
            node.aabb.expandBy(buildHierarchy(left));
        if (right)
            node.aabb.expandBy(buildHierarchy(right));
    }

    return node.aabb;
}

Spectrum BeamRadianceEstimator::query(const Ray &r, const Medium *medium) const {
    const Ray ray(r(r.mint), r.d, 0, r.maxt - r.mint, r.time);
    IndexType *stack = (IndexType *) alloca((m_depth+1) * sizeof(IndexType));
    IndexType index = 0, stackPos = 1;
    Spectrum result(0.0f);

    const Spectrum &sigmaT = medium->getSigmaT();
    const PhaseFunction *phase = medium->getPhaseFunction();
    MediumSamplingRecord mRec;

    while (stackPos > 0) {
        const BRENode &node = m_nodes[index];
        const Photon &photon = node.photon;

        /* Test against the node's bounding box */
        Float mint, maxt;
        if (!node.aabb.rayIntersect(ray, mint, maxt) || maxt < ray.mint || mint > ray.maxt) {
            index = stack[--stackPos];
            continue;
        }

        /* Recurse on inner photons */
        if (!photon.isLeaf()) {
            if (hasRightChild(index))
                stack[stackPos++] = photon.getRightIndex(index);
            index = photon.getLeftIndex(index);
        } else {
            index = stack[--stackPos];
        }

        Vector originToCenter = node.photon.getPosition() - ray.o;
        Float diskDistance = dot(originToCenter, ray.d), radSqr = node.radius * node.radius;
        Float distSqr = (ray(diskDistance) - node.photon.getPosition()).lengthSquared();

        if (diskDistance > 0 && distSqr < radSqr) {
            Float weight = K2(distSqr/radSqr)/radSqr;

            Vector wi = -node.photon.getDirection();

            Spectrum transmittance = Spectrum(-sigmaT * diskDistance).exp();
            result += transmittance * node.photon.getPower()
                * phase->eval(PhaseFunctionSamplingRecord(mRec, wi, -ray.d)) *
                (weight * m_scaleFactor);
        }
    }

    return result;
}

BeamRadianceEstimator::~BeamRadianceEstimator() {
    delete[] m_nodes;
}

MTS_IMPLEMENT_CLASS_S(BeamRadianceEstimator, false, Object)
MTS_NAMESPACE_END

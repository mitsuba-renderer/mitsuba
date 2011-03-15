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

#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include "bre.h"

MTS_NAMESPACE_BEGIN

BeamRadianceEstimator::BeamRadianceEstimator(const PhotonMap *pmap, size_t lookupSize) {
	size_t reducedLookupSize = std::sqrt(lookupSize);
	Float sizeFactor = (Float) lookupSize/ (Float) reducedLookupSize;

	PhotonMap::search_result *results = 
		new PhotonMap::search_result[reducedLookupSize+1];

	m_photonCount = pmap->getPhotonCount();
	m_scaleFactor = pmap->getScaleFactor();
	m_lastInnerNode = m_photonCount/2;
	m_lastRChildNode = (m_photonCount-1)/2;
	m_depth = log2i(m_photonCount)+1;

	m_nodes = new BRENode[m_photonCount+1];

	Log(EInfo, "Computing photon radii ..");
	for (size_t i=1; i<=m_photonCount; ++i) {
		const Photon &photon = pmap->getPhoton(i);
		BRENode &node = m_nodes[i];
		node.photon = photon;

		Float searchRadiusSqr = std::numeric_limits<Float>::infinity();
		pmap->nnSearch(photon.getPosition(), searchRadiusSqr, reducedLookupSize, results);

		/* Compute photon radius based on a locally uniform density assumption */
		node.radius = std::sqrt(searchRadiusSqr * sizeFactor);
	}

	Log(EInfo, "Generating a hierarchy for the beam radiance estimate");

	buildHierarchy(1);

	delete[] results;
}

BeamRadianceEstimator::BeamRadianceEstimator(Stream *stream, InstanceManager *manager) {
	m_photonCount = stream->readSize();
	m_scaleFactor = stream->readFloat();
	m_lastInnerNode = m_photonCount/2;
	m_lastRChildNode = (m_photonCount-1)/2;
	m_depth = log2i(m_photonCount)+1;
	m_nodes = new BRENode[m_photonCount+1];
	for (size_t i=1; i<=m_photonCount; ++i) {
		BRENode &node = m_nodes[i];
		node.aabb = AABB(stream);
		node.photon = Photon(stream);
		node.radius = stream->readFloat();
	}
}

void BeamRadianceEstimator::serialize(Stream *stream, InstanceManager *manager) const {
	stream->writeSize(m_photonCount);
	stream->writeFloat(m_scaleFactor);
	for (size_t i=1; i<=m_photonCount; ++i) {
		BRENode &node = m_nodes[i];
		node.aabb.serialize(stream);
		node.photon.serialize(stream);
		stream->writeFloat(node.radius);
	}
}

AABB BeamRadianceEstimator::buildHierarchy(size_t index) {
	BRENode &node = m_nodes[index];

	if (isInnerNode(index)) {
		node.aabb = buildHierarchy(leftChild(index));
		if (hasRightChild(index))
			node.aabb.expandBy(buildHierarchy(rightChild(index)));
	} else {
		Point center = node.photon.getPosition();
		Float radius = node.radius;
		node.aabb = AABB(
			center - Vector(radius, radius, radius),
			center + Vector(radius, radius, radius)
		);
	}

	return node.aabb;
}

inline Float K2(Float sqrParam) {
	Float tmp = 1-sqrParam;
	return (3/M_PI) * tmp * tmp;
}

Spectrum BeamRadianceEstimator::query(const Ray &ray, const Medium *medium) const {
	uint32_t *stack = (uint32_t *) alloca((m_depth+1) * sizeof(uint32_t));
	uint32_t index = 1, stackPos = 1;
	Spectrum result(0.0f);
	size_t nNodes = 0;

	const Spectrum &sigmaT = medium->getSigmaT();
	const PhaseFunction *phase = medium->getPhaseFunction();
	MediumSamplingRecord mRec;

	while (stackPos > 0) {
		const BRENode &node = m_nodes[index];

		/* Test against the node's bounding box */
		Float mint, maxt;
		if (!node.aabb.rayIntersect(ray, mint, maxt) || maxt < ray.mint || mint > ray.maxt) {
			index = stack[--stackPos];
			continue;
		}
		++nNodes;

		/* Recurse on inner nodes */
		if (isInnerNode(index)) {
			if (hasRightChild(index))
				stack[stackPos++] = leftChild(index);
			index = rightChild(index);
		} else {
			index = stack[--stackPos];
		}

		Vector originToCenter = node.photon.getPosition() - ray.o;
		Float diskDistance = dot(originToCenter, ray.d), radSqr = node.radius * node.radius;
		Float distSqr = (ray(diskDistance) - node.photon.getPosition()).lengthSquared();

		if (distSqr < radSqr) {
			Float weight = K2(distSqr/radSqr)/radSqr;

			Vector wi = -node.photon.getDirection();

			Spectrum transmittance = Spectrum(-sigmaT * diskDistance).exp();
			result += transmittance * node.photon.getPower()
				* phase->f(PhaseFunctionQueryRecord(mRec, wi, -ray.d)) *
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

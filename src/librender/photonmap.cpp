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

#include <mitsuba/render/photonmap.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/phase.h>
#include <fstream>

MTS_NAMESPACE_BEGIN

PhotonMap::PhotonMap(size_t maxPhotons) 
 : m_photonCount(0), m_maxPhotons(maxPhotons), m_balanced(false), m_scale(1.0f) {
	Assert(Photon::m_precompTableReady);

	/* For convenient heap addressing, the the photon list
	   entries start with number 1 */
	m_photons = (Photon *) allocAligned(sizeof(Photon) * (maxPhotons+1));
}
	
PhotonMap::PhotonMap(Stream *stream, InstanceManager *manager) { 
	m_aabb = AABB(stream);
	m_balanced = stream->readBool();
	m_maxPhotons = stream->readSize();
	m_lastInnerNode = stream->readSize();
	m_lastRChildNode = stream->readSize();
	m_scale = (Float) stream->readFloat();
	m_photonCount = stream->readSize();
	m_photons = new Photon[m_maxPhotons + 1];
	for (size_t i=1; i<=m_maxPhotons; ++i) 
		m_photons[i] = Photon(stream);
}

PhotonMap::~PhotonMap() {
	freeAligned(m_photons);
}

std::string PhotonMap::toString() const {
	std::ostringstream oss;
	oss << "PhotonMap[" << endl
		<< "  aabb = " << m_aabb.toString() << "," << endl
		<< "  photonCount = " << m_photonCount << "," << endl
		<< "  maxPhotons = " << m_maxPhotons << "," << endl
		<< "  balanced = " << m_balanced << "," << endl
		<< "  scale = " << m_scale << endl
		<< "]";
	return oss.str();
}

void PhotonMap::serialize(Stream *stream, InstanceManager *manager) const {
	Log(EDebug, "Serializing a photon map (%.2f KB)", 
		m_photonCount * 20.0f / 1024.0f);
	m_aabb.serialize(stream);
	stream->writeBool(m_balanced);
	stream->writeSize(m_maxPhotons);
	stream->writeSize(m_lastInnerNode);
	stream->writeSize(m_lastRChildNode);
	stream->writeFloat(m_scale);
	stream->writeSize(m_photonCount);
	for (size_t i=1; i<=m_maxPhotons; ++i)
		m_photons[i].serialize(stream);
}

void PhotonMap::dumpOBJ(const std::string &filename) {
	std::ofstream os(filename.c_str());
	os << "o Photons" << endl;
	for (size_t i=1; i<=getPhotonCount(); i++) {
		Point p = getPhoton(i).getPosition();
		os << "v " << p.x << " " << p.y << " " << p.z << endl;
	}

	/// Need to generate some fake geometry so that blender will import the points
	for (size_t i=3; i<=getPhotonCount(); i++) 
		os << "f " << i << " " << i-1 << " " << i-2 << endl;
	os.close();
}

MTS_IMPLEMENT_CLASS_S(PhotonMap, false, Object)
MTS_NAMESPACE_END

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

PhotonMap::PhotonMap(size_t photonCount) : m_kdtree(photonCount), m_scale(1.0f) {
	Assert(Photon::m_precompTableReady);
}

PhotonMap::PhotonMap(Stream *stream, InstanceManager *manager) { 
	m_aabb = AABB(stream);
	m_scale = (Float) stream->readFloat();
	m_kdtree.resize(stream->readSize());
	for (size_t i=0; i<m_kdtree.size(); ++i) 
		m_kdtree[i] = Photon(stream);
}

PhotonMap::~PhotonMap() {
}

std::string PhotonMap::toString() const {
	std::ostringstream oss;
	oss << "PhotonMap[" << endl
		<< "  aabb = " << m_aabb.toString() << "," << endl
		<< "  size = " << m_kdtree.size() << "," << endl
		<< "  capacity = " << m_kdtree.capacity() << "," << endl
		<< "  scale = " << m_scale << endl
		<< "]";
	return oss.str();
}

void PhotonMap::serialize(Stream *stream, InstanceManager *manager) const {
	Log(EDebug, "Serializing a photon map (%.2f KB)", 
		m_photonCount * sizeof(Photon) / 1024.0f);
	m_aabb.serialize(stream);
	stream->writeFloat(m_scale);
	stream->writeSize(m_tree.size());
	for (size_t i=0; i<m_tree.size(); ++i)
		m_photons[i].serialize(stream);
}

void PhotonMap::dumpOBJ(const std::string &filename) {
	std::ofstream os(filename.c_str());
	os << "o Photons" << endl;
	for (size_t i=0; i<m_tree.size(); ++i) {
		const Point &p = m_tree[i].getPosition();
		os << "v " << p.x << " " << p.y << " " << p.z << endl;
	}

	/// Need to generate some fake geometry so that blender will import the points
	for (size_t i=3; i<=m_tree.size(); i++) 
		os << "f " << i << " " << i-1 << " " << i-2 << endl;
	os.close();
}

MTS_IMPLEMENT_CLASS_S(PhotonMap, false, Object)
MTS_NAMESPACE_END

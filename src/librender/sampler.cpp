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

#include <mitsuba/render/sampler.h>

MTS_NAMESPACE_BEGIN

Sampler::Sampler(const Properties &props) 
 : ConfigurableObject(props), m_properties(props) {
}

Sampler::Sampler(Stream *stream, InstanceManager *manager) 
 : ConfigurableObject(stream, manager) {
	m_sampleCount = stream->readULong();
	unsigned int n1DArrays = stream->readUInt();
	for (unsigned int i=0; i<n1DArrays; ++i) 
		request1DArray(stream->readUInt());
	unsigned int n2DArrays = stream->readUInt();
	for (unsigned int i=0; i<n2DArrays; ++i) 
		request2DArray(stream->readUInt());
}

void Sampler::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);

	stream->writeULong(m_sampleCount);
	stream->writeUInt(m_req1D.size());
	for (size_t i=0; i<m_req1D.size(); ++i)
		stream->writeUInt(m_req1D[i]);
	stream->writeUInt(m_req2D.size());
	for (size_t i=0; i<m_req2D.size(); ++i)
		stream->writeUInt(m_req2D[i]);
}

void Sampler::generate() {
	m_sampleIndex = 0;
	m_sampleDepth1DArray = m_sampleDepth2DArray = 0;
}

void Sampler::advance() {
	m_sampleIndex++;
	m_sampleDepth1DArray = m_sampleDepth2DArray = 0;
}

void Sampler::setSampleIndex(uint64_t sampleIndex) {
	m_sampleIndex = sampleIndex;
	m_sampleDepth1DArray = m_sampleDepth2DArray = 0;
}

void Sampler::request1DArray(unsigned int size) {
	m_req1D.push_back(size);
	m_sampleArrays1D.push_back(new Float[m_sampleCount * size]);
}

void Sampler::request2DArray(unsigned int size) {
	m_req2D.push_back(size);
	m_sampleArrays2D.push_back(new Point2[m_sampleCount * size]);
}

Point2 *Sampler::next2DArray(unsigned int size) {
	Assert(m_sampleIndex < m_sampleCount);
	if (m_sampleDepth2DArray < (int) m_req2D.size()) {
		Assert(m_req2D[m_sampleDepth2DArray] == size);
		return m_sampleArrays2D[m_sampleDepth2DArray++] + m_sampleIndex * size;
	} else {
		Log(EError, "Tried to retrieve a size-%i 2D array,"
			" when this was not previously requested.", size);
		return NULL;
	}
}

Float *Sampler::next1DArray(unsigned int size) {
	Assert(m_sampleIndex < m_sampleCount);
	if (m_sampleDepth1DArray < (int) m_req1D.size()) {
		Assert(m_req1D[m_sampleDepth1DArray] == size);
		return m_sampleArrays1D[m_sampleDepth1DArray++] + m_sampleIndex * size;
	} else {
		Log(EError, "Tried to retrieve a size-%i 1D array,"
			" when this was not previously requested.", size);
		return NULL;
	}
}

Sampler::~Sampler() {
	for (size_t i=0; i<m_sampleArrays1D.size(); i++) {
		if (m_sampleArrays1D[i])
			delete[] m_sampleArrays1D[i];
	}
	for (size_t i=0; i<m_sampleArrays2D.size(); i++) {
		if (m_sampleArrays2D[i])
			delete[] m_sampleArrays2D[i];
	}
}

MTS_IMPLEMENT_CLASS(Sampler, true, ConfigurableObject)
MTS_NAMESPACE_END

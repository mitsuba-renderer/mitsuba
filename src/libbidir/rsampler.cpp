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

#include <mitsuba/bidir/rsampler.h>

MTS_NAMESPACE_BEGIN

ReplayableSampler::ReplayableSampler() : Sampler(Properties()) {
    m_initial = new Random();
    m_random = new Random();
    m_random->set(m_initial);
    m_sampleCount = 0;
    m_sampleIndex = 0;
}

ReplayableSampler::ReplayableSampler(Stream *stream, InstanceManager *manager)
    : Sampler(stream, manager) {
    m_initial = static_cast<Random *>(manager->getInstance(stream));
    m_random = new Random();
    m_random->set(m_initial);
    m_sampleCount = 0;
    m_sampleIndex = 0;
}

ReplayableSampler::~ReplayableSampler() {
}

void ReplayableSampler::serialize(Stream *stream, InstanceManager *manager) const {
    Sampler::serialize(stream, manager);
    manager->serialize(stream, m_initial.get());
}

ref<Sampler> ReplayableSampler::clone() {
    ref<ReplayableSampler> sampler = new ReplayableSampler();
    sampler->m_sampleCount = m_sampleCount;
    sampler->m_sampleIndex = m_sampleIndex;
    sampler->m_initial->set(m_initial);
    sampler->m_random->set(m_random);
    return sampler.get();
}

void ReplayableSampler::request1DArray(size_t size) {
    Log(EError, "ReplayableSampler::request2DArray() - unsupported!");
}

void ReplayableSampler::request2DArray(size_t size) {
    Log(EError, "ReplayableSampler::request2DArray() - unsupported!");
}

void ReplayableSampler::generate(const Point2i &) { }
void ReplayableSampler::advance() { }

void ReplayableSampler::setSampleIndex(size_t sampleIndex) {
    if (sampleIndex < m_sampleIndex) {
        m_sampleIndex = 0;
        m_random->set(m_initial);
    }

    while (m_sampleIndex != sampleIndex) {
        m_random->nextFloat();
        ++m_sampleIndex;
    }
}

Float ReplayableSampler::next1D() {
    ++m_sampleIndex;
    return m_random->nextFloat();
}

Point2 ReplayableSampler::next2D() {
    /// Enforce a specific order of evaluation
    Float value1 = m_random->nextFloat();
    Float value2 = m_random->nextFloat();
    m_sampleIndex += 2;
    return Point2(value1, value2);
}

std::string ReplayableSampler::toString() const {
    std::ostringstream oss;
    oss << "ReplayableSampler[" << endl
        << "  sampleCount = " << m_sampleCount << endl
        << "]";
    return oss.str();
}


MTS_IMPLEMENT_CLASS_S(ReplayableSampler, false, Sampler)
MTS_NAMESPACE_END

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

#include <mitsuba/core/rfilter.h>

MTS_NAMESPACE_BEGIN

ReconstructionFilter::ReconstructionFilter(const Properties &props)
 : ConfigurableObject(props) { }

ReconstructionFilter::~ReconstructionFilter() { }

ReconstructionFilter::ReconstructionFilter(Stream *stream, InstanceManager *manager)
 : ConfigurableObject(stream, manager) {
	 m_radius = stream->readFloat();
}

void ReconstructionFilter::serialize(Stream *stream, InstanceManager *manager) const {
	stream->writeFloat(m_radius);
}

void ReconstructionFilter::configure() {
	Assert(m_radius > 0);

	Float sum = 0.0f;
	/* Evaluate and normalize the filter */
	for (size_t i=0; i<MTS_FILTER_RESOLUTION; ++i) {
		Float value = eval((m_radius * i) / MTS_FILTER_RESOLUTION);
		m_values[i] = value;
		sum += value;
	}

	m_values[MTS_FILTER_RESOLUTION] = 0.0f;
	m_scaleFactor = MTS_FILTER_RESOLUTION / m_radius;
	m_borderSize = (int) std::ceil(m_radius - 0.5f);
	sum *= 2 * m_radius / MTS_FILTER_RESOLUTION;
	Float normalization = 1.0f / sum;
	for (size_t i=0; i<MTS_FILTER_RESOLUTION; ++i)
		m_values[i] *= normalization;
}

std::ostream &operator<<(std::ostream &os, const ReconstructionFilter::EBoundaryCondition &value) {
	switch (value) {
		case ReconstructionFilter::EClamp: os << "clamp"; break;
		case ReconstructionFilter::ERepeat: os << "repeat"; break;
		case ReconstructionFilter::EMirror: os << "mirror"; break;
		case ReconstructionFilter::EZero: os << "zero"; break;
		case ReconstructionFilter::EOne: os << "one"; break;
		default: os << "invalid"; break;
	}
	return os;
}

MTS_IMPLEMENT_CLASS(ReconstructionFilter, true, ConfigurableObject)
MTS_NAMESPACE_END

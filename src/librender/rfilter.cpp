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

#include <mitsuba/render/rfilter.h>

MTS_NAMESPACE_BEGIN

ReconstructionFilter::ReconstructionFilter(const Properties &props)
 : ConfigurableObject(props), m_properties(props) {
}
	
ReconstructionFilter::ReconstructionFilter(Stream *stream, InstanceManager *manager) 
 : ConfigurableObject(stream, manager) {
}

ReconstructionFilter::~ReconstructionFilter() {
}
	
Float ReconstructionFilter::sample(const Point2 &sample, Float &x, Float &y) const {
	x = (1 - sample.x * 2) * m_size.x;
	y = (1 - sample.y * 2) * m_size.y;
	return evaluate(x, y) / (4 * m_size.x * m_size.y);
}

TabulatedFilter::TabulatedFilter(const ReconstructionFilter *filter) {
	m_size = filter->getFilterSize();
	m_name = filter->getClass()->getName();
	m_factor = Vector2(
		FILTER_RESOLUTION / m_size.x,
		FILTER_RESOLUTION / m_size.y);
	Float sum = 0;
	/* Evaluate the filter and add a zero border (for performance
	   reasons during the evaluation) */
	for (int y=0; y<FILTER_RESOLUTION+1; ++y) {
		Float yPos = (y + 0.5f) / FILTER_RESOLUTION * m_size.y;
		for (int x=0; x<FILTER_RESOLUTION+1; ++x) {
			if (x == FILTER_RESOLUTION || y==FILTER_RESOLUTION)
				m_values[y][x] = 0;
			else
				m_values[y][x] = filter->evaluate(
					(x + 0.5f) / FILTER_RESOLUTION * m_size.x, yPos);
			sum += m_values[y][x];
		}
	}

	/* Normalize the filter (required for the particle tracer) */
	sum *= 4*filter->getFilterSize().x*filter->getFilterSize().y/
		(FILTER_RESOLUTION*FILTER_RESOLUTION);
	for (int y=0; y<FILTER_RESOLUTION+1; ++y) {
		for (int x=0; x<FILTER_RESOLUTION+1; ++x) {
			m_values[y][x] /= sum;
		}
	}
}

TabulatedFilter::TabulatedFilter(Stream *stream) {
	m_name = stream->readString();
	m_size = Vector2(stream);
	m_factor = Vector2(stream);
	for (int y=0; y<FILTER_RESOLUTION+1; ++y)
		stream->readFloatArray(m_values[y], FILTER_RESOLUTION+1);
}

void TabulatedFilter::serialize(Stream *stream) const {
	stream->writeString(m_name);
	m_size.serialize(stream);
	m_factor.serialize(stream);
	for (int y=0; y<FILTER_RESOLUTION+1; ++y)
		stream->writeFloatArray(m_values[y], FILTER_RESOLUTION+1);
}

TabulatedFilter::~TabulatedFilter() {
}

std::string TabulatedFilter::toString() const {
	std::ostringstream oss;
	oss << "TabulatedFilter[size=" << m_size.toString() << "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(ReconstructionFilter, true, ConfigurableObject)
MTS_IMPLEMENT_CLASS(TabulatedFilter, false, Object)
MTS_NAMESPACE_END

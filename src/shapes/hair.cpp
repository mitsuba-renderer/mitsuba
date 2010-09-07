#include <fstream>

#include <mitsuba/render/shape.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fresolver.h>

#include "hair.h"
#include "miterseg.h"

MTS_NAMESPACE_BEGIN

Hair::Hair(const Properties &props) : Shape(props) {
	std::string filename = props.getString("filename");
	m_radius = (Float) props.getFloat("radius", 0.05f);
	m_name = FileResolver::getInstance()->resolve(filename);

	Log(EInfo, "Loading hair geometry from \"%s\" ..", m_name.c_str());

	std::ifstream is(m_name.c_str());
	if (is.fail())
		Log(EError, "Could not open \"%s\"!", m_name.c_str());

	std::string line;
	bool newFiber = true;
	Point p;
	while (is.good()) {
		std::getline(is, line);
		if (line.length() > 0 && line[0] == '#')
			continue;
		if (line.length() == 0) {
			newFiber = true;
		} else {
			std::istringstream iss(line);
			iss >> p.x >> p.y >> p.z;
			if (iss.good()) {
				m_vertices.push_back(p);
				m_startFiber.push_back(newFiber);
				newFiber = false;
			}
		}
	}
	m_startFiber.push_back(true);

	buildSegIndex();

	Log(EDebug, "Read %i hair vertices, %i segments,", m_vertices.size(), m_segIndex.size());
}

Hair::Hair(Stream *stream, InstanceManager *manager) : Shape(stream, manager) {
	m_radius = stream->readFloat();
	size_t segmentCount = stream->readUInt();

	m_vertices.reserve(segmentCount);
	for (size_t i=0; i<segmentCount; ++i)
		m_vertices.push_back(Point(stream));

	m_startFiber.reserve(segmentCount+1);
	for (size_t i=0; i<segmentCount+1; ++i)
		m_startFiber.push_back(stream->readBool());

	buildSegIndex();
}

void Hair::serialize(Stream *stream, InstanceManager *manager) const {
	Shape::serialize(stream, manager);

	stream->writeFloat(m_radius);
	size_t segmentCount = m_vertices.size();
	stream->writeUInt(segmentCount);

	for (size_t i=0; i<segmentCount; ++i)
		m_vertices[i].serialize(stream);

	for (size_t i=0; i<segmentCount; ++i)
		stream->writeBool(m_startFiber[i]);
}

Shape *Hair::getElement(int index) {
	if ((size_t) index >= m_segIndex.size())
		return NULL;

	MiterHairSegment *segment = new MiterHairSegment(this, m_segIndex[index]);
	segment->addChild("bsdf", m_bsdf);
	segment->configure();
	return segment;
}


void Hair::buildSegIndex() {
	// Compute the index of the first vertex in each segment.
	m_segIndex.clear();
	for (size_t i=0; i<m_vertices.size(); i++)
		if (!m_startFiber[i+1])
			m_segIndex.push_back(i);
}



MTS_IMPLEMENT_CLASS_S(Hair, false, Shape)
MTS_EXPORT_PLUGIN(Hair, "Hair geometry");
MTS_NAMESPACE_END

#include <mitsuba/render/imageproc_wu.h>

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                          RectangularWorkUnit                         */
/* ==================================================================== */
	
void RectangularWorkUnit::set(const WorkUnit *wu) {
	const RectangularWorkUnit *rect = static_cast<const RectangularWorkUnit *>(wu);
	m_offset = rect->m_offset;
	m_size = rect->m_size;
}

void RectangularWorkUnit::load(Stream *stream) {
	int data[4];
	stream->readIntArray(data, 4);
	m_offset.x = data[0];
	m_offset.y = data[1];
	m_size.x   = data[2];
	m_size.y   = data[3];
}

void RectangularWorkUnit::save(Stream *stream) const {
	int data[4];
	data[0] = m_offset.x;
	data[1] = m_offset.y;
	data[2] = m_size.x;
	data[3] = m_size.y;
	stream->writeIntArray(data, 4);
}

std::string RectangularWorkUnit::toString() const {
	std::ostringstream oss;
	oss << "RectangularWorkUnit[offset=" << m_offset.toString() 
		<< ", size=" << m_size.toString() << "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(RectangularWorkUnit, false, WorkUnit)
MTS_NAMESPACE_END

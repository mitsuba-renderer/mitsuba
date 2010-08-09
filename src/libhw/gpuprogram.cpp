#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

GPUProgram::GPUProgram(const std::string &name)
 : m_name(name), m_maxVertices(0), m_bound(false) {
}

GPUProgram::~GPUProgram() {
}

std::string GPUProgram::toString() const {
	std::ostringstream oss;
	oss << "GPUProgram[name = '" << m_name<< "'";
	if (m_maxVertices != 0)
		oss << ", maxVertices=" << m_maxVertices;
	oss << "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(GPUProgram, true, Object)
MTS_NAMESPACE_END

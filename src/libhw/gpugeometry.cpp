#include <mitsuba/hw/gpugeometry.h>

MTS_NAMESPACE_BEGIN

GPUGeometry::GPUGeometry(const TriMesh *mesh)
 : m_mesh(mesh) {
}

GPUGeometry::~GPUGeometry() {
}

std::string GPUGeometry::toString() const {
	std::ostringstream oss;
	oss << "GPUGeometry[" << endl
		<< "  mesh = " << indent(m_mesh->toString()) << endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(GPUGeometry, true, Object)
MTS_NAMESPACE_END

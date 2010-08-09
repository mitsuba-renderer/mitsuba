#include <mitsuba/hw/gpusync.h>

MTS_NAMESPACE_BEGIN

GPUSync::GPUSync() { }
GPUSync::~GPUSync() { }
std::string GPUSync::toString() const {
	return "GPUSync[]";
}

MTS_IMPLEMENT_CLASS(GPUSync, true, Object)
MTS_NAMESPACE_END

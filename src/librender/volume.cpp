#include <mitsuba/render/volume.h>

#if defined(__LINUX__)
#include <sys/mman.h>
#include <fcntl.h>
#endif

MTS_NAMESPACE_BEGIN

VolumeDataSource::VolumeDataSource(const Properties &props) : ConfigurableObject(props) {
}

VolumeDataSource::VolumeDataSource(Stream *stream, InstanceManager *manager) :
	ConfigurableObject(stream, manager) {
	m_aabb = AABB(stream);
}

VolumeDataSource::~VolumeDataSource() {
}

void VolumeDataSource::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
	m_aabb.serialize(stream);
}

Float VolumeDataSource::lookupFloat(const Point &p) const {
	Log(EError, "'%s': does not implement lookupFloat()!", getClass()->getName().c_str());
	return 0;
}

Spectrum VolumeDataSource::lookupSpectrum(const Point &p) const {
	Log(EError, "'%s': does not implement lookupSpectrum()!", getClass()->getName().c_str());
	return Spectrum(0.0f);
}

Vector VolumeDataSource::lookupVector(const Point &p) const {
	Log(EError, "'%s': does not implement lookupVector()!", getClass()->getName().c_str());
	return Vector();
}

bool VolumeDataSource::supportsFloatLookups() const {
	return false;
}

bool VolumeDataSource::supportsSpectrumLookups() const {
	return false;
}

bool VolumeDataSource::supportsVectorLookups() const {
	return false;
}

MTS_IMPLEMENT_CLASS(VolumeDataSource, true, ConfigurableObject)
MTS_NAMESPACE_END


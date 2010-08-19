#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

BSDF::BSDF(const Properties &props)
 : ConfigurableObject(props), m_type(NULL), m_name(props.getID()) {
}

BSDF::BSDF(Stream *stream, InstanceManager *manager) 
 : ConfigurableObject(stream, manager), m_type(NULL) {
	m_name = stream->readString();
}

BSDF::~BSDF() {
}

/* Inefficient version in case this is not supported by the BSDF implementation */
Spectrum BSDF::sample(BSDFQueryRecord &bRec, Float &_pdf) const {
	if (sample(bRec).isBlack()) {
		_pdf = 0.0f;
		return Spectrum(0.0f);
	}
	/* Re-evaluation required because we want both the
		value and a matching probability density.
		(the previously sampled value may ony be wrt.
		one of the BSDF lobes) */
	_pdf = pdf(bRec);
	return f(bRec);
}

void BSDF::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
	stream->writeString(m_name);
}

void BSDF::addChild(const std::string &name, ConfigurableObject *obj) {
	ConfigurableObject::addChild(name, obj);
}

Float BSDF::pdfDelta(const BSDFQueryRecord &bRec) const {
	return 0.0f;
}

Spectrum BSDF::fDelta(const BSDFQueryRecord &bRec) const {
	return Spectrum(0.0f);
}

MTS_IMPLEMENT_CLASS(BSDF, true, ConfigurableObject)
MTS_NAMESPACE_END

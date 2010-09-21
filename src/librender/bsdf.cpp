/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/bsdf.h>
#include <mitsuba/core/frame.h>

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

void operator<<(const ETransportQuantity &quantity, std::ostream &os) {
	switch (quantity) {
		case EImportance: os << "importance"; break;
		case ERadiance:   os << "radiance"; break;
		default: os << "invalid"; break;
	};
}

std::string BSDFQueryRecord::toString() const {
	std::ostringstream oss;
	oss << "BSDFQueryRecord[" << std::endl
		<< "  wi = " << wi.toString() << "," << std::endl
		<< "  wo = " << wo.toString() << "," << std::endl
		<< "  sample = " << sample.toString() << "," << std::endl
		<< "  quantity = " << quantity << "," << std::endl
		<< "  typeMask = " << typeMask << "," << std::endl
		<< "  sampledType = " << sampledType << "," << std::endl
		<< "  component = " << component << "," << std::endl
		<< "  sampledComponent = " << sampledComponent << "," << std::endl
		<< "]";
	return oss.str();
}
	
MTS_IMPLEMENT_CLASS(BSDF, true, ConfigurableObject)
MTS_NAMESPACE_END

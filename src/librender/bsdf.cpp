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

#include <mitsuba/render/scene.h>
#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN

BSDF::BSDF(const Properties &props)
 : ConfigurableObject(props), m_name(props.getID()) {
}

BSDF::BSDF(Stream *stream, InstanceManager *manager) 
 : ConfigurableObject(stream, manager) {
	m_name = stream->readString();
}

BSDF::~BSDF() { }

void BSDF::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
	stream->writeString(m_name);
}

void BSDF::setParent(ConfigurableObject *parent) {
	/* BSDF's don't need to reference their parent -> do nothing */
}

void BSDF::addChild(const std::string &name, ConfigurableObject *obj) {
	ConfigurableObject::addChild(name, obj);
}

void BSDF::configure() {
	m_combinedType = 0;
	for (size_t i=0; i<m_components.size(); ++i)
		m_combinedType |= m_components[i];
}

Spectrum BSDF::getDiffuseReflectance(const Intersection &its) const {
	BSDFQueryRecord bRec(its);
	bRec.wi = bRec.wo = Vector(0, 0, 1);
	bRec.typeMask = EDiffuseReflection;
	return eval(bRec) * M_PI;
}

static std::string typeMaskToString(unsigned int typeMask) {
	std::ostringstream oss;
	oss << "{ ";
	if (typeMask & BSDF::EAll) { oss << "all "; typeMask &= ~BSDF::EAll; }
	if (typeMask & BSDF::ESmooth) { oss << "smooth "; typeMask &= ~BSDF::ESmooth; }
	if (typeMask & BSDF::EDiffuse) { oss << "diffuse "; typeMask &= ~BSDF::EDiffuse; }
	if (typeMask & BSDF::EGlossy) { oss << "glossy "; typeMask &= ~BSDF::EGlossy; }
	if (typeMask & BSDF::EDelta) { oss << "delta"; typeMask &= ~BSDF::EDelta; }
	if (typeMask & BSDF::EDelta1D) { oss << "delta1D "; typeMask &= ~BSDF::EDelta1D; }
	if (typeMask & BSDF::EDiffuseReflection) { oss << "diffuseReflection "; typeMask &= ~BSDF::EDiffuseReflection; }
	if (typeMask & BSDF::EDiffuseTransmission) { oss << "diffuseTransmission "; typeMask &= ~BSDF::EDiffuseTransmission; }
	if (typeMask & BSDF::EGlossyReflection) { oss << "glossyReflection "; typeMask &= ~BSDF::EGlossyReflection; }
	if (typeMask & BSDF::EGlossyTransmission) { oss << "glossyTransmission "; typeMask &= ~BSDF::EGlossyTransmission; }
	if (typeMask & BSDF::EDeltaReflection) { oss << "deltaReflection "; typeMask &= ~BSDF::EDeltaReflection; }
	if (typeMask & BSDF::EDeltaTransmission) { oss << "deltaTransmission "; typeMask &= ~BSDF::EDeltaTransmission; }
	if (typeMask & BSDF::EDelta1DReflection) { oss << "delta1DReflection "; typeMask &= ~BSDF::EDelta1DReflection; }
	if (typeMask & BSDF::EDelta1DTransmission) { oss << "delta1DTransmission "; typeMask &= ~BSDF::EDelta1DTransmission; }
	if (typeMask & BSDF::EAnisotropic) { oss << "anisotropic "; typeMask &= ~BSDF::EAnisotropic; }
	if (typeMask & BSDF::EFrontSide) { oss << "frontSide "; typeMask &= ~BSDF::EFrontSide; }
	if (typeMask & BSDF::EBackSide) { oss << "backSide "; typeMask &= ~BSDF::EBackSide; }
	if (typeMask & BSDF::ECanUseSampler) { oss << "canUseSampler "; typeMask &= ~BSDF::ECanUseSampler; }
	oss << "}";
	return oss.str();
}

std::string BSDFQueryRecord::toString() const {
	std::ostringstream oss;
	oss << "BSDFQueryRecord[" << std::endl
		<< "  wi = " << wi.toString() << "," << std::endl
		<< "  wo = " << wo.toString() << "," << std::endl
		<< "  quantity = " << quantity << "," << std::endl
		<< "  typeMask = " << typeMaskToString(typeMask) << "," << std::endl
		<< "  sampledType = " << typeMaskToString(sampledType) << "," << std::endl
		<< "  component = " << component << "," << std::endl
		<< "  sampledComponent = " << sampledComponent << "," << std::endl
		<< "]";
	return oss.str();
}
	
MTS_IMPLEMENT_CLASS(BSDF, true, ConfigurableObject)
MTS_NAMESPACE_END

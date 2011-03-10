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

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>

MTS_NAMESPACE_BEGIN

Medium::Medium(const Properties &props)
 : NetworkedObject(props) {
	Spectrum defaultSigmaS, defaultSigmaA;
	defaultSigmaA.fromLinearRGB(0.0014f, 0.0025f, 0.0142f);
	defaultSigmaS.fromLinearRGB(0.7f, 1.22f, 1.9f);

	if (props.hasProperty("sizeMultiplier"))
		Log(EError, "Deprecation error: the parameter sizeMultiplier"
			" has been renamed to densityMultiplier");

	m_densityMultiplier = props.getFloat("densityMultiplier", 1);
	m_sigmaA = props.getSpectrum("sigmaA", defaultSigmaA);
	m_sigmaS = props.getSpectrum("sigmaS", defaultSigmaS);
	m_sigmaA *= m_densityMultiplier;
	m_sigmaS *= m_densityMultiplier;
	m_sigmaT = m_sigmaA + m_sigmaS;
	m_albedo = (m_sigmaS/m_sigmaT).max();
}

Medium::Medium(Stream *stream, InstanceManager *manager)
 : NetworkedObject(stream, manager) {
	m_aabb = AABB(stream);
	m_densityMultiplier = stream->readFloat();
	m_sigmaA = Spectrum(stream);
	m_sigmaS = Spectrum(stream);
	m_sigmaT = m_sigmaA + m_sigmaS;
	m_albedo = (m_sigmaS/m_sigmaT).max();
	m_phaseFunction = static_cast<PhaseFunction *>(manager->getInstance(stream));
}
	
void Medium::addChild(const std::string &name, ConfigurableObject *child) {
	const Class *cClass = child->getClass();

	if (cClass->derivesFrom(MTS_CLASS(PhaseFunction))) {
		Assert(m_phaseFunction == NULL);
		m_phaseFunction = static_cast<PhaseFunction *>(child);
	} else {
		Log(EError, "Medium: Invalid child node! (\"%s\")",
			cClass->getName().c_str());
	}
}

void Medium::configure() {
	if (m_phaseFunction == NULL) {
		m_phaseFunction = static_cast<PhaseFunction *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(PhaseFunction), Properties("isotropic")));
	}
}

void Medium::serialize(Stream *stream, InstanceManager *manager) const {
	NetworkedObject::serialize(stream, manager);
	m_aabb.serialize(stream);
	stream->writeFloat(m_densityMultiplier);
	m_sigmaA.serialize(stream);
	m_sigmaS.serialize(stream);
	manager->serialize(stream, m_phaseFunction.get());
}

std::string MediumSamplingRecord::toString() const {
	std::ostringstream oss;
	oss << "MediumSamplingRecord[" << std::endl
		<< "  t = " << t << "," << std::endl
		<< "  p = " << p.toString() << "," << std::endl
		<< "  sigmaA = " << sigmaA.toString() << "," << std::endl
		<< "  sigmaS = " << sigmaS.toString() << "," << std::endl
		<< "  pdfFailure = " << pdfFailure << "," << std::endl
		<< "  pdfSuccess = " << pdfSuccess << "," << std::endl
		<< "  pdfSuccessRev = " << pdfSuccessRev << "," << std::endl
		<< "  transmittance = " << transmittance.toString()
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(Medium, true, NetworkedObject)
MTS_NAMESPACE_END

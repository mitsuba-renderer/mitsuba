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

MTS_NAMESPACE_BEGIN

/**
 * Perfect specular BRDF (i.e. an ideal mirror)
 */
class Mirror : public BSDF {
public:
	Mirror(const Properties &props) 
		: BSDF(props) {
		m_reflectance = props.getSpectrum("specularReflectance", Spectrum(0.8f));
		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDeltaReflection;
		m_combinedType = m_type[0];
		m_usesRayDifferentials = false;
	}

	Mirror(Stream *stream, InstanceManager *manager) 
		: BSDF(stream, manager) {
		m_reflectance = Spectrum(stream);
		m_componentCount = 1;
		m_type = new unsigned int[m_componentCount];
		m_type[0] = EDeltaReflection;
		m_combinedType = m_type[0];
		m_usesRayDifferentials = false;
	}

	virtual ~Mirror() {
		delete[] m_type;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		m_reflectance.serialize(stream);
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return Spectrum(0.0f);
	}
	
	Spectrum f(const BSDFQueryRecord &bRec) const {
		return Spectrum(0.0f);
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		return 0.0f;
	}

	inline void reflect(const Vector &wi, Vector &wo) const {
		wo = Vector(-wi.x, -wi.y, wi.z);
	}

	inline Spectrum sample(BSDFQueryRecord &bRec) const {
		Float pdf;
		return Mirror::sample(bRec, pdf);
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf) const {
		if (!(bRec.typeMask & m_combinedType))
			return Spectrum(0.0f);
		reflect(bRec.wi, bRec.wo);
		bRec.sampledComponent = 0;
		bRec.sampledType = EDeltaReflection;
		pdf = std::abs(Frame::cosTheta(bRec.wo));
		return m_reflectance;
	}

	Float pdfDelta(const BSDFQueryRecord &bRec) const {
		return std::abs(Frame::cosTheta(bRec.wo));
	}
	
	Spectrum fDelta(const BSDFQueryRecord &bRec) const {
		return m_reflectance;
	}


	MTS_DECLARE_CLASS()
protected:
	Spectrum m_reflectance;
};


MTS_IMPLEMENT_CLASS_S(Mirror, false, BSDF)
MTS_EXPORT_PLUGIN(Mirror, "Mirror BRDF");
MTS_NAMESPACE_END

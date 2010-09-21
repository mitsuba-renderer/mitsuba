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

#include <mitsuba/render/scene.h>
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/**
 * Portal luminaire -- can be used to turn a surface into a
 * portal that exposes luminaires behind it, such as an 
 * environment map. This is often necessary to make interior 
 * lighting work efficiently enough when using algorithms like 
 * path tracing or photon mapping.
 */
class PortalLuminaire : public Luminaire {
public:
	PortalLuminaire(const Properties &props) : Luminaire(props), m_shape(NULL) {
		AssertEx(m_luminaireToWorld.isIdentity(), "Error: non-identity transformation found. "
			"Portal luminaires inherit their transformation from the associated shape!");
		m_type = EDiffuseDirection | EOnSurface;
		m_intersectable = true;
	}

	virtual ~PortalLuminaire() {
		for (size_t i=0; i<m_luminaires.size(); ++i)
			m_luminaires[i]->decRef();
	}

	PortalLuminaire(Stream *stream, InstanceManager *manager) 
		: Luminaire(stream, manager) {
		m_shape = static_cast<Shape *>(manager->getInstance(stream));
		int luminaireCount = stream->readInt();
		for (int i=0; i<luminaireCount; ++i)
			addChild("", static_cast<Luminaire *>(manager->getInstance(stream)));
		configure();
	}
	
	void preprocess(const Scene *scene) {
		for (size_t i=0; i<m_luminaires.size(); ++i)
			m_luminaires[i]->preprocess(scene);
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Luminaire::serialize(stream, manager);

		manager->serialize(stream, m_shape);
		stream->writeInt((int) m_luminaires.size());
		for (size_t i=0; i<m_luminaires.size(); ++i)
			manager->serialize(stream, m_luminaires[i]);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();
		if (cClass->derivesFrom(Luminaire::m_theClass)) {
			Luminaire *luminaire = static_cast<Luminaire *>(child);
			m_luminaires.push_back(luminaire);
			luminaire->incRef();
		} else {
			Luminaire::addChild(name, child);
		}
	}

	void setParent(ConfigurableObject *parent) {
		ConfigurableObject::setParent(parent);

		if (parent->getClass()->derivesFrom(Shape::m_theClass)) {
			m_shape = static_cast<Shape *>(parent);
			parent->configure();
			m_surfaceArea = m_shape->getSurfaceArea();
		} else {
			Log(EError, "A portal light source must be child of a shape instance");
		}
	}

	void configure() {
		m_power = Spectrum(0.0f);
		if (m_luminaires.size() == 0)
			Log(EError, "Portal luminaire must have one or more child luminaires!");
		for (size_t i=0; i<m_luminaires.size(); ++i) {
			m_luminaires[i]->configure();
			m_power += m_luminaires[i]->getPower();
		}
	}

	Spectrum getPower() const {
		return m_power;
	}
	
	Spectrum Le(const EmissionRecord &eRec) const {
		Spectrum result(0.0f);

		if (dot(eRec.d, eRec.sRec.n) <= 0)
			return Spectrum(0.0f);

		for (size_t i=0; i<m_luminaires.size(); ++i)
			result += m_luminaires[i]->Le(Ray(eRec.sRec.p, -eRec.d));

		return result;
	}

	Spectrum Le(const LuminaireSamplingRecord &lRec) const {
		Spectrum result(0.0f);

		if (dot(lRec.d, lRec.sRec.n) <= 0)
			return Spectrum(0.0f);

		for (size_t i=0; i<m_luminaires.size(); ++i)
			result += m_luminaires[i]->Le(Ray(lRec.sRec.p, -lRec.d));

		return result;
	}

	inline void sample(const Point &p, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		lRec.pdf = m_shape->sampleSolidAngle(lRec.sRec, p, sample);
		lRec.d = p - lRec.sRec.p;

		if (EXPECT_TAKEN(lRec.pdf > 0 && dot(lRec.d, lRec.sRec.n) > 0)) {
			lRec.d = normalize(lRec.d);
			lRec.Le = Le(lRec);
		} else {
			lRec.pdf = 0;
		}
	}

	inline void sample(const Intersection &its, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		PortalLuminaire::sample(its.p, lRec, sample);
	}

	inline Float pdf(const Point &p, const LuminaireSamplingRecord &lRec) const {
		return m_shape->pdfSolidAngle(lRec.sRec, p);
	}

	Float pdf(const Intersection &its, const LuminaireSamplingRecord &lRec) const {
		return pdf(its.p, lRec);
	}

	void sampleEmission(EmissionRecord &eRec,
		const Point2 &sample1, const Point2 &sample2) const {
		eRec.pdfArea = m_shape->sampleArea(eRec.sRec, sample1);
		Vector wo = squareToHemispherePSA(sample2);
		eRec.pdfDir = Frame::cosTheta(wo) * INV_PI;
		eRec.d = Frame(eRec.sRec.n).toWorld(wo);
		eRec.P = Le(eRec);
	}

	void sampleEmissionArea(EmissionRecord &eRec, const Point2 &sample) const {
		eRec.pdfArea = m_shape->sampleArea(eRec.sRec, sample);
		eRec.P = Spectrum(1.0f);
	}

	Spectrum fArea(const EmissionRecord &eRec) const {
		return Spectrum(1.0f);
	}

	Spectrum sampleEmissionDirection(EmissionRecord &eRec, const Point2 &sample) const {
		Vector wo = squareToHemispherePSA(sample);
		eRec.d = Frame(eRec.sRec.n).toWorld(wo);
		eRec.pdfDir = Frame::cosTheta(wo) * INV_PI;
		return Le(eRec);
	}

	Spectrum f(const EmissionRecord &eRec) const {
		Float dp = dot(eRec.sRec.n, eRec.d);
		if (dp > 0)
			return Le(eRec);
		else
			return Spectrum(0.0f);
	}

	void pdfEmission(EmissionRecord &eRec) const {
		Float dp = dot(eRec.sRec.n, eRec.d);
		if (dp > 0)
			eRec.pdfDir = dp * INV_PI;
		else {
			eRec.pdfDir = 0;
		}
		eRec.pdfArea = m_shape->pdfArea(eRec.sRec);
	}

	std::string toString() const {
		return "PortalLuminaire[]";
	}
	
	MTS_DECLARE_CLASS()
private:
	const Shape *m_shape;
	Spectrum m_power;
	std::vector<Luminaire *> m_luminaires;
};

MTS_IMPLEMENT_CLASS_S(PortalLuminaire, false, Luminaire)
MTS_EXPORT_PLUGIN(PortalLuminaire, "Portal luminaire");
MTS_NAMESPACE_END

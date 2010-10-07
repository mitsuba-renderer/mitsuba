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
#include "maxexp.h"

MTS_NAMESPACE_BEGIN

/**
 * Homogeneous participating medium. An arbitrary (manifold) shape
 * must be specified as a child object.
 */
class HomogeneousMedium : public Medium {
public:
	enum ESamplingStrategy {
		ERandom,
		ESingle,
		EBalance,
		EManual,
		EMaximum
	};

	HomogeneousMedium(const Properties &props) 
		: Medium(props) {
		std::string strategy = props.getString("strategy", "random");

		if (strategy == "random")
			m_strategy = ERandom;
		else if (strategy == "single")
			m_strategy = ESingle;
		else if (strategy == "balance")
			m_strategy = EBalance;
		else if (strategy == "maximum")
			m_strategy = EMaximum;
		else if (strategy == "manual")
			m_strategy = EManual;
		else
			Log(EError, "Specified an unknown sampling strategy");

		m_channel = props.getInteger("channel", 0);
		m_sigma = props.getFloat("sigma", -1); // for the manual strategy
		Assert(m_channel >= 0 && m_channel < SPECTRUM_SAMPLES);
		if (props.getBoolean("mono", false)) {
			m_sigmaA = Spectrum(m_sigmaA[m_channel]);
			m_sigmaS = Spectrum(m_sigmaS[m_channel]);
			m_sigmaT = m_sigmaA + m_sigmaS;
		}

		std::vector<Float> coeffs(SPECTRUM_SAMPLES);
		for (int i=0; i<SPECTRUM_SAMPLES; ++i)
			coeffs[i] = m_sigmaT[i];

		m_maxExpDist = new MaxExpDist(coeffs);
		m_kdTree = new KDTree();
	}

	HomogeneousMedium(Stream *stream, InstanceManager *manager)
		: Medium(stream, manager) {
		m_kdTree = new KDTree();
		size_t shapeCount = stream->readUInt();
		for (size_t i=0; i<shapeCount; ++i) 
			addChild("", static_cast<Shape *>(manager->getInstance(stream)));
		m_channel = stream->readInt();
		m_strategy = (ESamplingStrategy) stream->readInt();
		m_sigma = stream->readFloat();

		std::vector<Float> coeffs(SPECTRUM_SAMPLES);
		for (int i=0; i<SPECTRUM_SAMPLES; ++i)
			coeffs[i] = m_sigmaT[i];

		configure();
		m_maxExpDist = new MaxExpDist(coeffs);
	}

	virtual ~HomogeneousMedium() {
		for (size_t i=0; i<m_shapes.size(); ++i)
			m_shapes[i]->decRef();
		delete m_maxExpDist;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Medium::serialize(stream, manager);
		stream->writeUInt((uint32_t) m_shapes.size());
		for (size_t i=0; i<m_shapes.size(); ++i)
			manager->serialize(stream, m_shapes[i]);
		stream->writeInt(m_channel);
		stream->writeInt(m_strategy);
		stream->writeFloat(m_sigma);
	}
	
	void configure() {
		Medium::configure();
		if (m_shapes.size() == 0)
			Log(EError, "This medium requires one or more Shape instance as a child");
		m_kdTree->build();
		m_aabb = m_kdTree->getAABB();
	}

	bool isInside(const Ray &r) const {
		Ray ray(r(r.mint + Epsilon), r.d);
		Intersection its;
		if (!m_kdTree->rayIntersect(ray, its))
			return false;
		return dot(ray.d, its.geoFrame.n) > 0;
	}

	Spectrum tau(const Ray &r) const {
		Float dLength = r.d.length();
		Ray ray(r(r.mint), r.d / dLength, Epsilon, std::numeric_limits<Float>::infinity());
		Float coveredLength = 0, remaining = (r.maxt - r.mint) * dLength;
		bool inside = isInside(r);
		Intersection its;
		int iterations = 0;

		while (remaining > 0 && m_kdTree->rayIntersect(ray, its)) {
			if (inside)
				coveredLength += std::min(remaining, its.t);
			remaining -= its.t;
			inside = !inside;
			ray.o = its.p;

			if (++iterations > 20) {
				/// Just a precaution..
				Log(EWarn, "tau(): round-off error issues?");
				break;
			}
		}

		return m_sigmaT * coveredLength;
	}

	bool sampleDistance(const Ray &theRay, Float distSurf, 
			MediumSamplingRecord &mRec, Sampler *sampler) const {
		Intersection its;
		Ray ray(theRay.o, theRay.d);
		int iterations = 0;

		/* Check if the start of the ray is already inside the medium */
		bool inside = isInside(ray);
		Point orig(theRay.o);

		/* Remaining distance in the medium */
		Float distMed;
		Float sigmaT;

		switch (m_strategy) {
			case EManual: {
					sigmaT = m_sigma;
					Float desiredAttenuation = 1 - sampler->next1D();
					distMed = -std::log(desiredAttenuation) / sigmaT;
					mRec.pdf = sigmaT * std::exp(-sigmaT * distMed);
					mRec.miWeight = 1;
				}
				break;
			case ESingle: {
					sigmaT = m_sigmaT[m_channel];
					Float desiredAttenuation = 1 - sampler->next1D();
					distMed = -std::log(desiredAttenuation) / sigmaT;
					mRec.pdf = sigmaT * std::exp(-sigmaT * distMed);
					mRec.miWeight = 1;
				}
				break;
			case ERandom: {
					Float rv = sampler->next1D();
					int channel = std::min(
						(int) (rv * SPECTRUM_SAMPLES), SPECTRUM_SAMPLES-1);
					sigmaT = m_sigmaT[channel];
					Float desiredAttenuation = 1 - sampler->next1D();
					distMed = -std::log(desiredAttenuation) / sigmaT;
					mRec.pdf = sigmaT * std::exp(-sigmaT * distMed);
					mRec.miWeight = 1;
				}
				break;
			case EBalance: {
					Float rv = sampler->next1D();
					int channel = std::min(
						(int) (rv * SPECTRUM_SAMPLES), SPECTRUM_SAMPLES-1);
					sigmaT = m_sigmaT[channel];
					Float desiredAttenuation = 1 - sampler->next1D();
					distMed = -std::log(desiredAttenuation) / sigmaT;
					Spectrum pdf = m_sigmaT * (Spectrum(-distMed) * m_sigmaT).exp();
					mRec.pdf = pdf[channel] / SPECTRUM_SAMPLES;

					Float sum = 0;
					for (int i=0; i<SPECTRUM_SAMPLES; ++i)
						sum += pdf[i];
					mRec.miWeight = pdf[channel]/sum;
				}
				break;
			case EMaximum: {
					Float rv = sampler->next1D();
					distMed = m_maxExpDist->sample(rv, mRec.pdf);
					mRec.miWeight = 1;
					sigmaT = -1; // make the compiler happy
				};
				break;
			default:
				Log(EError, "Unknown sampling strategy!");
				return false;
		}

		Float traveled = 0,  // Traveled ray distance
			  covered = 0;   // Distance covered by the medium

		while (m_kdTree->rayIntersect(ray, its)) {
			if (inside) {
				/* Moving through the medium */
				if (its.t > distMed && distMed < distSurf) {
					/* A medium interaction occurred */
					mRec.p = ray(distMed);
					mRec.sigmaA = m_sigmaA;
					mRec.sigmaS = m_sigmaS;
					mRec.albedo = m_albedo;
					mRec.medium = this;
					mRec.attenuation = (m_sigmaT * (-covered-distMed)).exp();
					mRec.t = traveled + distMed;
					return true;
				} else if (its.t > distSurf) {
					/* A surface interaction occurred */
					covered += distSurf;
					traveled += distSurf;
					break;
				} else {
					/* Still moving through the medium */
					covered += its.t;
					distMed -= its.t;
				}
			} else {
				/* Moving through space outside of the medium */
				if (its.t > distSurf) {
					/* A surface interaction occurred */
					traveled += distSurf;
					break;
				}
			}
			traveled += its.t;
			distSurf -= its.t;
			inside = !inside;
			ray.o = its.p;

			if (++iterations > 10) {
				/// Just a precaution..
				Log(EWarn, "selectDistance(): round-off error issues?");
				break;
			}
		}

		/* There was no medium interaction inside the permitted ray interval.
		   This occurred with the probability tau[channel](0, maxDist) */
		if (m_strategy != EMaximum)
			mRec.pdf = std::exp(-sigmaT * covered);
		else
			mRec.pdf = 1-m_maxExpDist->cdf(covered);
		mRec.t = traveled;
		mRec.attenuation = (m_sigmaT * (-covered)).exp();
		mRec.miWeight = 1;

		return false;
	}

	void setParent(ConfigurableObject *parent) {
		if (parent->getClass()->derivesFrom(Shape::m_theClass))
			Log(EError, "Medium cannot be a parent of a shape");
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(Shape::m_theClass)) {
			Shape *shape = static_cast<Shape *>(child);
			if (shape->isCompound()) {
				int ctr = 0;
				while (true) {
					ref<Shape> childShape = shape->getElement(ctr++);
					if (!childShape)
						break;
					addChild("", childShape);
				}
			} else {
				m_kdTree->addShape(shape);
				shape->incRef();
				m_shapes.push_back(shape);
			}
		} else {
			Medium::addChild(name, child);
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HomogeneousMedium[" << endl
			<< "  sigmaA = " << m_sigmaA.toString() << "," << std::endl
			<< "  sigmaS = " << m_sigmaS.toString() << "," << std::endl
			<< "  sigmaT = " << m_sigmaT.toString() << "," << std::endl
			<< "  phase = " << indent(m_phaseFunction->toString()) << "," << std::endl 
			<< "  shapes = " << indent(listToString(m_shapes)) << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<KDTree> m_kdTree;
	std::vector<Shape *> m_shapes;
	int m_channel;
	Float m_sigma;
	ESamplingStrategy m_strategy;
	MaxExpDist *m_maxExpDist;
};

MTS_IMPLEMENT_CLASS_S(HomogeneousMedium, false, Medium)
MTS_EXPORT_PLUGIN(HomogeneousMedium, "Homogeneous medium");
MTS_NAMESPACE_END

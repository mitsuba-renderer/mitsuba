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
	/**
	 * The following sampling strategies for choosing a 
	 * suitable distribution in the in-scattering line integral
	 * are available.
	 */
	enum ESamplingStrategy {
		EBalance,  /// Exp. distrib.; pick a random channel each time
		ESingle,  /// Exp. distrib.; pick a specified channel
		EManual,  /// Exp. distrib.; manually specify a value
		EMaximum  /// Max-of-exp. distribs
	};

	HomogeneousMedium(const Properties &props) 
		: Medium(props), m_maxExpDist(NULL) {
		std::string strategy = props.getString("strategy", "balance");

		/**
		 * The goal of the medium sampling weight is to be able to
		 * sample medium intarctions according to
		 *    sigma_s(t) * tau(0 <-> t)
		 * as opposed to
		 *    sigma_t(t) * tau(0 <-> t)
		 * See the separate writeup for more details.
		 */
		m_mediumSamplingWeight = props.getFloat("mediumSamplingWeight", -1);
		if (m_mediumSamplingWeight == -1) {
			for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
				/// Record the highest albedo values across channels
				Float albedo = m_sigmaS[i] / m_sigmaT[i];
				if (albedo > m_mediumSamplingWeight)
					m_mediumSamplingWeight = albedo;
			}
			if (m_mediumSamplingWeight > 0) {
				/* The medium scatters some light -> place at least half 
				   of the samples in it, otherwise we will render lots
				   of spatially varying noise where one pixel has a
				   medium interaction while the neighbors don't */
				m_mediumSamplingWeight = std::max(m_mediumSamplingWeight, 
					(Float) 0.5f);
			}
		}

		if (strategy == "balance") {
			m_strategy = EBalance;
		} else if (strategy == "single") {
			m_strategy = ESingle;

			int channel = 0;
			Float smallest = std::numeric_limits<Float>::infinity();
			for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
				if (m_sigmaT[i] < smallest) {
					smallest = m_sigmaT[i];
					channel = i;
				}
			}

			channel = props.getInteger("channel", channel);
			Assert(channel >= 0 && channel < SPECTRUM_SAMPLES);
			m_samplingDensity = m_sigmaT[channel];

			if (props.getBoolean("monochromatic", false)) {
				/* Optionally configure as a monochromatic medium */
				m_sigmaA = Spectrum(m_sigmaA[channel]);
				m_sigmaS = Spectrum(m_sigmaS[channel]);
				m_sigmaT = m_sigmaA + m_sigmaS;
			}
		} else if (strategy == "maximum") {
			m_strategy = EMaximum;
			std::vector<Float> coeffs(SPECTRUM_SAMPLES);
			for (int i=0; i<SPECTRUM_SAMPLES; ++i)
				coeffs[i] = m_sigmaT[i];
			m_maxExpDist = new MaxExpDist(coeffs);
		} else if (strategy == "manual") {
			m_strategy = EManual;
			m_samplingDensity= props.getFloat("sigma");
		} else {
			Log(EError, "Specified an unknown sampling strategy");
		}

		m_kdTree = new ShapeKDTree();
	}

	HomogeneousMedium(Stream *stream, InstanceManager *manager)
		: Medium(stream, manager), m_maxExpDist(NULL) {
		m_kdTree = new ShapeKDTree();
		size_t shapeCount = stream->readUInt();
		for (size_t i=0; i<shapeCount; ++i) 
			addChild("", static_cast<Shape *>(manager->getInstance(stream)));
		m_strategy = (ESamplingStrategy) stream->readInt();
		m_samplingDensity= stream->readFloat();
		m_mediumSamplingWeight = stream->readFloat();

		if (m_strategy == EMaximum) {
			std::vector<Float> coeffs(SPECTRUM_SAMPLES);
			for (int i=0; i<SPECTRUM_SAMPLES; ++i)
				coeffs[i] = m_sigmaT[i];
			m_maxExpDist = new MaxExpDist(coeffs);
		}

		configure();
	}

	virtual ~HomogeneousMedium() {
		for (size_t i=0; i<m_shapes.size(); ++i)
			m_shapes[i]->decRef();
		if (m_maxExpDist)
			delete m_maxExpDist;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Medium::serialize(stream, manager);
		stream->writeUInt((uint32_t) m_shapes.size());
		for (size_t i=0; i<m_shapes.size(); ++i)
			manager->serialize(stream, m_shapes[i]);
		stream->writeInt(m_strategy);
		stream->writeFloat(m_samplingDensity);
		stream->writeFloat(m_mediumSamplingWeight);
	}
	
	void configure() {
		Medium::configure();
		if (m_shapes.size() == 0)
			Log(EError, "This medium requires one or more Shape instance as a child");
		m_kdTree->build();
		m_aabb = m_kdTree->getAABB();
	}

	Float getCoveredLength(const Ray &r) const {
		Ray ray(r(r.mint), r.d,
			0, std::numeric_limits<Float>::infinity(), 0.0f);
		Float coveredLength = 0, remaining = r.maxt - r.mint;
		Intersection its;
		int iterations = 0;

		while (remaining > 0 && m_kdTree->rayIntersect(ray, its)) {
			bool inside = dot(its.geoFrame.n, ray.d) > 0;
			if (inside)
				coveredLength += std::min(remaining, its.t);
			remaining -= its.t;
			ray.o = its.p;
			ray.mint = Epsilon;

			if (++iterations > 10) {
				/// Just a precaution..
				Log(EWarn, "getCoveredLength(): round-off error issues?");
				break;
			}
		}
		return coveredLength;
	}

	Spectrum tau(const Ray &ray) const {
		return (m_sigmaT * (-getCoveredLength(ray))).exp();
	}

	bool sampleDistance(const Ray &r, MediumSamplingRecord &mRec, 
			Sampler *sampler) const {
		Intersection its;
		Ray ray(r(r.mint), r.d,
			0, std::numeric_limits<Float>::infinity(), 0.0f);
		Float distSurf = r.maxt - r.mint;
		int iterations = 0;

		/* Remaining distance in the medium */
		Float distMed, samplingDensity = m_samplingDensity;
		Float rand = sampler->next1D();
		if (m_strategy == EBalance) {
			int channel = std::min((int) (sampler->next1D() 
				* SPECTRUM_SAMPLES), SPECTRUM_SAMPLES-1);
			samplingDensity = m_sigmaT[channel];
		}

		if (rand > m_mediumSamplingWeight) {
			/* Don't generate a medium interaction, but still
			   compute the attenuation etc */
			distMed = std::numeric_limits<Float>::infinity();
		} else {
			rand /= m_mediumSamplingWeight;
			if (m_strategy != EMaximum) {
				distMed = -std::log(1-rand) / samplingDensity;
			} else {
				distMed = m_maxExpDist->sample(1-rand, mRec.pdfSuccess);
			}
		}

		Float traveled = 0,  // Traveled ray distance
			  coveredLength = 0;   // Distance covered by the medium

		bool success = false;

		while (m_kdTree->rayIntersect(ray, its)) {
			bool inside = dot(its.geoFrame.n, ray.d) > 0;
			if (inside) {
				/* Moving through the medium */
				if (its.t > distMed && distMed < distSurf) {
					/* A medium interaction occurred */
					mRec.p = ray(distMed);
					mRec.sigmaA = m_sigmaA;
					mRec.sigmaS = m_sigmaS;
					mRec.albedo = m_mediumSamplingWeight;
					mRec.medium = this;
					coveredLength += distMed;
					traveled += distMed;
					success = true;
					break;
				} else if (its.t > distSurf) {
					/* A surface interaction occurred */
					coveredLength += distSurf;
					traveled += distSurf;
					break;
				} else {
					/* Still moving through the medium */
					coveredLength += its.t;
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
			ray.o = its.p;
			ray.mint = Epsilon;

			if (++iterations > 10) {
				/// Just a precaution..
				Log(EWarn, "sampleDistance(): round-off error issues?");
				break;
			}
		}
		mRec.attenuation = (m_sigmaT * (-coveredLength)).exp();
		mRec.t = traveled;

		switch (m_strategy) {
			case EMaximum:
				mRec.pdfFailure = 1-m_maxExpDist->cdf(coveredLength);
				break;

			case EBalance: 
				mRec.pdfFailure = 0;
				mRec.pdfSuccess = 0;
				for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
					Float tmp = std::exp(-m_sigmaT[i] * coveredLength);
					mRec.pdfFailure += tmp;
					mRec.pdfSuccess += m_sigmaT[i] * tmp;
				}
				mRec.pdfFailure /= SPECTRUM_SAMPLES;
				mRec.pdfSuccess /= SPECTRUM_SAMPLES;
				break;

			case ESingle:
			case EManual:
				mRec.pdfFailure = std::exp(-samplingDensity * coveredLength);
				mRec.pdfSuccess = samplingDensity * mRec.pdfFailure;
				break;

			default:
				Log(EError, "Unknown sampling strategy!");
		}

		mRec.pdfSuccessRev = mRec.pdfSuccess = mRec.pdfSuccess * m_mediumSamplingWeight;
		mRec.pdfFailure = m_mediumSamplingWeight * mRec.pdfFailure + (1-m_mediumSamplingWeight);

		return success;
	}

	void pdfDistance(const Ray &ray, Float t, MediumSamplingRecord &mRec) const {
		Float coveredLength = getCoveredLength(Ray(ray, 0, t));
		switch (m_strategy) {
			case EManual: 
			case ESingle: {
					Float temp = std::exp(-m_samplingDensity * coveredLength);
					mRec.pdfSuccess = m_samplingDensity * temp;
					mRec.pdfFailure = temp;
				}
				break;
			case EBalance: {
					mRec.pdfSuccess = 0;
					mRec.pdfFailure = 0;
					for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
						Float temp = std::exp(-m_sigmaT[i] * coveredLength);
						mRec.pdfSuccess += m_sigmaT[i] * temp;
						mRec.pdfFailure += temp;
					}
					mRec.pdfSuccess /= SPECTRUM_SAMPLES;
					mRec.pdfFailure /= SPECTRUM_SAMPLES;
				}
				break;
			case EMaximum:
				mRec.pdfSuccess = m_maxExpDist->pdf(coveredLength);
				mRec.pdfFailure = 1-m_maxExpDist->cdf(coveredLength);
				break;
			default:
				Log(EError, "Unknown sampling strategy!");
		}
		mRec.attenuation = (m_sigmaT * (-coveredLength)).exp();
		mRec.pdfSuccess = mRec.pdfSuccessRev = mRec.pdfSuccess * m_mediumSamplingWeight;
		mRec.pdfFailure = mRec.pdfFailure * m_mediumSamplingWeight + (1-m_mediumSamplingWeight);
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
			<< "  sigmaA = " << m_sigmaA.toString() << "," << endl
			<< "  sigmaS = " << m_sigmaS.toString() << "," << endl
			<< "  sigmaT = " << m_sigmaT.toString() << "," << endl
			<< "  mediumSamplingWeight = " << m_mediumSamplingWeight << "," << endl
			<< "  samplingDensity = " << m_samplingDensity << "," << endl
			<< "  strategy = ";
		switch (m_strategy) {
			case ESingle: oss << "single," << endl; break;
			case EManual: oss << "manual," << endl; break;
			case EBalance: oss << "balance," << endl; break;
			case EMaximum: oss << "maximum," << endl; break;
		}

		oss << "  phase = " << indent(m_phaseFunction.toString()) << "," << endl 
			<< "  shapes = " << indent(listToString(m_shapes)) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	ref<ShapeKDTree> m_kdTree;
	std::vector<Shape *> m_shapes;
	Float m_samplingDensity, m_mediumSamplingWeight;
	ESamplingStrategy m_strategy;
	MaxExpDist *m_maxExpDist;
};

MTS_IMPLEMENT_CLASS_S(HomogeneousMedium, false, Medium)
MTS_EXPORT_PLUGIN(HomogeneousMedium, "Homogeneous medium");
MTS_NAMESPACE_END

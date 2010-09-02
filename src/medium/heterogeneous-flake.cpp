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
#include <mitsuba/render/volume.h>
#include <mitsuba/core/shvector4d.h>
#include <fstream>

/// Header identification record for phase function caches
#define PHASE_FUNCTION_HEADER 0x40044004

MTS_NAMESPACE_BEGIN

struct AbsCos {
	Float operator()(const Vector &w) const { return std::abs(w.z); }
};

struct FlakeDistr {
	Float exp;
	bool fiber;

	FlakeDistr(Float exp, bool fiber) : exp(exp), fiber(fiber) {
	}

	inline Float operator()(const Vector &w) const {
		if (fiber)
			return std::pow(1-std::pow(w.z, 2), exp/2.0f);
		else
			return std::pow(std::abs(w.z), exp);
	}
};

template <typename Distr> struct FlakePhaseFunctor {
	Distr d;

	FlakePhaseFunctor(Distr d) : d(d) { }

	inline Float operator()(const Vector &wi, const Vector &wo) const {
		Vector H = normalize(wi + wo);
		return d(H) + d(-H);
	}
};

class HeterogeneousFlakeMedium;

class FlakePhaseFunction : public PhaseFunction {
public:
	FlakePhaseFunction(const SHVector &sigmaT, int samplingRecursions, const SHVector4D *phaseExpansion, 
		Float exponent, Float normalization, bool fiber, HeterogeneousFlakeMedium *medium) : PhaseFunction(Properties()), 
		m_sigmaT(sigmaT), m_samplingRecursions(samplingRecursions), m_phaseExpansion(phaseExpansion), 
		m_exponent(exponent), m_normalization(normalization), m_fiber(fiber), m_medium(medium) {
		initialize();
	}

	FlakePhaseFunction(Stream *stream, InstanceManager *manager) :
			PhaseFunction(stream, manager) {
		/* Not implemented */
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);
		/* Not implemented */
	}

	void initialize() {
		m_shSampler = new SHSampler(m_phaseExpansion->getBands(), m_samplingRecursions);
		Log(EInfo, "Constructing a SH sampler: %s", m_shSampler->toString().c_str());
	}

	virtual ~FlakePhaseFunction() { }


	Spectrum sample(const MediumSamplingRecord &mRec, const Vector &wi, Vector &wo, 
		ESampledType &sampledType, const Point2 &_sample) const {
		Float pdf;
		Spectrum value = sample(mRec, wi, wo, sampledType, pdf, _sample);
		return value/pdf;
	}

	/* Implementations are below */
	Spectrum f(const MediumSamplingRecord &mRec, const Vector &_wi, const Vector &_wo) const;

	Spectrum sample(const MediumSamplingRecord &mRec, const Vector &_wi, Vector &_wo, 
			ESampledType &sampledType, Float &pdf, 
				const Point2 &_sample) const;

	Float pdf(const MediumSamplingRecord &mRec, const Vector &wi, const Vector &wo) const;

	std::string toString() const {
		return "FlakePhaseFunction[]";
	}

	MTS_DECLARE_CLASS()
private:
	SHVector m_sigmaT;
	ref<SHSampler> m_shSampler;
	int m_samplingRecursions;
	const SHVector4D *m_phaseExpansion;
	mutable PrimitiveThreadLocal<SHVector> m_temp;
	Float m_exponent;
	Float m_normalization;
	bool m_fiber;
	const HeterogeneousFlakeMedium *m_medium; 
};

/*
 * When inverting an integral of the form f(t):=\int_0^t [...] dt
 * using composite Simpson's rule quadrature, the following constant
 * specified how many times the step size will be reduced when we
 * have almost reached the desired value.
 */
#define NUM_REFINES 10

class HeterogeneousFlakeMedium : public Medium {
public:
	HeterogeneousFlakeMedium(const Properties &props) 
		: Medium(props) {
		m_stepSize = props.getFloat("stepSize", 0);

		/* Flake model configuration */
		int bands = props.getInteger("bands", 7);
		int phaseBands = props.getInteger("phaseBands", bands);
		Assert(bands > 0 && phaseBands > 0);

		m_samplingRecursions = props.getInteger("samplingRecursions", 10);
		m_fiber = props.getBoolean("fiber", true);
		const int numSamples = 80;
		ref<FileStream> stream;

		SHVector absCos = SHVector(bands);
		absCos.project(AbsCos(), numSamples);
		Float absCosError = absCos.l2Error(AbsCos(), numSamples);

		m_exponent = props.getFloat("exponent", 4);
		FlakeDistr distrFunc(m_exponent, m_fiber);
		m_D = SHVector(bands);
		m_D.project(distrFunc, numSamples);
		m_normalization = 1 / (m_D(0, 0) * 2 * std::sqrt(M_PI));
		Float dError = m_D.l2Error(distrFunc, numSamples);
		m_D.normalize();

		Log(EInfo, "=======  Flake medium properties =======");
		Log(EInfo, "       distribution type = %s", m_fiber ? " sin^p" : "cos^p");
		Log(EInfo, "   distribution exponent = %f", m_exponent);
		Log(EInfo, "    SH bands (D, sigmaT) = %i", bands);
		Log(EInfo, "    SH bands (phase fct) = %i", phaseBands);
		Log(EInfo, "  |cos| projection error = %f", absCosError);
		Log(EInfo, "     D  projection error = %f", dError);
		Log(EInfo, "     sampling recursions = %i", m_samplingRecursions);

		bool computePhaseProjection = true;
		if (FileStream::exists("flake-phase.dat")) {
			/* Avoid recomputing this every time */
			stream = new FileStream("flake-phase.dat", FileStream::EReadOnly);
			unsigned int header = stream->readUInt();
			int nBands = stream->readInt();
			bool isFiber = stream->readBool();
			Float exponent = stream->readFloat();
			if (header == PHASE_FUNCTION_HEADER && nBands == phaseBands 
				&& exponent == m_exponent && isFiber == m_fiber) {
				m_phaseExpansion = SHVector4D(stream);
				computePhaseProjection = false;
			}
			stream->close();
		}

		if (computePhaseProjection) {
			int res = 50;
			FlakePhaseFunctor<FlakeDistr> phaseFunctor(distrFunc);
			m_phaseExpansion = SHVector4D(res, 2*res, phaseBands);
			Log(EDebug, "Phase function parameterization: %s", m_phaseExpansion.toString().c_str());
			m_phaseExpansion.project(phaseFunctor, numSamples);
			m_phaseExpansion.normalize();
			stream = new FileStream("flake-phase.dat", FileStream::ETruncReadWrite);
			stream->writeInt(PHASE_FUNCTION_HEADER);
			stream->writeInt(phaseBands);
			stream->writeBool(m_fiber);
			stream->writeFloat(m_exponent);
			m_phaseExpansion.serialize(stream);
			stream->close();
		}
		
		for (int i=0; i<phaseBands; ++i)
			Log(EInfo, "Energy in phase function band %i: %f", i, m_phaseExpansion.energy(i));

		m_sigmaT = SHVector(m_D);
		m_sigmaT.convolve(absCos);
		m_sigmaT *= 2;
		Assert(m_sigmaT.isAzimuthallyInvariant());
	}

	/* Unserialize from a binary data stream */
	HeterogeneousFlakeMedium(Stream *stream, InstanceManager *manager) 
		: Medium(stream, manager) {
		m_shape = static_cast<Shape *>(manager->getInstance(stream));
		m_densities = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_albedo = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_orientations = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_stepSize = stream->readFloat();
		m_exponent = stream->readFloat();
		m_normalization = stream->readFloat();
		m_fiber = stream->readBool();
		m_D = SHVector(stream);
		m_sigmaT = SHVector(stream);
		m_phaseExpansion = SHVector4D(stream);
		m_samplingRecursions = stream->readInt();
		m_kdTree = new KDTree();
		if (m_shape != NULL) {
			m_kdTree->addShape(m_shape.get());
			m_kdTree->build();
		}
		configure();
	}
	
	virtual ~HeterogeneousFlakeMedium() {
	}

	/* Serialize the volume to a binary data stream */
	void serialize(Stream *stream, InstanceManager *manager) const {
		Medium::serialize(stream, manager);
		manager->serialize(stream, m_shape.get());
		manager->serialize(stream, m_densities.get());
		manager->serialize(stream, m_albedo.get());
		manager->serialize(stream, m_orientations.get());
		stream->writeFloat(m_stepSize);
		stream->writeFloat(m_exponent);
		stream->writeFloat(m_normalization);
		stream->writeBool(m_fiber);
		m_D.serialize(stream);
		m_sigmaT.serialize(stream);
		m_phaseExpansion.serialize(stream);
		stream->writeInt(m_samplingRecursions);
	}

	void configure() {
		Medium::configure();
		if (m_densities.get() == NULL)
			Log(EError, "No densities specified!");
		if (m_albedo.get() == NULL)
			Log(EError, "No albedo specified!");
		if (m_orientations.get() == NULL)
			Log(EError, "No orientations specified!");

		if (m_stepSize == 0) {
			m_stepSize = std::min(std::min(
				m_densities->getStepSize(),
				m_albedo->getStepSize()),
				m_orientations->getStepSize());

			if (m_stepSize == std::numeric_limits<Float>::infinity()) 
				Log(EError, "Unable to infer step size, please specify!");
		}

		m_phaseFunction = new FlakePhaseFunction(m_sigmaT, m_samplingRecursions,
			&m_phaseExpansion, m_exponent, m_normalization, m_fiber, this);

		if (m_shape != NULL)
			m_aabb = m_kdTree->getAABB();
		else
			m_aabb = m_densities->getAABB();
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(Shape::m_theClass)) {
			Assert(m_shape == NULL);
			m_shape = static_cast<Shape *>(child);
			m_kdTree = new KDTree();
			m_kdTree->addShape(m_shape.get());
			m_kdTree->build();
		} else if (child->getClass()->derivesFrom(VolumeDataSource::m_theClass)) {
			VolumeDataSource *volume = static_cast<VolumeDataSource *>(child);

			if (name == "albedo") {
				Assert(volume->supportsSpectrumLookups());
				m_albedo = volume;
			} else if (name == "density") {
				Assert(volume->supportsFloatLookups());
				m_densities = volume;
			} else if (name == "orientations") {
				Assert(volume->supportsVectorLookups());
				m_orientations = volume;
			}
		} else {
			Medium::addChild(name, child);
		}
	}

	Float distanceToMediumEntry(const Ray &ray) const {
		if (m_shape != NULL) {
			Ray r(ray, Epsilon, std::numeric_limits<Float>::infinity());
			Intersection its;
			if (!m_kdTree->rayIntersect(r, its)) 
				return std::numeric_limits<Float>::infinity(); 
			return dot(ray.d, its.geoFrame.n) > 0 ? 0 : its.t;
		} else {
			Float mint, maxt;
			if (!m_aabb.rayIntersect(ray, mint, maxt))
				return std::numeric_limits<Float>::infinity(); 
			if (mint <= Epsilon && maxt <= Epsilon)
				return std::numeric_limits<Float>::infinity(); 
			else
				return (mint < 0) ? 0 : mint;
		}
	}

	Float distanceToMediumExit(const Ray &ray) const {
		if (m_shape != NULL) {
			Ray r(ray, Epsilon, std::numeric_limits<Float>::infinity());
			Intersection its;
			if (!m_kdTree->rayIntersect(r, its)) 
				return 0;
			return dot(ray.d, its.geoFrame.n) < 0 ? 0 : its.t;
		} else {
			Float mint, maxt;
			if (!m_aabb.rayIntersect(ray, mint, maxt)) 
				return 0;
			if (mint < Epsilon && maxt > 0)
				return maxt;
			else
				return 0;
		}
	}

	Spectrum tau(const Ray &r) const {
		Float dLength = r.d.length();
		Ray ray(r(r.mint), r.d / dLength);
		Float remaining = (r.maxt - r.mint) * dLength;
		Float integral = 0.0f;
		int iterations = 0;

		Float entry = distanceToMediumEntry(ray);
		while (entry < remaining && entry != std::numeric_limits<Float>::infinity()) {
			ray.o = ray(entry);
			remaining -= entry;
			Float exit = distanceToMediumExit(ray);

			if (exit != std::numeric_limits<Float>::infinity())
				integral += integrateDensities(ray, std::min(remaining, exit));

			remaining -= exit;
			if (remaining <= 0) 
				break;

			ray.o = ray(exit);
			entry = distanceToMediumEntry(ray);
			if (++iterations > 10) {
				/// Just a precaution..
				Log(EWarn, "tau(): round-off error issues?");
				break;
			}
		}

		return Spectrum(integral * m_sizeMultiplier);
	}
	
	Float integrateDensities(Ray ray, Float length) const {
		int nParts = (int) std::ceil(length/m_stepSize);
		nParts += nParts % 2;
		const Float stepSize = length/nParts;
		const Vector increment = ray.d * stepSize;

		Float m=4;
		Float accumulatedTau = 
			sigmaT(ray.o, ray.d) + sigmaT(ray(length), ray.d);

		ray.o += increment;
		for (int i=1; i<nParts; ++i) {
			Float value = sigmaT(ray.o, ray.d);
			accumulatedTau += value*m;
			ray.o += increment;
			m = 6-m;
		}
		accumulatedTau *= stepSize/3;

		return accumulatedTau;
	}

	Float invertDensityIntegral(Ray ray, Float maxDist, 
		Float &accumulatedTau, 
		Float &currentSigmaT, 
		Spectrum &currentAlbedo, 
		Float desiredTau) const {

		int nParts = (int) std::ceil(maxDist/m_stepSize);
		Float stepSize = maxDist/nParts;
		Vector fullIncrement = ray.d * stepSize,
			   halfIncrement = fullIncrement * .5f;
	
		Float lastNode = sigmaT(ray.o, ray.d);
		Float t = 0;

		int numRefines = 0;
		while (t <= maxDist) {
			const Float node1 = lastNode,
						node2 = sigmaT(ray.o + halfIncrement, ray.d),
						node3 = sigmaT(ray.o + fullIncrement, ray.d);
			const Float newDensity = accumulatedTau 
				+ (node1+node2*4+node3) * (stepSize/6);

			if (newDensity > desiredTau) {
				if (++numRefines > NUM_REFINES) {
					currentAlbedo = m_albedo->lookupSpectrum(ray.o+fullIncrement);
					currentSigmaT = node3;
					return t;
				}

				stepSize *= .5f;
				fullIncrement *= .5f;
				halfIncrement *= .5f;
				continue;
			}

			accumulatedTau = newDensity;
			lastNode = node3;
			t += stepSize;
			ray.o += fullIncrement;
		}
		return std::numeric_limits<Float>::infinity();
	}

	bool sampleDistance(const Ray &r, Float maxDist, 
			MediumSamplingRecord &mRec,  Sampler *sampler) const {
		Float dLength = r.d.length();
		Ray ray(r(r.mint), r.d / dLength);
		Float remaining      = (maxDist - r.mint) * dLength,
			  desiredTau     = -std::log(1-sampler->next1D())/m_sizeMultiplier,
			  accumulatedTau = 0.0f,
			  currentSigmaT  = 0.0f;
		Spectrum currentAlbedo(0.0f), integral(0.0f);
		int iterations = 0;
		bool success = false;

		mRec.miWeight = 1;
		Float entry = distanceToMediumEntry(ray);
		while (entry < remaining && entry != std::numeric_limits<Float>::infinity()) {
			ray.o = ray(entry);
			remaining -= entry;
			Float exit = distanceToMediumExit(ray);

			if (exit != std::numeric_limits<Float>::infinity()) {
				Float t = invertDensityIntegral(ray, std::min(remaining, exit), 
						accumulatedTau, currentSigmaT, currentAlbedo, desiredTau);
				if (t != std::numeric_limits<Float>::infinity()) {
					mRec.p = ray(t);
					mRec.t = (mRec.p-r.o).length();
					success = true;
					break;
				}
			} else {
				break;
			}

			remaining -= exit;

			if (remaining < 0) {
				break;
			}

			ray.o = ray(exit);
			entry = distanceToMediumEntry(ray);
			if (++iterations > 10) {
				/// Just a precaution..
				Log(EWarn, "sampleDistance(): round-off error issues?");
				break;
			}
		}
		accumulatedTau *= m_sizeMultiplier;
		currentSigmaT *= m_sizeMultiplier;

		if (!success) {
			/* Could not achieve the desired density - hit the next surface instead */
			mRec.pdf = std::max((Float) Epsilon, std::exp(-accumulatedTau));
			mRec.attenuation = (-Spectrum(accumulatedTau)).exp();
			return false;
		}

		mRec.pdf = std::exp(-accumulatedTau) * std::max((Float) Epsilon, currentSigmaT);
		mRec.attenuation = (-Spectrum(accumulatedTau)).exp();

		mRec.sigmaS = currentAlbedo * currentSigmaT;
		mRec.sigmaA = Spectrum(currentSigmaT) - mRec.sigmaS;
		mRec.albedo = currentAlbedo.max();
		mRec.orientation = m_orientations->lookupVector(mRec.p);
		mRec.medium = this;
		return true;
	}


	inline Float sigmaT(const Point &p, const Vector &d) const {
		Float density = m_densities->lookupFloat(p);
		if (density == 0)
			return 0.0f;
		Vector orientation = m_orientations->lookupVector(p);
		Float length = orientation.length();
		if (length == 0)
			return 0.0f;
		else
			orientation /= length;
		Vector localD = normalize(Frame(orientation).toLocal(d));

		return density * m_sigmaT.evalAzimuthallyInvariant(localD);
	}

	inline Float lookupDensity(const Point &p, const Vector &d) const {
		return m_densities->lookupFloat(p);
	}

	inline Vector lookupOrientation(const Point &p) const {
		return m_orientations->lookupVector(p);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HeterogeneousFlakeMedium[" << endl
			<< "  albedo=" << indent(m_albedo.toString()) << endl
			<< "  orientations=" << indent(m_orientations.toString()) << endl
			<< "  densities=" << indent(m_densities.toString()) << endl
			<< "  shape =" << indent(m_shape.toString()) << endl
			<< "]";
		return oss.str();
	}
	MTS_DECLARE_CLASS()
private:
	ref<VolumeDataSource> m_densities;
	ref<VolumeDataSource> m_albedo;
	ref<VolumeDataSource> m_orientations;
	ref<Shape> m_shape;
	ref<KDTree> m_kdTree;
	Float m_stepSize;
	/* Flake model-related */
	SHVector m_D, m_sigmaT;
	SHVector4D m_phaseExpansion;
	int m_samplingRecursions;
	bool m_fiber;
	Float m_exponent, m_normalization;
};

Spectrum FlakePhaseFunction::f(const MediumSamplingRecord &mRec, const Vector &_wi, const Vector &_wo) const {
	Vector orientation = m_medium->lookupOrientation(mRec.p);
	if (orientation.lengthSquared() == 0)
		return Spectrum(0.0f);

	Frame frame(orientation);
	Vector wi = frame.toLocal(_wi);
	Vector wo = frame.toLocal(_wo);
#if 0 
	/* Evaluate spherical harmonics representation - might be inaccurate */
	SHVector temp(m_phaseExpansion->getBands());
	m_phaseExpansion->lookup(wi, temp);
	Float result = std::max((Float) 0, temp.eval(wo));
	return Spectrum(result);
#else
	/* Evaluate the real phase function */
	FlakeDistr D(m_exponent, m_fiber);
	Float sigmaT = m_sigmaT.evalAzimuthallyInvariant(wi);
	Vector H = wi + wo;
	Float length = H.length();
	if (length == 0)
		return Spectrum(0.0f);
	else
		H /= length;
	return (D(H) + D(-H)) * m_normalization / (2*sigmaT);
#endif
}

Spectrum FlakePhaseFunction::sample(const MediumSamplingRecord &mRec, const Vector &_wi, Vector &_wo, 
		ESampledType &sampledType, Float &pdf, const Point2 &_sample) const {
	sampledType = PhaseFunction::ENormal;
#if 1
	/* Importance sampling using the interpolated phase function */
	Frame frame(m_medium->lookupOrientation(mRec.p));
	Point2 sample(_sample);
	Vector wi = frame.toLocal(_wi);

	SHVector &temp = m_temp.get();
	if (EXPECT_NOT_TAKEN(temp.getBands() == 0))
		temp = SHVector(m_phaseExpansion->getBands());
	m_phaseExpansion->lookup(wi, temp);
	pdf = m_shSampler->warp(temp, sample);

	Vector wo = sphericalDirection(sample.x, sample.y);
	_wo = frame.toWorld(wo);
	return Spectrum(std::max((Float) 0, temp.eval(sample.x, sample.y)));
//	return f(mRec, _wi, _wo);
#else
	/* Uniform sampling */
	_wo = squareToSphere(_sample);
	pdf = 1/(4*M_PI);
	return f(mRec, _wi, _wo);
#endif
}

Float FlakePhaseFunction::pdf(const MediumSamplingRecord &mRec, const Vector &_wi, const Vector &_wo) const {
	Frame frame(m_medium->lookupOrientation(mRec.p));
	SHVector temp(m_phaseExpansion->getBands());
	Vector wi = frame.toLocal(_wi);
	Vector wo = frame.toLocal(_wo);
	m_phaseExpansion->lookup(wi, temp);

	return std::max((Float) 0, temp.eval(wo));
}

MTS_IMPLEMENT_CLASS_S(HeterogeneousFlakeMedium, false, Medium)
MTS_IMPLEMENT_CLASS_S(FlakePhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(HeterogeneousFlakeMedium, "Stenciled heterogeneous medium");
MTS_NAMESPACE_END


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
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/brent.h>
#include <mitsuba/core/shvector4d.h>
#include <mitsuba/core/statistics.h>
#include <fstream>

/*
 * When inverting an integral of the form f(t):=\int_0^t [...] dt
 * using composite Simpson's rule quadrature, the following constant
 * specified how many times the step size will be reduced when we
 * have almost reached the desired value.
 */
#define NUM_REFINES 8 

/**
 * Allow to stop integrating densities when the resulting segment
 * has a throughput of less than 'Epsilon'
 */
#define HET_EARLY_EXIT 1

/// Header identification record for phase function caches
#define PHASE_FUNCTION_HEADER 0x40044004

#define USE_REJECTION_SAMPLING 1
#define USE_COSTHETA_DISTR 1

MTS_NAMESPACE_BEGIN

static StatsCounter totalSamples("Micro-flake model", "Sample generations");
static StatsCounter avgSampleIterations("Micro-flake model", "Average rejection sampling iterations", EAverage);
static StatsCounter brentSolves("Micro-flake model", "Brent solver calls");
static StatsCounter avgBrentFunEvals("Micro-flake model", "Average Brent solver function evaluations", EAverage);

struct AbsCos {
	Float operator()(const Vector &w) const { return std::abs(w.z); }
};

#define USE_COSTHETA_DISTR 1

enum EMode {
	ESinP = 0,
	ECosP
};

struct FlakeDistr {
	Float exp;
	EMode mode;

	FlakeDistr() { }
	FlakeDistr(EMode mode, Float exp) : exp(exp), mode(mode) { }

	inline Float operator()(const Vector &w) const {
		switch (mode) {
			case ESinP:
				return std::pow(1-std::pow(w.z, 2), exp/2.0f);
			case ECosP:
				return std::pow(std::abs(w.z), exp);
			default:
				SLog(EError, "Invalid flake distribution!");
				return 0.0f;
		}
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
		Float exponent, Float normalization, EMode mode, HeterogeneousFlakeMedium *medium) : PhaseFunction(Properties()), 
		m_sigmaT(sigmaT), m_samplingRecursions(samplingRecursions), m_distr(mode, exponent), m_phaseExpansion(phaseExpansion), 
		m_exponent(exponent), m_normalization(normalization), m_mode(mode), m_medium(medium) {
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

	/* Implementations are below */
	Spectrum f(const PhaseFunctionQueryRecord &pRec) const;

	Spectrum sample(PhaseFunctionQueryRecord &pRec,
		Float &pdf, Sampler *sampler) const;
	
	Spectrum sample(PhaseFunctionQueryRecord &pRec,
		Sampler *sampler) const;

	Float pdf(const PhaseFunctionQueryRecord &pRec) const;

	std::string toString() const {
		return "FlakePhaseFunction[]";
	}

	MTS_DECLARE_CLASS()
private:
	SHVector m_sigmaT;
	ref<SHSampler> m_shSampler;
	int m_samplingRecursions;
	FlakeDistr m_distr;
	const SHVector4D *m_phaseExpansion;
	mutable PrimitiveThreadLocal<SHVector> m_temp;
	Float m_exponent;
	Float m_normalization;
	EMode m_mode;
	const HeterogeneousFlakeMedium *m_medium; 
};

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
		std::string mode = props.getString("mode", "sinp");

		if (mode == "sinp") {
			m_mode = ESinP;
		} else if (mode == "cosp") {
			m_mode = ECosP;
		} else {
			Log(EError, "Unknown mode: %s!", mode.c_str());
		}	

		const int numSamples = 80;
		ref<FileStream> stream;

		SHVector absCos = SHVector(bands);
		absCos.project(AbsCos(), numSamples);
		Float absCosError = absCos.l2Error(AbsCos(), numSamples);
		
		m_exponent = props.getFloat("exponent", 4);
		FlakeDistr distrFunc(m_mode, m_exponent);

		m_D = SHVector(bands);
		m_D.project(distrFunc, numSamples);
		m_normalization = 1 / (m_D(0, 0) * 2 * std::sqrt(M_PI));
		Float dError = m_D.l2Error(distrFunc, numSamples);
		m_D.normalize();

		Log(EInfo, "=======  Flake medium properties =======");
		Log(EInfo, "       distribution type = %s", mode.c_str());
		Log(EInfo, "   distribution exponent = %f", m_exponent);
		Log(EInfo, "    SH bands (D, sigmaT) = %i", bands);
		Log(EInfo, "    SH bands (phase fct) = %i", phaseBands);
		Log(EInfo, "  |cos| projection error = %f", absCosError);
		Log(EInfo, "     D  projection error = %f", dError);
		Log(EInfo, "     sampling recursions = %i", m_samplingRecursions);

		bool computePhaseProjection = true;
		if (fs::exists("flake-phase.dat")) {
			/* Avoid recomputing this every time */
			stream = new FileStream("flake-phase.dat", FileStream::EReadOnly);
			unsigned int header = stream->readUInt();
			int nBands = stream->readInt();
			int mode = stream->readInt();
			Float exponent = stream->readFloat();
			if (header == PHASE_FUNCTION_HEADER && nBands == phaseBands 
				&& exponent == m_exponent && mode == m_mode) {
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
			stream->writeInt(m_mode);
			stream->writeFloat(m_exponent);
			m_phaseExpansion.serialize(stream);
			stream->close();
		}

		for (int i=0; i<phaseBands; ++i)
			Log(EInfo, "Energy in phase function band %i: %f", i, m_phaseExpansion.energy(i));

		/* Compute sigma_t except for the a*rho factor */
		m_sigmaT = SHVector(m_D);
		m_sigmaT.convolve(absCos);
		Assert(m_sigmaT.isAzimuthallyInvariant());
		m_kdTree = new ShapeKDTree();
	}
	
	/* Unserialize from a binary data stream */
	HeterogeneousFlakeMedium(Stream *stream, InstanceManager *manager) 
		: Medium(stream, manager) {
		m_kdTree = new ShapeKDTree();
		size_t shapeCount = stream->readUInt();
		for (size_t i=0; i<shapeCount; ++i) 
			addChild("", static_cast<Shape *>(manager->getInstance(stream)));
		m_densities = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_albedo = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_orientations = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_stepSize = stream->readFloat();
		m_exponent = stream->readFloat();
		m_normalization = stream->readFloat();
		m_mode = (EMode) stream->readInt();
		m_D = SHVector(stream);
		m_sigmaT = SHVector(stream);
		m_phaseExpansion = SHVector4D(stream);
		m_samplingRecursions = stream->readInt();
		configure();
	}

	virtual ~HeterogeneousFlakeMedium() {
		for (size_t i=0; i<m_shapes.size(); ++i)
			m_shapes[i]->decRef();
	}

	/* Serialize the volume to a binary data stream */
	void serialize(Stream *stream, InstanceManager *manager) const {
		Medium::serialize(stream, manager);
		stream->writeUInt((uint32_t) m_shapes.size());
		for (size_t i=0; i<m_shapes.size(); ++i)
			manager->serialize(stream, m_shapes[i]);
		manager->serialize(stream, m_densities.get());
		manager->serialize(stream, m_albedo.get());
		manager->serialize(stream, m_orientations.get());
		stream->writeFloat(m_stepSize);
		stream->writeFloat(m_exponent);
		stream->writeFloat(m_normalization);
		stream->writeInt(m_mode);
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
			&m_phaseExpansion, m_exponent, m_normalization, m_mode, this);

		if (m_shapes.size() > 0) {
			m_kdTree->build();
			m_aabb = m_kdTree->getAABB();
		} else {
			m_aabb = m_densities->getAABB();
		}
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

	/**
	 * \brief Intersect a ray with the stencil and 
	 * determine whether it started on the inside
	 */
	inline bool rayIntersect(const Ray &ray, Float &t, bool &inside) const {
		if (m_shapes.size() != 0) {
			Intersection its;
			if (!m_kdTree->rayIntersect(ray, its)) 
				return false;
			inside = dot(its.geoFrame.n, ray.d) > 0;
			t = its.t;
		} else {
			Float mint, maxt;
			if (!m_aabb.rayIntersect(ray, mint, maxt))
				return false;
			if (mint <= ray.mint && maxt <= ray.mint)
				return false;
			inside = mint <= 0 && maxt > 0;
			t = (mint <= 0) ? maxt : mint;
		}
		return true;
	}

	inline Vector lookupOrientation(const Point &p) const {
		Vector orientation = m_orientations->lookupVector(p);
		Float lengthSqr = orientation.lengthSquared();
		if (lengthSqr != 0)
			return orientation / std::sqrt(lengthSqr);
		else
			return Vector(0.0f);
	}

	inline Float sigmaT(const Point &p, const Vector &d) const {
		/**
		 * Double the densities so that a uniform flake distribution
		 * reproduces isotropic scattering
		 */
		Float density = 2 * m_densities->lookupFloat(p) * m_densityMultiplier;
		if (density == 0)
			return 0.0f;
		Vector orientation = lookupOrientation(p);
		if (orientation == Vector(0.0f))
			return 0.0f;

		Vector localD = Frame(orientation).toLocal(d);

		return density * m_sigmaT.evalAzimuthallyInvariant(localD);
	}

	/// Integrate densities using the composite Simpson's rule
	inline Float integrateDensity(Ray ray, Float length) const {
		if (length == 0)
			return 0.0f;
		int nParts = (int) std::ceil(length/m_stepSize);
		nParts += nParts % 2;
		const Float stepSize = length/nParts;
		const Vector increment = ray.d * stepSize;

		Float m = 4;
		Float accumulatedDensity = 
			sigmaT(ray.o, ray.d) + sigmaT(ray(length), ray.d);
#ifdef HET_EARLY_EXIT
		const Float stopAfterDensity = -std::log(Epsilon);
		const Float stopValue = stopAfterDensity*3/stepSize;
#endif

		ray.o += increment;
		for (int i=1; i<nParts; ++i) {
			Float value = sigmaT(ray.o, ray.d);
			accumulatedDensity += value*m;
			Point tmp = ray.o;
			ray.o += increment;

			if (ray.o == tmp) {
				Log(EWarn, "integrateDensity(): failed to make forward progress -- "
						"round-off error issues? The step size was %.18f", stepSize);
				break;
			}

#ifdef HET_EARLY_EXIT
			if (accumulatedDensity > stopValue) // Stop early
				return std::numeric_limits<Float>::infinity();
#endif
			m = 6 - m;
		}
		accumulatedDensity *= stepSize/3;

		return accumulatedDensity;
	}

	/**
	 * \brief Integrate the density covered by the specified ray
	 * segment, while clipping the volume to the underlying
	 * stencil volume.
	 */
	inline Float integrateDensity(const Ray &r) const {
		Ray ray(r(r.mint), r.d,
			0, std::numeric_limits<Float>::infinity(), 0.0f);
		Float integral = 0, remaining = r.maxt - r.mint;
		int iterations = 0;
		bool inside = false;
		Float t = 0;

		while (remaining > 0 && rayIntersect(ray, t, inside)) {
			if (inside) 
				integral += integrateDensity(ray, std::min(remaining, t));

			remaining -= t;
			ray.o = ray(t);
			ray.mint = Epsilon;

			if (++iterations > 100) {
				/// Just a precaution..
				Log(EWarn, "integrateDensity(): round-off error issues?");
				break;
			}
		}
		return integral;
	}
	
	/**
	 * Attempts to solve the following 1D integral equation for 't'
	 * \int_0^t density(ray.o + x * ray.d) * dx + accumulatedDensity == desiredDensity.
	 * When no solution can be found in [0, maxDist] the function returns
	 * false. For convenience, the function returns the current values of sigmaT 
	 * and the albedo, as well as the 3D position 'ray.o+t*ray.d' upon
	 * success.
	 */
	bool invertDensityIntegral(Ray ray, Float maxDist, 
			Float &accumulatedDensity, Float desiredDensity,
			Float &currentSigmaT, Spectrum &currentAlbedo, 
			Point &currentPoint) const {
		if (maxDist == 0)
			return std::numeric_limits<Float>::infinity();
		int nParts = (int) std::ceil(maxDist/m_stepSize);
		Float stepSize = maxDist/nParts;
		Vector fullIncrement = ray.d * stepSize,
			   halfIncrement = fullIncrement * .5f;
		Float t = 0;
		int numRefines = 0;
		bool success = false;
		Float node1 = sigmaT(ray.o, ray.d);

		while (t <= maxDist) {
			Float node2 = sigmaT(ray.o + halfIncrement, ray.d);
			Float node3 = sigmaT(ray.o + fullIncrement, ray.d);
			const Float newDensity = accumulatedDensity + 
				(node1+node2*4+node3) * (stepSize/6);

			if (newDensity > desiredDensity) {
				if (t+stepSize/2 <= maxDist) {
					/* Record the last "good" scattering event */
					success = true;
					if (EXPECT_TAKEN(node2 != 0)) {
						currentPoint = ray.o + halfIncrement;
						currentSigmaT = node2;
					} else if (node3 != 0 && t+stepSize <= maxDist) {
						currentPoint = ray.o + fullIncrement;
						currentSigmaT = node3;
					} else if (node1 != 0) {
						currentPoint = ray.o;
						currentSigmaT = node1;
					}
				}

				if (++numRefines > NUM_REFINES)
					break;

				stepSize *= .5f;
				fullIncrement = halfIncrement;
				halfIncrement *= .5f;
				continue;
			}

			if (ray.o+fullIncrement == ray.o) /* Unable to make forward progress - stop */
				break;

			accumulatedDensity = newDensity;
			node1 = node3;
			t += stepSize;
			ray.o += fullIncrement;
		}

		if (success) 
			currentAlbedo = m_albedo->lookupSpectrum(currentPoint);

		return success;
	}

	/// Evaluate the attenuation along a ray segment
	Spectrum tau(const Ray &ray) const {
		return Spectrum(std::exp(-integrateDensity(ray)));
	}

	bool sampleDistance(const Ray &r, MediumSamplingRecord &mRec,
			Sampler *sampler) const {
		Ray ray(r(r.mint), r.d, 0.0f);
		Float distSurf       = r.maxt - r.mint,
			  desiredDensity     = -std::log(1-sampler->next1D()),
			  accumulatedDensity = 0.0f,
			  currentSigmaT  = 0.0f;
		Spectrum currentAlbedo(0.0f);
		int iterations = 0;
		bool inside = false, success = false;
		Float t = 0;
		Float sigmaTOrigin = sigmaT(ray.o, ray.d);

		while (rayIntersect(ray, t, inside)) {
			if (inside) {
				Point currentPoint(0.0f);
				success = invertDensityIntegral(ray, std::min(t, distSurf), 
						accumulatedDensity, desiredDensity, currentSigmaT, currentAlbedo,
						currentPoint);
				if (success) {
					/* A medium interaction occurred */
					mRec.p = currentPoint;
					mRec.t = (mRec.p-r.o).length();
					success = true;
					break;
				}
			}

			distSurf -= t;

			if (distSurf < 0) {
				/* A surface interaction occurred */
				break;
			}
			ray.o = ray(t);
			ray.mint = Epsilon;

			if (++iterations > 100) {
				/// Just a precaution..
				Log(EWarn, "sampleDistance(): round-off error issues?");
				break;
			}
		}
		Float expVal = std::max(Epsilon, std::exp(-accumulatedDensity));

		if (success)
			mRec.orientation = lookupOrientation(mRec.p);

		mRec.pdfFailure = expVal;
		mRec.pdfSuccess = expVal * std::max((Float) Epsilon, currentSigmaT);
		mRec.pdfSuccessRev = expVal * std::max((Float) Epsilon, sigmaTOrigin);
		mRec.attenuation = Spectrum(expVal);

		if (!success)
			return false;

		mRec.sigmaS = currentAlbedo * currentSigmaT;
		mRec.sigmaA = Spectrum(currentSigmaT) - mRec.sigmaS;
		mRec.albedo = currentAlbedo.max();
		mRec.medium = this;
		return true;
	}

	void pdfDistance(const Ray &ray, Float t, MediumSamplingRecord &mRec) const {
		Float expVal = std::exp(-integrateDensity(Ray(ray, 0, t)));

		mRec.attenuation = Spectrum(expVal);
		mRec.pdfFailure = expVal;
		mRec.pdfSuccess = expVal * std::max((Float) Epsilon, 
			sigmaT(ray(t), ray.d));
		mRec.pdfSuccessRev = expVal * std::max((Float) Epsilon, 
			sigmaT(ray(ray.mint), ray.d));
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HeterogeneousFlakeMedium[" << endl
			<< "  albedo=" << indent(m_albedo.toString()) << "," << endl
			<< "  orientations=" << indent(m_orientations.toString()) << "," << endl
			<< "  densities=" << indent(m_densities.toString()) << "," << endl
			<< "  stepSize=" << m_stepSize << "," << endl
			<< "  densityMultiplier=" << m_densityMultiplier << endl
			<< "]";
		return oss.str();
	}
	MTS_DECLARE_CLASS()
private:
	ref<VolumeDataSource> m_densities;
	ref<VolumeDataSource> m_albedo;
	ref<VolumeDataSource> m_orientations;
	ref<ShapeKDTree> m_kdTree;
	std::vector<Shape *> m_shapes;
	Float m_stepSize;
	/* Flake model-related */
	SHVector m_D, m_sigmaT;
	SHVector4D m_phaseExpansion;
	int m_samplingRecursions;
	EMode m_mode;
	Float m_exponent, m_normalization;
};

Spectrum FlakePhaseFunction::f(const PhaseFunctionQueryRecord &pRec) const {
	Frame frame(pRec.mRec.orientation);
	Vector wi = frame.toLocal(pRec.wi);
	Vector wo = frame.toLocal(pRec.wo);
#if defined(FLAKE_EVAL_APPROXIMATION) 
	/* Evaluate spherical harmonics representation - might be inaccurate */
	SHVector temp(m_phaseExpansion->getBands());
	m_phaseExpansion->lookup(wi, temp);
	Float result = std::max((Float) 0, temp.eval(wo));
	return Spectrum(result);
#else
	/* Evaluate the real phase function */
	Float sigmaT = m_sigmaT.evalAzimuthallyInvariant(wi);
	Vector H = wi + wo;
	Float length = H.length();

	if (length == 0)
		return Spectrum(0.0f);
	else
		H /= length;

	Float result = m_distr(H) / (2*sigmaT);
	result *= m_normalization;
	return Spectrum(result);
#endif
}

Spectrum FlakePhaseFunction::sample(PhaseFunctionQueryRecord &pRec,
		Float &_pdf, Sampler *sampler) const {
#if defined(FLAKE_SAMPLE_UNIFORM)
	/* Uniform sampling */
	pRec.wo = squareToSphere(sampler->next2D());
	_pdf = 1/(4*M_PI);
	return f(pRec);
#else
	Point2 sample(sampler->next2D());
	Frame frame(pRec.mRec.orientation);
	Vector wi = frame.toLocal(pRec.wi);
	/* Importance sampling using the interpolated phase function */
	SHVector &temp = m_temp.get();
	if (EXPECT_NOT_TAKEN(temp.getBands() == 0))
		temp = SHVector(m_phaseExpansion->getBands());
	m_phaseExpansion->lookup(wi, temp);
	_pdf = m_shSampler->warp(temp, sample);

	Vector wo = sphericalDirection(sample.x, sample.y);
	pRec.wo = frame.toWorld(wo);
	#if defined(FLAKE_EVAL_APPROXIMATION) 
		return Spectrum(std::max((Float) 0, temp.eval(sample.x, sample.y)));
	#else
		return f(pRec);
	#endif
#endif
}

Spectrum FlakePhaseFunction::sample(PhaseFunctionQueryRecord &pRec,
		Sampler *sampler) const {
	Float pdf;
	Spectrum result = FlakePhaseFunction::sample(pRec, pdf, sampler);
	if (!result.isZero())
		return result / pdf;
	else
		return Spectrum(0.0f);
}

Float FlakePhaseFunction::pdf(const PhaseFunctionQueryRecord &pRec) const {
#if defined(FLAKE_SAMPLE_UNIFORM)
	return 1/(4*M_PI);
#else
	Frame frame(pRec.mRec.orientation);
	Vector wi = frame.toLocal(pRec.wi);
	Vector wo = frame.toLocal(pRec.wo);
	SHVector temp(m_phaseExpansion->getBands());
	m_phaseExpansion->lookup(wi, temp);

	return std::max((Float) 0, temp.eval(wo));
#endif
}

MTS_IMPLEMENT_CLASS_S(HeterogeneousFlakeMedium, false, Medium)
MTS_IMPLEMENT_CLASS_S(FlakePhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(HeterogeneousFlakeMedium, "Heterogeneous micro-flake medium");
MTS_NAMESPACE_END


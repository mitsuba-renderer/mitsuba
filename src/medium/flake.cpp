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
#include <mitsuba/core/shvector4d.h>
#include <mitsuba/core/fstream.h>
#include <boost/numeric/ublas/symmetric.hpp>
#ifdef MTS_HAVE_LAPACK
#include <boost/numeric/bindings/lapack/syev.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#endif
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iomanip>

MTS_NAMESPACE_BEGIN

#ifdef MTS_HAVE_LAPACK
namespace lapack = boost::numeric::bindings::lapack;
#endif

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

/* Generate MATLAB-style printouts */
std::string mtxToString(const ublas::matrix<Float> &mtx) {
	std::ostringstream oss;
	oss << std::setprecision(5);
	oss << "[";
	for (size_t i=0; i<mtx.size1(); ++i) {
		for (size_t j=0; j<mtx.size2(); ++j) {
			oss << mtx(i, j);
			if (j+1<mtx.size1())
				oss << " ";
		}
		if (i+1<mtx.size2())
			oss << ";";
	}
	oss << "]";
	return oss.str();
}

std::string vecToString(const ublas::vector<Float> &vec) {
	std::ostringstream oss;
	oss << std::setprecision(5);
	oss << "[";
	for (size_t i=0; i<vec.size(); ++i) {
		oss << vec(i);
		if (i+1<vec.size())
			oss << " ";
	}
	oss << "]";
	return oss.str();
}

/**
 * Phase function for use in a flake medium
 */
class FlakePhaseFunction : public PhaseFunction {
public:
	FlakePhaseFunction(const SHVector &distr, const SHVector &sigmaS,
		const SHVector4D &phaseExpansion, Float area, Float rho, const Frame &frame,
		int samplingRecursions) : PhaseFunction(Properties()), 
		m_distr(distr), m_sigmaS(sigmaS), m_phaseExpansion(phaseExpansion), m_area(area), 
		m_rho(rho), m_frame(frame), m_samplingRecursions(samplingRecursions) {
		initialize();
	}

	FlakePhaseFunction(Stream *stream, InstanceManager *manager) :
			PhaseFunction(stream, manager) {
		m_distr  = SHVector(stream);
		m_sigmaS = SHVector(stream);
		m_phaseExpansion = SHVector4D(stream);
		m_area = stream->readFloat();
		m_rho = stream->readFloat();
		m_frame = Frame(stream);
		m_samplingRecursions = stream->readInt();
		initialize();
	}

	void initialize() {
		m_shSampler = new SHSampler(m_phaseExpansion.getBands(), m_samplingRecursions);
		Log(EInfo, "Constructing a SH sampler: %s", m_shSampler->toString().c_str());
	}

	virtual ~FlakePhaseFunction() { }

	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);

		m_distr.serialize(stream);
		m_sigmaS.serialize(stream);
		m_phaseExpansion.serialize(stream);
		stream->writeFloat(m_area);
		stream->writeFloat(m_rho);
		m_frame.serialize(stream);
		stream->writeInt(m_samplingRecursions);
	}

	Spectrum f(const MediumSamplingRecord &mRec, const Vector &_wi, const Vector &_wo) const {
		/*
		Vector h = normalize(wi+wo);
		return m_area*m_rho*(
			std::max((Float) 0, m_distr.eval(h)) + 
			std::max((Float) 0, m_distr.eval(-h)))
			/ (8*m_sigmaS.eval(wi));
		*/

		Vector wi = m_frame.toLocal(_wi);
		Vector wo = m_frame.toLocal(_wo);
		SHVector temp(m_phaseExpansion.getBands());
		m_phaseExpansion.lookup(wi, temp);

		return Spectrum(std::max((Float) 0, temp.eval(wo)));
	}

	Spectrum sample(const MediumSamplingRecord &mRec, const Vector &wi, Vector &wo, 
		ESampledType &sampledType, const Point2 &_sample) const {
		Float pdf;
		Spectrum value = sample(mRec, wi, wo, sampledType, pdf, _sample);
		return value/pdf;
	}

	Spectrum sample(const MediumSamplingRecord &mRec, const Vector &_wi, Vector &_wo, 
			ESampledType &sampledType, Float &pdf, const Point2 &_sample) const {
		sampledType = PhaseFunction::ENormal;
#if 1
		Point2 sample(_sample);
		Vector wi = m_frame.toLocal(_wi);
		SHVector temp(m_phaseExpansion.getBands());
		m_phaseExpansion.lookup(wi, temp);

		pdf = m_shSampler->warp(temp, sample);
		Vector wo = sphericalDirection(sample.x, sample.y);
		_wo = m_frame.toWorld(wo);
		return Spectrum(std::max((Float) 0, temp.eval(sample.x, sample.y)));
#else
		wo = squareToSphere(_sample);
		pdf = 1/(4*M_PI);
		return f(wi, wo);
#endif
	}

	Float pdf(const MediumSamplingRecord &mRec, const Vector &wi, const Vector &wo) const {
		Log(EError, "Not implemented!");
		return 0.0f;
	}

	std::string toString() const {
		return "FlakePhaseFunction[]";
	}

	MTS_DECLARE_CLASS()
private:
	SHVector m_distr, m_sigmaS;
	ref<SHSampler> m_shSampler;
	SHVector4D m_phaseExpansion;
	Float m_area, m_rho;
	Frame m_frame;
	int m_samplingRecursions;
};

/**
 * Homogeneous participating medium based on a flake distribution. 
 * An arbitrary (manifold) shape must be specified as a child object.
 */
class FlakeMedium : public Medium {
public:
	FlakeMedium(const Properties &props) 
		: Medium(props) {
		int bands = props.getInteger("bands", 8);
		Assert(bands>0);

		m_albedo = props.getFloat("albedo", 0.95);
		int samplingRecursions = props.getInteger("samplingRecursions", 10);
		m_rho = props.getFloat("rho", 2);
		m_area = props.getFloat("area", 1.3);
		bool fiber = props.getBoolean("fiber", true);
		const int numSamples = 80;
		ref<FileStream> stream;

		SHVector absCos = SHVector(bands);
		absCos.project(AbsCos(), numSamples);
		Float absCosError = absCos.l2Error(AbsCos(), numSamples);

		Vector centralAxis = normalize(props.getVector("axis", Vector(1, 2, 3)));
		Float exponent = props.getFloat("exponent", 10);
		FlakeDistr distrFunc(exponent, fiber);

		m_D = SHVector(bands);
		m_D.project(distrFunc, numSamples);
		Float dError = m_D.l2Error(distrFunc, numSamples);
		m_D.normalize();

		SHVector4D phaseExpansion;
		if (FileStream::exists("flake-phase.dat")) {
			stream = new FileStream("flake-phase.dat", FileStream::EReadOnly);
			phaseExpansion = SHVector4D(stream);
			stream->close();
		} else {
			int res = 50;
			FlakePhaseFunctor<FlakeDistr> phaseFunctor(distrFunc);
			phaseExpansion = SHVector4D(res, 2*res, bands);
			Log(EDebug, "Phase function parameterization: %s", phaseExpansion.toString().c_str());
			phaseExpansion.project(phaseFunctor, numSamples);
			phaseExpansion.normalize();
			stream = new FileStream("flake-phase.dat", FileStream::ETruncReadWrite);
			phaseExpansion.serialize(stream);
			stream->close();
		}

		m_sigmaT = SHVector(m_D);
		m_sigmaT.convolve(absCos);
		m_sigmaT *= m_area*m_rho;
		Assert(m_sigmaT.isAzimuthallyInvariant());

		m_sigmaS = SHVector(m_sigmaT);
		m_sigmaS *= m_albedo;

		ublas::identity_matrix<Float> I(3);
		ublas::symmetric_matrix<Float> m_M = 9*m_area*m_rho/8.0f * (
			(1+3*m_albedo)*m_D.mu2() + 
			(1-m_albedo)*I
		);

		ublas::matrix<Float, ublas::column_major> Q(m_M), sqrtW(3,3), w(3,3);
		ublas::vector<Float> eigs(3);

#ifdef MTS_HAVE_LAPACK
		int result = lapack::syev('V', 'U', Q, eigs);
		if (result != 0)
			SLog(EError, "Unable to diagonalize the diffusion tensor!");

		sqrtW.clear(); w.clear();
		for (int i=0; i<3; ++i) {
			w(i, i) = eigs(i);
			sqrtW(i, i) = std::sqrt(eigs(i));
		}

		m_frame = Frame(centralAxis);
		ublas::matrix<Float, ublas::column_major> rot(3,3);
		rot(0, 0) = m_frame.s.x; rot(1, 0) = m_frame.s.y; rot(2, 0) = m_frame.s.z;
		rot(0, 1) = m_frame.t.x; rot(1, 1) = m_frame.t.y; rot(2, 1) = m_frame.t.z;
		rot(0, 2) = m_frame.n.x; rot(1, 2) = m_frame.n.y; rot(2, 2) = m_frame.n.z;
		Q = prod(rot, Q);

		m_M = ublas::prod(Q, ublas::matrix<Float>(ublas::prod(w, ublas::trans(Q))));
		ublas::symmetric_matrix<Float> m_P = 
			ublas::prod(Q, ublas::matrix<Float>(ublas::prod(sqrtW, ublas::trans(Q))));
#endif

		stream = new FileStream("flake-sigmaT.dat", FileStream::ETruncReadWrite);
		m_sigmaT.serialize(stream);
		stream->close();
		stream = new FileStream("flake-sigmaS.dat", FileStream::ETruncReadWrite);
		m_sigmaS.serialize(stream);
		stream->close();
		stream = new FileStream("flake-distr.dat", FileStream::ETruncReadWrite);
		m_D.serialize(stream);
		stream->close();

		m_phaseFunction = new FlakePhaseFunction(m_D, m_sigmaS, 
			phaseExpansion, m_area, m_rho, m_frame, samplingRecursions);

		Log(EInfo, "===  Flake medium - Material properties ===");
		Log(EInfo, "   particle surface area = %f", m_area);
		Log(EInfo, "        particle density = %f", m_rho);
		Log(EInfo, "                  albedo = %f", m_albedo);
		if (fiber)
			Log(EInfo, "                    D(w) = (1-dot(w, %s)^2)^%f", centralAxis.toString().c_str(), exponent/2);
		else
			Log(EInfo, "                    D(w) = abs(dot(w, %s))^%f", centralAxis.toString().c_str(), exponent);
		Log(EInfo, "");
		Log(EInfo, "      number of SH bands = %i", bands);
		Log(EInfo, "  |cos| projection error = %f", absCosError);
		Log(EInfo, "     D  projection error = %f", dError);
		Log(EInfo, "     sampling recursions = %i", samplingRecursions);
		Log(EInfo, "");
		Log(EInfo, "                coeff(D) = %s", m_D.toString().c_str());
		Log(EInfo, "                  mu2(D) = %s", mtxToString(m_D.mu2()).c_str());
		Log(EInfo, "           coeff(sigmaT) = %s", m_sigmaT.toString().c_str());
		Log(EInfo, "             mu2(sigmaT) = %s", mtxToString(m_sigmaT.mu2()).c_str());
		Log(EInfo, "             mu0(sigmaA) = %f", 2*m_area*m_rho*M_PI*(1-m_albedo));
		Log(EInfo, "               sigmaT(n) = %f", m_sigmaT.eval(m_frame.toLocal(Vector(0,0,1))));
#ifdef MTS_HAVE_LAPACK
		Log(EInfo, "                       M = %s", mtxToString(m_M).c_str());
		Log(EInfo, "                 eigs(M) = %s", vecToString(eigs).c_str());
		Log(EInfo, "                       P = %s", mtxToString(m_P).c_str());
#endif
	}

	FlakeMedium(Stream *stream, InstanceManager *manager)
		: Medium(stream, manager) {
		m_shape = static_cast<Shape *>(manager->getInstance(stream));
		m_kdTree = new KDTree();
		m_kdTree->addShape(m_shape.get());
		m_kdTree->build();
		
		m_D = SHVector(stream);
		m_sigmaS = SHVector(stream);
		m_sigmaT = SHVector(stream);
		m_area = stream->readFloat();
		m_rho = stream->readFloat();
		m_frame = Frame(stream);
	}

	virtual ~FlakeMedium() {
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Medium::serialize(stream, manager);
		manager->serialize(stream, m_shape.get());
		m_D.serialize(stream);
		m_sigmaS.serialize(stream);
		m_sigmaT.serialize(stream);
		stream->writeFloat(m_area);
		stream->writeFloat(m_rho);
		m_frame.serialize(stream);
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
		Ray ray(r(r.mint), r.d / dLength);
		Float coveredLength = 0, remaining = (r.maxt - r.mint) * dLength;
		bool inside = isInside(r);
		Intersection its;
		int iterations = 0;
		ray.mint = remaining * Epsilon;

		Float sigmaT = m_sigmaT.eval(m_frame.toLocal(ray.d));
		Assert(sigmaT > 0);

		while (remaining > 0 && m_kdTree->rayIntersect(ray, its)) {
			if (inside)
				coveredLength += std::min(remaining, its.t);
			remaining -= its.t;
			inside = !inside;
			ray.o = its.p;

			if (++iterations > 10) {
				/// Just a precaution..
				Log(EWarn, "tau(): round-off error issues?");
				break;
			}
		}

		return sigmaT * coveredLength;
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
		Float sigmaT = m_sigmaT.eval(m_frame.toLocal(theRay.d));
		Float sigmaS = m_albedo * sigmaT;
		Float sigmaA = (1-m_albedo) * sigmaT;

		Float desiredAttenuation = 1 - sampler->next1D();
		Float distMed = -std::log(desiredAttenuation) / sigmaT;
		mRec.pdf = sigmaT * std::exp(-sigmaT * distMed);
		mRec.miWeight = 1;

		Float traveled = 0,  // Traveled ray distance
			  covered = 0;   // Distance covered by the medium

		while (m_kdTree->rayIntersect(ray, its)) {
			if (inside) {
				/* Moving through the medium */
				if (its.t > distMed && distMed < distSurf) {
					/* A medium interaction occurred */
					mRec.p = ray(distMed);
					mRec.sigmaA = sigmaA;
					mRec.sigmaS = sigmaS;
					mRec.albedo = m_albedo;
					mRec.medium = this;
					mRec.attenuation = (Spectrum(sigmaT) * (-distMed)).exp();
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
		mRec.pdf = std::exp(-sigmaT * covered);
		mRec.t = traveled;
		mRec.attenuation = (Spectrum(sigmaT) * (-covered)).exp();
		mRec.miWeight = 1;

		return false;
	}

	void setParent(ConfigurableObject *parent) {
		if (parent->getClass()->derivesFrom(Shape::m_theClass))
			Log(EError, "Medium shape cannot be part of the scene");
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(Shape::m_theClass)) {
			Assert(m_shape == NULL);
			m_shape = static_cast<Shape *>(child);
			m_kdTree = new KDTree();
			m_kdTree->addShape(m_shape.get());
			m_kdTree->build();
			m_aabb = m_kdTree->getAABB();
		} else {
			Medium::addChild(name, child);
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "FlakeMedium[" << endl
			<< "  area = " << m_area << "," << std::endl
			<< "  rho = " << m_rho << "," << std::endl
			<< "  shape = " << indent(m_shape->toString()) << std::endl
			<< "]";
		return oss.str();
	}
	
	MTS_DECLARE_CLASS()
private:
	ref<Shape> m_shape;
	ref<KDTree> m_kdTree;
	SHVector m_D, m_sigmaS, m_sigmaT;
	Float m_area, m_rho;
	Frame m_frame;
};

MTS_IMPLEMENT_CLASS_S(FlakePhaseFunction, false, PhaseFunction)
MTS_IMPLEMENT_CLASS_S(FlakeMedium, false, Medium)
MTS_EXPORT_PLUGIN(FlakeMedium, "Flake medium");
MTS_NAMESPACE_END

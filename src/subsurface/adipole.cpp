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
#include <mitsuba/core/plugin.h>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/bindings/lapack/syev.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "irrtree.h"

namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

MTS_NAMESPACE_BEGIN

/**
 * Computes the combined diffuse radiant exitance 
 * caused by a number of dipole sources
 */
struct AnisotropicDipoleQuery {
	inline AnisotropicDipoleQuery(const ublas::symmetric_matrix<Float> *P,
		const Point *xr, const Point *xv, const Spectrum &detP, const Spectrum &beta,
			const Frame &frame, Float Fdt, const Point &p) 
			: P(P), xr(xr), xv(xv), detP(detP), beta(beta), frame(frame), Fdt(Fdt), p(p) {
		count = 0;
	}

	inline void operator()(const IrradianceSample &sample) {
		/* Project onto slab */
		Vector toP = frame.toLocal(sample.p - p);
		Float length = toP.length();
		toP.z = 0;
		toP = toP / toP.length() * length;
		Spectrum dMo;
		Float temp[3];

		for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
			Vector x = toP - xr[i], xp;
			for (int k=0; k<3; ++k)
				for (int l=0; l<3; ++l)
						xp[k] += P[i](k, l) * x[l];

			Float dr = xp.length(); 
			x = toP - xv[i]; xp = Vector();
			for (int k=0; k<3; ++k)
				for (int l=0; l<3; ++l)
						xp[k] += P[i](k, l) * x[l];

			Float dv = xp.length();

			Float zr = -xr[i].z, zv = xv[i].z;

			temp[i] = detP[i]/(4*M_PI) * (
				zr*(beta[i]*dr+1)*std::fastexp(-beta[i]*dr)/(dr*dr*dr) + 
				zv*(beta[i]*dv+1)*std::fastexp(-beta[i]*dv)/(dv*dv*dv)
			);
		}
		dMo.fromLinearRGB(temp[0], temp[1], temp[2]);

		result += dMo * sample.E * (sample.area * Fdt);
		count++;
	}

	inline const Spectrum &getResult() const {
		return result;
	}

	const ublas::symmetric_matrix<Float> *P;
	const Point *xr, *xv;
	const Spectrum &detP, &beta;
	Frame frame;
	Spectrum result;

	int count;
	Float Fdt;
	Point p;
};

static ref<Mutex> irrOctreeMutex = new Mutex();
static int irrOctreeIndex = 0;

static Float parseFloat(const std::string &name, const std::string &str, Float defVal = -1) {
	char *end_ptr = NULL;
	if (str == "") {
		if (defVal == -1)
			SLog(EError, "Missing floating point value (in <%s>)", name.c_str());
		return defVal;
	}
	Float result = (Float) std::strtod(str.c_str(), &end_ptr);
	if (*end_ptr != '\0')
		SLog(EError, "Invalid floating point value specified (in <%s>)", name.c_str());
	return result;
}

class AnisotropicDipole : public Subsurface {
public:
	AnisotropicDipole(const Properties &props) 
		: Subsurface(props) {
		irrOctreeMutex->lock();
		m_octreeIndex = irrOctreeIndex++;
		irrOctreeMutex->unlock();

		/* Multiplicative factor, which can be used to adjust the number of
		   irradiance samples */
		m_sampleMultiplier = props.getFloat("sampleMultiplier", 2.0f);
		/* Error threshold - lower means better quality */
		m_minDelta= props.getFloat("quality", 0.1f);
		/* Max. depth of the created octree */
		m_maxDepth = props.getInteger("maxDepth", 40);
		/* Multiplicative factor for the subsurface term - can be used to remove
		   this contribution completely, making it possible to use this integrator
		   for other interesting things.. */
		m_ssFactor = props.getSpectrum("ssFactor", Spectrum(1.0f));
		m_maxDepth = props.getInteger("maxDepth", 40);

		std::vector<std::string> tokens = tokenize(props.getString("D"), ",;[] ");
		if (tokens.size() != 9)
			Log(EError, "Invalid parameter 'mu2D'!");

		m_D = ublas::symmetric_matrix<Float>(3, 3);
		int index = 0;
		for (int i=0; i<3; ++i)
			for (int j=0; j<3; ++j)
				m_D(i, j) = parseFloat("D", tokens[index++]);

		m_sigmaTn = m_sigmaT * props.getFloat("sigmaTn");

		m_ready = false;
		m_octreeResID = -1;
	}
	
	AnisotropicDipole(Stream *stream, InstanceManager *manager) 
	 : Subsurface(stream, manager) {
		m_ssFactor = Spectrum(stream);
		m_sampleMultiplier = stream->readFloat();
		m_minDelta = stream->readFloat();
		m_maxDepth = stream->readInt();
		m_octreeIndex = stream->readInt();
		m_D = ublas::symmetric_matrix<Float>(3, 3);
		for (int i=0; i<3; ++i) 
			for (int j=0; j<3; ++j) 
				m_D(i, j) = stream->readFloat();
		m_sigmaTn = Spectrum(stream);
		m_ready = false;
		m_octreeResID = -1;
		configure();
	}

	virtual ~AnisotropicDipole() {
		if (m_octreeResID != -1)
			Scheduler::getInstance()->unregisterResource(m_octreeResID);
	}

	void bindUsedResources(ParallelProcess *proc) const {
		if (m_octreeResID != -1)
			proc->bindResource(formatString("irrOctree%i", m_octreeIndex), m_octreeResID);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Subsurface::serialize(stream, manager);
		m_ssFactor.serialize(stream);
		stream->writeFloat(m_sampleMultiplier);
		stream->writeFloat(m_minDelta);
		stream->writeInt(m_maxDepth);
		stream->writeInt(m_octreeIndex);

		for (int i=0; i<3; ++i) 
			for (int j=0; j<3; ++j) 
				stream->writeFloat(m_D(i, j));
		m_sigmaTn.serialize(stream);
	}

	inline ublas::vector<Float> cross(const ublas::vector<Float> &v1, const ublas::vector<Float> &v2) {
		/* Left-handed vector cross product */
		ublas::vector<Float> result(3);
		result(0) = v1(1)*v2(2) - v1(2)*v2(1);
		result(1) = v1(2)*v2(0) - v1(0)*v2(2);
		result(2) = v1(0)*v2(1) - v1(1)*v2(0);
		return result;
	}

	Spectrum Lo(const Intersection &its, const Vector &d) const {
		if (!m_ready || m_ssFactor.isZero())
			return Spectrum(0.0f);
		AnisotropicDipoleQuery query(m_P, m_xr, m_xv, m_detP, m_beta, its.shFrame, m_Fdt, its.p);
	
		const Normal &n = its.shFrame.n;
		m_octree->execute(query);

		if (m_eta == 1.0f) {
			return query.getResult() * m_ssFactor * INV_PI;
		} else {
			Float Ft = 1.0f - fresnel(absDot(n, d));
			return query.getResult() * m_ssFactor * INV_PI * (Ft / m_Fdr);
		}
	}

	Spectrum Li(const Ray &ray, const Normal &n) const {
		return Spectrum(0.0f);
	}

	void configure() {
		/* Average reflectance due to mismatched indices of refraction
		   at the boundary - [Groenhuis et al. 1983]*/
		m_Fdr = -1.440f / (m_eta * m_eta) + 0.710f / m_eta 
			+ 0.668f + 0.0636f * m_eta;

		/* Average transmittance at the boundary */
		m_Fdt = 1.0f - m_Fdr;

		if (m_eta == 1.0f) {
			m_Fdr = (Float) 0.0f;
			m_Fdt = (Float) 1.0f;
		}

		Spectrum albedo = m_sigmaS/m_sigmaT;
	

		/* Approximate dipole boundary condition term */
		m_A = (1 + m_Fdr) / m_Fdt;


		Vector centralAxis = normalize(Vector(1,1,0));
		for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
			ublas::identity_matrix<Float> I(3);
			ublas::symmetric_matrix<Float> M = 9*4*m_sigmaT[i]/16.0f * (
				(1+3*albedo[i])*m_D + 
				(1-albedo[i])*I
			);
			ublas::matrix<Float, ublas::column_major> Q(M), sqrtW(3,3), invSqrtW(3,3);
			ublas::vector<Float> eigs(3);

			int result = lapack::syev('V', 'U', Q, eigs);
			if (result != 0)
				SLog(EError, "Unable to diagonalize the diffusion tensor!");

			sqrtW.clear(); invSqrtW.clear();
			Float det = 1;
			for (int j=0; j<3; ++j) {
				sqrtW(j, j) = std::sqrt(eigs(j));
				invSqrtW(j, j) = 1/sqrtW(j,j);
				det *= sqrtW(j, j);
			}

			Frame frame(centralAxis);
			ublas::matrix<Float, ublas::column_major> rot(3,3);
			rot(0, 0) = frame.s.x; rot(1, 0) = frame.s.y; rot(2, 0) = frame.s.z;
			rot(0, 1) = frame.t.x; rot(1, 1) = frame.t.y; rot(2, 1) = frame.t.z;
			rot(0, 2) = frame.n.x; rot(1, 2) = frame.n.y; rot(2, 2) = frame.n.z;
			Q = prod(rot, Q);

			m_P[i] = ublas::prod(Q, ublas::matrix<Float>(ublas::prod(sqrtW, ublas::trans(Q))));
			ublas::matrix<Float> PInv = 
				ublas::prod(Q, ublas::matrix<Float>(ublas::prod(invSqrtW, ublas::trans(Q))));
			ublas::matrix<Float> MInv = prod(PInv, PInv);

			m_beta[i] = std::sqrt(m_sigmaA[i]);
				m_dp[i] = 2*m_A*MInv(2,2);
			ublas::vector<Float> n = prod(PInv, cross(column(m_P[i], 0), column(m_P[i], 1)));
			m_xr[i] = Point(0, 0, -1/m_sigmaTn[i]);
			m_xv[i] = m_xr[i] + Vector(n(0) / n(2), n(1) / n(2), 1) * 2 * (1/m_sigmaTn[i] + m_dp[i]);
			m_detP[i] = det;
		}
	}

	/// Unpolarized fresnel reflection term for dielectric materials
	Float fresnel(Float cosThetaI) const {
		Float g = std::sqrt(m_eta*m_eta - 1.0f + cosThetaI * cosThetaI);
		Float temp1 = (g - cosThetaI)/(g + cosThetaI);
		Float temp2 = (cosThetaI * (g + cosThetaI) - 1) / 
			(cosThetaI * (g - cosThetaI) + 1.0f);
		return 0.5f * temp1 * temp1 * (1.0f + temp2 * temp2);
	}

	void preprocess(const Scene *scene, int sceneResID, 
		int cameraResID, int samplerResID) {
		if (m_ready)
			return;

		if (!scene->getIntegrator()->getClass()
				->derivesFrom(MTS_CLASS(SampleIntegrator))) {
			Log(EError, "The dipole subsurface integrator requires "
				"a sampling-based surface integrator!");
		}

		m_octree = new IrradianceOctree(m_maxDepth, m_minDelta, 
			scene->getKDTree()->getAABB());

		Float sa = 0;
		for (std::vector<Shape *>::iterator it = m_shapes.begin(); 
			it != m_shapes.end(); ++it)
			sa += (*it)->getSurfaceArea();

		Float minMFP = (Spectrum(1)/m_sigmaTn).min();
		size_t sampleCount = (size_t) std::ceil(sa / (M_PI * minMFP * minMFP)
			* m_sampleMultiplier);

		ref<Scheduler> sched = Scheduler::getInstance();

		/* This could be a bit more elegant.. - inform the irradiance
		   sampler about the index of this subsurface integrator */
		std::vector<Subsurface *> ssIntegrators
			= scene->getSubsurfaceIntegrators();
		int index = -1;
		for (size_t i=0; i<ssIntegrators.size(); ++i) {
			if (ssIntegrators[i] == this) {
				index = i;
				break;
			}
		}
		Assert(index != -1);

		ref<IrradianceSamplingProcess> proc = new IrradianceSamplingProcess(
			sampleCount, (size_t) std::ceil(sampleCount/100.0f), index);

		proc->bindResource("scene", sceneResID);
		scene->bindUsedResources(proc);
		sched->schedule(proc);
		sched->wait(proc);

		const IrradianceRecordVector &results = *proc->getSamples();
		for (size_t i=0; i<results.size(); ++i) 
			m_octree->addSample(results[i]);

		m_octree->preprocess();
		m_octreeResID = Scheduler::getInstance()->registerResource(m_octree);

		m_ready = true;
	}

	void wakeup(std::map<std::string, SerializableObject *> &params) {
		std::string octreeName = formatString("irrOctree%i", m_octreeIndex);
		if (!m_octree.get() && params.find(octreeName) != params.end()) {
			m_octree = static_cast<IrradianceOctree *>(params[octreeName]);
			m_ready = true;
		}
	}

	MTS_DECLARE_CLASS()
private:
	Float m_sampleMultiplier;
	Float m_Fdr, m_Fdt, m_A, m_minDelta;
	Spectrum m_ssFactor;

	ref<IrradianceOctree> m_octree;
	int m_octreeResID, m_octreeIndex;
	int m_maxDepth;
	bool m_ready, m_requireSample;

	ublas::symmetric_matrix<Float> m_D, m_P[SPECTRUM_SAMPLES];
	Point m_xr[SPECTRUM_SAMPLES];
	Point m_xv[SPECTRUM_SAMPLES];
	Spectrum m_detP, m_dp, m_beta;
	Spectrum m_mu0SigmaA, m_sigmaTn;
};

MTS_IMPLEMENT_CLASS_S(AnisotropicDipole, false, Subsurface)
MTS_EXPORT_PLUGIN(AnisotropicDipole, "Anisotropic dipole model");
MTS_NAMESPACE_END

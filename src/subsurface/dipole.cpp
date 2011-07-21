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
#include "../medium/materials.h"
#include "irrtree.h"

MTS_NAMESPACE_BEGIN

/**
 * Computes the combined diffuse radiant exitance 
 * caused by a number of dipole sources
 */
struct IsotropicDipoleQuery {
#if !defined(MTS_SSE) || (SPECTRUM_SAMPLES != 3)
	inline IsotropicDipoleQuery(const Spectrum &zr, const Spectrum &zv, 
		const Spectrum &sigmaTr, Float Fdt, const Point &p) 
		: zr(zr), zv(zv), sigmaTr(sigmaTr), result(0.0f), Fdt(Fdt), p(p) {
			count = 0;
	}

	inline void operator()(const IrradianceSample &sample) {
		Spectrum rSqr = Spectrum((p - sample.p).lengthSquared());
		/* Distance to the real source */
		Spectrum dr = (rSqr + zr*zr).sqrt();
		/* Distance to the image point source */
		Spectrum dv = (rSqr + zv*zv).sqrt();
		Spectrum C1 = zr * (sigmaTr + Spectrum(1.0f) / dr);
		Spectrum C2 = zv * (sigmaTr + Spectrum(1.0f) / dv);

		/* Do not include the reduced albedo - will be canceled out later */
		Spectrum dMo = Spectrum(0.25f * INV_PI) *
			 (C1 * ((-sigmaTr * dr).exp()) / (dr * dr)
			+ C2 * ((-sigmaTr * dv).exp()) / (dv * dv));
		result += dMo * sample.E * (sample.area * Fdt);
		count++;
	}

	inline const Spectrum &getResult() const {
		return result;
	}

	Spectrum zr, zv, sigmaTr, result;
#else
	inline IsotropicDipoleQuery(const Spectrum &_zr, const Spectrum &_zv, 
		const Spectrum &_sigmaTr, Float Fdt, const Point &p) : Fdt(Fdt), p(p) {
		zr = _mm_set_ps(_zr[0], _zr[1], _zr[2], 0);
		zv = _mm_set_ps(_zv[0], _zv[1], _zv[2], 0);
		sigmaTr = _mm_set_ps(_sigmaTr[0], _sigmaTr[1], _sigmaTr[2], 0);
		zrSqr = _mm_mul_ps(zr, zr);
		zvSqr = _mm_mul_ps(zv, zv);
		result.ps = _mm_setzero_ps();
		count = 0;
	}

	inline void operator()(const IrradianceSample &sample) {
		/* Distance to the positive point source of the dipole */
		const __m128 lengthSquared = _mm_set1_ps((p - sample.p).lengthSquared()),
			drSqr = _mm_add_ps(zrSqr, lengthSquared), 
			dvSqr = _mm_add_ps(zvSqr, lengthSquared),
			dr = _mm_sqrt_ps(drSqr), dv = _mm_sqrt_ps(dvSqr), 
			one = _mm_set1_ps(1.0f),
			factor = _mm_mul_ps(_mm_set1_ps(0.25f*INV_PI*sample.area * Fdt),
				_mm_set_ps(sample.E[0], sample.E[1], sample.E[2], 0)),
			C1fac = _mm_div_ps(_mm_mul_ps(zr, _mm_add_ps(sigmaTr, _mm_div_ps(one, dr))), drSqr),
			C2fac = _mm_div_ps(_mm_mul_ps(zv, _mm_add_ps(sigmaTr, _mm_div_ps(one, dv))), dvSqr);
		SSEVector temp1(_mm_mul_ps(dr, sigmaTr)), temp2(_mm_mul_ps(dv, sigmaTr));
		const __m128
			exp1 = _mm_set_ps(expf(-temp1.f[3]), expf(-temp1.f[2]), expf(-temp1.f[1]), 0),
			exp2 = _mm_set_ps(expf(-temp2.f[3]), expf(-temp2.f[2]), expf(-temp2.f[1]), 0);
		result.ps = _mm_add_ps(result.ps, _mm_mul_ps(factor, _mm_add_ps(
			_mm_mul_ps(C1fac, exp1), _mm_mul_ps(C2fac, exp2))));
	}

	Spectrum getResult() {
		Spectrum value;
		for (int i=0; i<3; ++i)
			value[i] = result.f[3-i];
		return value;
	}

	__m128 zr, zv, zrSqr, zvSqr, sigmaTr;
	SSEVector result;
#endif

	int count;
	Float Fdt;
	Point p;
};

static ref<Mutex> irrOctreeMutex = new Mutex();
static int irrOctreeIndex = 0;

/**
 * Subsurface scattering integrator using Jensen's fast hierarchical 
 * dipole approximation scheme.
 *
 * ("A Rapid Hierarhical Rendering Technique for Translucent 
 *   Materials" by Herik Wann Jensen and Juan Buhler, in SIGGRAPH 02)
 */
class IsotropicDipole : public Subsurface {
public:
	IsotropicDipole(const Properties &props) 
		: Subsurface(props) {
		irrOctreeMutex->lock();
		m_octreeIndex = irrOctreeIndex++;
		irrOctreeMutex->unlock();
		
		/* How many samples should be taken when estimating 
		   the irradiance at a given point in the scene? */
		m_irrSamples = props.getInteger("irrSamples", 32);
				
		/* When estimating the irradiance at a given point, should indirect illumination be included
		   in the final estimate? */
		m_irrIndirect = props.getBoolean("irrIndirect", true);

		/* Multiplicative factor, which can be used to adjust the number of
		   irradiance samples */
		m_sampleMultiplier = props.getFloat("sampleMultiplier", 2.0f);
		/* Error threshold - lower means better quality */
		m_minDelta= props.getFloat("quality", 0.1f);
		/* Max. depth of the created octree */
		m_maxDepth = props.getInteger("maxDepth", 40);
		/* Multiplicative factor for the subsurface term */
		m_ssFactor = props.getSpectrum("ssFactor", Spectrum(1.0f));
		m_maxDepth = props.getInteger("maxDepth", 40);
	
		/* Asymmetry parameter of the phase function */
		m_g = props.getFloat("g", 0);
		m_ready = false;
		m_octreeResID = -1;

		lookupMaterial(props, m_sigmaS, m_sigmaA, &m_eta);
	}

	IsotropicDipole(Stream *stream, InstanceManager *manager) 
	 : Subsurface(stream, manager) {
		m_sigmaS = Spectrum(stream);
		m_sigmaA = Spectrum(stream);
		m_ssFactor = Spectrum(stream);
		m_g = stream->readFloat();
		m_eta = stream->readFloat();
		m_sampleMultiplier = stream->readFloat();
		m_minDelta = stream->readFloat();
		m_maxDepth = stream->readInt();
		m_octreeIndex = stream->readInt();
		m_irrSamples = stream->readInt();
		m_irrIndirect = stream->readBool();
		m_ready = false;
		m_octreeResID = -1;
		configure();
	}

	virtual ~IsotropicDipole() {
		if (m_octreeResID != -1)
			Scheduler::getInstance()->unregisterResource(m_octreeResID);
	}

	void bindUsedResources(ParallelProcess *proc) const {
		if (m_octreeResID != -1)
			proc->bindResource(formatString("irrOctree%i", m_octreeIndex), m_octreeResID);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Subsurface::serialize(stream, manager);
		m_sigmaS.serialize(stream);
		m_sigmaA.serialize(stream);
		m_ssFactor.serialize(stream);
		stream->writeFloat(m_g);
		stream->writeFloat(m_eta);
		stream->writeFloat(m_sampleMultiplier);
		stream->writeFloat(m_minDelta);
		stream->writeInt(m_maxDepth);
		stream->writeInt(m_octreeIndex);
		stream->writeInt(m_irrSamples);
		stream->writeBool(m_irrIndirect);
	}

	Spectrum Lo(const Scene *scene, Sampler *sampler,
			const Intersection &its, const Vector &d, int depth) const {
		if (!m_ready || m_ssFactor.isZero())
			return Spectrum(0.0f);
		IsotropicDipoleQuery query(m_zr, m_zv, m_sigmaTr, m_Fdt, its.p);
	
		const Normal &n = its.shFrame.n;
		m_octree->execute(query);

		if (m_eta == 1.0f) {
			return query.getResult() * m_ssFactor * INV_PI;
		} else {
			Float Ft = 1.0f - fresnel(absDot(n, d));
			return query.getResult() * m_ssFactor * INV_PI * (Ft / m_Fdr);
		}
	}

	void configure() {
		m_sigmaSPrime = m_sigmaS * (1-m_g);
		m_sigmaTPrime = m_sigmaSPrime + m_sigmaA;

		/* Mean-free path (avg. distance traveled through the medium) */
		m_mfp = Spectrum(1.0f) / m_sigmaTPrime;

		/* Also find the smallest mean-free path for all wavelengths */
		m_minMFP = std::numeric_limits<Float>::max();
		for (int lambda=0; lambda<SPECTRUM_SAMPLES; lambda++)
			m_minMFP = std::min(m_minMFP, m_mfp[lambda]);

		/* Average reflectance due to mismatched indices of refraction
		   at the boundary - [Groenhuis et al. 1983]*/
		m_Fdr = -1.440f / (m_eta * m_eta) + 0.710f / m_eta 
			+ 0.668f + 0.0636f * m_eta;

		/* Reduced albedo */
		m_alphaPrime = m_sigmaSPrime / m_sigmaTPrime;

		/* Average transmittance at the boundary */
		m_Fdt = 1.0f - m_Fdr;

		if (m_eta == 1.0f) {
			m_Fdr = (Float) 0.0f;
			m_Fdt = (Float) 1.0f;
		}

		/* Approximate dipole boundary condition term */
		m_A = (1 + m_Fdr) / m_Fdt;

		/* Effective transport extinction coefficient */
		m_sigmaTr = (m_sigmaA * m_sigmaTPrime * 3.0f).sqrt();

		/* Diffusion coefficient */
		m_D = Spectrum(1.0f) / (m_sigmaTPrime * 3.0f);

		/* Distance of the dipole point sources to the surface */
		m_zr = m_mfp; 
		m_zv = m_mfp * (1.0f + 4.0f/3.0f * m_A);
	}

	/// Unpolarized fresnel reflection term for dielectric materials
	Float fresnel(Float cosThetaI) const {
		Float g = std::sqrt(m_eta*m_eta - 1.0f + cosThetaI * cosThetaI);
		Float temp1 = (g - cosThetaI)/(g + cosThetaI);
		Float temp2 = (cosThetaI * (g + cosThetaI) - 1) / 
			(cosThetaI * (g - cosThetaI) + 1.0f);
		return 0.5f * temp1 * temp1 * (1.0f + temp2 * temp2);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int cameraResID, int samplerResID) {
		if (m_ready)
			return true;

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
		size_t sampleCount = (size_t) std::ceil(sa / (M_PI * m_minMFP * m_minMFP)
			* m_sampleMultiplier);
		Log(EInfo, "Generating " SIZE_T_FMT " irradiance samples..", sampleCount);

		ref<Scheduler> sched = Scheduler::getInstance();

		/* This could be a bit more elegant.. - inform the irradiance
		   sampler about the index of this subsurface integrator */
		std::vector<Subsurface *> ssIntegrators
			= scene->getSubsurfaceIntegrators();
		int index = -1;
		for (size_t i=0; i<ssIntegrators.size(); ++i) {
			if (ssIntegrators[i] == this) {
				index = (int) i;
				break;
			}
		}
		Assert(index != -1);

		ref<IrradianceSamplingProcess> proc = new IrradianceSamplingProcess(
			sampleCount, (size_t) std::ceil(sampleCount/100.0f), index, 
			m_irrSamples, m_irrIndirect, job);

		proc->bindResource("scene", sceneResID);
		scene->bindUsedResources(proc);
		m_proc = proc;
		sched->schedule(proc);
		sched->wait(proc);
		m_proc = NULL;
		if (proc->getReturnStatus() != ParallelProcess::ESuccess)
			return false;

		const IrradianceRecordVector &results = *proc->getSamples();
		for (size_t i=0; i<results.size(); ++i) 
			m_octree->addSample(results[i]);

		m_octree->preprocess();
		m_octreeResID = Scheduler::getInstance()->registerResource(m_octree);

		m_ready = true;
		return true;
	}

	void wakeup(std::map<std::string, SerializableObject *> &params) {
		std::string octreeName = formatString("irrOctree%i", m_octreeIndex);
		if (!m_octree.get() && params.find(octreeName) != params.end()) {
			m_octree = static_cast<IrradianceOctree *>(params[octreeName]);
			m_ready = true;
		}
	}

	void cancel() {
		Scheduler::getInstance()->cancel(m_proc);
	}

	MTS_DECLARE_CLASS()
private:
	Float m_minMFP, m_sampleMultiplier;
	Float m_Fdr, m_Fdt, m_A, m_minDelta, m_g, m_eta;
	Spectrum m_sigmaS, m_sigmaA;
	Spectrum m_mfp, m_sigmaTr, m_zr, m_zv, m_alphaPrime;
	Spectrum m_sigmaSPrime, m_sigmaTPrime, m_D, m_ssFactor;
	ref<IrradianceOctree> m_octree;
	ref<ParallelProcess> m_proc;
	int m_octreeResID, m_octreeIndex;
	int m_maxDepth;
	int m_irrSamples;
	bool m_irrIndirect;
	bool m_ready, m_requireSample;
};

MTS_IMPLEMENT_CLASS_S(IsotropicDipole, false, Subsurface)
MTS_EXPORT_PLUGIN(IsotropicDipole, "Isotropic dipole model");
MTS_NAMESPACE_END

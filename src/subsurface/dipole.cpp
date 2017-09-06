/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

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
#include <mitsuba/core/sse.h>
#include <mitsuba/core/ssemath.h>
#include "../medium/materials.h"
#include "irrtree.h"
#include "bluenoise.h"

MTS_NAMESPACE_BEGIN

/**
 * Computes the combined diffuse radiant exitance
 * caused by a number of dipole sources
 */
struct IsotropicDipoleQuery {
#if !defined(MTS_SSE) || SPECTRUM_SAMPLES != 3
    inline IsotropicDipoleQuery(const Spectrum &zr, const Spectrum &zv,
        const Spectrum &sigmaTr, const Point &p)
        : zr(zr), zv(zv), sigmaTr(sigmaTr), result(0.0f), p(p) {
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
        Spectrum dMo = Spectrum(INV_FOURPI) *
             (C1 * ((-sigmaTr * dr).exp()) / (dr * dr)
            + C2 * ((-sigmaTr * dv).exp()) / (dv * dv));

        result += dMo * sample.E * sample.area;
    }

    inline const Spectrum &getResult() const {
        return result;
    }

    const Spectrum &zr, &zv, &sigmaTr;
    Spectrum result;
#else
    inline IsotropicDipoleQuery(const Spectrum &_zr, const Spectrum &_zv,
        const Spectrum &_sigmaTr, const Point &p) : p(p) {
        zr = _mm_set_ps(_zr[0], _zr[1], _zr[2], 0);
        zv = _mm_set_ps(_zv[0], _zv[1], _zv[2], 0);
        sigmaTr = _mm_set_ps(_sigmaTr[0], _sigmaTr[1], _sigmaTr[2], 0);
        zrSqr = _mm_mul_ps(zr, zr);
        zvSqr = _mm_mul_ps(zv, zv);
        result.ps = _mm_setzero_ps();
    }

    inline void operator()(const IrradianceSample &sample) {
        /* Distance to the positive point source of the dipole */
        const __m128 lengthSquared = _mm_set1_ps((p - sample.p).lengthSquared()),
            drSqr = _mm_add_ps(zrSqr, lengthSquared),
            dvSqr = _mm_add_ps(zvSqr, lengthSquared),
            dr = _mm_sqrt_ps(drSqr), dv = _mm_sqrt_ps(dvSqr),
            one = _mm_set1_ps(1.0f),
            factor = _mm_mul_ps(_mm_set1_ps(INV_FOURPI*sample.area),
                _mm_set_ps(sample.E[0], sample.E[1], sample.E[2], 0)),
            C1fac = _mm_div_ps(_mm_mul_ps(zr, _mm_add_ps(sigmaTr, _mm_div_ps(one, dr))), drSqr),
            C2fac = _mm_div_ps(_mm_mul_ps(zv, _mm_add_ps(sigmaTr, _mm_div_ps(one, dv))), dvSqr),
            sigmaTrNeg = negate_ps(sigmaTr),
            exp1 = math::exp_ps(_mm_mul_ps(dr, sigmaTrNeg)),
            exp2 = math::exp_ps(_mm_mul_ps(dv, sigmaTrNeg));

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

    Point p;
};

static ref<Mutex> irrOctreeMutex = new Mutex();
static int irrOctreeIndex = 0;

/*!\plugin{dipole}{Dipole-based subsurface scattering model}
 * \parameters{
 *     \parameter{material}{\String}{
 *         Name of a material preset, see
 *         \tblref{medium-coefficients}. \default{\texttt{skin1}}
 *     }
 *     \parameter{sigmaA, sigmaS}{\Spectrum}{
 *         Absorption and scattering
 *         coefficients of the medium in inverse scene units.
 *         These parameters are mutually exclusive with \code{sigmaT} and \code{albedo}
 *         \default{configured based on \code{material}}
 *     }
 *     \parameter{sigmaT, albedo}{\Spectrum}{
 *         Extinction coefficient in inverse scene units
 *         and a (unitless) single-scattering albedo.
 *         These parameters are mutually exclusive with \code{sigmaA} and \code{sigmaS}
 *         \default{configured based on \code{material}}
 *     }
 *     \parameter{scale}{\Float}{
 *         Optional scale factor that will be applied to the \code{sigma*} parameters.
 *         It is provided for convenience when accomodating data based on different units,
 *         or to simply tweak the density of the medium. \default{1}}
 *     \parameter{intIOR}{\Float\Or\String}{Interior index of refraction specified
 *      numerically or using a known material name. \default{based on \code{material}}}
 *     \parameter{extIOR}{\Float\Or\String}{Exterior index of refraction specified
 *      numerically or using a known material name. \default{based on \code{material}}}
 *     \parameter{irrSamples}{\Integer}{
 *         Number of samples to use when estimating the
 *         irradiance at a point on the surface \default{16}
 *     }
 * }
 *
 * \renderings{
 *    \rendering{The material test ball rendered with the \code{skimmilk}
 *    material preset}{subsurface_dipole.jpg}
 *    \rendering{The material test ball rendered with the \code{skin1}
 *    material preset}{subsurface_dipole_2.jpg}
 * }
 * \renderings{
 *    \rendering{\code{scale=1}}{subsurface_dipole_dragon.jpg}
 *    \rendering{\code{scale=0.2}}{subsurface_dipole_dragon2.jpg}
 *    \caption{The dragon model rendered with the \code{skin2}
 *    material preset (model courtesy of XYZ RGB). The \code{scale}
 *    parameter is useful to communicate the relative size of
 *    an object to the viewer.}
 * }

 * This plugin implements the classic dipole subsurface scattering model
 * from radiative transport and medical physics \cite{Eason1978Theory,
 * Farrell1992Diffusion} in the form proposed by Jensen et al.
 * \cite{Jensen2001Practical}. It relies on the assumption that light entering
 * a material will undergo many (i.e. hundreds) of internal scattering
 * events, such that diffusion theory becomes applicable. In this
 * case\footnote{and after making several fairly strong simplifications:
 * the geometry is assumed to be a planar half-space, and the internal
 * scattering from the material boundary is only considered approximately.}
 * a simple analytic solution of the subsurface scattering profile is available
 * that enables simulating this effect without having to account for the vast
 * numbers of internal scattering events individually.
 *
 * For each \code{dipole} instance in the scene, the plugin adds a pre-process pass
 * to the rendering that computes the irradiance on a large set of sample positions
 * spread uniformly over the surface in question. The locations of these
 * points are chosen using a technique by Bowers et al. \cite{Bowers2010Parallel}
 * that creates particularly well-distributed (blue noise) samples. Later during
 * rendering, these  illumination samples are convolved with the diffusion profile
 * using a fast hierarchical technique proposed by Jensen and Buhler \cite{Jensen2005Rapid}.
 *
 * There are two different ways of configuring the medium properties.
 * One possibility is to load a material preset
 * using the \code{material} parameter---see \tblref{medium-coefficients}
 * for details. Alternatively, when specifying parameters by hand, they
 * can either be provided using the scattering and absorption coefficients,
 * or by declaring the extinction coefficient and single scattering albedo
 * (whichever is more convenient). Mixing these parameter initialization
 * methods is not allowed.
 *
 * All scattering parameters (named \code{sigma*}) should
 * be provided in inverse scene units. For instance, when a world-space
 * distance of 1 unit corresponds to a meter, the scattering coefficents must
 * be in units of inverse meters. For convenience, the \code{scale}
 * parameter can be used to correct this. For instance, when the scene is
 * in meters and the coefficients are in inverse millimeters, set
 * \code{scale=1000}.
 *
 * Note that a subsurface integrator can be associated with an \code{id}
 * and shared by several shapes using the reference mechanism introduced in
 * \secref{format}. This can be useful when an object is made up of many
 * separate sub-shapes.
 *
 * \renderings{
 *    \medrendering{Rendered using \pluginref{dipole}}{subsurface_dipole_bad1.jpg}
 *    \medrendering{Rendered using \pluginref{homogeneous}}{subsurface_dipole_bad2.jpg}
 *    \medrendering{\code{irrSamples} set too low}{subsurface_dipole_bad3.jpg}
 *    \caption{Two problem cases that may occur when rendering with the \pluginref{dipole}:
 *     \textbf{(a)-(b)}: These two renderings show a glass ball filled with diluted milk
 *     rendered using diffusion theory and radiative transport, respectively.
 *     The former produces an incorrect result, since the assumption of
 *     many scattering events breaks down.
 *     \textbf{(c)}: When the number of irradiance samples is too low when rendering
 *     with the dipole model, the resulting noise becomes visible as ``blotchy'' artifacts
 *     in the rendering.}
 * }
 *
 * \subsubsection*{Typical material setup}
 * To create a realistic material with subsurface scattering, it is necessary
 * to associate the underlying shape with an appropriately configured surface
 * and subsurface scattering model. Both should be aware of the material's
 * index of refraction.
 *
 * Because the \pluginref{dipole} plugin is responsible for all internal
 * scattering, the surface scattering model should only account for specular
 * reflection due to the index of refraction change. There are two models
 * in Mitsuba that can do this: \pluginref{plastic} and
 * \pluginref{roughplastic} (for smooth and rough interfaces, respectively).
 * An example is given on the next page.
 * \pagebreak
 * \begin{xml}
 * <shape type="...">
 *     <subsurface type="dipole">
 *         <string name="intIOR" value="water"/>
 *         <string name="extIOR" value="air"/>
 *         <rgb name="sigmaS" value="87.2, 127.2, 143.2"/>
 *         <rgb name="sigmaA" value="1.04, 5.6, 11.6"/>
 *         <integer name="irrSamples" value="64"/>
 *     </subsurface>
 *
 *     <bsdf type="plastic">
 *         <string name="intIOR" value="water"/>
 *         <string name="extIOR" value="air"/>
 *         <!-- Note: the diffuse component must be disabled! -->
 *         <spectrum name="diffuseReflectance" value="0"/>
 *     </bsdf>
 * <shape>
 * \end{xml}
 *
 * \remarks{
 *    \item This plugin only implements the multiple scattering component of
 *    the dipole model, i.e. single scattering is omitted. Furthermore, the
 *    numerous assumptions built into the underlying theory can cause severe
 *    inaccuracies.
 *
 *    For this reason, this plugin is the right choice for making pictures
 *    that ``look nice'', but it should be avoided when the output must hold
 *    up to real-world measurements. In this case, please use participating media
 *    (\secref{media}).
 *
 *   \item It is quite important that the \code{sigma*} parameters have the right units.
 *   For instance: if the \code{sigmaT} parameter is accidentally set to a value that
 *   is too small by a factor of 1000, the plugin will attempt to create
 *   one million times as many irradiance samples, which will likely cause
 *   the rendering process to crash with an ``out of memory'' failure.
 * }
 */

class IsotropicDipole : public Subsurface {
public:
    IsotropicDipole(const Properties &props)
        : Subsurface(props) {
        {
            LockGuard lock(irrOctreeMutex);
            m_octreeIndex = irrOctreeIndex++;
        }

        /* How many samples should be taken when estimating
           the irradiance at a given point in the scene? */
        m_irrSamples = props.getInteger("irrSamples", 16);

        /* When estimating the irradiance at a given point,
           should indirect illumination be included in the final estimate? */
        m_irrIndirect = props.getBoolean("irrIndirect", true);

        /* Multiplicative factor, which can be used to adjust the number of
           irradiance samples */
        m_sampleMultiplier = props.getFloat("sampleMultiplier", 1.0f);

        /* Error threshold - lower means better quality */
        m_quality = props.getFloat("quality", 0.2f);

        /* Asymmetry parameter of the phase function */
        m_octreeResID = -1;

        lookupMaterial(props, m_sigmaS, m_sigmaA, m_g, &m_eta);
    }

    IsotropicDipole(Stream *stream, InstanceManager *manager)
     : Subsurface(stream, manager) {
        m_sigmaS = Spectrum(stream);
        m_sigmaA = Spectrum(stream);
        m_g = Spectrum(stream);
        m_eta = stream->readFloat();
        m_sampleMultiplier = stream->readFloat();
        m_quality = stream->readFloat();
        m_octreeIndex = stream->readInt();
        m_irrSamples = stream->readInt();
        m_irrIndirect = stream->readBool();
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
        m_g.serialize(stream);
        stream->writeFloat(m_eta);
        stream->writeFloat(m_sampleMultiplier);
        stream->writeFloat(m_quality);
        stream->writeInt(m_octreeIndex);
        stream->writeInt(m_irrSamples);
        stream->writeBool(m_irrIndirect);
    }

    Spectrum Lo(const Scene *scene, Sampler *sampler,
            const Intersection &its, const Vector &d, int depth) const {
        if (!m_active || dot(its.shFrame.n, d) < 0)
            return Spectrum(0.0f);
        IsotropicDipoleQuery query(m_zr, m_zv, m_sigmaTr, its.p);

        m_octree->performQuery(query);
        Spectrum result(query.getResult() * INV_PI);

        if (m_eta != 1.0f)
            result *= 1.0f - fresnelDielectricExt(dot(its.shFrame.n, d), m_eta);

        return result;
    }

    void configure() {
        m_sigmaSPrime = m_sigmaS * (Spectrum(1.0f) - m_g);
        m_sigmaTPrime = m_sigmaSPrime + m_sigmaA;

        /* Find the smallest mean-free path over all wavelengths */
        Spectrum mfp = Spectrum(1.0f) / m_sigmaTPrime;
        m_radius = std::numeric_limits<Float>::max();
        for (int lambda=0; lambda<SPECTRUM_SAMPLES; lambda++)
            m_radius = std::min(m_radius, mfp[lambda]);

        /* Average diffuse reflectance due to mismatched indices of refraction */
        m_Fdr = fresnelDiffuseReflectance(1 / m_eta);

        /* Dipole boundary condition distance term */
        Float A = (1 + m_Fdr) / (1 - m_Fdr);

        /* Effective transport extinction coefficient */
        m_sigmaTr = (m_sigmaA * m_sigmaTPrime * 3.0f).sqrt();

        /* Distance of the two dipole point sources to the surface */
        m_zr = mfp;
        m_zv = mfp * (1.0f + 4.0f/3.0f * A);
    }

    bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
        int sceneResID, int cameraResID, int _samplerResID) {
        if (m_octree)
            return true;

        if (!scene->getIntegrator()->getClass()
                ->derivesFrom(MTS_CLASS(SamplingIntegrator)))
            Log(EError, "The dipole subsurface scattering model requires "
                "a sampling-based surface integrator!");

        ref<Scheduler> sched = Scheduler::getInstance();
        ref<Timer> timer = new Timer();

        AABB aabb;
        Float sa;

        ref<PositionSampleVector> points = new PositionSampleVector();
        /* It is necessary to increase the sampling resolution to
           prevent low-frequency noise in the output */
        Float actualRadius = m_radius / std::sqrt(m_sampleMultiplier * 20);
        blueNoisePointSet(scene, m_shapes, actualRadius, points, sa, aabb, job);

        /* 2. Gather irradiance in parallel */
        const Sensor *sensor = scene->getSensor();
        ref<IrradianceSamplingProcess> proc = new IrradianceSamplingProcess(
            points, 1024, m_irrSamples, m_irrIndirect,
            sensor->getShutterOpen() + 0.5f * sensor->getShutterOpenTime(), job);

        /* Create a sampler instance for every core */
        ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
            createObject(MTS_CLASS(Sampler), Properties("independent")));
        std::vector<SerializableObject *> samplers(sched->getCoreCount());
        for (size_t i=0; i<sched->getCoreCount(); ++i) {
            ref<Sampler> clonedSampler = sampler->clone();
            clonedSampler->incRef();
            samplers[i] = clonedSampler.get();
        }

        int samplerResID = sched->registerMultiResource(samplers);
        int integratorResID = sched->registerResource(
            const_cast<Integrator *>(scene->getIntegrator()));

        proc->bindResource("scene", sceneResID);
        proc->bindResource("integrator", integratorResID);
        proc->bindResource("sampler", samplerResID);
        scene->bindUsedResources(proc);
        m_proc = proc;
        sched->schedule(proc);
        sched->wait(proc);
        m_proc = NULL;
        for (size_t i=0; i<samplers.size(); ++i)
            samplers[i]->decRef();

        sched->unregisterResource(samplerResID);
        sched->unregisterResource(integratorResID);
        if (proc->getReturnStatus() != ParallelProcess::ESuccess)
            return false;

        Log(EDebug, "Done gathering (took %i ms), clustering ..", timer->getMilliseconds());
        timer->reset();

        std::vector<IrradianceSample> &samples = proc->getIrradianceSampleVector()->get();
        sa /= samples.size();

        for (size_t i=0; i<samples.size(); ++i)
            samples[i].area = sa;

        m_octree = new IrradianceOctree(aabb, m_quality, samples);

        Log(EDebug, "Done clustering (took %i ms).", timer->getMilliseconds());
        m_octreeResID = Scheduler::getInstance()->registerResource(m_octree);

        return true;
    }

    void wakeup(ConfigurableObject *parent,
        std::map<std::string, SerializableObject *> &params) {
        std::string octreeName = formatString("irrOctree%i", m_octreeIndex);
        if (!m_octree.get() && params.find(octreeName) != params.end()) {
            m_octree = static_cast<IrradianceOctree *>(params[octreeName]);
            m_active = true;
        }
    }

    void cancel() {
        Scheduler::getInstance()->cancel(m_proc);
    }

    MTS_DECLARE_CLASS()
private:
    Float m_radius, m_sampleMultiplier;
    Float m_Fdr, m_quality, m_eta;
    Spectrum m_sigmaS, m_sigmaA, m_g;
    Spectrum m_sigmaTr, m_zr, m_zv;
    Spectrum m_sigmaSPrime, m_sigmaTPrime;
    ref<IrradianceOctree> m_octree;
    ref<ParallelProcess> m_proc;
    int m_octreeResID, m_octreeIndex;
    int m_irrSamples;
    bool m_irrIndirect;
};

MTS_IMPLEMENT_CLASS_S(IsotropicDipole, false, Subsurface)
MTS_EXPORT_PLUGIN(IsotropicDipole, "Isotropic dipole model");
MTS_NAMESPACE_END

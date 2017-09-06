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
#include <mitsuba/core/statistics.h>
#include <boost/math/distributions/normal.hpp>

MTS_NAMESPACE_BEGIN

/*!\plugin{adaptive}{Adaptive integrator}
 * \order{13}
 * \parameters{
 *     \parameter{maxError}{\Float}{Maximum relative error
 *         threshold\default{0.05}}
 *     \parameter{pValue}{\Float}{
 *         Required p-value to accept a sample \default{0.05}
 *     }
 *     \parameter{maxSampleFactor}{\Integer}{
 *         Maximum number of samples to be generated \emph{relative} to the
 *         number of configured pixel samples. The adaptive integrator
 *         will stop after this many samples, regardless of whether
 *         or not the error criterion was satisfied.
 *         A negative value will be interpreted as $\infty$.
 *         \default{32---for instance, when 64 pixel samples are configured in
 *         the \code{sampler}, this means that the adaptive integrator
 *         will give up after 32*64=2048 samples}
 *     }
 * }
 *
 * This ``meta-integrator'' repeatedly invokes a provided sub-integrator
 * until the computed radiance values satisfy a specified relative error bound
 * (5% by default) with a certain probability (95% by default). Internally,
 * it uses a Z-test to decide when to stop collecting samples. While repeatedly
 * applying a Z-test in this manner is not good practice in terms of
 * a rigorous statistical analysis, it provides a useful mathematically
 * motivated stopping criterion.
 *
 * \begin{xml}[caption={An example how to make the \pluginref{path} integrator adaptive}]
 * <integrator type="adaptive">
 *     <integrator type="path"/>
 * </integrator>
 * \end{xml}
 *
 * \remarks{
 *    \item The adaptive integrator needs a variance estimate to work
 *     correctly. Hence, the underlying sample generator should be set to a reasonably
 *     large number of pixel samples (e.g. 64 or higher) so that this estimate can be obtained.
 *    \item This plugin uses a relatively simplistic error heuristic that does not
 *    share information between pixels and only reasons about variance in image space.
 *    In the future, it will likely be replaced with something more robust.
 * }
 */
class AdaptiveIntegrator : public SamplingIntegrator {
public:
    AdaptiveIntegrator(const Properties &props) : SamplingIntegrator(props) {
        /* Maximum relative error threshold */
        m_maxError = props.getFloat("maxError", 0.05f);
        /* Maximum number of samples to take (relative to the number of pixel samples
           that were configured in the sampler). The sample collection
           will stop after this many samples even if the variance is still
           too high. A negative value will be interpreted as infinity. */
        m_maxSampleFactor = props.getInteger("maxSampleFactor", 32);
        /* Required P-value to accept a sample. */
        m_pValue = props.getFloat("pValue", 0.05f);
        m_verbose = props.getBoolean("verbose", false);
    }

    AdaptiveIntegrator(Stream *stream, InstanceManager *manager)
     : SamplingIntegrator(stream, manager) {
        m_subIntegrator = static_cast<SamplingIntegrator *>(manager->getInstance(stream));
        m_maxSampleFactor = stream->readInt();
        m_maxError = stream->readFloat();
        m_quantile = stream->readFloat();
        m_averageLuminance = stream->readFloat();
        m_pValue = stream->readFloat();
        m_verbose = false;
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        const Class *cClass = child->getClass();

        if (cClass->derivesFrom(MTS_CLASS(Integrator))) {
            if (!cClass->derivesFrom(MTS_CLASS(SamplingIntegrator)))
                Log(EError, "The sub-integrator must be derived from the class SamplingIntegrator");
            m_subIntegrator = static_cast<SamplingIntegrator *>(child);
        } else {
            Integrator::addChild(name, child);
        }
    }

    void configureSampler(const Scene *scene, Sampler *sampler) {
        SamplingIntegrator::configureSampler(scene, sampler);
        m_subIntegrator->configureSampler(scene, sampler);
    }

    bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID) {
        if (!SamplingIntegrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID))
            return false;
        if (m_subIntegrator == NULL)
            Log(EError, "No sub-integrator was specified!");
        Sampler *sampler = static_cast<Sampler *>(Scheduler::getInstance()->getResource(samplerResID, 0));
        Sensor *sensor = static_cast<Sensor *>(Scheduler::getInstance()->getResource(sensorResID));
        if (sampler->getClass()->getName() != "IndependentSampler")
            Log(EError, "The error-controlling integrator should only be "
                "used in conjunction with the independent sampler");
        if (!m_subIntegrator->preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID))
            return false;

        Vector2i filmSize = sensor->getFilm()->getSize();
        bool needsApertureSample = sensor->needsApertureSample();
        bool needsTimeSample = sensor->needsTimeSample();
        const int nSamples = 10000;
        Float luminance = 0;

        Point2 apertureSample(0.5f);
        Float timeSample = 0.5f;
        RadianceQueryRecord rRec(scene, sampler);

        /* Estimate the overall luminance on the image plane */
        for (int i=0; i<nSamples; ++i) {
            sampler->generate(Point2i(0));

            rRec.newQuery(RadianceQueryRecord::ERadiance, sensor->getMedium());
            rRec.extra = RadianceQueryRecord::EAdaptiveQuery;

            Point2 samplePos(rRec.nextSample2D());
            samplePos.x *= filmSize.x;
            samplePos.y *= filmSize.y;

            if (needsApertureSample)
                apertureSample = rRec.nextSample2D();
            if (needsTimeSample)
                timeSample = rRec.nextSample1D();

            RayDifferential eyeRay;
            Spectrum sampleValue = sensor->sampleRay(
                eyeRay, samplePos, apertureSample, timeSample);

            sampleValue *= m_subIntegrator->Li(eyeRay, rRec);
            luminance += sampleValue.getLuminance();
        }

        m_averageLuminance = luminance / nSamples;

        boost::math::normal dist(0, 1);
        m_quantile = (Float) boost::math::quantile(dist, 1-m_pValue/2);
        Log(EInfo, "Configuring for a %.1f%% confidence interval, quantile=%f, avg. luminance=%f",
            (1-m_pValue)*100, m_quantile, m_averageLuminance);
        return true;
    }

    void renderBlock(const Scene *scene, const Sensor *sensor,
            Sampler *sampler, ImageBlock *block, const bool &stop,
            const std::vector< TPoint2<uint8_t> > &points) const {
        typedef TSpectrum<Float, SPECTRUM_SAMPLES + 2> SpectrumAlphaWeight;

        bool needsApertureSample = sensor->needsApertureSample();
        bool needsTimeSample = sensor->needsTimeSample();

        if (sampler->getSampleCount() < 8)
            Log(EError, "Starting the adaptive integrator with less than 8 "
                "samples per pixel does not make much sense -- giving up.");

        RayDifferential eyeRay;
        RadianceQueryRecord rRec(scene, sampler);

        Float diffScaleFactor = 1.0f /
            std::sqrt((Float) sampler->getSampleCount());

        Point2 apertureSample(0.5f);
        Float timeSample = 0.5f;
        int borderSize = sensor->getFilm()->getReconstructionFilter()->getBorderSize();

        size_t sampleCount;
        block->clear();

        SpectrumAlphaWeight *target = (SpectrumAlphaWeight *) block->getBitmap()->getUInt8Data();
        SpectrumAlphaWeight *snapshot = (SpectrumAlphaWeight *) alloca(sizeof(SpectrumAlphaWeight)
            * (2*borderSize+1)*(2*borderSize+1));

        for (size_t i=0; i<points.size(); ++i) {
            Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
            sampler->generate(offset);

            /* Before starting to place samples within the area of a single pixel, the
               following takes a snapshot of all surrounding spectrum+weight+alpha
               values. Those are then used later to ensure that adjacent pixels will
               not be disproportionately biased by this pixel's contributions. */
            for (int y=0; y<2*borderSize+1; ++y) {
                SpectrumAlphaWeight *src = target + ((y+points[i].y)
                    * block->getBitmap()->getWidth() + points[i].x);
                SpectrumAlphaWeight *dst = snapshot + y*(2*borderSize+1);
                memcpy(dst, src, sizeof(SpectrumAlphaWeight) * (2*borderSize+1));
            }

            Float mean = 0, meanSqr = 0.0f;
            sampleCount = 0;

            while (true) {
                if (stop)
                    return;

                rRec.newQuery(RadianceQueryRecord::ESensorRay, sensor->getMedium());
                rRec.extra = RadianceQueryRecord::EAdaptiveQuery;

                Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));
                if (needsApertureSample)
                    apertureSample = rRec.nextSample2D();
                if (needsTimeSample)
                    timeSample = rRec.nextSample1D();

                Spectrum sampleValue = sensor->sampleRayDifferential(
                    eyeRay, samplePos, apertureSample, timeSample);
                eyeRay.scaleDifferential(diffScaleFactor);

                sampleValue *= m_subIntegrator->Li(eyeRay, rRec);

                Float sampleLuminance;
                if (block->put(samplePos, sampleValue, rRec.alpha)) {
                    /* Check for problems with the sample */
                    sampleLuminance = sampleValue.getLuminance();
                } else {
                    sampleLuminance = 0.0f;
                }
                ++sampleCount;
                sampler->advance();

                /* Numerically robust online variance estimation using an
                   algorithm proposed by Donald Knuth (TAOCP vol.2, 3rd ed., p.232) */
                const Float delta = sampleLuminance - mean;
                mean += delta / sampleCount;
                meanSqr += delta * (sampleLuminance - mean);

                if (m_maxSampleFactor >= 0 && sampleCount >= m_maxSampleFactor * sampler->getSampleCount()) {
                    break;
                } else if (sampleCount >= sampler->getSampleCount()) {
                    /* Variance of the primary estimator */
                    const Float variance = meanSqr / (sampleCount-1);

                    Float stdError = std::sqrt(variance/sampleCount);

                    /* Half width of the confidence interval */
                    Float ciWidth = stdError * m_quantile;

                    /* Relative error heuristic */
                    Float base = std::max(mean, m_averageLuminance * 0.01f);

                    if (m_verbose && (sampleCount % 100) == 0)
                        Log(EDebug, "%i samples, mean=%f, stddev=%f, std error=%f, ci width=%f, max allowed=%f", sampleCount, mean,
                            std::sqrt(variance), stdError, ciWidth, base * m_maxError);

                    if (ciWidth <= m_maxError * base)
                        break;
                }
            }

            /* Ensure that a large amounts of samples in one pixel do not
               bias neighboring pixels (due to the reconstruction filter) */
            Float factor = 1.0f / sampleCount;
            for (int y=0; y<2*borderSize+1; ++y) {
                SpectrumAlphaWeight *dst = target + ((y+points[i].y)
                    * block->getBitmap()->getWidth() + points[i].x);
                SpectrumAlphaWeight *backup = snapshot + y*(2*borderSize+1);

                for (int x=0; x<2*borderSize+1; ++x)
                    dst[x] = backup[x] * (1-factor) + dst[x] * factor;
            }
        }
    }

    Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
        return m_subIntegrator->Li(ray, rRec);
    }

    Spectrum E(const Scene *scene, const Intersection &its, const Medium *medium,
            Sampler *sampler, int nSamples, bool includeIndirect) const {
        return m_subIntegrator->E(scene, its, medium,
            sampler, nSamples, includeIndirect);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);
        manager->serialize(stream, m_subIntegrator.get());

        stream->writeInt(m_maxSampleFactor);
        stream->writeFloat(m_maxError);
        stream->writeFloat(m_quantile);
        stream->writeFloat(m_averageLuminance);
        stream->writeFloat(m_pValue);
    }

    void bindUsedResources(ParallelProcess *proc) const {
        m_subIntegrator->bindUsedResources(proc);
    }

    void wakeup(ConfigurableObject *parent,
            std::map<std::string, SerializableObject *> &params) {
        m_subIntegrator->wakeup(this, params);
    }

    void cancel() {
        SamplingIntegrator::cancel();
        m_subIntegrator->cancel();
    }

    const Integrator *getSubIntegrator(int idx) const {
        if (idx != 0)
            return NULL;
        return m_subIntegrator.get();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "AdaptiveIntegrator[" << endl
            << "  maxSamples = " << m_maxSampleFactor << "," << endl
            << "  maxError = " << m_maxError << "," << endl
            << "  quantile = " << m_quantile << "," << endl
            << "  pvalue = " << m_pValue << "," << endl
            << "  subIntegrator = " << indent(m_subIntegrator->toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    ref<SamplingIntegrator> m_subIntegrator;
    Float m_maxError, m_quantile, m_pValue, m_averageLuminance;
    int m_maxSampleFactor;
    bool m_verbose;
};

MTS_IMPLEMENT_CLASS_S(AdaptiveIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(AdaptiveIntegrator, "Adaptive integrator");
MTS_NAMESPACE_END

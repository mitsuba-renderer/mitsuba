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

#include <mitsuba/bidir/util.h>
#include <mitsuba/bidir/mempool.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/serialization.h>
#include <mitsuba/core/cobject.h>

MTS_NAMESPACE_BEGIN

ref<Bitmap> BidirectionalUtils::renderDirectComponent(Scene *scene, int sceneResID,
        int sensorResID, RenderQueue *queue, const RenderJob *job, size_t directSamples) {
    ref<PluginManager> pluginMgr = PluginManager::getInstance();
    ref<Scheduler> scheduler = Scheduler::getInstance();
    const Film *film = scene->getFilm();
    Integrator *integrator = scene->getIntegrator();
    /* Render the direct illumination component separately */
    ref<Bitmap> directImage = new Bitmap(Bitmap::ERGBA, Bitmap::EFloat32, film->getCropSize());
    bool hasMedia = scene->getMedia().size() > 0;
    bool hasDOF = scene->getSensor()->needsApertureSample();
    size_t pixelSamples = directSamples;
    Properties integratorProps(hasMedia ? "volpath" : "direct");

    if (hasMedia || hasDOF) {
        integratorProps.setInteger("maxDepth", 2);
    } else {
        /* No participating media / DoF -> we can more carefully
           distribute samples between shading and visibility */
        int shadingSamples = 1;
        while (pixelSamples > 8) {
            pixelSamples /= 2;
            shadingSamples *= 2;
        }
        integratorProps.setSize("shadingSamples", shadingSamples);
    }

    ref<Integrator> directIntegrator = static_cast<Integrator *> (pluginMgr->
            createObject(Integrator::m_theClass, integratorProps));
    /* Create a low discrepancy sampler instance for every core */
    Properties samplerProps("ldsampler");
    samplerProps.setSize("sampleCount", pixelSamples);
    ref<Sampler> ldSampler = static_cast<Sampler *> (pluginMgr->
            createObject(Sampler::m_theClass, samplerProps));
    ldSampler->configure();
    directIntegrator->configure();
    directIntegrator->configureSampler(scene, ldSampler);
    std::vector<SerializableObject *> samplers(scheduler->getCoreCount());
    for (size_t i=0; i<scheduler->getCoreCount(); ++i) {
        ref<Sampler> clonedSampler = ldSampler->clone();
        clonedSampler->incRef();
        samplers[i] = clonedSampler.get();
    }
    int ldSamplerResID = scheduler->registerMultiResource(samplers);
    for (size_t i=0; i<scheduler->getCoreCount(); ++i)
        samplers[i]->decRef();

    integrator->incRef();
    scene->setIntegrator(directIntegrator);
    bool success = directIntegrator->render(scene, queue, job,
        sceneResID, sensorResID, ldSamplerResID);
    scene->setIntegrator(integrator);
    integrator->decRef();
    scheduler->unregisterResource(ldSamplerResID);

    if (success) {
        ref<Bitmap> bitmap = new Bitmap(
            Bitmap::ESpectrum, Bitmap::EFloat,
            film->getCropSize());
        film->develop(Point2i(0, 0),
            film->getCropSize(), Point2i(0, 0), bitmap);
        return bitmap;
    } else {
        return NULL;
    }
}

ref<Bitmap> BidirectionalUtils::mltLuminancePass(Scene *scene, int sceneResID,
        RenderQueue *queue, int sizeFactor, ref<RenderJob> &nestedJob) {
    ref<PluginManager> pluginMgr = PluginManager::getInstance();
    ref<Scheduler> scheduler = Scheduler::getInstance();
    Properties integratorProps = scene->getIntegrator()->getProperties();

    Vector2i origCropSize   = scene->getFilm()->getCropSize();
    Vector2i origSize       = scene->getFilm()->getSize();

    Vector2i reducedSize = Vector2i(
        std::max(1, origSize.x / sizeFactor),
        std::max(1, origSize.y / sizeFactor));

    Vector2i reducedCropSize = Vector2i(
        std::max(1, origCropSize.x / sizeFactor),
        std::max(1, origCropSize.y / sizeFactor));

    Point2i reducedCropOffset =
        scene->getFilm()->getCropOffset()/sizeFactor;

    size_t sampleCount = scene->getSampler()->getSampleCount();
    const Sensor *sensor = scene->getSensor();

    Properties filmProps("hdrfilm");
    filmProps.setInteger("width", reducedSize.x, false);
    filmProps.setInteger("height", reducedSize.y, false);
    filmProps.setInteger("cropWidth", reducedCropSize.x, false);
    filmProps.setInteger("cropHeight", reducedCropSize.y, false);
    filmProps.setInteger("cropOffsetX", reducedCropOffset.x, false);
    filmProps.setInteger("cropOffsetY", reducedCropOffset.x, false);
    ref<Film> nestedFilm = static_cast<Film *>(
        pluginMgr->createObject(Film::m_theClass, filmProps));
    nestedFilm->configure();

    /* Use a higher number of mutations/pixel compared to the second stage */
    Properties samplerProps("independent");
    samplerProps.setSize("sampleCount", sampleCount * sizeFactor);
    ref<Sampler> nestedSampler = static_cast<Sampler *>(
        pluginMgr->createObject(Sampler::m_theClass, samplerProps));
    nestedSampler->configure();
    std::vector<SerializableObject *> samplers(scheduler->getCoreCount());
    for (size_t i=0; i<scheduler->getCoreCount(); ++i) {
        ref<Sampler> clonedSampler = nestedSampler->clone();
        clonedSampler->incRef();
        samplers[i] = clonedSampler.get();
    }
    int nestedSamplerResID = scheduler->registerMultiResource(samplers);
    for (size_t i=0; i<scheduler->getCoreCount(); ++i)
        samplers[i]->decRef();

    /* Configure the sensor */
    Properties sensorProps = sensor->getProperties();
    ref<Sensor> nestedSensor = static_cast<Sensor *>
        (pluginMgr->createObject(Sensor::m_theClass, sensorProps));
    nestedSensor->addChild(nestedSampler);
    nestedSensor->addChild(nestedFilm);
    nestedSensor->configure();
    int nestedSensorResID = scheduler->registerResource(nestedSensor);

    integratorProps.setBoolean("firstStage", true, false);
    ref<Integrator> nestedIntegrator = static_cast<Integrator *> (pluginMgr->
            createObject(Integrator::m_theClass, integratorProps));

    ref<Scene> nestedScene = new Scene(scene);
    nestedScene->setSensor(nestedSensor);
    nestedScene->setIntegrator(nestedIntegrator);
    nestedScene->configure();
    nestedScene->initialize();

    nestedJob = new RenderJob("mlti", nestedScene, queue,
        sceneResID, nestedSensorResID, nestedSamplerResID);

    nestedJob->start();
    if (!nestedJob->wait()) {
        nestedJob = NULL;
        scheduler->unregisterResource(nestedSensorResID);
        scheduler->unregisterResource(nestedSamplerResID);
        return NULL;
    }
    nestedJob = NULL;

    scheduler->unregisterResource(nestedSensorResID);
    scheduler->unregisterResource(nestedSamplerResID);

    /* Instantiate a Gaussian reconstruction filter */
    ref<ReconstructionFilter> rfilter = static_cast<ReconstructionFilter *> (
        PluginManager::getInstance()->createObject(
        MTS_CLASS(ReconstructionFilter), Properties("gaussian")));
    rfilter->configure();

    /* Develop the rendered image into a luminance bitmap */
    ref<Bitmap> luminanceMap = new Bitmap(Bitmap::ELuminance,
        Bitmap::EFloat, reducedCropSize);
    nestedFilm->develop(Point2i(0, 0), reducedCropSize,
        Point2i(0, 0), luminanceMap);

    /* Up-sample the low resolution luminance map */
    luminanceMap = luminanceMap->resample(rfilter,
        ReconstructionFilter::EClamp,
        ReconstructionFilter::EClamp, origCropSize,
        0.0f, std::numeric_limits<Float>::infinity());

    return luminanceMap;
}

MTS_NAMESPACE_END

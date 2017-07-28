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

#include <mitsuba/core/statistics.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/qmc.h>
#include <mitsuba/hw/vpl.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/hw/gpugeometry.h>
#include <mitsuba/hw/gputexture.h>
#include "../shapes/instance.h"
#include "data/shaders.h"

/**
 * \brief Use paraboloid maps instead of hemicube shadow maps?
 *
 * On current hardware, and with the implemented optimizations,
 * this technique sadly does not yet pull its own weight.
 */
//#define MTS_VPL_USE_PARABOLOID_MAPS

/**
 * \brief Render cube/hemicube shadow maps in a single pass?
 * (using multiple render targets and geometry shaders)
 *
 * Some benchmarking on different ATI/NVidia hardware shows that
 * this is *always* slower! So it is disabled for now.
 */
//#define MTS_VPL_USE_SINGLE_PASS

#define SUPPLEMENTAL_CODE_MARKER "{{ SUPPLEMENTAL CODE }}"

MTS_NAMESPACE_BEGIN

static StatsCounter statsMaxResidentShaders("VPL renderer", "Max. # of resident shaders", EMaximumValue);
static StatsCounter statsMaxTriangles("VPL renderer", "Max. # of drawn triangles", EMaximumValue);

VPLShaderManager::VPLShaderManager(Renderer *renderer)
        : m_renderer(renderer) {
    m_shadowMapResolution = 0;
    m_clamping = 0.1f;
    m_diffuseSources = true;
    m_diffuseReceivers = false;
    m_vplIndex = 0;
}

VPLShaderManager::~VPLShaderManager() {
    cleanup();
}

void VPLShaderManager::init() {
    if (!m_shadowGen) {
        m_shadowGen = new ShadowMapGenerator(m_renderer);
        m_shadowGen->init();
    }
}

void VPLShaderManager::cleanup() {
    setScene(NULL);
    if (m_shadowMapCube) {
        m_shadowMapCube->cleanup();
        m_shadowMapCube = NULL;
    }
    if (m_shadowMap2D) {
        m_shadowMap2D->cleanup();
        m_shadowMap2D = NULL;
    }
    if (m_shadowGen) {
        m_shadowGen->cleanup();
        m_shadowGen = NULL;
    }
}

void VPLShaderManager::setScene(const Scene *scene) {
    if (scene == m_scene || (scene && m_scene && scene->getShapes() == m_scene->getShapes())) {
        m_scene = scene;
        return; /* This is still the same scene! Do nothing.. */
    }

    if (m_scene) {
        /* Unregister any content from the previous scene (which was uploaded to the GPU) */
        const ref_vector<Shape> &oldShapes = m_scene->getShapes();
        for (size_t i=0; i<oldShapes.size(); ++i) {
            const Shape *shape = oldShapes[i].get();

            if (shape->getClass()->getName() == "Instance") {
                const Instance *instance = static_cast<const Instance *>(shape);
                const std::vector<const Shape *> &instantiatedShapes =
                        instance->getShapeGroup()->getKDTree()->getShapes();

                for (size_t j=0; j<instantiatedShapes.size(); ++j) {
                    shape = instantiatedShapes[j];
                    if (!m_renderer->unregisterGeometry(shape))
                        continue;
                    const BSDF *bsdf = shape->getBSDF();
                    if (!bsdf)
                        bsdf = const_cast<Shape *>(shape)->createTriMesh()->getBSDF();
                    m_renderer->unregisterShaderForResource(bsdf);
                }
            } else {
                const BSDF *bsdf = shape->getBSDF();
                if (!bsdf)
                    bsdf = const_cast<Shape *>(shape)->createTriMesh()->getBSDF();

                if (!m_renderer->unregisterGeometry(shape))
                    continue;
                m_renderer->unregisterShaderForResource(bsdf);
            }
        }

        const ref_vector<Emitter> &emitters = m_scene->getEmitters();
        for (size_t i=0; i<emitters.size(); ++i)
            m_renderer->unregisterShaderForResource(emitters[i].get());

        /* Release all active GPU programs */
        std::map<std::string, VPLConfiguration>::iterator it = m_configurations.begin();
        for (; it != m_configurations.end(); ++it) {
            GPUProgram *program = (*it).second.program;
            program->cleanup();
            program->decRef();
        }
        m_configurations.clear();

        if (m_backgroundProgram) {
            m_backgroundProgram->cleanup();
            m_backgroundProgram = NULL;
        }
    }

    m_scene = scene;
    m_vplIndex = 0;

    if (!scene)
        return;

    const ref_vector<Shape> &shapes = scene->getShapes();
    m_geometry.clear();
    m_geometry.reserve(shapes.size());
    m_opaqueGeometry.clear();
    m_opaqueGeometry.reserve(shapes.size());
    m_animatedGeometry.clear();

    Matrix4x4 identityTrafo;
    identityTrafo.setIdentity();

    /* Upload all geometry to the GPU, create shaders for scattering models */
    for (size_t i=0; i<shapes.size(); ++i) {
        const Shape *shape = shapes[i].get();

        if (shape->getClass()->getName() == "Instance") {
            const Instance *instance = static_cast<const Instance *>(shape);
            const std::vector<const Shape *> &instantiatedShapes =
                    instance->getShapeGroup()->getKDTree()->getShapes();
            const AnimatedTransform *atrafo = instance->getWorldTransform();
            const Matrix4x4 &trafo = atrafo->eval(0).getMatrix();

            for (size_t j=0; j<instantiatedShapes.size(); ++j) {
                shape = instantiatedShapes[j];
                GPUGeometry *gpuGeo = m_renderer->registerGeometry(shape);
                if (!gpuGeo)
                    continue;

                const BSDF *bsdf = shape->getBSDF();
                if (!bsdf)
                    bsdf = gpuGeo->getTriMesh()->getBSDF();

                Shader *shader = m_renderer->registerShaderForResource(bsdf);
                if (shader && !shader->isComplete()) {
                    m_renderer->unregisterShaderForResource(bsdf);
                    shader = NULL;
                }

                gpuGeo->setShader(shader);
                ssize_t geometryIndex = (ssize_t) m_geometry.size(), opaqueGeometryIndex = -1;
                m_geometry.push_back(std::make_pair(gpuGeo, trafo));

                if (shader && !(shader->getFlags() & Shader::ETransparent)) {
                    opaqueGeometryIndex = (ssize_t) m_opaqueGeometry.size();
                    m_opaqueGeometry.push_back(std::make_pair(gpuGeo, trafo));
                }

                if (!atrafo->isStatic()) {
                    m_animatedGeometry.push_back(AnimatedGeometryRecord(atrafo,
                        geometryIndex, opaqueGeometryIndex));
                }
            }
        } else {
            GPUGeometry *gpuGeo = m_renderer->registerGeometry(shape);
            if (!gpuGeo)
                continue;
            const BSDF *bsdf = shape->getBSDF();
            if (!bsdf)
                bsdf = gpuGeo->getTriMesh()->getBSDF();
            Shader *shader = m_renderer->registerShaderForResource(bsdf);
            if (shader && !shader->isComplete()) {
                m_renderer->unregisterShaderForResource(bsdf);
                shader = NULL;
            }
            gpuGeo->setShader(shader);
            m_geometry.push_back(std::make_pair(gpuGeo, identityTrafo));

            if (shader && !(shader->getFlags() & Shader::ETransparent))
                m_opaqueGeometry.push_back(std::make_pair(gpuGeo, identityTrafo));
        }
    }

    const ref_vector<Emitter> &emitters = scene->getEmitters();
    for (size_t i=0; i<emitters.size(); ++i)
        m_renderer->registerShaderForResource(emitters[i].get());

    if (scene->hasEnvironmentEmitter()) {
        Shader *shader = m_renderer->getShaderForResource(m_scene->getEnvironmentEmitter());
        m_backgroundDependencies = DependencyNode(shader);

        int id = 0;
        std::ostringstream oss;
        std::string evalName = m_backgroundDependencies.generateCode(oss, id);
        std::string marker = SUPPLEMENTAL_CODE_MARKER;
        std::string code = sh_background_frag;
        size_t insertionPos = code.find(marker);
        Assert(insertionPos != std::string::npos);
        code.replace(insertionPos, marker.length(), oss.str());

        ref<GPUProgram> prog = m_renderer->createGPUProgram("Background program");
        prog->setSource(GPUProgram::EVertexProgram, sh_background_vert);
        prog->setSource(GPUProgram::EFragmentProgram, code.c_str());
        prog->define("BACKGROUND_EVAL_NAME", evalName + "_background");

        if (scene->getSensor()->getType() & Sensor::EOrthographicCamera)
            prog->define("DIRECTIONAL_CAMERA");

        prog->init();

        id = 0;
        m_backgroundDependencies.resolve(prog, id);
        m_backgroundProgram = prog;
        m_backgroundParam_camPosition = prog->getParameterID("camPosition", false);
        m_backgroundParam_camDirection = prog->getParameterID("camDirection", false);
        m_backgroundParam_clipToWorld = prog->getParameterID("clipToWorld", false);
        m_backgroundParam_emitterScale = prog->getParameterID("emitterScale", false);
    }

    std::vector<size_t> geometryPermutation(m_geometry.size()),
                        opaqueGeometryPermutation(m_opaqueGeometry.size());

    for (size_t i=0; i<m_geometry.size(); ++i)
        geometryPermutation[i] = i;
    for (size_t i=0; i<m_opaqueGeometry.size(); ++i)
        opaqueGeometryPermutation[i] = i;

    /// Sort using the MaterialOrder to reduce material changes/pipeline flushes
    std::sort(geometryPermutation.begin(),
              geometryPermutation.end(), MaterialOrder(m_geometry));
    std::sort(opaqueGeometryPermutation.begin(),
              opaqueGeometryPermutation.end(), MaterialOrder(m_opaqueGeometry));

    if (!m_animatedGeometry.empty()) {
        std::vector<size_t> geometryPermutationInv(m_geometry.size()),
                            opaqueGeometryPermutationInv(m_opaqueGeometry.size());

        for (size_t i=0; i<m_geometry.size(); ++i)
            geometryPermutationInv[geometryPermutation[i]] = i;
        for (size_t i=0; i<m_opaqueGeometry.size(); ++i)
            opaqueGeometryPermutationInv[opaqueGeometryPermutation[i]] = i;

        for (size_t i=0; i<m_animatedGeometry.size(); ++i) {
            AnimatedGeometryRecord &agRec = m_animatedGeometry[i];
            if (agRec.geometryIndex >= 0)
                agRec.geometryIndex = geometryPermutationInv[agRec.geometryIndex];
            if (agRec.opaqueGeometryIndex >= 0)
                agRec.opaqueGeometryIndex = opaqueGeometryPermutationInv[agRec.opaqueGeometryIndex];
        }
    }

    permute_inplace(&m_geometry[0], geometryPermutation);
    permute_inplace(&m_opaqueGeometry[0], opaqueGeometryPermutation);
}

void VPLShaderManager::setVPL(const VPL &vpl) {
    const size_t sampleCount = 250;

    /* Estimate good near and far plane locations by tracing some rays */
    m_nearClip =  std::numeric_limits<Float>::infinity();
    m_farClip  = -std::numeric_limits<Float>::infinity();

    /* Update animations */
    for (size_t i=0; i<m_animatedGeometry.size(); ++i) {
        const AnimatedGeometryRecord &agRec = m_animatedGeometry[i];
        const Matrix4x4 &matrix = agRec.trafo->eval(vpl.its.time).getMatrix();

        if (agRec.geometryIndex >= 0)
            m_geometry[agRec.geometryIndex].second = matrix;

        if (agRec.opaqueGeometryIndex >= 0)
            m_opaqueGeometry[agRec.opaqueGeometryIndex].second = matrix;
    }

    if (vpl.type != EDirectionalEmitterVPL) {
        /* Trace a few rays from the VPL to estimate a suitable depth range */
        for (size_t i=0; i<sampleCount; ++i) {
            Point2 sample = sample02(i);

            Ray ray;
            if (vpl.type == ESurfaceVPL || (vpl.type == EPointEmitterVPL &&
                        vpl.emitter->getType() & Emitter::EOnSurface)) {
                ray = Ray(vpl.its.p, vpl.its.shFrame.toWorld(
                            warp::squareToCosineHemisphere(sample)), 0);

                #if defined(MTS_VPL_USE_PARABOLOID_MAPS)
                    m_shadowMapType = ShadowMapGenerator::EParaboloid;
                #elif defined(MTS_VPL_USE_SINGLE_PASS)
                    m_shadowMapType = ShadowMapGenerator::EHemicubeSinglePass;
                #else
                    m_shadowMapType = ShadowMapGenerator::EHemicube;
                #endif
            } else if (vpl.type == EPointEmitterVPL) {
                ray = Ray(vpl.its.p, warp::squareToUniformSphere(sample), 0);
                #if defined(MTS_VPL_USE_SINGLE_PASS)
                    m_shadowMapType = ShadowMapGenerator::ECubeSinglePass;
                #else
                    m_shadowMapType = ShadowMapGenerator::ECube;
                #endif
            } else {
                Log(EError, "Unsupported VPL type!");
            }

            ConstShapePtr shape; Normal n; Point2 uv; /* unused */
            Float t, accum = 0;
            bool hit = false;

            for (int it=0; it<5; ++it) {
                if (!m_scene->rayIntersect(ray, t, shape, n, uv))
                    break;

                accum += t;
                if (!(shape->getBSDF()->getType() & BSDF::ETransmission)) {
                    hit = true;
                    break;
                }

                ray.o = ray(t);
            }

            if (hit && accum > 0) {
                m_nearClip = std::min(m_nearClip, accum);
                m_farClip = std::max(m_farClip, accum);
            }
        }

        if (m_nearClip >= m_farClip) {
            BSphere bsphere(m_scene->getKDTree()->getAABB().getBSphere());
            Float minDist = 0;

            if ((vpl.type == ESurfaceVPL || vpl.type == EPointEmitterVPL) &&
                    !bsphere.contains(vpl.its.p))
                minDist = (bsphere.center - vpl.its.p).length() - bsphere.radius;

            /* If everything fails, use really conservative bounds */
            m_farClip = minDist + bsphere.radius * 2.25f;
            m_nearClip = std::max(minDist - 0.25f * bsphere.radius, m_farClip * 1e-5f);
        } else {
            m_nearClip = std::max(m_nearClip / 1.5f, m_farClip * 1e-5f);
            m_farClip *= 1.5f;
        }

        m_shadowMapTransform =
            (Transform::translate(Vector(vpl.its.p)) *
            Transform::fromFrame(vpl.its.shFrame)).inverse();
    } else {
        m_shadowMapType = ShadowMapGenerator::EDirectional;
        m_shadowMapTransform = m_shadowGen->directionalFindGoodFrame(
            m_scene->getKDTree()->getAABB(), vpl.its.shFrame.n);
    }

    bool is2D =
        m_shadowMapType == ShadowMapGenerator::EParaboloid ||
        m_shadowMapType == ShadowMapGenerator::EDirectional;

    if (is2D) {
        if (!m_shadowMap2D || m_shadowMap2D->getSize().x != m_shadowMapResolution) {
            if (m_shadowMap2D)
                m_shadowMap2D->cleanup();
            m_shadowMap2D = m_shadowGen->allocate(m_renderer,
                m_shadowMapType, m_shadowMapResolution);
        }
        m_shadowMap = m_shadowMap2D;
    } else {
        if (!m_shadowMapCube || m_shadowMapCube->getSize().x != m_shadowMapResolution) {
            if (m_shadowMapCube)
                m_shadowMapCube->cleanup();
            m_shadowMapCube = m_shadowGen->allocate(m_renderer,
                m_shadowMapType, m_shadowMapResolution);
        }
        m_shadowMap = m_shadowMapCube;
    }

    Float sample = sampleTEAFloat(m_vplIndex++, 0x12345);
    m_shadowGen->render(m_renderer, m_shadowMap, m_shadowMapType,
            m_shadowMapTransform, m_nearClip, m_farClip,
            sample > 0.3 ? m_opaqueGeometry : m_geometry);

    /* Convert between the Mitsuba and OpenGL matrix conventions */
    m_shadowMapTransform = Transform::scale(
        Vector(-1.0f, 1.0f, -1.0f)) * m_shadowMapTransform;
}

void VPLShaderManager::bind(const VPL &vpl, const BSDF *bsdf, const Sensor *sensor,
        const Emitter *emitter, const Matrix4x4 &instanceTransform, bool faceNormals) {
    Shader *bsdfShader = m_renderer->getShaderForResource(bsdf);
    Shader *vplShader = (vpl.type == EPointEmitterVPL || vpl.type == EDirectionalEmitterVPL)
        ? m_renderer->getShaderForResource(vpl.emitter)
        : m_renderer->getShaderForResource(vpl.its.getBSDF());
    Shader *emitterShader = (emitter == NULL) ? NULL
        : m_renderer->getShaderForResource(emitter);

    /* Find situations in which we won't be able to render the
       object properly ...  Case 1: one of the shaders is missing */
    bool unsupported = bsdfShader == NULL || vplShader == NULL ||
        (emitter != NULL && emitterShader == NULL);

    /* The material uses face normals (which have to be generated on the fly), and:
         case 2: anisotropy is active as well
         case 3: the graphics card does not support geometry shaders
    */
    if (faceNormals)
        unsupported |= (bsdf->getType() & BSDF::EAnisotropic)
            || !m_renderer->getCapabilities()->isSupported(
                RendererCapabilities::EGeometryShaders);

    if (unsupported) {
        std::map<std::string, VPLConfiguration>::iterator it
            = m_configurations.find("unsupported");
        m_targetConfiguration = VPLConfiguration();
        if (it != m_configurations.end()) {
            /* A program for this configuration has been created previously */
            m_currentProgram = (*it).second;
        } else {
            m_currentProgram = VPLConfiguration();
            ref<GPUProgram> prog = m_renderer->createGPUProgram("Unsupported material program");
            prog->setSource(GPUProgram::EVertexProgram, sh_unsupported_vert);
            prog->setSource(GPUProgram::EFragmentProgram, sh_unsupported_frag);
            prog->init();
            prog->incRef();
            m_currentProgram.program = prog;

            m_currentProgram.param_instanceTransform = prog->getParameterID("instanceTransform", false);
            m_configurations["unsupported"] = m_currentProgram;
        }
        GPUProgram *prog = m_currentProgram.program;
        prog->bind();
        prog->setParameter(m_currentProgram.param_instanceTransform, instanceTransform);
        return;
    }

    m_targetConfiguration = VPLConfiguration(vplShader,
        bsdfShader, emitterShader, faceNormals);

    m_alpha = bsdfShader->getAlpha();

    /* Generate a fingerprint of this shader chain, and check if it is known */
    std::string fingerprint = m_targetConfiguration.toString();

    if (m_diffuseSources)
        fingerprint += ", ds";
    if (m_diffuseReceivers)
        fingerprint += ", dr";

    std::map<std::string, VPLConfiguration>::iterator it
        = m_configurations.find(fingerprint);

    bool directionalCamera = false;
    if (sensor->getType() & Sensor::EOrthographicCamera) {
        directionalCamera = true;
    } else if (!(sensor->getType() & Sensor::EPerspectiveCamera)) {
        /* Sensor is neither a perspective nor orthographic camera */
        Log(EError, "Unsupported sensor type!");
    }

    if (it != m_configurations.end()) {
        /* A program for this configuration has been created previously */
        m_currentProgram = (*it).second;
    } else {
        std::ostringstream oss;
        std::string vplEvalName, bsdfEvalName, emitterEvalName;
        m_targetConfiguration.generateCode(oss, vplEvalName, bsdfEvalName, emitterEvalName);
        m_currentProgram = m_targetConfiguration;

        ref<GPUProgram> prog = m_renderer->createGPUProgram(fingerprint);
        prog->setSource(GPUProgram::EVertexProgram, sh_render_vert);

        if (faceNormals) {
            prog->setSource(GPUProgram::EGeometryProgram, sh_render_geom);
            prog->setInputGeometryType(GPUProgram::ETriangles);
            prog->setOutputGeometryType(GPUProgram::ETriangleStrips);
            prog->setMaxVertices(3);
        }

        std::string code(sh_render_frag);
        std::string marker(SUPPLEMENTAL_CODE_MARKER);
        size_t insertionPos = code.find(marker);
        Assert(insertionPos != std::string::npos);
        code.replace(insertionPos, marker.length(), oss.str());

        prog->setSource(GPUProgram::EFragmentProgram, code);

        if (bsdf->getType() & BSDF::EAnisotropic)
            prog->define("ANISOTROPIC");

        if (fingerprint.find("VertexColor") != std::string::npos)
            prog->define("VERTEX_COLORS");

        if (m_shadowMapType == ShadowMapGenerator::EDirectional)
            prog->define("DIRECTIONAL_VPL");
        else if (m_shadowMapType == ShadowMapGenerator::EParaboloid)
            prog->define("PARABOLOIDAL_VPL");
        else
            prog->define("CUBEMAP_VPL");

        if (directionalCamera)
            prog->define("DIRECTIONAL_CAMERA");

        if (faceNormals)
            prog->define("FACE_NORMALS");

        if (vpl.type == EDirectionalEmitterVPL || vpl.type == EPointEmitterVPL) {
            prog->define("EMITTER_VPL");
            prog->define("VPL_EVAL_NAME", vplEvalName + "_dir");
        } else {
            if (m_diffuseSources)
                prog->define("VPL_EVAL_NAME", vplEvalName + "_diffuse");
            else
                prog->define("VPL_EVAL_NAME", vplEvalName);
        }

        if (vpl.type == ESurfaceVPL || (vpl.type == EPointEmitterVPL &&
                vpl.emitter->getType() & Emitter::EOnSurface))
            prog->define("VPL_ON_SURFACE");

        if (emitterShader) {
            prog->define("EMITTER_AREA_EVAL_NAME", emitterEvalName + "_area");
            prog->define("EMITTER_DIR_EVAL_NAME", emitterEvalName + "_dir");
        }

        if (m_diffuseReceivers)
            prog->define("BSDF_EVAL_NAME", bsdfEvalName + "_diffuse");
        else
            prog->define("BSDF_EVAL_NAME", bsdfEvalName);

        prog->init();
        prog->incRef();
        m_currentProgram.program = prog;

        m_currentProgram.param_camPosition       = prog->getParameterID("camPosition", false);
        m_currentProgram.param_camDirection      = prog->getParameterID("camDirection", false);
        m_currentProgram.param_vplPosition       = prog->getParameterID("vplPosition", false);
        m_currentProgram.param_vplDirection      = prog->getParameterID("vplDirection", false);
        m_currentProgram.param_vplPower          = prog->getParameterID("vplPower", false);
        m_currentProgram.param_vplTransform      = prog->getParameterID("vplTransform", false);
        m_currentProgram.param_vplFrame          = prog->getParameterID("vplFrame", false);
        m_currentProgram.param_vplUV             = prog->getParameterID("vplUV", false);
        m_currentProgram.param_vplWi             = prog->getParameterID("vplWi", false);
        m_currentProgram.param_minDistSqr        = prog->getParameterID("minDistSqr", false);
        m_currentProgram.param_emitterScale      = prog->getParameterID("emitterScale", false);
        m_currentProgram.param_depthRange        = prog->getParameterID("depthRange", false);
        m_currentProgram.param_instanceTransform = prog->getParameterID("instanceTransform", false);
        m_currentProgram.param_shadowMap         = prog->getParameterID("shadowMap", false);
        m_currentProgram.resolve(prog);
        m_configurations[fingerprint] = m_currentProgram;

        statsMaxResidentShaders.recordMaximum(m_configurations.size() + m_shadowGen->getShaderCount());
    }

    GPUProgram *prog = m_currentProgram.program;
    prog->bind();

    Float minDist = m_nearClip + (m_farClip - m_nearClip) * m_clamping;
    prog->setParameter(m_currentProgram.param_instanceTransform, instanceTransform);
    prog->setParameter(m_currentProgram.param_vplTransform, m_shadowMapTransform);
    prog->setParameter(m_currentProgram.param_depthRange, Vector2(m_nearClip, m_farClip));
    prog->setParameter(m_currentProgram.param_vplPower, vpl.P);
    prog->setParameter(m_currentProgram.param_vplFrame,
        Matrix3x3(vpl.its.shFrame.s, vpl.its.shFrame.t, vpl.its.shFrame.n));
    prog->setParameter(m_currentProgram.param_vplUV, vpl.its.uv);
    prog->setParameter(m_currentProgram.param_vplWi, vpl.its.wi);
    prog->setParameter(m_currentProgram.param_shadowMap, m_shadowMap);
    prog->setParameter(m_currentProgram.param_minDistSqr, minDist*minDist);
    prog->setParameter(m_currentProgram.param_emitterScale, vpl.emitterScale);

    if (directionalCamera) {
        Vector d = sensor->getWorldTransform()->eval(0)(Vector(0.0f, 0.0f, 1.0f));
        prog->setParameter(m_currentProgram.param_camDirection, d);
    } else {
        Point p = sensor->getWorldTransform()->eval(0).transformAffine(Point(0.0f));
        prog->setParameter(m_currentProgram.param_camPosition, p);
    }

    if (m_shadowMapType == ShadowMapGenerator::EDirectional)
        prog->setParameter(m_currentProgram.param_vplDirection, vpl.its.shFrame.n);
    else
        prog->setParameter(m_currentProgram.param_vplPosition, vpl.its.p);

    m_targetConfiguration.bind(m_currentProgram, 1);
}

void VPLShaderManager::unbind() {
    if (!m_currentProgram.program)
        return;
    m_shadowMap->unbind();
    m_currentProgram.program->unbind();
    m_currentProgram.program = NULL;
    m_targetConfiguration.unbind();
}

void VPLShaderManager::drawAllGeometryForVPL(const VPL &vpl, const Sensor *sensor) {
    const Emitter *currentEmitter = NULL;
    const BSDF *currentBSDF = NULL;
    bool currentHasNormals = false;

    m_renderer->setDepthTest(true);

    Matrix4x4 currentObjTrafo;
    currentObjTrafo.setIdentity();
    m_shadowMap->bind(0);

    m_renderer->beginDrawingMeshes();

    size_t nTriangles = 0;
    for (std::vector<Renderer::TransformedGPUGeometry>::const_iterator it = m_geometry.begin();
            it != m_geometry.end(); ++it) {
        const GPUGeometry *geo = (*it).first;
        const Matrix4x4 &trafo = (*it).second;
        const BSDF *bsdf = geo->getTriMesh()->getBSDF();
        const Emitter *emitter = geo->getTriMesh()->getEmitter();
        bool hasNormals = !geo->getTriMesh()->hasVertexNormals();

        nTriangles += geo->getTriMesh()->getTriangleCount();

        if (emitter != currentEmitter || bsdf != currentBSDF || hasNormals != currentHasNormals) {
            currentBSDF = bsdf; currentEmitter = emitter; currentHasNormals = hasNormals;

            if (m_currentProgram.program) {
                m_currentProgram.program->unbind();
                m_currentProgram.program = NULL;
                m_targetConfiguration.unbind();
            }

            bind(vpl, bsdf, sensor, emitter, trafo, hasNormals);
            currentObjTrafo = trafo;
        } else if (trafo != currentObjTrafo) {
            if (m_currentProgram.program)
                m_currentProgram.program->setParameter(
                    m_currentProgram.param_instanceTransform, trafo);
            currentObjTrafo = trafo;
        }

        if (m_alpha != 1.0f && sampleTEAFloat((uint32_t)
                (it - m_geometry.begin()), m_vplIndex, 8) > m_alpha)
            continue;

        m_renderer->drawMesh(geo);
    }

    statsMaxTriangles.recordMaximum(nTriangles);

    m_renderer->endDrawingMeshes();
    unbind();
    m_renderer->checkError();
}

void VPLShaderManager::drawBackground(const Sensor *sensor,
        const Transform &projectionTransform, Float scaleFactor) {
    if (m_backgroundProgram == NULL)
        return;

    const Transform &trafo = sensor->getWorldTransform()->eval(0);

    Transform clipToWorld = trafo
        * Transform::scale(Vector(-1, 1, -1)) * projectionTransform.inverse();

    GPUProgram *prog = m_backgroundProgram;
    int tuOffset = 0;
    prog->bind();
    m_backgroundDependencies.bind(prog, m_backgroundDependencies, tuOffset);

    if (sensor->getType() & Sensor::EOrthographicCamera) {
        Vector d = trafo(Vector(0.0f, 0.0f, 1.0f));
        prog->setParameter(m_backgroundParam_camDirection, d);
    } else {
        Point p = trafo(Point(0.0f));
        prog->setParameter(m_backgroundParam_camPosition, p);
    }

    prog->setParameter(m_backgroundParam_emitterScale, scaleFactor);
    prog->setParameter(m_backgroundParam_clipToWorld, clipToWorld);
    m_renderer->blitQuad(false);
    prog->unbind();
    m_backgroundDependencies.unbind();
}

MTS_IMPLEMENT_CLASS(VPLShaderManager, false, Object)
MTS_NAMESPACE_END

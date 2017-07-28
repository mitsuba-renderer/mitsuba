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

#include <mitsuba/hw/shadow.h>
#include <mitsuba/hw/gputexture.h>
#include "data/shaders.h"

MTS_NAMESPACE_BEGIN

ShadowMapGenerator::ShadowMapGenerator(Renderer *renderer) {
    /* Standard directional (i.e. orthographic) shadow map generator */
    GPUProgram *prog = renderer->createGPUProgram("Directional shadow map generator");
    prog->setSource(GPUProgram::EVertexProgram, sh_directional_vert);
    prog->setSource(GPUProgram::EFragmentProgram, sh_directional_frag);
    m_program[EDirectional] = prog;

    /* Cube (i.e. omnidirectional) shadow map generator */
    prog = renderer->createGPUProgram("Cube/Hemicube shadow map generator (5/6 pass version)");
    prog->setSource(GPUProgram::EVertexProgram, sh_cube_6pass_vert);
    prog->setSource(GPUProgram::EFragmentProgram, sh_cube_6pass_frag);
    m_program[ECube] = prog;

    if (renderer->getCapabilities()->isSupported(RendererCapabilities::EGeometryShaders)) {
        /* Paraboloid shadow map generator based on "Fast Non-Linear Projections
           using Graphics Hardware" by Jean-Dominique Gascuel, Nicolas Holzschuch,
           Gabrier Fournier, and Bernhard Peroche (I3D, 2008) */
        prog = renderer->createGPUProgram("Paraboloid shadow map generator");
        prog->setSource(GPUProgram::EVertexProgram, sh_paraboloid_vert);
        prog->setSource(GPUProgram::EGeometryProgram, sh_paraboloid_geom);
        prog->setSource(GPUProgram::EFragmentProgram, sh_paraboloid_frag);
        prog->setInputGeometryType(GPUProgram::ETriangles);
        prog->setOutputGeometryType(GPUProgram::ETriangleStrips);
        prog->setMaxVertices(4);
        m_program[EParaboloid] = prog;

        /* Cube shadow map generator (single pass via a geometry shader) */
        prog = renderer->createGPUProgram("Cube shadow map generator (1 pass version)");
        prog->setSource(GPUProgram::EVertexProgram, sh_cube_1pass_vert);
        prog->setSource(GPUProgram::EGeometryProgram, sh_cube_1pass_geom);
        prog->setSource(GPUProgram::EFragmentProgram, sh_cube_1pass_frag);
        prog->setInputGeometryType(GPUProgram::ETriangles);
        prog->setOutputGeometryType(GPUProgram::ETriangleStrips);
        prog->setMaxVertices(6*3); /* One triangle per cube map face */
        m_program[ECubeSinglePass] = prog;

        /* Hemicube shadow map generator (single pass via a geometry shader) */
        prog = renderer->createGPUProgram("Hemicube shadow map generator (1 pass version)");
        prog->setSource(GPUProgram::EVertexProgram, sh_hemicube_1pass_vert);
        prog->setSource(GPUProgram::EGeometryProgram, sh_hemicube_1pass_geom);
        prog->setSource(GPUProgram::EFragmentProgram, sh_hemicube_1pass_frag);
        prog->setInputGeometryType(GPUProgram::ETriangles);
        prog->setOutputGeometryType(GPUProgram::ETriangleStrips);
        prog->setMaxVertices(5*3); /* One triangle per hemicube map face */
        m_program[EHemicubeSinglePass] = prog;

    }
    m_program[EHemicube] = m_program[ECube];

    /* By default, assume that cube depth maps are supported by the driver */
    m_cubeDepthMapsSupported = true;
}

size_t ShadowMapGenerator::getShaderCount() const {
    size_t result = 0;
    for (int i=0; i<ETypeCount; ++i) {
        if (m_program[i].get())
            ++result;
    }
    return result;
}

void ShadowMapGenerator::init() {
    for (int i=0; i<ETypeCount; ++i) {
        if (!m_program[i] || i == EHemicube)
            continue;
        m_program[i]->init();

        switch (i) {
            case ECube:
                m_cubeTransform = m_program[i]->getParameterID("transform");
                m_cubeProjDir = m_program[i]->getParameterID("projDir", false);
                break;

            case ECubeSinglePass:
                for (int j=0; j<6; ++j) {
                    m_cubeSinglePassTransform[j] = m_program[i]->getParameterID(formatString("transform[%i]", j));
                    m_cubeSinglePassProjDir[j] = m_program[i]->getParameterID(formatString("projDir[%i]", j));
                }
                break;

            case EHemicubeSinglePass:
                for (int j=0; j<5; ++j) {
                    m_hemicubeSinglePassTransform[j] = m_program[i]->getParameterID(formatString("transform[%i]", j));
                    m_hemicubeSinglePassProjDir[j] = m_program[i]->getParameterID(formatString("projDir[%i]", j));
                }
                break;

            case EParaboloid:
                m_paraboloidMinDepth = m_program[i]->getParameterID("minDepth");
                m_paraboloidInvDepthRange = m_program[i]->getParameterID("invDepthRange");
                break;

            default:
                break;
        }
    }
}

void ShadowMapGenerator::cleanup() {
    for (int i=0; i<ETypeCount; ++i) {
        if (m_program[i] && i != EHemicube)
            m_program[i]->cleanup();
    }
}

ref<GPUTexture> ShadowMapGenerator::allocate(Renderer *renderer,
        EShadowMapType type, int res) {
    ref<GPUTexture> result = renderer->createGPUTexture("Shadow map");

    result->setSize(Point3i(res, res, 1));
    result->setComponentFormat(GPUTexture::EFloat32);
    result->setPixelFormat(GPUTexture::EDepth);
    result->setFrameBufferType(GPUTexture::EDepthBuffer);
    result->setDepthMode(GPUTexture::ENormal);
    result->setMipMapped(false);

retry:
    if (type == ECube || type == ECubeSinglePass ||
            type == EHemicube || type == EHemicubeSinglePass) {
        result->setWrapType(GPUTexture::EClampToEdge);
        result->setType(GPUTexture::ETextureCubeMap);

        if (!m_cubeDepthMapsSupported) {
            /* A few graphics cards can't deal with depth cube maps although
               they claim to provide all of the necessary extensions. Use
               a rather wasteful workaround in this case: create a RGB cube map
               (luminance is not widely supported ..) with a single depth buffer */
            result->setFrameBufferType(GPUTexture::EColorBuffer);
            result->setPixelFormat(GPUTexture::ERGB);
        }
    } else {
        result->setType(GPUTexture::ETexture2D);
        result->setBorderColor(Color3(0.0f));
        result->setWrapType(GPUTexture::EClampToBorder);
    }

    try {
        result->init();
    } catch (const std::exception &ex) {
        if (result->getType() == GPUTexture::ETextureCubeMap && m_cubeDepthMapsSupported) {
            m_cubeDepthMapsSupported = false;
            Log(EWarn, "Graphics driver refused to create a cube depth map! [%s] "
                "Attempting to use a workaround ...", ex.what());
            result->cleanup();

            /* Generate new cube map shader with workaround */
            ref<GPUProgram> prog = renderer->createGPUProgram("Alternate Cube/Hemicube "
                "shadow map generator (5/6 pass version)");
            prog->define("DEPTH_CUBEMAPS_UNSUPPORTED");
            prog->setSource(GPUProgram::EVertexProgram, sh_cube_6pass_vert);
            prog->setSource(GPUProgram::EFragmentProgram, sh_cube_6pass_frag);
            prog->init();
            m_program[EHemicube] = m_program[ECube] = prog;

            goto retry;
        } else {
            throw;
        }
    }

    return result;
}

void ShadowMapGenerator::render(Renderer *renderer, GPUTexture *shadowMap,
        EShadowMapType type, const Transform &trafo, Float minDepth, Float maxDepth,
        const std::vector<Renderer::TransformedGPUGeometry> &geo) {
    GPUProgram *prog = m_program[type];

    if (!prog)
        Log(EError, "Cannot render shadow map (the "
            "graphics card has insufficient capabilities)");

    float invDepthRange = (float) (1.0f / (maxDepth - minDepth));
    shadowMap->activateTarget();
    renderer->setDepthTest(true);
    prog->bind();

    switch (type) {
        case EDirectional: {
                Matrix4x4 identity;
                identity.setIdentity();
                shadowMap->clear();
                renderer->setCamera(identity, trafo.getMatrix());
                renderer->drawAll(geo);
            }
            break;

        case EParaboloid: {
                Matrix4x4 identity;
                identity.setIdentity();
                shadowMap->clear();
                renderer->setCamera(identity, trafo.getMatrix());
                prog->setParameter(m_paraboloidMinDepth, minDepth);
                prog->setParameter(m_paraboloidInvDepthRange, invDepthRange);
                renderer->drawAll(geo);
            }
            break;

        case EHemicube:
        case ECube : {
                renderer->clearTransforms();
                Transform projTrafo = Transform::glPerspective(90.0f, minDepth, maxDepth);
                renderer->setClearColor(Color3(0));

                /* Render each cube face separately */
                for (int i=0; i<6; ++i) {
                    Transform viewTrafo;

                    Point o(0.0f);
                    switch (i) {
                        /* These come from the OpenGL spec, see section 3.8.10 ("Texture application") */
                        case 0: viewTrafo = Transform::lookAt(o, Point( 1,  0,  0), Vector(0, -1,  0)); break;
                        case 1: viewTrafo = Transform::lookAt(o, Point(-1,  0,  0), Vector(0, -1,  0)); break;
                        case 2: viewTrafo = Transform::lookAt(o, Point( 0,  1,  0), Vector(0,  0,  1)); break;
                        case 3: viewTrafo = Transform::lookAt(o, Point( 0, -1,  0), Vector(0,  0, -1)); break;
                        case 4: viewTrafo = Transform::lookAt(o, Point( 0,  0,  1), Vector(0, -1,  0)); break;
                        case 5: viewTrafo = Transform::lookAt(o, Point( 0,  0, -1), Vector(0, -1,  0)); break;
                    }

                    /* Need a left-handed view matrix to render the faces (see the ARB_texture_cube_map spec) */
                    viewTrafo = Transform::scale(Vector(1, -1, 1)) * viewTrafo.inverse() * trafo;
                    prog->setParameter(m_cubeTransform, projTrafo * viewTrafo);
                    prog->setParameter(m_cubeProjDir,
                        (-viewTrafo.getMatrix().row(2) - Vector4(0, 0, 0, minDepth)) * invDepthRange);
                    shadowMap->activateSide(i);
                    if (type == ECube || i != 4) {
                        shadowMap->clear();
                        renderer->drawAll(geo);
                    } else {
                        renderer->setClearDepth(0);
                        shadowMap->clear();
                        renderer->setClearDepth(1);
                    }
                }
            }
            break;

        case ECubeSinglePass: {
                renderer->clearTransforms();
                shadowMap->activateSide(-1);
                shadowMap->clear();

                Transform projTrafo = Transform::glPerspective(90.0f, minDepth, maxDepth);

                /* Set the parameters for each cube map face */
                for (int i=0; i<6; ++i) {
                    Transform viewTrafo;
                    Point o(0.0f);
                    switch (i) {
                        /* These come from the OpenGL spec, see section 3.8.10 ("Texture application") */
                        case 0: viewTrafo = Transform::lookAt(o, Point( 1,  0,  0), Vector(0, -1,  0)); break;
                        case 1: viewTrafo = Transform::lookAt(o, Point(-1,  0,  0), Vector(0, -1,  0)); break;
                        case 2: viewTrafo = Transform::lookAt(o, Point( 0,  1,  0), Vector(0,  0,  1)); break;
                        case 3: viewTrafo = Transform::lookAt(o, Point( 0, -1,  0), Vector(0,  0, -1)); break;
                        case 4: viewTrafo = Transform::lookAt(o, Point( 0,  0,  1), Vector(0, -1,  0)); break;
                        case 5: viewTrafo = Transform::lookAt(o, Point( 0,  0, -1), Vector(0, -1,  0)); break;
                    }

                    /* Need a left-handed view matrix to render the faces (see the ARB_texture_cube_map spec) */
                    viewTrafo = Transform::scale(Vector(1, -1, 1)) * viewTrafo.inverse() * trafo;
                    prog->setParameter(m_cubeSinglePassTransform[i], projTrafo * viewTrafo);
                    prog->setParameter(m_cubeSinglePassProjDir[i],
                        (-viewTrafo.getMatrix().row(2) - Vector4(0, 0, 0, minDepth)) * invDepthRange);
                }

                /* Send the geometry to the GPU once */
                renderer->drawAll(geo);
            }
            break;

        case EHemicubeSinglePass: {
                renderer->clearTransforms();
                for (int i=0; i<6; ++i) {
                    shadowMap->activateSide(i);

                    if (i != 4) {
                        shadowMap->clear();
                    } else {
                        renderer->setClearDepth(0);
                        shadowMap->clear();
                        renderer->setClearDepth(1);
                    }
                }
                shadowMap->activateSide(-1);

                Transform projTrafo = Transform::glPerspective(90.0f, minDepth, maxDepth);

                /* Set the parameters for each cube map face */
                for (int i=0; i<5; ++i) {
                    Transform viewTrafo;
                    Point o(0.0f);
                    switch (i) {
                        /* These come from the OpenGL spec, see section 3.8.10 ("Texture application") */
                        case 0: viewTrafo = Transform::lookAt(o, Point( 1,  0,  0), Vector(0, -1,  0)); break;
                        case 1: viewTrafo = Transform::lookAt(o, Point(-1,  0,  0), Vector(0, -1,  0)); break;
                        case 2: viewTrafo = Transform::lookAt(o, Point( 0,  1,  0), Vector(0,  0,  1)); break;
                        case 3: viewTrafo = Transform::lookAt(o, Point( 0, -1,  0), Vector(0,  0, -1)); break;
                        case 4: viewTrafo = Transform::lookAt(o, Point( 0,  0, -1), Vector(0, -1,  0)); break;
                    }

                    /* Need a left-handed view matrix to render the faces (see the ARB_texture_cube_map spec) */
                    viewTrafo = Transform::scale(Vector(1, -1, 1)) * viewTrafo.inverse() * trafo;
                    prog->setParameter(m_hemicubeSinglePassTransform[i], projTrafo * viewTrafo);
                    prog->setParameter(m_hemicubeSinglePassProjDir[i],
                        (-viewTrafo.getMatrix().row(2) - Vector4(0, 0, 0, minDepth)) * invDepthRange);
                }

                /* Send the geometry to the GPU once */
                renderer->drawAll(geo);
            }
            break;

        default:
            Log(EError, "Invalid shadow map type!");
    }

    prog->unbind();
    shadowMap->releaseTarget();
}

Transform ShadowMapGenerator::directionalFindGoodFrame(const AABB &aabb, const Vector &d) const {
    /* Start with an arbitrary frame */
    Frame frame(d);
    Vector samples[8], center(0.0f);
    for (int i=0; i<8; ++i) {
        Vector v = frame.toLocal(Vector(aabb.getCorner(i)));
        samples[i] = v;
        center += v;
    }
    center *= 1.0f / 8.0f;
    for (int i=0; i<8; ++i)
        samples[i] -= center;

    /* Heuristic: do a 3x3 principal components analysis to try
       to align the projection with the coordinate axis. This is to
       make best use of the shadow map resolution */
    Matrix2x2 M;
    M.setZero();
    for (int i=0; i<2; ++i) {
        for (int j=0; j<=i; ++j) {
            Float sum = 0.0f;
            for (int k=0; k<8; ++k)
                sum += samples[k][i] * samples[k][j];
            M(i, j) = sum;
        }
    }
    M(0, 1) = M(1, 0);

    Float eig[2];
    Matrix2x2 Q;
    M.symEig(Q, eig);

    Frame newFrame(
        frame.s*Q(1, 0) + frame.t*Q(1, 1),
        frame.s*Q(0, 0) + frame.t*Q(0, 1), d);

    Transform trafo = Transform::fromFrame(newFrame).inverse() *
        Transform::translate(-frame.s*center.x-frame.t*center.y);

    AABB aabb2;
    for (int i=0; i<8; ++i)
        aabb2.expandBy(trafo(aabb.getCorner(i)));

    Vector extents = aabb2.getExtents();
    aabb2.min -= extents * 1e-3f;
    aabb2.max += extents * 1e-3f;

    return Transform::glOrthographic(
         aabb2.min.x,  aabb2.max.x,
         aabb2.min.y,  aabb2.max.y,
        -aabb2.max.z, -aabb2.min.z) * trafo;
}


MTS_IMPLEMENT_CLASS(ShadowMapGenerator, false, Object)
MTS_NAMESPACE_END

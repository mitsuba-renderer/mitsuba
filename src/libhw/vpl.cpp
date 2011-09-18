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

#include <mitsuba/core/plugin.h>
#include <mitsuba/hw/vpl.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/hw/gputexture.h>
#include "../shapes/instance.h"

MTS_NAMESPACE_BEGIN

VPLShaderManager::VPLShaderManager(const Scene *scene, Renderer *renderer)
	 : m_scene(scene), m_renderer(renderer), m_clamping(0.1f),
 	   m_maxClipDist(std::numeric_limits<Float>::infinity()), m_initialized(false), 
	   m_shadowMapResolution(512), m_singlePass(false), 
	   m_diffuseSources(false), m_diffuseReceivers(false) {
}

VPLShaderManager::~VPLShaderManager() {
	if (m_initialized)
		cleanup();
}

void VPLShaderManager::init() {
    if (m_renderer->getCapabilities()->isSupported(RendererCapabilities::EGeometryShaders)) {
		m_shadowProgram = m_renderer->createGPUProgram("Shadow Program");
	
		m_shadowProgram->setSource(GPUProgram::EVertexProgram,
			"void main() {\n"
			"	gl_Position = gl_ModelViewMatrix * gl_Vertex;\n"
			"}"
		);

		m_shadowProgram->setSource(GPUProgram::EGeometryProgram,
			"#version 120\n"
			"#extension GL_EXT_geometry_shader4 : enable\n"
			"\n"
			"uniform mat4 cubeMapTransform[6];\n"
			"uniform vec4 depthVec[6];\n"
			"varying float depth;\n"
			"\n"
			"void main() {\n"
			"	depth = 0;\n" // avoid an (incorrect?) warning
			"   for (int side = 0; side < 6; side++) {\n"
			"	    gl_Layer = side;\n"
			"       for (int i = 0; i < gl_VerticesIn; i++) {\n"
			"           gl_Position = cubeMapTransform[side] * gl_PositionIn[i];\n"
			"           depth = dot(depthVec[side], gl_PositionIn[i]);\n"
			"           EmitVertex();\n"
			"       }\n"
			"       EndPrimitive();\n"
			"   }\n"
			"}\n"
		);

		m_shadowProgram->setSource(GPUProgram::EFragmentProgram,
			"#version 120\n"
			"varying float depth;\n"
			"void main() {\n"
			"   float dx = dFdx(depth), dy = dFdy(depth);"
			"   gl_FragDepth = depth + sqrt(dx*dx + dy*dy);"
			"}\n"
		);

		/* Six output triangles per input triangle */
		m_shadowProgram->setMaxVertices(18); 
		m_shadowProgram->init();

		for (int i=0; i<6; ++i)	{
			m_shadowProgramParam_cubeMapTransform[i] = 
				m_shadowProgram->getParameterID(formatString("cubeMapTransform[%i]", i));
			m_shadowProgramParam_depthVec[i] = 
				m_shadowProgram->getParameterID(formatString("depthVec[%i]", i));
		}
	}

	m_altShadowProgram = m_renderer->createGPUProgram("Alternative Shadow Program");
	m_altShadowProgram->setSource(GPUProgram::EVertexProgram,
		"uniform mat4 cubeMapTransform;\n"
		"uniform vec4 depthVec;\n"
		"varying float depth;\n"
		"void main() {\n"
		"	gl_Position = gl_ModelViewMatrix * (cubeMapTransform * gl_Vertex);\n"
		"	depth = dot(depthVec, gl_Vertex);\n"
		"}\n"
	);

	m_altShadowProgram->setSource(GPUProgram::EFragmentProgram,
		"#version 120\n"
		"varying float depth;\n"
		"void main() {\n"
		"	float dx = dFdx(depth), dy = dFdy(depth);"
		"   gl_FragDepth = depth + sqrt(dx*dx + dy*dy);"
		"}\n"
	);

	m_altShadowProgram->init();

	m_altShadowProgramParam_cubeMapTransform = 
		m_altShadowProgram->getParameterID("cubeMapTransform");
	m_altShadowProgramParam_depthVec = 
		m_altShadowProgram->getParameterID("depthVec");

	const std::vector<Shape *> shapes = m_scene->getShapes();
	const std::vector<Luminaire *> luminaires = m_scene->getLuminaires();
	
	for (size_t i=0; i<shapes.size(); ++i) {
		ref<TriMesh> triMesh = shapes[i]->createTriMesh();
		if (!triMesh) {
			std::string shapeClass = shapes[i]->getClass()->getName();
			if (shapeClass == "Instance") {
				const Instance *instance = static_cast<const Instance *>(shapes[i]);
				const std::vector<const Shape *> &subShapes = 
						instance->getShapeGroup()->getKDTree()->getShapes(); 
				for (size_t j=0; j<subShapes.size(); ++j) {
					triMesh = const_cast<Shape *>(subShapes[j])->createTriMesh();
					if (!triMesh)
						continue;
					GPUGeometry *gpuGeo = m_renderer->registerGeometry(triMesh);
					Shader *shader = triMesh->hasBSDF() ? 
						m_renderer->registerShaderForResource(triMesh->getBSDF()) : NULL;
					if (shader != NULL && !shader->isComplete()) {
						m_renderer->unregisterShaderForResource(triMesh->getBSDF());
					} else if (shader != NULL && shader->getFlags() & Shader::ETransparent) {
						m_transparentMeshes.push_back(std::make_pair(triMesh.get(), instance->getWorldTransform()));
						continue;
					}
					m_meshes.push_back(std::make_pair(triMesh.get(), instance->getWorldTransform()));
					if (gpuGeo)
						m_drawList.push_back(std::make_pair(gpuGeo, instance->getWorldTransform()));
				}
			}
			continue;
		}
		GPUGeometry *gpuGeo = m_renderer->registerGeometry(triMesh);
		Shader *shader = triMesh->hasBSDF() ? 
			m_renderer->registerShaderForResource(triMesh->getBSDF()) : NULL;
		if (shader != NULL && !shader->isComplete()) {
			m_renderer->unregisterShaderForResource(triMesh->getBSDF());
		} else if (shader != NULL && shader->getFlags() & Shader::ETransparent) {
			m_transparentMeshes.push_back(std::make_pair(triMesh.get(), Transform()));
			continue;
		}
		m_meshes.push_back(std::make_pair(triMesh.get(), Transform()));
		if (gpuGeo)
			m_drawList.push_back(std::make_pair(gpuGeo, Transform()));
	}

	for (size_t i=0; i<luminaires.size(); ++i)
		m_renderer->registerShaderForResource(luminaires[i]);

	if (m_scene->hasBackgroundLuminaire() && 
		m_renderer->getShaderForResource(m_scene->getBackgroundLuminaire()) != NULL) {
		Shader *shader = m_renderer->getShaderForResource(m_scene->getBackgroundLuminaire());
		m_backgroundDependencies = VPLDependencyNode(shader);
		int id = 0;
		std::ostringstream oss;
		std::string evalName = m_backgroundDependencies.recursiveGenerateCode(oss, id);

		m_backgroundProgram = m_renderer->createGPUProgram("Background program");
		m_backgroundProgram->setSource(GPUProgram::EVertexProgram,
			"uniform mat4 clipToWorld;\n"
			"varying vec3 d;\n"
			"void main() {\n"
			"	gl_Position = ftransform();\n"
			"	vec4 tmp = clipToWorld * (gl_ModelViewProjectionMatrix * gl_Vertex);\n"
			"   d = tmp.xyz/tmp.w;"
			"}\n"
		);

		oss << "varying vec3 d;" << endl
			<< "uniform vec3 camPos;" << endl
			<< "uniform float scale;" << endl
			<< "void main() {" << endl
			<< "  gl_FragColor.rgb = scale * " << evalName << "_background(normalize(d - camPos));" << endl
			<< "  gl_FragColor.a = 1.0;" << endl
			<< "}" << endl;

		m_backgroundProgram->setSource(GPUProgram::EFragmentProgram, oss.str());
		m_backgroundProgram->init();

		id = 0;
		m_backgroundDependencies.recursiveResolve(m_backgroundProgram, id);
	}

	m_initialized = true;
}

void VPLShaderManager::cleanup() {
	for (std::map<std::string, ProgramAndConfiguration>::iterator it = m_programs.begin();
		it != m_programs.end(); ++it) {
		(*it).second.program->cleanup();
		(*it).second.program->decRef();
	}

	if (m_shadowMap)
		m_shadowMap->cleanup();

	if (m_backgroundProgram) {
		m_backgroundProgram->cleanup();
		m_backgroundProgram = NULL;
	}

	const std::vector<Luminaire *> luminaires = m_scene->getLuminaires();

	for (size_t i=0; i<m_meshes.size(); ++i) {
		m_renderer->unregisterGeometry(m_meshes[i].first);
		m_renderer->unregisterShaderForResource(m_meshes[i].first->getBSDF());
	}

	m_meshes.clear();
	m_drawList.clear();

	for (size_t i=0; i<luminaires.size(); ++i)
		m_renderer->unregisterShaderForResource(luminaires[i]);
	m_initialized = false;
}


void VPLShaderManager::setVPL(const VPL &vpl) {
	Point p = vpl.its.p + vpl.its.shFrame.n * 0.01;
	Intersection its;

	/* Estimate good near and far plane locations by tracing some rays */
	Float nearClip =  std::numeric_limits<Float>::infinity(),
		  farClip  = -std::numeric_limits<Float>::infinity();
	Ray ray;
	ray.o = p;

	if (m_shadowMap == NULL || m_shadowMapResolution != m_shadowMap->getSize().x) {
		m_shadowMap = m_renderer->createGPUTexture("Shadow cube map", NULL);
		m_shadowMap->setSize(Point3i(m_shadowMapResolution, m_shadowMapResolution, 1));
		m_shadowMap->setFrameBufferType(GPUTexture::EDepthBuffer);
		m_shadowMap->setType(GPUTexture::ETextureCubeMap);
		m_shadowMap->setWrapType(GPUTexture::EClampToEdge);
		m_shadowMap->setFilterType(GPUTexture::ENearest);
		m_shadowMap->setDepthMode(GPUTexture::ENormal);
		m_shadowMap->init();
	}

	const int sampleCount = 200;
	const Float invSampleCount = 1.0f/sampleCount;

	for (int i=1; i<=sampleCount; ++i) {
		Vector dir;
		Point2 seed(i*invSampleCount, radicalInverse(2, i)); // Hammersley seq.
		if (vpl.type == ESurfaceVPL || vpl.luminaire->getType() & Luminaire::EOnSurface)
			dir = vpl.its.shFrame.toWorld(squareToHemispherePSA(seed));
		else
			dir = squareToSphere(seed);
		ray.setDirection(dir);
		if (m_scene->rayIntersect(ray, its)) {
			nearClip = std::min(nearClip, its.t);
			farClip = std::max(farClip, its.t);
		}
	}

	m_minDist = nearClip + (farClip - nearClip) * m_clamping;

	nearClip = std::min(nearClip, (Float) 0.001f);
	farClip = std::min(farClip * 1.5f, m_maxClipDist);

	if (farClip < 0 || nearClip >= farClip) {
		/* Unable to find any surface - just default values based on the scene size */
		nearClip = 1e-3f * m_scene->getBSphere().radius;
		farClip = 2 * m_scene->getBSphere().radius;
		m_minDist = 0;
	}
	farClip = std::min(farClip, 5.0f*m_scene->getBSphere().radius);

	m_nearClip = nearClip;
	m_invClipRange = 1/(farClip-nearClip);
	Transform lightViewTrafo, lightProjTrafo = Transform::glPerspective(90.0f, nearClip, farClip);
	Matrix4x4 identity;
	identity.setIdentity();
	m_renderer->setCamera(identity, identity);

	m_shadowMap->activateTarget();
	if (m_singlePass && m_shadowProgram != NULL) {
		/* "Fancy": render the whole cube map in a single pass using 
		   a geometry program. On anything but brand-new hardware, this 
		   is actually slower. */

		m_shadowMap->activateSide(-1);
		m_shadowMap->clear();
		m_shadowProgram->bind();
		try {
			for (int i=0; i<6; ++i) {
				switch (i) {
					case 0: lightViewTrafo = Transform::lookAt(p, p + Vector(1, 0, 0), Vector(0, 1, 0)).inverse(); break;
					case 1: lightViewTrafo = Transform::lookAt(p, p + Vector(-1, 0, 0), Vector(0, 1, 0)).inverse(); break;
					case 2: lightViewTrafo = Transform::lookAt(p, p + Vector(0, 1, 0), Vector(0, 0, -1)).inverse(); break;
					case 3: lightViewTrafo = Transform::lookAt(p, p + Vector(0, -1, 0), Vector(0, 0, 1)).inverse(); break;
					case 4: lightViewTrafo = Transform::lookAt(p, p + Vector(0, 0, 1), Vector(0, 1, 0)).inverse(); break;
					case 5: lightViewTrafo = Transform::lookAt(p, p + Vector(0, 0, -1), Vector(0, 1, 0)).inverse(); break;
				}
				lightViewTrafo = Transform::scale(Vector(-1, 1, 1)) * lightViewTrafo;
				const Matrix4x4 &viewMatrix = lightViewTrafo.getMatrix();
				m_shadowProgram->setParameter(m_shadowProgramParam_cubeMapTransform[i], lightProjTrafo * lightViewTrafo);
				m_shadowProgram->setParameter(m_shadowProgramParam_depthVec[i], Vector4(
					-viewMatrix.m[2][0] * m_invClipRange,
					-viewMatrix.m[2][1] * m_invClipRange,
					-viewMatrix.m[2][2] * m_invClipRange,
					(-viewMatrix.m[2][3] - m_nearClip) * m_invClipRange
				));
			}
			m_renderer->drawAll(m_drawList);
		} catch (const std::exception &ex) {
			m_shadowProgram->unbind();
			throw ex;
		}
		m_shadowProgram->unbind();
	} else {
		/* Old-fashioned: render 6 times, once for each cube map face */
		m_altShadowProgram->bind();
		try {
			for (int i=0; i<6; ++i) {
				switch (i) {
					case 0: lightViewTrafo = Transform::lookAt(p, p + Vector(1, 0, 0), Vector(0, 1, 0)).inverse(); break;
					case 1: lightViewTrafo = Transform::lookAt(p, p + Vector(-1, 0, 0), Vector(0, 1, 0)).inverse(); break;
					case 2: lightViewTrafo = Transform::lookAt(p, p + Vector(0, 1, 0), Vector(0, 0, -1)).inverse(); break;
					case 3: lightViewTrafo = Transform::lookAt(p, p + Vector(0, -1, 0), Vector(0, 0, 1)).inverse(); break;
					case 4: lightViewTrafo = Transform::lookAt(p, p + Vector(0, 0, 1), Vector(0, 1, 0)).inverse(); break;
					case 5: lightViewTrafo = Transform::lookAt(p, p + Vector(0, 0, -1), Vector(0, 1, 0)).inverse(); break;
				}
				lightViewTrafo = Transform::scale(Vector(-1, 1, 1)) * lightViewTrafo;
				const Matrix4x4 &viewMatrix = lightViewTrafo.getMatrix();

				m_altShadowProgram->setParameter(m_altShadowProgramParam_cubeMapTransform, lightProjTrafo * lightViewTrafo);
				m_altShadowProgram->setParameter(m_altShadowProgramParam_depthVec, Vector4(
					-viewMatrix.m[2][0] * m_invClipRange,
					-viewMatrix.m[2][1] * m_invClipRange,
					-viewMatrix.m[2][2] * m_invClipRange,
					(-viewMatrix.m[2][3] - m_nearClip) * m_invClipRange
				));

				m_shadowMap->activateSide(i);
				m_shadowMap->clear();
				m_renderer->drawAll(m_drawList);
			}
		} catch (std::exception &ex) {
			m_altShadowProgram->unbind();
			throw ex;
		}

		m_altShadowProgram->unbind();
	}
	m_shadowMap->releaseTarget();
}

void VPLShaderManager::configure(const VPL &vpl, const BSDF *bsdf, 
			const Luminaire *luminaire, const Point &camPos, bool faceNormals) {
	Shader *bsdfShader = m_renderer->getShaderForResource(bsdf);
	Shader *vplShader = (vpl.type == ELuminaireVPL)
		? m_renderer->getShaderForResource(vpl.luminaire)
		: m_renderer->getShaderForResource(vpl.its.shape->getBSDF());
	Shader *lumShader = (luminaire == NULL) ? NULL 
		: m_renderer->getShaderForResource(luminaire);
	std::ostringstream oss;

	if (bsdfShader == NULL || vplShader == NULL ||
		(luminaire != NULL && lumShader == NULL)) {
		/* Unsupported! */
		m_renderer->setColor(Spectrum(0.0f));
		return;
	}

	bool anisotropic = bsdf->getType() & BSDF::EAnisotropic;

	m_targetConfig = VPLProgramConfiguration(vplShader, bsdfShader, 
			lumShader, faceNormals);
	m_targetConfig.toString(oss);
	std::string configName = oss.str();
	std::map<std::string, ProgramAndConfiguration>::iterator it =
		m_programs.find(configName);
	GPUProgram *program = NULL;

	if (it != m_programs.end()) {
		/* A program for this configuration has been created previously */
		m_current = (*it).second;
		program = m_current.program;
	} else {
		/* No program for this particular combination exists -- create one */
		program = m_renderer->createGPUProgram(configName);

		if (faceNormals) {
			/* Generate face normals in a geometry shader */
	
    		if (!m_renderer->getCapabilities()->isSupported(
					RendererCapabilities::EGeometryShaders))
				Log(EError, "Face normals require geometry shader support!");
			if (anisotropic)
				Log(EError, "Anisotropy and face normals can't be combined at the moment");
	
			oss.str("");
			oss << "#version 120" << endl
				<< "#extension GL_EXT_geometry_shader4 : enable" << endl
				<< "varying in vec3 lightVec_vertex[3], camVec_vertex[3];" << endl
				<< "varying in vec2 uv_vertex[3];" << endl
				<< "varying in vec3 vertexColor_vertex[3];" << endl
				<< "varying out vec3 normal;" << endl
				<< "varying out vec3 lightVec, camVec;" << endl
				<< "varying out vec2 uv;" << endl
				<< "varying out vec3 vertexColor;" << endl
				<< endl
				<< "void main() {" << endl
				<< "   vec3 edge1 = camVec_vertex[0]-camVec_vertex[1];" << endl
				<< "   vec3 edge2 = camVec_vertex[0]-camVec_vertex[2];" << endl
				<< "   normal = normalize(cross(edge1, edge2));" << endl
				<< "   gl_Position = vec4(0.0);" << endl
				<< "   lightVec = camVec = vec3(0.0);" << endl
				<< "   for (int i=0; i<gl_VerticesIn; ++i) {" << endl
				<< "      gl_Position = gl_PositionIn[i];" << endl
				<< "      uv = uv_vertex[i];" << endl
				<< "      vertexColor = vertexColor_vertex[i];" << endl
				<< "      lightVec = lightVec_vertex[i];" << endl
				<< "      camVec = camVec_vertex[i];" << endl
				<< "      EmitVertex();" << endl
				<< "   }" << endl
				<< "   EndPrimitive();" << endl
				<< "}" << endl;

			program->setMaxVertices(3); 
			program->setSource(GPUProgram::EGeometryProgram, oss.str());
		}

		/* Vertex program */
		oss.str("");
        oss << "#version 120" << endl;
		if (anisotropic)
			oss << "varying vec3 tangent;" << endl;
		oss << "uniform vec3 vplPos, camPos;" << endl;

		if (!faceNormals) {
			oss << "varying vec3 lightVec, camVec;" << endl
				<< "varying vec2 uv;" << endl
				<< "varying vec3 normal;" << endl
				<< "varying vec3 vertexColor;" << endl
				<< endl
				<< "void main() {" << endl
				<< "   uv = gl_MultiTexCoord0.xy;" << endl
				<< "   camVec = camPos - gl_Vertex.xyz;" << endl
				<< "   lightVec = vplPos - gl_Vertex.xyz;" << endl
				<< "   gl_Position = ftransform();" << endl
				<< "   vertexColor = gl_Color.rgb;" << endl
				<< "   normal = gl_Normal;" << endl;
		} else {
			oss << "varying vec3 lightVec_vertex, camVec_vertex;" << endl
				<< "varying vec2 uv_vertex;" << endl
				<< "varying vec3 vertexColor_vertex;" << endl
				<< endl
				<< "void main() {" << endl
				<< "   uv_vertex = gl_MultiTexCoord0.xy;" << endl
				<< "   camVec_vertex = camPos - gl_Vertex.xyz;" << endl
				<< "   lightVec_vertex = vplPos - gl_Vertex.xyz;" << endl
				<< "   gl_Position = ftransform();" << endl
				<< "   vertexColor_vertex = gl_Color.rgb;" << endl;
		}
		if (anisotropic)
			oss << "   tangent = gl_MultiTexCoord1.xyz;" << endl;
		oss << "}" << endl;

		program->setSource(GPUProgram::EVertexProgram, oss.str());
		oss.str("");

		oss << "#version 120" << endl
			<< endl
			<< "/* Uniform inputs */" << endl
			<< "uniform samplerCube shadowMap;" << endl
			<< "uniform vec3 vplPower, vplS, vplT, vplN, vplWi;" << endl
			<< "uniform float nearClip, invClipRange, minDist, alpha;" << endl
			<< "uniform vec2 vplUV;" << endl
			<< "uniform bool diffuseSources, diffuseReceivers;" << endl
			<< "varying vec3 vertexColor;" << endl
			<< endl
			<< "/* Inputs <- Vertex program */" << endl
			<< "varying vec3 normal, lightVec, camVec;" << endl
			<< "varying vec2 uv;" << endl;
		if (anisotropic)
			oss << "varying vec3 tangent;" << endl;

		oss << endl
			<< "/* Some helper functions for BSDF implementations */" << endl
			<< "float cosTheta(vec3 v) { return v.z; }" << endl
			<< "float sinTheta2(vec3 v) { return 1.0-v.z*v.z; }" << endl
			<< "float sinTheta(vec3 v) { float st2 = sinTheta2(v); if (st2 <= 0) return 0.0; else return sqrt(sinTheta2(v)); }" << endl
			<< "float tanTheta(vec3 v) { return sinTheta(v)/cosTheta(v); }" << endl
			<< "float sinPhi(vec3 v) { return v.y/sinTheta(v); }" << endl
			<< "float cosPhi(vec3 v) { return v.x/sinTheta(v); }" << endl
			<< "const float pi = 3.141592653589;" << endl
			<< "const float inv_pi = 0.318309886183791;" << endl
			<< endl;

		std::string vplEvalName, bsdfEvalName, lumEvalName;
		m_targetConfig.generateCode(oss, vplEvalName, bsdfEvalName, lumEvalName);

		oss << "void main() {" << endl
			<< "   /* Set up an ONB */" << endl
			<< "   vec3 N = normalize(normal);" << endl;
		if (anisotropic) {
			oss << "   vec3 S = normalize(tangent - dot(tangent, N)*N);" << endl;
		} else {
			oss << "   vec3 S;" << endl
				<< "   if (abs(N.x) > abs(N.y)) {" << endl
				<< "        float invLen = 1.0 / sqrt(N.x*N.x + N.z*N.z);" << endl
				<< "        S = vec3(-N.z * invLen, 0.0, N.x * invLen);" << endl
				<< "   } else {" << endl
				<< "        float invLen = 1.0 / sqrt(N.y*N.y + N.z*N.z);" << endl
				<< "        S = vec3(0.0, -N.z * invLen, N.y * invLen);" << endl
				<< "   }" << endl;
		}
		oss << "   vec3 T = cross(N, S);" << endl
			<< endl
			<< "   /* Compute shadows */" << endl
			<< "   float d = length(lightVec);" << endl
			<< "   vec3 nLightVec = lightVec/d, absLightVec = abs(lightVec);" << endl
			<< "   float depth = max(max(absLightVec.x, absLightVec.y), absLightVec.z);" << endl
			<< "   depth = (depth-nearClip) * invClipRange - 0.005;" << endl
			<< "   float shadow = textureCube(shadowMap, nLightVec).r > depth ? 1.0 : 0.0;" << endl
			<< endl
			<< "   /* Shading */" << endl
			<< "   vec3 nCamVec = normalize(camVec);" << endl
			<< "   vec3 wo = vec3(dot(S, nLightVec)," << endl
			<< "                  dot(T, nLightVec)," << endl
			<< "                  dot(N, nLightVec));" << endl
			<< "   vec3 wi = vec3(dot(S, nCamVec)," << endl
			<< "                  dot(T, nCamVec)," << endl
			<< "                  dot(N, nCamVec));" << endl
			<< "   vec3 vplWo = -vec3(dot(vplS, nLightVec)," << endl
			<< "                      dot(vplT, nLightVec)," << endl
			<< "                      dot(vplN, nLightVec));" << endl
			<< "   vec3 contrib = vplPower;" << endl
			<< "   if (!diffuseSources)" << endl 
			<< "      contrib *= " << vplEvalName;
			if (vpl.type == ESurfaceVPL)
				oss << "(vplUV, vplWi, vplWo);" << endl;
			else
				oss << "_dir(vplWo);" << endl;
			if (vpl.type == ESurfaceVPL)
				oss << "   else contrib *= max(0, cosTheta(vplWo));" << endl;
		oss << "   if (d < minDist) d = minDist;" << endl
			<< "   if (!diffuseReceivers)" << endl
			<< "      contrib *= "<< bsdfEvalName << "(uv, wi, wo);" << endl
			<< "   else" << endl
			<< "      contrib *= " << bsdfEvalName << "_diffuse(uv, wi, wo);" << endl
			<< "   gl_FragColor.rgb = contrib";
		if (vpl.type == ELuminaireVPL 
				&& (vpl.luminaire->getType() & Luminaire::EOnSurface))
			oss << " * (shadow * abs(cosTheta(vplWo)) / (d*d))";
		else 
			oss << " * (shadow / (d*d))";
		if (luminaire != NULL) {
			oss << endl;
			oss << "                      + " << lumEvalName << "_area(uv)"
				<< " * " << lumEvalName << "_dir(wi);" << endl;
		} else {
			oss << ";" << endl;
		}
		oss << "   gl_FragColor.a = alpha;" << endl
			<< "}" << endl;

		program->setSource(GPUProgram::EFragmentProgram, oss.str());
		try {
			program->init();
		} catch (const std::exception &) {
			Log(EWarn, "Unable to compile the following VPL program:\n%s", oss.str().c_str());
			throw;
		}

		m_targetConfig.resolve(program);
		m_targetConfig.param_shadowMap = program->getParameterID("shadowMap", false);
		m_targetConfig.param_vplPos = program->getParameterID("vplPos", false);
		m_targetConfig.param_camPos = program->getParameterID("camPos", false);
		m_targetConfig.param_vplPower = program->getParameterID("vplPower", false);
		m_targetConfig.param_vplN = program->getParameterID("vplN", false);
		m_targetConfig.param_vplS = program->getParameterID("vplS", false);
		m_targetConfig.param_vplT = program->getParameterID("vplT", false);
		m_targetConfig.param_vplWi = program->getParameterID("vplWi", false);
		m_targetConfig.param_vplUV = program->getParameterID("vplUV", false);
		m_targetConfig.param_nearClip = program->getParameterID("nearClip", false);
		m_targetConfig.param_invClipRange = program->getParameterID("invClipRange", false);
		m_targetConfig.param_minDist = program->getParameterID("minDist", false);
		m_targetConfig.param_diffuseSources = program->getParameterID("diffuseSources", false);
		m_targetConfig.param_diffuseReceivers = program->getParameterID("diffuseReceivers", false);
		m_targetConfig.param_alpha = program->getParameterID("alpha", false);
		m_current.program = program;
		m_current.config = m_targetConfig;
		m_programs[configName] = m_current;
		program->incRef();
	}

	program->bind();
	m_shadowMap->bind(0);

	const VPLProgramConfiguration &config = m_current.config;

	program->setParameter(config.param_shadowMap, m_shadowMap);
	program->setParameter(config.param_vplPos, vpl.its.p);
	program->setParameter(config.param_camPos, camPos);
	program->setParameter(config.param_vplN, vpl.its.shFrame.n);
	program->setParameter(config.param_vplS, vpl.its.shFrame.s);
	program->setParameter(config.param_vplT, vpl.its.shFrame.t);
	
	program->setParameter(config.param_alpha, 
		bsdfShader->getFlags() & Shader::ETransparent ? 0.5f : 1.0f);

	if (vpl.type == ESurfaceVPL) {
		program->setParameter(config.param_vplWi, vpl.its.wi);
		program->setParameter(config.param_vplUV, vpl.its.uv);
		program->setParameter(config.param_diffuseSources, m_diffuseSources);
	}

	Spectrum power = vpl.P;
	if (m_diffuseSources && vpl.type == ESurfaceVPL)
		power *= vpl.its.shape->getBSDF()->getDiffuseReflectance(vpl.its) * INV_PI;
	program->setParameter(config.param_vplPower, power);
	program->setParameter(config.param_diffuseReceivers, m_diffuseReceivers);
	program->setParameter(config.param_nearClip, m_nearClip);
	program->setParameter(config.param_invClipRange, m_invClipRange);
	program->setParameter(config.param_minDist, m_minDist);

	int textureUnitOffset = 1;
	m_targetConfig.bind(program, config, textureUnitOffset);
}

void VPLShaderManager::drawBackground(const Transform &clipToWorld, const Point &camPos, Float scaleFactor) {
	if (m_backgroundProgram == NULL)
		return;
	int textureUnitOffset = 0;	
	m_backgroundProgram->bind();
	m_backgroundDependencies.recursiveBind(m_backgroundProgram, 
		m_backgroundDependencies, textureUnitOffset);
	m_backgroundProgram->setParameter("clipToWorld", clipToWorld, false);
	m_backgroundProgram->setParameter("camPos", camPos, false);
	m_backgroundProgram->setParameter("scale", scaleFactor);
	m_renderer->blitQuad(false);
	m_backgroundProgram->unbind();
	m_backgroundDependencies.recursiveUnbind();
}

void VPLShaderManager::unbind() {
	if (m_current.program && m_current.program->isBound()) {
		m_targetConfig.unbind();
		m_current.program->unbind();
		m_shadowMap->unbind();
	}
}

MTS_IMPLEMENT_CLASS(VPLShaderManager, false, Object)
MTS_NAMESPACE_END

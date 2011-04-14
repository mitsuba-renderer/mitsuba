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

#include <mitsuba/hw/renderer.h>
#if defined(WIN32)
#include <mitsuba/hw/wglrenderer.h>
#elif defined(__OSX__)
#include <mitsuba/hw/nsglrenderer.h>
#else
#include <mitsuba/hw/glxrenderer.h>
#endif
#include <mitsuba/render/shader.h>
#include <mitsuba/hw/gpugeometry.h>

MTS_NAMESPACE_BEGIN

Renderer::Renderer(Session *session)
 : m_session(session), m_initialized(false), m_logLevel(EDebug),
   m_warnLogLevel(EWarn) {
	m_capabilities = new RendererCapabilities();
}

Renderer::~Renderer() {
}

void Renderer::init(Device *device, Renderer *other) {
	Assert(!m_initialized);
	m_device = device;
}

void Renderer::shutdown() {
	Assert(m_initialized);
}

Renderer *Renderer::create(Session *session) {
#if defined(WIN32)
	return new WGLRenderer(static_cast<WGLSession *>(session));
#elif defined(__OSX__)
	return new NSGLRenderer(static_cast<NSGLSession *>(session));
#else
	return new GLXRenderer(static_cast<X11Session *>(session));
#endif
}

Shader *Renderer::registerShaderForResource(const HWResource *resource) {
	if (m_shaders.find(resource) != m_shaders.end()) {
		ShaderRecord &sRec = m_shaders[resource];
		sRec.refCount++;
		return sRec.shader;
	}

	ShaderRecord sRec;
	sRec.shader = resource->createShader(this);
	if (sRec.shader == NULL) {
		Log(EWarn, "Resource does not have a hardware shader implementation: %s", 
			dynamic_cast<const Object *>(resource)->toString().c_str());
		return NULL;
	}

	sRec.shader->incRef();
	sRec.refCount = 1;
	m_shaders[resource] = sRec;
	return sRec.shader;
}

Shader *Renderer::getShaderForResource(const HWResource *resource) {
	if (m_shaders.find(resource) == m_shaders.end())
		return NULL;
	return m_shaders[resource].shader;
}

void Renderer::unregisterShaderForResource(const HWResource *resource) {
	if (resource == NULL)
		return;
	if (m_shaders.find(resource) == m_shaders.end())
		return;
	ShaderRecord &sRec = m_shaders[resource];
	if (--sRec.refCount == 0) {
		sRec.shader->cleanup(this);
		sRec.shader->decRef();
		m_shaders.erase(resource);
	}
}

GPUGeometry *Renderer::registerGeometry(const TriMesh *mesh) {
	if (!m_capabilities->isSupported(RendererCapabilities::EVertexBufferObjects))
		return NULL;

	GPUGeometry *geometry;
	if (m_geometry.find(mesh) == m_geometry.end()) {
		geometry = createGPUGeometry(mesh);
		m_geometry[mesh] = geometry;
		geometry->init();
	} else {
		geometry = m_geometry[mesh];
	}
	geometry->incRef();
	return geometry;
}

void Renderer::unregisterGeometry(const TriMesh *mesh) {
	if (!m_capabilities->isSupported(RendererCapabilities::EVertexBufferObjects))
		return;
	if (m_geometry.find(mesh) == m_geometry.end())
		Log(EError, "Cannot unregister geometry!");
	GPUGeometry *geometry = m_geometry[mesh];
	int refCount = geometry->getRefCount();
	if (refCount == 1) {
		geometry->cleanup();
		m_geometry.erase(mesh);
	}
	geometry->decRef();
}

std::string RendererCapabilities::toString() const {
	std::ostringstream oss;
	oss << "RenderCapabilities[\n";
	if (m_capabilities[EShadingLanguage])
		oss << "\tEShadingLanguage,\n";
	if (m_capabilities[ERenderToTexture])
		oss << "\tERenderToTexture,\n";
	if (m_capabilities[EBufferBlit])
		oss << "\tEBufferBlit,\n";
	if (m_capabilities[EFloatingPointBuffer])
		oss << "\tEFloatingPointBuffer,\n";
	if (m_capabilities[EFloatingPointTextures])
		oss << "\tEFloatingPointTextures,\n";
	if (m_capabilities[EMultisampleRenderToTexture])
		oss << "\tEEMultisampleRenderToTexture,\n";
	if (m_capabilities[EVertexBufferObjects])
		oss << "\tEVertexBufferObjects,\n";
	if (m_capabilities[EGeometryShaders])
		oss << "\tEGeometryShaders,\n";
	if (m_capabilities[ESyncObjects])
		oss << "\tESyncObjects,\n";
	if (m_capabilities[EBindless])
		oss << "\tEBindless,\n";
	oss << "]";
	std::string result = oss.str();
	return result.substr(0, result.length()-3) + "\n]";
}

MTS_IMPLEMENT_CLASS(Renderer, true, Object)
MTS_IMPLEMENT_CLASS(RendererCapabilities, false, Object)
MTS_NAMESPACE_END

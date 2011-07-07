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

#include <mitsuba/render/texture.h>
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

Texture::Texture(const Properties &props)
 : ConfigurableObject(props) {
}

Texture::Texture(Stream *stream, InstanceManager *manager) 
 : ConfigurableObject(stream, manager) {
}
	
Vector3i Texture::getResolution() const {
	return Vector3i(0);
}

Texture::~Texture() {
}

void Texture::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
}

Texture2D::Texture2D(const Properties &props) : Texture(props) {
	if (props.getString("coordinates", "uv") == "uv") {
		m_uvOffset = Point2(
			props.getFloat("uoffset", 0.0f),
			props.getFloat("voffset", 0.0f)
		);
		m_uvScale = Vector2(
			props.getFloat("uscale", 1.0f),
			props.getFloat("vscale", 1.0f)
		);
	} else {
		Log(EError, "Only UV coordinates are supported at the moment!");
	}
}

Texture2D::Texture2D(Stream *stream, InstanceManager *manager) 
 : Texture(stream, manager) {
	m_uvOffset = Point2(stream);
	m_uvScale = Vector2(stream);
}

Texture2D::~Texture2D() {
}

void Texture2D::serialize(Stream *stream, InstanceManager *manager) const {
	Texture::serialize(stream, manager);
	m_uvOffset.serialize(stream);
	m_uvScale.serialize(stream);
}

Spectrum Texture2D::getValue(const Intersection &its) const {
	Point2 uv = Point2(its.uv.x * m_uvScale.x, its.uv.y * m_uvScale.y) + m_uvOffset;
	if (its.hasUVPartials) {
		return getValue(uv, 
			its.dudx * m_uvScale.x, its.dudy * m_uvScale.x,
			its.dvdx * m_uvScale.y, its.dvdy * m_uvScale.y);
	} else {
		return getValue(uv);
	}
}

ConstantSpectrumTexture::ConstantSpectrumTexture(Stream *stream, InstanceManager *manager) 
 : Texture(stream, manager) {
	m_value = Spectrum(stream);
}

void ConstantSpectrumTexture::serialize(Stream *stream, InstanceManager *manager) const {
	Texture::serialize(stream, manager);

	m_value.serialize(stream);
}

ConstantFloatTexture::ConstantFloatTexture(Stream *stream, InstanceManager *manager) 
 : Texture(stream, manager) {
	m_value = stream->readFloat();
}

void ConstantFloatTexture::serialize(Stream *stream, InstanceManager *manager) const {
	Texture::serialize(stream, manager);
	stream->writeFloat(m_value);
}

ScaleTexture::ScaleTexture(Stream *stream, InstanceManager *manager) 
 : Texture(stream, manager) {
	m_nested = static_cast<Texture *>(manager->getInstance(stream));
	m_scale = stream->readFloat();
}

void ScaleTexture::serialize(Stream *stream, InstanceManager *manager) const {
	Texture::serialize(stream, manager);
	manager->serialize(stream, m_nested.get());
	stream->writeFloat(m_scale);
}

class ConstantSpectrumTextureShader : public Shader {
public:
	ConstantSpectrumTextureShader(Renderer *renderer, const Spectrum &value) 
		: Shader(renderer, ETextureShader), m_value(value) {
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform vec3 " << evalName << "_value;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv) {" << endl
			<< "    return " << evalName << "_value;" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_value"));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_value);
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum m_value;
};

class ConstantFloatTextureShader : public Shader {
public:
	ConstantFloatTextureShader(Renderer *renderer, const Float &value) 
		: Shader(renderer, ETextureShader), m_value(value) {
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform float " << evalName << "_value;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv) {" << endl
			<< "    return vec3(" << evalName << "_value);" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_value"));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_value);
	}

	MTS_DECLARE_CLASS()
private:
	Float m_value;
};

class ScaleTextureShader : public Shader {
public:
	ScaleTextureShader(Renderer *renderer, const Texture *nested, const Float &scale) 
		: Shader(renderer, ETextureShader), m_nested(nested), m_scale(scale) {
		m_nestedShader = renderer->registerShaderForResource(m_nested.get());
	}

	bool isComplete() const {
		return m_nestedShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_nested.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_nestedShader.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform float " << evalName << "_scale;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv) {" << endl
			<< "    return " << depNames[0] << "(uv) * " << evalName << "_scale;" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_scale"));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &nestedUnitOffset) const {
		program->setParameter(parameterIDs[0], m_scale);
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_nested;
	ref<Shader> m_nestedShader;
	Float m_scale;
};

Shader *ConstantSpectrumTexture::createShader(Renderer *renderer) const { 
	return new ConstantSpectrumTextureShader(renderer, m_value);
}

Shader *ConstantFloatTexture::createShader(Renderer *renderer) const { 
	return new ConstantFloatTextureShader(renderer, m_value);
}

Shader *ScaleTexture::createShader(Renderer *renderer) const { 
	return new ScaleTextureShader(renderer, m_nested.get(), m_scale);
}

MTS_IMPLEMENT_CLASS(Texture, true, ConfigurableObject)
MTS_IMPLEMENT_CLASS(Texture2D, true, Texture)
MTS_IMPLEMENT_CLASS_S(ConstantSpectrumTexture, false, Texture)
MTS_IMPLEMENT_CLASS(ConstantSpectrumTextureShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(ConstantFloatTexture, false, Texture)
MTS_IMPLEMENT_CLASS(ConstantFloatTextureShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(ScaleTexture, false, Texture)
MTS_IMPLEMENT_CLASS(ScaleTextureShader, false, Shader)
MTS_NAMESPACE_END

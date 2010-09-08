/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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

Texture::~Texture() {
}
	
void Texture::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
}

ConstantTexture::ConstantTexture(Stream *stream, InstanceManager *manager) 
 : Texture(stream, manager) {
	m_value = Spectrum(stream);
}

void ConstantTexture::serialize(Stream *stream, InstanceManager *manager) const {
	Texture::serialize(stream, manager);

	m_value.serialize(stream);
}

class ConstantTextureShader : public Shader {
public:
	ConstantTextureShader(Renderer *renderer, const Spectrum &value) 
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

Shader *ConstantTexture::createShader(Renderer *renderer) const { 
	return new ConstantTextureShader(renderer, m_value);
}

MTS_IMPLEMENT_CLASS(Texture, true, ConfigurableObject)
MTS_IMPLEMENT_CLASS_S(ConstantTexture, false, Texture)
MTS_IMPLEMENT_CLASS(ConstantTextureShader, false, Shader)
MTS_NAMESPACE_END

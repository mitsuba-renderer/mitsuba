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

#include <mitsuba/render/shader.h>

MTS_NAMESPACE_BEGIN

Shader::Shader(Renderer *renderer, EShaderType type)
  : m_type(type), m_flags(0) {
}

Shader::~Shader() {
}

/* These do nothing by default */
void Shader::putDependencies(std::vector<Shader *> &deps) { }
void Shader::bind(GPUProgram *program, const std::vector<int> &parameterIDs,
	int &textureUnitOffset) const { }
void Shader::resolve(const GPUProgram *program, const std::string &evalName,
	std::vector<int> &parameterIDs) const { }
void Shader::unbind() const { }
void Shader::cleanup(Renderer *renderer) { }
bool Shader::isComplete() const { return true; }
Float Shader::getAlpha() const { return 1.0f; }

Shader *HWResource::createShader(Renderer *renderer) const {
	return NULL;
}

MTS_IMPLEMENT_CLASS(Shader, true, Object)
MTS_NAMESPACE_END

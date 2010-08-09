#include <mitsuba/render/shader.h>

MTS_NAMESPACE_BEGIN
	
Shader::Shader(Renderer *renderer, EShaderType type)
  : m_type(type) {
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

Shader *HWResource::createShader(Renderer *renderer) const {
	return NULL;
}

MTS_IMPLEMENT_CLASS(Shader, true, Object)
MTS_NAMESPACE_END

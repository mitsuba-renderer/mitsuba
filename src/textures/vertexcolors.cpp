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
#include <mitsuba/render/shape.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

/**
 * Vertex colors passthrough texture
 */
class VertexColors : public Texture {
public:
	VertexColors(const Properties &props) : Texture(props) {
	}

	VertexColors(Stream *stream, InstanceManager *manager) 
	 : Texture(stream, manager) {
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture::serialize(stream, manager);
	}

	Spectrum getValue(const Intersection &its) const {
		return its.color;
	}

	bool usesRayDifferentials() const {
		return false;
	}

	Spectrum getAverage() const {
		return Spectrum(1.0f);
	}
	
	Spectrum getMaximum() const {
		return Spectrum(1.0f);
	}

	std::string toString() const {
		return "VertexColors[]";
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
protected:
	Spectrum m_reflectance;
};

// ================ Hardware shader implementation ================ 

class VertexColorShader : public Shader {
public:
	VertexColorShader(Renderer *renderer) : Shader(renderer, ETextureShader) {
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv) {" << endl
			<< "    return vertexColor;" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum m_brightReflectance;
	Spectrum m_darkReflectance;
	Float m_width;
};

Shader *VertexColors::createShader(Renderer *renderer) const {
	return new VertexColorShader(renderer);
}
	
MTS_IMPLEMENT_CLASS(VertexColorShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(VertexColors, false, Texture)
MTS_EXPORT_PLUGIN(VertexColors, "VertexColors texture");
MTS_NAMESPACE_END

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
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/**
 * Checkerboard texture
 */
class Checkerboard : public Texture2D {
public:
	Checkerboard(const Properties &props) : Texture2D(props) {
		m_brightColor = props.getSpectrum("brightColor", Spectrum(.4f));
		m_darkColor = props.getSpectrum("darkColor", Spectrum(.2f));
	}

	Checkerboard(Stream *stream, InstanceManager *manager) 
	 : Texture2D(stream, manager) {
		m_brightColor = Spectrum(stream);
		m_darkColor = Spectrum(stream);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture2D::serialize(stream, manager);
		m_brightColor.serialize(stream);
		m_darkColor.serialize(stream);
	}

	inline Spectrum getValue(const Point2 &uv) const {
		int x = 2*modulo((int) (uv.x * 2), 2) - 1, 
			y = 2*modulo((int) (uv.y * 2), 2) - 1;

		if (x*y == 1)
			return m_brightColor;
		else
			return m_darkColor;
	}

	Spectrum getValue(const Point2 &uv, Float dudx, Float dudy, Float dvdx, Float dvdy) const {
		return Checkerboard::getValue(uv);
	}

	bool usesRayDifferentials() const {
		return false;
	}

	Spectrum getAverage() const {
		return m_darkColor * .5f;
	}
	
	Spectrum getMaximum() const {
		return m_brightColor;
	}

	std::string toString() const {
		return "Checkerboard[]";
	}
	
	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
protected:
	Spectrum m_darkColor;
	Spectrum m_brightColor;
};

// ================ Hardware shader implementation ================ 

class CheckerboardShader : public Shader {
public:
	CheckerboardShader(Renderer *renderer, const Spectrum &brightColor, 
		const Spectrum &darkColor, const Point2 &uvOffset,
		const Vector2 &uvScale) : Shader(renderer, ETextureShader),
		m_brightColor(brightColor), m_darkColor(darkColor), 
		m_uvOffset(uvOffset), m_uvScale(uvScale) {
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform vec3 " << evalName << "_brightColor;" << endl
			<< "uniform vec3 " << evalName << "_darkColor;" << endl
			<< "uniform vec2 " << evalName << "_uvOffset;" << endl
			<< "uniform vec2 " << evalName << "_uvScale;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv) {" << endl
			<< "    uv = vec2(" << endl
			<< "        uv.x * " << evalName << "_uvScale.x + " << evalName << "_uvOffset.x," << endl
			<< "        uv.y * " << evalName << "_uvScale.y + " << evalName << "_uvOffset.y);" << endl
			<< "    float x = 2*(mod(int(uv.x*2), 2)) - 1, y = 2*(mod(int(uv.y*2), 2)) - 1;" << endl
			<< "    if (x*y == 1)" << endl
			<< "        return " << evalName << "_brightColor;" << endl
			<< "    else" << endl
			<< "        return " << evalName << "_darkColor;" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_brightColor", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_darkColor", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_uvOffset", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_uvScale", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, 
		int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_brightColor);
		program->setParameter(parameterIDs[1], m_darkColor);
		program->setParameter(parameterIDs[2], m_uvOffset);
		program->setParameter(parameterIDs[3], m_uvScale);
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum m_brightColor;
	Spectrum m_darkColor;
	Point2 m_uvOffset;
	Vector2 m_uvScale;
};

Shader *Checkerboard::createShader(Renderer *renderer) const {
	return new CheckerboardShader(renderer, m_brightColor, m_darkColor, 
		m_uvOffset, m_uvScale);
}
	
MTS_IMPLEMENT_CLASS(CheckerboardShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Checkerboard, false, Texture2D)
MTS_EXPORT_PLUGIN(Checkerboard, "Checkerboard texture");
MTS_NAMESPACE_END

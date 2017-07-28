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

#include <mitsuba/render/texture.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{gridtexture}{Procedural grid texture}
 * \order{2}
 * \parameters{
 *     \parameter{color0}{\Spectrum}{
 *       Color values of the background
 *       \default{0.2}
 *     }
 *     \parameter{color1}{\Spectrum}{
 *       Color value of the lines
 *       \default{0.4}
 *     }
 *     \parameter{lineWidth}{\Float}{Width of the grid lines in UV space
 *        \default{0.01}
 *     }
 *     \parameter{uscale, vscale}{\Float}{
 *       Multiplicative factors that should be applied to UV values before a lookup
 *     }
 *     \parameter{uoffset, voffset}{\Float}{
 *       Numerical offset that should be applied to UV values before a lookup
 *     }
 * }
 * \renderings{
 *     \rendering{Grid texture applied to the material test object}{tex_gridtexture}
 * }
 * This plugin implements a simple procedural grid texture with customizable
 * colors and line width.
 */
class GridTexture : public Texture2D {
public:
    GridTexture(const Properties &props) : Texture2D(props) {
        m_color0 = props.getSpectrum("color0", Spectrum(.2f));
        m_color1 = props.getSpectrum("color1", Spectrum(.4f));
        m_lineWidth = props.getFloat("lineWidth", .01f);
    }

    GridTexture(Stream *stream, InstanceManager *manager)
     : Texture2D(stream, manager) {
        m_color0 = Spectrum(stream);
        m_color1 = Spectrum(stream);
        m_lineWidth = stream->readFloat();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Texture2D::serialize(stream, manager);
        m_color0.serialize(stream);
        m_color1.serialize(stream);
        stream->writeFloat(m_lineWidth);
    }

    inline Spectrum eval(const Point2 &uv) const {
        Float x = uv.x - math::floorToInt(uv.x);
        Float y = uv.y - math::floorToInt(uv.y);

        if (x > .5)
            x -= 1;
        if (y > .5)
            y -= 1;

        if (std::abs(x) < m_lineWidth || std::abs(y) < m_lineWidth)
            return m_color1;
        else
            return m_color0;
    }

    Spectrum eval(const Point2 &uv,
            const Vector2 &d0, const Vector2 &d1) const {
        /* Filtering is currently not supported */
        return GridTexture::eval(uv);
    }

    bool usesRayDifferentials() const {
        return false;
    }

    Spectrum getMaximum() const {
        Spectrum max;
        for (int i=0; i<SPECTRUM_SAMPLES; ++i)
            max[i] = std::max(m_color0[i], m_color1[i]);
        return max;
    }

    Spectrum getMinimum() const {
        Spectrum min;
        for (int i=0; i<SPECTRUM_SAMPLES; ++i)
            min[i] = std::min(m_color0[i], m_color1[i]);
        return min;
    }

    Spectrum getAverage() const {
        Float interiorWidth = std::max((Float) 0.0f, 1-2*m_lineWidth),
              interiorArea = interiorWidth * interiorWidth,
              lineArea = 1 - interiorArea;
        return m_color1 * lineArea + m_color0 * interiorArea;
    }

    bool isConstant() const {
        return false;
    }

    bool isMonochromatic() const {
        return Spectrum(m_color0[0]) == m_color0
            && Spectrum(m_color1[0]) == m_color1;
    }

    std::string toString() const {
        return "GridTexture[]";
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
protected:
    Spectrum m_color0;
    Spectrum m_color1;
    Float m_lineWidth;
};

// ================ Hardware shader implementation ================

class GridTextureShader : public Shader {
public:
    GridTextureShader(Renderer *renderer, const Spectrum &color0,
        const Spectrum &color1, Float lineWidth, const Point2 &uvOffset,
        const Vector2 &uvScale) : Shader(renderer, ETextureShader),
        m_color0(color0), m_color1(color1),
        m_lineWidth(lineWidth), m_uvOffset(uvOffset), m_uvScale(uvScale) {
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "uniform vec3 " << evalName << "_color0;" << endl
            << "uniform vec3 " << evalName << "_color1;" << endl
            << "uniform float " << evalName << "_lineWidth;" << endl
            << "uniform vec2 " << evalName << "_uvOffset;" << endl
            << "uniform vec2 " << evalName << "_uvScale;" << endl
            << endl
            << "vec3 " << evalName << "(vec2 uv) {" << endl
            << "    uv = vec2(" << endl
            << "        uv.x * " << evalName << "_uvScale.x + " << evalName << "_uvOffset.x," << endl
            << "        uv.y * " << evalName << "_uvScale.y + " << evalName << "_uvOffset.y);" << endl
            << "    float x = uv.x - floor(uv.x);" << endl
            << "    float y = uv.y - floor(uv.y);" << endl
            << "    if (x > .5) x -= 1.0;" << endl
            << "    if (y > .5) y -= 1.0;" << endl
            << "    if (abs(x) < " << evalName << "_lineWidth || abs(y) < " << evalName << "_lineWidth)" << endl
            << "        return " << evalName << "_color1;" << endl
            << "    else" << endl
            << "        return " << evalName << "_color0;" << endl
            << "}" << endl;
    }

    void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
        parameterIDs.push_back(program->getParameterID(evalName + "_color0", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_color1", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_lineWidth", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_uvOffset", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_uvScale", false));
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs,
        int &textureUnitOffset) const {
        program->setParameter(parameterIDs[0], m_color0);
        program->setParameter(parameterIDs[1], m_color1);
        program->setParameter(parameterIDs[2], m_lineWidth);
        program->setParameter(parameterIDs[3], m_uvOffset);
        program->setParameter(parameterIDs[4], m_uvScale);
    }

    MTS_DECLARE_CLASS()
private:
    Spectrum m_color0;
    Spectrum m_color1;
    Float m_lineWidth;
    Point2 m_uvOffset;
    Vector2 m_uvScale;
};

Shader *GridTexture::createShader(Renderer *renderer) const {
    return new GridTextureShader(renderer, m_color0, m_color1,
            m_lineWidth, m_uvOffset, m_uvScale);
}

MTS_IMPLEMENT_CLASS(GridTextureShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(GridTexture, false, Texture2D)
MTS_EXPORT_PLUGIN(GridTexture, "Grid texture");
MTS_NAMESPACE_END

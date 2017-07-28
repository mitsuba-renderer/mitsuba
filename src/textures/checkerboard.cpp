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

/*!\plugin{checkerboard}{Checkerboard}
 * \order{2}
 * \parameters{
 *     \parameter{color0, color1}{\Spectrum}{
 *       Color values for the two differently-colored patches
 *       \default{0.4 and 0.2}
 *     }
 *     \parameter{uoffset, voffset}{\Float}{
 *       Numerical offset that should be applied to UV values before a lookup
 *     }
 *     \parameter{uscale, vscale}{\Float}{
 *       Multiplicative factors that should be applied to UV values before a lookup
 *     }
 * }
 * \renderings{
 *     \rendering{Checkerboard applied to the material test object
 *                as well as the ground plane}{tex_checkerboard}
 * }
 * This plugin implements a simple procedural checkerboard texture with
 * customizable colors.
 */
class Checkerboard : public Texture2D {
public:
    Checkerboard(const Properties &props) : Texture2D(props) {
        m_color0 = props.getSpectrum("color0", Spectrum(.4f));
        m_color1 = props.getSpectrum("color1", Spectrum(.2f));
    }

    Checkerboard(Stream *stream, InstanceManager *manager)
     : Texture2D(stream, manager) {
        m_color0 = Spectrum(stream);
        m_color1 = Spectrum(stream);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Texture2D::serialize(stream, manager);
        m_color0.serialize(stream);
        m_color1.serialize(stream);
    }

    inline Spectrum eval(const Point2 &uv) const {
        int x = 2*math::modulo((int) (uv.x * 2), 2) - 1,
            y = 2*math::modulo((int) (uv.y * 2), 2) - 1;

        if (x*y == 1)
            return m_color0;
        else
            return m_color1;
    }

    Spectrum eval(const Point2 &uv,
            const Vector2 &d0, const Vector2 &d1) const {
        /* Filtering is currently not supported */
        return Checkerboard::eval(uv);
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
        return (m_color0 + m_color1) * 0.5f;
    }

    bool isConstant() const {
        return false;
    }

    bool isMonochromatic() const {
        return Spectrum(m_color0[0]) == m_color0
            && Spectrum(m_color1[0]) == m_color1;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "Checkerboard[" << endl
            << "    color1 = " << m_color1.toString() << "," << endl
            << "    color0 = " << m_color0.toString() << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
protected:
    Spectrum m_color1;
    Spectrum m_color0;
};

// ================ Hardware shader implementation ================

class CheckerboardShader : public Shader {
public:
    CheckerboardShader(Renderer *renderer, const Spectrum &color0,
        const Spectrum &color1, const Point2 &uvOffset,
        const Vector2 &uvScale) : Shader(renderer, ETextureShader),
        m_color0(color0), m_color1(color1),
        m_uvOffset(uvOffset), m_uvScale(uvScale) {
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "uniform vec3 " << evalName << "_color0;" << endl
            << "uniform vec3 " << evalName << "_color1;" << endl
            << "uniform vec2 " << evalName << "_uvOffset;" << endl
            << "uniform vec2 " << evalName << "_uvScale;" << endl
            << endl
            << "vec3 " << evalName << "(vec2 uv) {" << endl
            << "    uv = vec2(" << endl
            << "        uv.x * " << evalName << "_uvScale.x + " << evalName << "_uvOffset.x," << endl
            << "        uv.y * " << evalName << "_uvScale.y + " << evalName << "_uvOffset.y);" << endl
            << "    float x = 2*(mod(int(uv.x*2), 2)) - 1, y = 2*(mod(int(uv.y*2), 2)) - 1;" << endl
            << "    if (x*y == 1)" << endl
            << "        return " << evalName << "_color0;" << endl
            << "    else" << endl
            << "        return " << evalName << "_color1;" << endl
            << "}" << endl;
    }

    void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
        parameterIDs.push_back(program->getParameterID(evalName + "_color0", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_color1", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_uvOffset", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_uvScale", false));
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs,
        int &textureUnitOffset) const {
        program->setParameter(parameterIDs[0], m_color0);
        program->setParameter(parameterIDs[1], m_color1);
        program->setParameter(parameterIDs[2], m_uvOffset);
        program->setParameter(parameterIDs[3], m_uvScale);
    }

    MTS_DECLARE_CLASS()
private:
    Spectrum m_color0;
    Spectrum m_color1;
    Point2 m_uvOffset;
    Vector2 m_uvScale;
};

Shader *Checkerboard::createShader(Renderer *renderer) const {
    return new CheckerboardShader(renderer, m_color0, m_color1,
        m_uvOffset, m_uvScale);
}

MTS_IMPLEMENT_CLASS(CheckerboardShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Checkerboard, false, Texture2D)
MTS_EXPORT_PLUGIN(Checkerboard, "Checkerboard texture");
MTS_NAMESPACE_END

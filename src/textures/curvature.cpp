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
#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{curvature}{Curvature texture}
 * \order{7}
 * \parameters{
 *     \parameter{curvature}{\String}{
 *       Specifies what should be shown -- must be equal to
 *       \code{mean} or \code{gaussian}.
 *     }
 *     \parameter{scale}{\Float}{
 *        A scale factor to bring curvature values into the
 *        displayable range [-1, 1]. Everything outside of this range
 *        will be clamped.
 *     }
 * }
 *
 * \renderings{
 *     \rendering{Mean curvature}{tex_curvature_mean}
 *     \rendering{Gaussian curvature}{tex_curvature_gaussian}
 * }
 *
 * This texture can visualize the mean and Gaussian curvature of the underlying
 * shape for inspection. Red and blue denote positive and negative values,
 * respectively.
 */
class Curvature : public Texture {
public:
    Curvature(const Properties &props) : Texture(props) {
        m_scale = props.getFloat("scale");
        std::string curvature = props.getString("curvature", "gaussian");
        if (curvature == "gaussian")
            m_showK = true;
        else if (curvature == "mean")
            m_showK = false;
        else
            Log(EError, "Invalid 'curvature' parameter: must be set to 'gaussian' or ' mean'!");
    }

    Curvature(Stream *stream, InstanceManager *manager)
     : Texture(stream, manager) {
         m_scale = stream->readFloat();
         m_showK = stream->readBool();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Texture::serialize(stream, manager);
        stream->writeFloat(m_scale);
        stream->writeBool(m_showK);
    }

    Spectrum lookupGradient(Float value) const {
        Spectrum result(0.0f);
        if (value < 0)
            result.fromLinearRGB(0, 0, std::min(-value*m_scale, (Float) 1.0f));
        if (value > 0)
            result.fromLinearRGB(std::min(value*m_scale, (Float) 1.0f), 0.0f, 0.0f);
        return result;
    }

    Spectrum eval(const Intersection &its, bool /* unused */) const {
        Float H, K;
        its.shape->getCurvature(its, H, K);
        return lookupGradient(m_showK ? K : H);
    }

    bool usesRayDifferentials() const {
        /// Intentionally return \c true here so that required tangent
        /// space information is provided for triangle meshes
        return true;
    }

    Spectrum getAverage() const {
        /// No idea, really
        return Spectrum(0.5f);
    }

    Spectrum getMinimum() const {
        return Spectrum(0.0f);
    }

    Spectrum getMaximum() const {
        return Spectrum(1.0f);
    }

    bool isMonochromatic() const {
        return false;
    }

    bool isConstant() const {
        return false;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "Curvature[" << endl
            << "   scale = " << m_scale << "," << endl
            << "   showK = " << m_showK << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    Float m_scale;
    bool m_showK;
};

// ================ Hardware shader implementation ================

class CurvatureShader : public Shader {
public:
    CurvatureShader(Renderer *renderer, const Spectrum &value)
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
        parameterIDs.push_back(program->getParameterID(evalName + "_value", false));
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
        program->setParameter(parameterIDs[0], m_value);
    }

    MTS_DECLARE_CLASS()
private:
    Spectrum m_value;
};

Shader *Curvature::createShader(Renderer *renderer) const {
    return new CurvatureShader(renderer, Spectrum(0.5f));
}

MTS_IMPLEMENT_CLASS(CurvatureShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Curvature, false, Texture)
MTS_EXPORT_PLUGIN(Curvature, "Curvature texture");
MTS_NAMESPACE_END

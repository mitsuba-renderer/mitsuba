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

/*!\plugin{wireframe}{Wireframe texture}
 * \order{6}
 * \parameters{
 *     \parameter{interiorColor}{\Spectrum}{
 *       Color value of the interior of triangles
 *       \default{0.5}
 *     }
 *     \parameter{edgeColor}{\Spectrum}{
 *       Edge color value
 *       \default{0.1}
 *     }
 *     \parameter{lineWidth}{\Float}{
 *        World-space width of the mesh edges
 *        \default{automatic}
 *     }
 *     \parameter{stepWidth}{\Float}{
 *        Controls the width of the step function used for the
 *        color transition. It is specified as a value between zero
 *        and one (relative to the \code{lineWidth} parameter)
 *        \default{0.5}
 *     }
 * }
 * \renderings{
 *     \rendering{Wireframe texture applied to the material test object}{tex_wireframe}
 * }
 *
 * This plugin implements a simple two-color wireframe texture map
 * that reveals the structure of a triangular mesh.
 */
class WireFrame : public Texture {
public:
    WireFrame(const Properties &props) : Texture(props) {
        m_lineWidth = props.getFloat("lineWidth", 0.0f);
        m_stepWidth = props.getFloat("stepWidth", 0.5f);
        m_edgeColor = props.getSpectrum("edgeColor", Spectrum(0.1f));
        m_interiorColor = props.getSpectrum("interiorColor", Spectrum(.5f));
        m_stepWidth = std::max((Float) 0.0f, std::min(m_stepWidth, (Float) 1.0f));
        m_mutex = new Mutex();
    }

    WireFrame(Stream *stream, InstanceManager *manager)
     : Texture(stream, manager) {
        m_mutex = new Mutex();
        m_edgeColor = Spectrum(stream);
        m_interiorColor = Spectrum(stream);
        m_lineWidth = stream->readFloat();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Texture::serialize(stream, manager);
        m_edgeColor.serialize(stream);
        m_interiorColor.serialize(stream);
        stream->writeFloat(m_lineWidth);
    }

    Spectrum eval(const Intersection &its, bool /* unused */) const {
        if (!its.shape->getClass()->derivesFrom(MTS_CLASS(TriMesh)))
            return m_interiorColor;

        const TriMesh *triMesh = static_cast<const TriMesh *>(its.shape);
        const Point *positions = triMesh->getVertexPositions();
        if (its.primIndex >= triMesh->getTriangleCount())
            return m_interiorColor;

        if (m_lineWidth == 0) {
            /* Somewhat hacky but probably helpful in many cases.
               This tries to find a suitable line width, which is set
               to 10% of the average average edge length */
            LockGuard lock(m_mutex);
            if (m_lineWidth == 0) {
                Float lineWidth = 0;
                for (size_t i=0; i<triMesh->getTriangleCount(); ++i) {
                    const Triangle &tri = triMesh->getTriangles()[i];
                    for (int j=0; j<3; ++j)
                        lineWidth += (positions[tri.idx[j]]
                            - positions[tri.idx[(j+1)%3]]).length();
                }

                m_lineWidth = 0.1f * lineWidth / (3 * triMesh->getTriangleCount());
            }
        }

        const Triangle &tri = triMesh->getTriangles()[its.primIndex];

        Float minDist = std::numeric_limits<Float>::infinity();
        for (int i=0; i<3; ++i) {
            const Point& cur  = positions[tri.idx[i]];
            const Point& next = positions[tri.idx[(i+1)%3]];

            Vector d1 = normalize(next - cur),
                   d2 = its.p - cur;

            minDist = std::min(minDist, (cur + d1 * dot(d1, d2) - its.p).lengthSquared());
        }

        Float a = math::smoothStep(m_lineWidth*(1.f-m_stepWidth), m_lineWidth, std::sqrt(minDist));
        return m_edgeColor*(1-a) + m_interiorColor*a;
    }

    bool usesRayDifferentials() const {
        return false;
    }

    Spectrum getAverage() const {
        Spectrum value;
        /* Approximate ... */
        for (int i=0; i<SPECTRUM_SAMPLES; ++i)
            value[i] = 0.5f * (m_edgeColor[i] + m_interiorColor[i]);
        return value;
    }

    Spectrum getMinimum() const {
        Spectrum value;
        for (int i=0; i<SPECTRUM_SAMPLES; ++i)
            value[i] = std::min(m_edgeColor[i], m_interiorColor[i]);
        return value;
    }

    Spectrum getMaximum() const {
        Spectrum value;
        for (int i=0; i<SPECTRUM_SAMPLES; ++i)
            value[i] = std::max(m_edgeColor[i], m_interiorColor[i]);
        return value;
    }

    bool isConstant() const {
        return false;
    }

    bool isMonochromatic() const {
        return Spectrum(m_edgeColor[0]) == m_edgeColor
            && Spectrum(m_interiorColor[0]) == m_interiorColor;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "WireFrame[" << endl
            << "  edgeColor = " << m_edgeColor.toString() << "," << endl
            << "  interiorColor = " << m_interiorColor.toString() << "," << endl
            << "  lineWidth = " << m_lineWidth << endl
            << "  stepWidth = " << m_stepWidth << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
protected:
    mutable Float m_lineWidth;
    mutable ref<Mutex> m_mutex;
    Float m_stepWidth;
    Spectrum m_edgeColor;
    Spectrum m_interiorColor;
};

// ================ Hardware shader implementation ================

class WireFrameShader : public Shader {
public:
    WireFrameShader(Renderer *renderer, const Spectrum &value)
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

Shader *WireFrame::createShader(Renderer *renderer) const {
    return new WireFrameShader(renderer, m_interiorColor);
}

MTS_IMPLEMENT_CLASS(WireFrameShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(WireFrame, false, Texture)
MTS_EXPORT_PLUGIN(WireFrame, "Wireframe texture");
MTS_NAMESPACE_END

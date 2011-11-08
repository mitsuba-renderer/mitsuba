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
#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{wireframe}{Wireframe texture}
 */
class WireFrame : public Texture {
public:
	WireFrame(const Properties &props) : Texture(props) {
		m_edgeColor = props.getSpectrum("edgeColor", Spectrum(0.0f));
		m_interiorColor = props.getSpectrum("interiorColor", Spectrum(.5f));
		m_lineWidth = props.getFloat("lineWidth", .01f);
	}

	WireFrame(Stream *stream, InstanceManager *manager) 
	 : Texture(stream, manager) {

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

	Spectrum getValue(const Intersection &its) const {
		if (!its.shape->getClass()->derivesFrom(MTS_CLASS(TriMesh)))
			return m_interiorColor;

		const TriMesh *trimesh = static_cast<const TriMesh *>(its.shape);
		if (its.primIndex >= trimesh->getTriangleCount())
			return m_interiorColor;

		const Triangle &tri = trimesh->getTriangles()[its.primIndex];
		const Point *positions = trimesh->getVertexPositions();

		Point pos[] = { positions[tri.idx[0]],
			positions[tri.idx[1]],
			positions[tri.idx[2]] };

		Float minDist = std::numeric_limits<Float>::infinity();
		for (int i=0; i<3; ++i) {
			Point cur = pos[i], next = pos[(i+1)%3];

			Vector d1 = normalize(next - cur),
				   d2 = its.p - cur;

			minDist = std::min(minDist, (cur + d1 * dot(d1, d2) - its.p).length());
		}

		Float value = 1-smoothStep(0, m_lineWidth, minDist);

		return m_edgeColor * value + m_interiorColor * (1-value);
	}

	bool usesRayDifferentials() const {
		return false;
	}

	Spectrum getAverage() const {
		Spectrum value;
		/* Approximate ... */
		for (size_t i=0; i<SPECTRUM_SAMPLES; ++i)
			value[i] = 0.5f * (m_edgeColor[i] + m_interiorColor[i]);
		return value;
	}

	Spectrum getMinimum() const {
		Spectrum value;
		for (size_t i=0; i<SPECTRUM_SAMPLES; ++i)
			value[i] = std::min(m_edgeColor[i], m_interiorColor[i]);
		return value;
	}

	Spectrum getMaximum() const {
		Spectrum value;
		for (size_t i=0; i<SPECTRUM_SAMPLES; ++i)
			value[i] = std::max(m_edgeColor[i], m_interiorColor[i]);
		return value;
	}

	bool isConstant() const {
		return false;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "WireFrame[" << endl
			<< "  edgeColor = " << m_edgeColor.toString() << ","
			<< "  interiorColor = " << m_interiorColor.toString() << ","
			<< "  lineWidth = " << m_lineWidth << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
protected:
	Spectrum m_edgeColor;
	Spectrum m_interiorColor;
	Float m_lineWidth;
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
MTS_EXPORT_PLUGIN(WireFrame, "Vertex color texture");
MTS_NAMESPACE_END

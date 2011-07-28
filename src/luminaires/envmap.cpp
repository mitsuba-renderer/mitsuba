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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/mipmap.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/hw/gputexture.h>
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/**
 * Basic environment map implementation without importance sampling.
 * Uses the scene's bounding sphere to simulate an infinitely far-away
 * light source. Expects an EXR image in latitude-longitude 
 * (equirectangular) format.
 */
class EnvMapLuminaire : public Luminaire {
public:
	EnvMapLuminaire(const Properties &props) : Luminaire(props) {
		m_intensityScale = props.getFloat("intensityScale", 1);
		ref<Bitmap> bitmap;

		if (props.hasProperty("bitmap")) {
			bitmap = reinterpret_cast<Bitmap *>(props.getData("bitmap").ptr);
			m_path = "<unknown>";
		} else {
			m_path = Thread::getThread()->getFileResolver()->resolve(props.getString("filename"));
			Log(EInfo, "Loading environment map \"%s\"", m_path.leaf().c_str());
			ref<Stream> is = new FileStream(m_path, FileStream::EReadOnly);
			bitmap = new Bitmap(Bitmap::EEXR, is);
		}

		m_mipmap = MIPMap::fromBitmap(bitmap, MIPMap::ETrilinear,
				MIPMap::ERepeat, 0.0f, Spectrum::EIlluminant);
		m_average = m_mipmap->triangle(m_mipmap->getLevels()-1, 0, 0) * m_intensityScale;
		m_type = EOnSurface;
	}

	EnvMapLuminaire(Stream *stream, InstanceManager *manager) 
		: Luminaire(stream, manager) {
		m_intensityScale = stream->readFloat();
		m_path = stream->readString();
		m_bsphere = BSphere(stream);
		Log(EInfo, "Unserializing environment map \"%s\"", m_path.leaf().c_str());
		uint32_t size = stream->readUInt();
		ref<MemoryStream> mStream = new MemoryStream(size);
		stream->copyTo(mStream, size);
		mStream->setPos(0);
		ref<Bitmap> bitmap = new Bitmap(Bitmap::EEXR, mStream);
		m_mipmap = MIPMap::fromBitmap(bitmap, MIPMap::ETrilinear,
				MIPMap::ERepeat, 0.0f, Spectrum::EIlluminant);
		m_average = m_mipmap->triangle(m_mipmap->getLevels()-1, 0, 0) * m_intensityScale;
		m_surfaceArea = 4 * m_bsphere.radius * m_bsphere.radius * M_PI;
		m_invSurfaceArea = 1/m_surfaceArea;

		if (Scheduler::getInstance()->hasRemoteWorkers()
			&& !fs::exists(m_path)) {
			/* This code is running on a machine different from
			   the one that created the stream. Because we might
			   later have to handle a call to serialize(), the
			   whole bitmap must be kept in memory */
			m_stream = mStream;
		}
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Luminaire::serialize(stream, manager);
		stream->writeFloat(m_intensityScale);
		stream->writeString(m_path.file_string());
		m_bsphere.serialize(stream);

		if (m_stream.get()) {
			stream->writeUInt((unsigned int) m_stream->getSize());
			stream->write(m_stream->getData(), m_stream->getSize());
		} else {
			ref<Stream> is = new FileStream(m_path, FileStream::EReadOnly);
			stream->writeUInt((uint32_t) is->getSize());
			is->copyTo(stream);
		}
	}

	void configure() {
		int mipMapLevel = std::min(3, m_mipmap->getLevels()-1);
		m_pdfResolution = m_mipmap->getLevelResolution(mipMapLevel);
		m_pdfInvResolution = Vector2(1.0f / m_pdfResolution.x, 
				1.0f / m_pdfResolution.y);

		Log(EDebug, "Creating a %ix%i sampling density", 
				m_pdfResolution.x, m_pdfResolution.y);
		const Spectrum *coarseImage = m_mipmap->getImageData(mipMapLevel);
		int index = 0;
		m_pdf = DiscretePDF(m_pdfResolution.x * m_pdfResolution.y);
		for (int y=0; y<m_pdfResolution.y; ++y) {
			float sinFactor = std::sin(M_PI * (y + .5f) / m_pdfResolution.y);

			for (int x=0; x<m_pdfResolution.x; ++x)
				m_pdf[index++] = coarseImage[x + y * m_pdfResolution.x].getLuminance() * sinFactor;
		}
		m_pdfPixelSize = Vector2(2 * M_PI / m_pdfResolution.x, M_PI / m_pdfResolution.y);
		m_pdf.build();
	}

	void preprocess(const Scene *scene) {
		if (m_bsphere.isEmpty()) {
			/* Get the scene's bounding sphere and slightly enlarge it */
			m_bsphere = scene->getBSphere();
			m_bsphere.radius *= 1.01f;
		}
		if (scene->getCamera()) {
			BSphere old = m_bsphere;
			m_bsphere.expandBy(scene->getCamera()->getPosition());
			if (old != m_bsphere)
				m_bsphere.radius *= 1.01f;
		}
		m_surfaceArea = 4 * m_bsphere.radius * m_bsphere.radius * M_PI;
		m_invSurfaceArea = 1/m_surfaceArea;
	}

	Spectrum getPower() const {
		return m_average * m_surfaceArea * M_PI;
	}

	/// Sample an emission direction
	Vector sampleDirection(Point2 sample, Float &pdf, Spectrum &value) const {
		int idx = m_pdf.sampleReuse(sample.x, pdf);
		int row = idx / m_pdfResolution.x;
		int col = idx - m_pdfResolution.x * row;
		Float x = col + sample.x, y = row + sample.y;
		value = m_mipmap->triangle(0, x * m_pdfInvResolution.x, 
			y * m_pdfInvResolution.y) * m_intensityScale;
		Float theta = m_pdfPixelSize.y * y,
			  phi   = m_pdfPixelSize.x * x;

		/* Spherical-to-cartesian coordinate mapping with 
		   theta=0 => Y=1 */
		Float cosTheta = std::cos(theta),
			  sinTheta = std::sqrt(1-cosTheta*cosTheta),
			  cosPhi = std::cos(phi),
			  sinPhi = std::sin(phi);

		Vector sampledDirection(sinTheta * sinPhi, 
			cosTheta, -sinTheta*cosPhi);
		pdf = pdf / (m_pdfPixelSize.x * m_pdfPixelSize.y * sinTheta);

		return m_luminaireToWorld(-sampledDirection);
	}

	Point2 fromSphere(const Vector &d) const {
		Float u = std::atan2(d.x,-d.z) * (0.5f * INV_PI),
			  v = std::acos(std::max((Float) -1.0f, 
				  std::min((Float) 1.0f, d.y))) * INV_PI;
		if (u < 0)
			u += 1;
		return Point2(u, v);
	}

	inline Spectrum Le(const Vector &direction) const {
		Point2 uv = fromSphere(m_worldToLuminaire(direction));

		return m_mipmap->triangle(0, uv.x, uv.y)
			* m_intensityScale;
	}

	inline Spectrum Le(const Ray &ray) const {
		return Le(normalize(ray.d));
	}

	void sample(const Point &p, LuminaireSamplingRecord &lRec,
		const Point2 &sample) const {
		Vector d = sampleDirection(sample, lRec.pdf, lRec.value);

		Float nearHit, farHit;
		if (m_bsphere.contains(p) && m_bsphere.rayIntersect(Ray(p, -d, 0.0f), nearHit, farHit)) {
			lRec.sRec.p = p - d * nearHit;
			lRec.sRec.n = normalize(m_bsphere.center - lRec.sRec.p);
			lRec.d = d;
		} else {
			lRec.pdf = 0.0f;
		}
	}

	Float pdf(const Point &p, const LuminaireSamplingRecord &lRec, bool delta) const {
		const Vector d = m_worldToLuminaire(-lRec.d);	
		Point2 xy = fromSphere(d);
		xy.x *= m_pdfResolution.x;
		xy.y *= m_pdfResolution.y;
		int xPos = std::min(std::max((int) std::floor(xy.x), 0), m_pdfResolution.x-1);
		int yPos = std::min(std::max((int) std::floor(xy.y), 0), m_pdfResolution.y-1);

		Float pdf = m_pdf[xPos + yPos * m_pdfResolution.x];
		Float sinTheta = std::sqrt(std::max((Float) Epsilon, 1-d.y*d.y));

		return pdf / (m_pdfPixelSize.x * m_pdfPixelSize.y * sinTheta);
	}

	/**
	 * This is the tricky bit - we want to sample a ray that
	 * has uniform density over the set of all rays passing
	 * through the scene.
	 * For more detail, see "Using low-discrepancy sequences and 
	 * the Crofton formula to compute surface areas of geometric models"
	 * by Li, X. and Wang, W. and Martin, R.R. and Bowyer, A. 
	 * (Computer-Aided Design vol 35, #9, pp. 771--782)
	 */
	void sampleEmission(EmissionRecord &eRec, 
		const Point2 &sample1, const Point2 &sample2) const {
		Assert(eRec.type == EmissionRecord::ENormal);
		/* Chord model - generate the ray passing through two uniformly
		   distributed points on a sphere containing the scene */
		Vector d = squareToSphere(sample1);
		eRec.sRec.p = m_bsphere.center + d * m_bsphere.radius;
		eRec.sRec.n = Normal(-d);
		Point p2 = m_bsphere.center + squareToSphere(sample2) * m_bsphere.radius;
		eRec.d = p2 - eRec.sRec.p;
		Float length = eRec.d.length();

		if (length == 0) {
			eRec.value = Spectrum(0.0f);
			eRec.pdfArea = eRec.pdfDir = 1.0f;
			return;
		}

		eRec.d /= length;
		eRec.pdfArea = m_invSurfaceArea;
		eRec.pdfDir = INV_PI * dot(eRec.sRec.n, eRec.d);
		eRec.value = Le(-eRec.d);
	}

	void sampleEmissionArea(EmissionRecord &eRec, const Point2 &sample) const {
		if (eRec.type == EmissionRecord::ENormal) {
			Vector d = squareToSphere(sample);
			eRec.sRec.p = m_bsphere.center + d * m_bsphere.radius;
			eRec.sRec.n = Normal(-d);
			eRec.pdfArea = 1.0f / (4 * M_PI * m_bsphere.radius * m_bsphere.radius);
			eRec.value = Spectrum(M_PI);
		} else {
			/* Preview mode, which is more suitable for VPL-based rendering: approximate 
			   the infinitely far-away source with set of diffuse point sources */
			const Float radius = m_bsphere.radius * 1.5f;
			Vector d = squareToSphere(sample);
			eRec.sRec.p = m_bsphere.center + d * radius;
			eRec.sRec.n = Normal(-d);
			eRec.pdfArea = 1.0f / (4 * M_PI * radius * radius);
			eRec.value = Le(d) * M_PI;
		}
	}

	Spectrum sampleEmissionDirection(EmissionRecord &eRec, const Point2 &sample) const {
		Float radius = m_bsphere.radius;
		if (eRec.type == EmissionRecord::EPreview) 
			radius *= 1.5f;
		Point p2 = m_bsphere.center + squareToSphere(sample) * radius;
		eRec.d = p2 - eRec.sRec.p;
		Float length = eRec.d.length();

		if (length == 0.0f) {
			eRec.pdfDir = 1.0f;
			return Spectrum(0.0f);
		}

		eRec.d /= length;
		eRec.pdfDir = INV_PI * dot(eRec.sRec.n, eRec.d);
		if (eRec.type == EmissionRecord::ENormal)
			return Le(-eRec.d) * INV_PI;
		else
			return Spectrum(INV_PI);
	}

	void pdfEmission(EmissionRecord &eRec, bool delta) const {
		Assert(eRec.type == EmissionRecord::ENormal);
		Float dp = dot(eRec.sRec.n, eRec.d);
		if (dp > 0)
			eRec.pdfDir = delta ? 0.0f : INV_PI * dp;
		else
			eRec.pdfDir = 0;
		eRec.pdfArea = delta ? 0.0f : m_invSurfaceArea;
	}

	Spectrum evalDirection(const EmissionRecord &eRec) const {
		if (eRec.type == EmissionRecord::ENormal)
			return Le(-eRec.d) * INV_PI;
		else
			return Spectrum(INV_PI);
	}

	Spectrum evalArea(const EmissionRecord &eRec) const {
		Assert(eRec.type == EmissionRecord::ENormal);
		return Spectrum(M_PI);
	}

	bool createEmissionRecord(EmissionRecord &eRec, const Ray &ray) const {
		Float nearHit, farHit;
		if (!m_bsphere.contains(ray.o) || !m_bsphere.rayIntersect(ray, nearHit, farHit)) {
			Log(EWarn, "Could not create an emission record -- the ray "
				"in question appears to be outside of the scene bounds!");
			return false;
		}

		eRec.type = EmissionRecord::ENormal;
		eRec.sRec.p = ray(nearHit);
		eRec.sRec.n = normalize(m_bsphere.center - eRec.sRec.p);
		eRec.pdfArea = m_invSurfaceArea;
		eRec.pdfDir = INV_PI * dot(eRec.sRec.n, eRec.d);
		eRec.d = -ray.d;
		eRec.value = Le(ray.d);
		eRec.luminaire = this;
		return true;
	}

	bool isBackgroundLuminaire() const {
		return true;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "EnvMapLuminaire[" << std::endl
			<< "  name = \"" << m_name << "\"," << std::endl
			<< "  path = \"" << m_path << "\"," << std::endl
			<< "  intensityScale = " << m_intensityScale << "," << std::endl
			<< "  power = " << getPower().toString() << "," << std::endl
			<< "  bsphere = " << m_bsphere.toString() << std::endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	Spectrum m_average;
	BSphere m_bsphere;
	Float m_intensityScale;
	Float m_surfaceArea;
	Float m_invSurfaceArea;
	fs::path m_path;
	ref<MIPMap> m_mipmap;
	ref<MemoryStream> m_stream;
	DiscretePDF m_pdf;
	Vector2i m_pdfResolution;
	Vector2 m_pdfInvResolution;
	Vector2 m_pdfPixelSize;
};

// ================ Hardware shader implementation ================ 

class EnvMapLuminaireShader : public Shader {
public:
	EnvMapLuminaireShader(Renderer *renderer, const fs::path &filename, ref<Bitmap> bitmap, 
			Float intensityScale, const Transform &worldToLuminaire) : Shader(renderer, ELuminaireShader) {
		m_gpuTexture = renderer->createGPUTexture(filename.leaf(), bitmap);
		m_gpuTexture->setWrapType(GPUTexture::ERepeat);
		m_gpuTexture->setMaxAnisotropy(8);
		m_gpuTexture->init();
		/* Release the memory on the host side */
		m_gpuTexture->setBitmap(0, NULL);
		m_intensityScale = intensityScale;
		m_worldToLuminaire = worldToLuminaire;
	}

	void cleanup(Renderer *renderer) {
		m_gpuTexture->cleanup();
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_texture", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_intensityScale", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_worldToLuminaire", false));
	}

	void generateCode(std::ostringstream &oss, const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform sampler2D " << evalName << "_texture;" << endl
			<< "uniform float " << evalName << "_intensityScale;" << endl
			<< "uniform mat4 " << evalName << "_worldToLuminaire;" << endl
			<< endl
			<< "vec3 " << evalName << "_dir(vec3 wo) {" << endl
			<< "   return vec3(0.318309);" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_background(vec3 wo) {" << endl
			<< "   vec3 d = normalize((" << evalName << "_worldToLuminaire * vec4(wo, 0.0)).xyz);" << endl
			<< "   float u = atan(d.x, -d.z) * 0.15915;" << endl
			<< "   if (u < 0.0)" << endl
			<< "       u += 1.0;" << endl
			<< "   float v = acos(max(-1.0, min(1.0, d.y))) * 0.318309;" << endl
			// The following is not very elegant, but necessary to trick GLSL
			// into doing correct texture filtering across the u=0 to u=1 seam.
			<< "   if (u < 0.1)" << endl 
			<< "       return texture2D(" << evalName << "_texture, vec2(u+1.0, v)).rgb * " << evalName << "_intensityScale;" << endl
			<< "   else" << endl
			<< "       return texture2D(" << evalName << "_texture, vec2(u, v)).rgb * " << evalName << "_intensityScale;" << endl
			<< "}" << endl;
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, 
		int &textureUnitOffset) const {
		m_gpuTexture->bind(textureUnitOffset++);
		program->setParameter(parameterIDs[0], m_gpuTexture.get());
		program->setParameter(parameterIDs[1], m_intensityScale);
		program->setParameter(parameterIDs[2], m_worldToLuminaire);
	}

	void unbind() const {
		m_gpuTexture->unbind();
	}

	MTS_DECLARE_CLASS()
private:
	ref<GPUTexture> m_gpuTexture;
	Float m_intensityScale;
	Transform m_worldToLuminaire;
};

Shader *EnvMapLuminaire::createShader(Renderer *renderer) const { 
	return new EnvMapLuminaireShader(renderer, m_path, m_mipmap->getBitmap(), 
			m_intensityScale, m_worldToLuminaire);
}

MTS_IMPLEMENT_CLASS_S(EnvMapLuminaire, false, Luminaire)
MTS_IMPLEMENT_CLASS(EnvMapLuminaireShader, false, Shader)
MTS_EXPORT_PLUGIN(EnvMapLuminaire, "Environment map luminaire");
MTS_NAMESPACE_END

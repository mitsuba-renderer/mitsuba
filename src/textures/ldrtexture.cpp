/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/mipmap.h>
#include <mitsuba/hw/renderer.h>
#include <mitsuba/hw/gputexture.h>
#include <mitsuba/hw/gpuprogram.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

/**
 * Gamma-corrected bitmap texture using the JPG, PNG, TGA or BMP
 */
class LDRTexture : public Texture2D {
public:
	enum EFilterType {
		EEWAFilter = 0,
		EIsotropicFilter
	};

	LDRTexture(const Properties &props) : Texture2D(props) {
		m_filename = Thread::getThread()->getFileResolver()->resolve(
			props.getString("filename"));
		m_gamma = props.getFloat("gamma", -1); /* -1 means sRGB */
		Log(EInfo, "Loading texture \"%s\"", m_filename.leaf().c_str());

		ref<FileStream> fs = new FileStream(m_filename, FileStream::EReadOnly);
		std::string extension = boost::to_lower_copy(m_filename.extension());

		std::string filterType = props.getString("filterType", "ewa");
		std::string wrapMode = props.getString("wrapMode", "repeat");

		if (filterType == "ewa")
			m_anisotropic = true;
		else if (filterType == "isotropic")
			m_anisotropic = false;
		else
			Log(EError, "Unknown filter type '%s' -- must be "
				"'ewa' or 'isotropic'!", filterType.c_str());

		if (wrapMode == "repeat")
			m_wrapMode = MIPMap::ERepeat;
		else if (wrapMode == "clamp")
			m_wrapMode = MIPMap::EClamp;
		else if (wrapMode == "black")
			m_wrapMode = MIPMap::EBlack;
		else if (wrapMode == "white")
			m_wrapMode = MIPMap::EWhite;
		else
			Log(EError, "Unknown wrap mode '%s' -- must be "
				"'repeat', 'clamp', 'black', or 'white'!", filterType.c_str());
	
		m_maxAnisotropy = props.getFloat("maxAnisotropy", 8);

		if (extension == ".jpg" || extension == ".jpeg")
			m_format = Bitmap::EJPEG;
		else if (extension == ".png")
			m_format = Bitmap::EPNG;
		else if (extension == ".tga")
			m_format = Bitmap::ETGA;
		else if (extension == ".bmp")
			m_format = Bitmap::EBMP;
		else
			Log(EError, "Cannot deduce the file type of '%s'!", m_filename.file_string().c_str());

		ref<Bitmap> bitmap = new Bitmap(m_format, fs);
		initializeFrom(bitmap);
	}

	LDRTexture(Stream *stream, InstanceManager *manager) 
	 : Texture2D(stream, manager) {
		m_filename = stream->readString();
		Log(EInfo, "Unserializing texture \"%s\"", m_filename.leaf().c_str());
		m_gamma = stream->readFloat();
		m_format = static_cast<Bitmap::EFileFormat>(stream->readInt());
		m_anisotropic = stream->readBool();
		m_wrapMode = (MIPMap::EWrapMode) stream->readBool();
		m_maxAnisotropy = stream->readFloat();
		unsigned int size = stream->readUInt();
		ref<MemoryStream> mStream = new MemoryStream(size);
		stream->copyTo(mStream, size);
		mStream->setPos(0);
		ref<Bitmap> bitmap = new Bitmap(m_format, mStream);
		initializeFrom(bitmap);

		if (Scheduler::getInstance()->hasRemoteWorkers()
			&& !fs::exists(m_filename)) {
			/* This code is running on a machine different from
			   the one that created the stream. Because we might
			   later have to handle a call to serialize(), the
			   whole bitmap must be kept in memory */
			m_stream = mStream;
		}
	}

	inline Float fromSRGBComponent(Float value) {
		if (value <= (Float) 0.04045)
			return value / (Float) 12.92;
		return std::pow((value + (Float) 0.055)
			/ (Float) (1.0 + 0.055), (Float) 2.4);
	}

	void initializeFrom(Bitmap *bitmap) {
		ref<Bitmap> corrected = new Bitmap(bitmap->getWidth(), bitmap->getHeight(), 128);

		float tbl[256];
		if (m_gamma == -1) {
			for (int i=0; i<256; ++i) 
				tbl[i] = fromSRGBComponent(i/255.0f);
		} else {
			for (int i=0; i<256; ++i)
				tbl[i] = std::pow(i/255.0f, m_gamma);
		}

		uint8_t *data = bitmap->getData();
		float *flData = corrected->getFloatData();
		if (bitmap->getBitsPerPixel() == 32) {
			for (int y=0; y<bitmap->getHeight(); ++y) {
				for (int x=0; x<bitmap->getWidth(); ++x) {
					float
						r = tbl[*data++],
						g = tbl[*data++],
						b = tbl[*data++],
						a = *data++ / 255.0f;
					*flData++ = r;
					*flData++ = g;
					*flData++ = b;
					*flData++ = a;
				}
			}
		} else if (bitmap->getBitsPerPixel() == 24) {
			for (int y=0; y<bitmap->getHeight(); ++y) {
				for (int x=0; x<bitmap->getWidth(); ++x) {
					float
						r = tbl[*data++],
						g = tbl[*data++],
						b = tbl[*data++];
					*flData++ = r;
					*flData++ = g;
					*flData++ = b;
					*flData++ = 1.0f;
				}
			}
		} else if (bitmap->getBitsPerPixel() == 16) {
			for (int y=0; y<bitmap->getHeight(); ++y) {
				for (int x=0; x<bitmap->getWidth(); ++x) {
					float col = tbl[*data++],
						a = *data++ / 255.0f;
					*flData++ = col;
					*flData++ = col;
					*flData++ = col;
					*flData++ = a;
				}
			}
		} else if (bitmap->getBitsPerPixel() == 8) {
			for (int y=0; y<bitmap->getHeight(); ++y) {
				for (int x=0; x<bitmap->getWidth(); ++x) {
					float col = tbl[*data++];
					*flData++ = col;
					*flData++ = col;
					*flData++ = col;
					*flData++ = 1.0f;
				}
			}
		} else {
			Log(EError, "%i bpp images are currently not supported!", bitmap->getBitsPerPixel());
		}

		m_mipmap = MIPMap::fromBitmap(corrected, !m_anisotropic,
				m_wrapMode, m_maxAnisotropy);
		m_average = m_mipmap->triangle(m_mipmap->getLevels()-1, 0, 0);
		m_maximum = m_mipmap->getMaximum();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture2D::serialize(stream, manager);
		stream->writeString(m_filename.file_string());
		stream->writeFloat(m_gamma);
		stream->writeInt(m_format);
		stream->writeBool(m_anisotropic);
		stream->writeUInt(m_wrapMode);
		stream->writeFloat(m_maxAnisotropy);
		if (m_stream.get()) {
			stream->writeUInt((unsigned int) m_stream->getSize());
			stream->write(m_stream->getData(), m_stream->getSize());
		} else {
			ref<Stream> mStream = new MemoryStream();
			ref<Stream> is = new FileStream(m_filename, FileStream::EReadOnly);
			stream->writeUInt((unsigned int) is->getSize());
			is->copyTo(stream);
		}
	}

	Spectrum getValue(const Point2 &uv) const {
		return m_mipmap->triangle(0, uv.x, uv.y);
	}

	Spectrum getValue(const Point2 &uv, Float dudx, 
			Float dudy, Float dvdx, Float dvdy) const {
		return m_mipmap->getValue(uv.x, uv.y, dudx, dudy, dvdx, dvdy);
	}

	Spectrum getAverage() const {
		return m_average;
	}

	Spectrum getMaximum() const {
		return m_maximum;
	}

	bool usesRayDifferentials() const {
		return true;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "LDRTexture[filename=\"" << m_filename << "\", gamma=" << m_gamma << "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
protected:
	ref<MIPMap> m_mipmap;
	ref<MemoryStream> m_stream;
	fs::path m_filename;
	Bitmap::EFileFormat m_format;
	Spectrum m_average, m_maximum;
	Float m_gamma;
	bool m_anisotropic;
	MIPMap::EWrapMode m_wrapMode;
	Float m_maxAnisotropy;
};

// ================ Hardware shader implementation ================ 

class LDRTextureShader : public Shader {
public:
	LDRTextureShader(Renderer *renderer, std::string filename, ref<Bitmap> bitmap,
			const Point2 &uvOffset, const Vector2 &uvScale, MIPMap::EWrapMode wrapMode, 
			Float maxAnisotropy) 
		: Shader(renderer, ETextureShader), m_uvOffset(uvOffset), m_uvScale(uvScale) {
		m_gpuTexture = renderer->createGPUTexture(filename, bitmap);
		if (wrapMode == MIPMap::ERepeat)
			m_gpuTexture->setWrapType(GPUTexture::ERepeat);
		else
			m_gpuTexture->setWrapType(GPUTexture::EClampToEdge);
		m_gpuTexture->setMaxAnisotropy((int) maxAnisotropy);
		m_gpuTexture->init();
		/* Release the memory on the host side */
		m_gpuTexture->setBitmap(0, NULL);
	}

	void cleanup(Renderer *renderer) {
		m_gpuTexture->cleanup();
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform sampler2D " << evalName << "_texture;" << endl
			<< "uniform vec2 " << evalName << "_uvOffset;" << endl
			<< "uniform vec2 " << evalName << "_uvScale;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv) {" << endl
			<< "    return texture2D(" << evalName << "_texture, vec2(" << endl
			<< "          uv.x * " << evalName << "_uvScale.x + " << evalName << "_uvOffset.x," << endl 
			<< "          uv.y * " << evalName << "_uvScale.y + " << evalName << "_uvOffset.y)).rgb;" << endl 
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_texture", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_uvOffset", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_uvScale", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, 
		int &textureUnitOffset) const {
		m_gpuTexture->bind(textureUnitOffset++);
		program->setParameter(parameterIDs[0], m_gpuTexture.get());
		program->setParameter(parameterIDs[1], m_uvOffset);
		program->setParameter(parameterIDs[2], m_uvScale);
	}

	void unbind() const {
		m_gpuTexture->unbind();
	}

	MTS_DECLARE_CLASS()
private:
	ref<GPUTexture> m_gpuTexture;
	Point2 m_uvOffset;
	Vector2 m_uvScale;
};

Shader *LDRTexture::createShader(Renderer *renderer) const {
	return new LDRTextureShader(renderer, m_filename.leaf(), 
			m_mipmap->getLDRBitmap(), m_uvOffset, m_uvScale,
			m_wrapMode, m_anisotropic ? m_maxAnisotropy : 1.0f);
}

MTS_IMPLEMENT_CLASS_S(LDRTexture, false, Texture2D)
MTS_IMPLEMENT_CLASS(LDRTextureShader, false, Shader)
MTS_EXPORT_PLUGIN(LDRTexture, "LDR texture (JPG/PNG/TGA/BMP)");
MTS_NAMESPACE_END

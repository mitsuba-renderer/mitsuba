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

/*!\plugin{bitmap}{Bitmap texture}
 * \parameters{
 *     \parameter{filename}{\String}{
 *       Filename of the bitmap to be loaded
 *     }
 *     \parameter{gamma}{\Float}{
 *       Gamma value of the source bitmap file 
 *       \default{\emph{automatic}, i.e. linear for EXR input, 
 *       and sRGB for everything else.}
 *     }
 *     \parameter{filterType}{\String}{
 *       Specifies the texture filturing that should be used for lookups
 *       \begin{enumerate}[(i)]
 *           \item \code{ewa}: Elliptically weighted average (a.k.a.
 *           anisotropic filtering). This produces the best quality
 *           \item \code{trilinear}: Simple trilinear (isotropic) filtering.
 *           \item \code{none}: No filtering, do nearest neighbor lookups.
 *       \end{enumerate}
 *       Default: \code{ewa}.
 *     }
 *     \parameter{wrapMode}{\String}{
 *       This parameter defines the behavior of the texture outside of the $[0,1]$ $uv$ range.
 *       \begin{enumerate}[(i)]
 *           \item \code{repeat}: Repeat the texture (i.e. $uv$ coordinates
 *           are taken modulo 2)
 *           \item \code{clamp}: Clamp $uv$ coordinates to $[0,1]$
 *           \item \code{black}: Switch to a zero-valued texture
 *           \item \code{white}: Switch to a one-valued texture
 *       \end{enumerate}
 *       Default: \code{repeat}.
 *     }
 *     \parameter{maxAnisotropy}{\Float}{
 *        Specifies an upper limit on the amount of anisotropy 
 *        of \code{ewa} lookups\default{8}
 *     }
 *     \parameter{uscale, vscale}{\Float}{
 *       Multiplicative factors that should be applied to UV values before a lookup
 *     }
 *     \parameter{uoffset, voffset}{\Float}{
 *       Numerical offset that should be applied to UV values before a lookup
 *     }
 * }
 * This plugin implements a bitmap-based texture, which supports the following
 * file formats:
 * \begin{itemize}
 *     \item OpenEXR
 *     \item JPEG
 *     \item PNG (Portable Network Graphics)
 *     \item TGA (Targa)
 *     \item BMP (Windows bitmaps)
 * \end{itemize}
 *
 * The plugin internally converts all bitmap data into a \emph{linear} space to ensure
 * a proper workflow.
 */

class BitmapTexture : public Texture2D {
public:
	BitmapTexture(const Properties &props) : Texture2D(props) {
		m_filename = Thread::getThread()->getFileResolver()->resolve(
			props.getString("filename"));

		/* -1 means sRGB. Gamma is ignored when loading EXR files */
		m_gamma = props.getFloat("gamma", -1);
		Log(EInfo, "Loading texture \"%s\"", m_filename.leaf().c_str());

		ref<FileStream> fs = new FileStream(m_filename, FileStream::EReadOnly);
		std::string extension = boost::to_lower_copy(m_filename.extension());
		std::string filterType = boost::to_lower_copy(props.getString("filterType", "ewa"));
		std::string wrapMode = boost::to_lower_copy(props.getString("wrapMode", "repeat"));

		if (filterType == "ewa")
			m_filterType = MIPMap::EEWA;
		else if (filterType == "trilinear")
			m_filterType = MIPMap::ETrilinear;
		else if (filterType == "none")
			m_filterType = MIPMap::ENone;
		else
			Log(EError, "Unknown filter type '%s' -- must be "
				"'ewa', 'isotropic', or 'none'!", filterType.c_str());

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

		if (extension == ".exr")
			m_format = Bitmap::EEXR;
		else if (extension == ".jpg" || extension == ".jpeg")
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

	BitmapTexture(Stream *stream, InstanceManager *manager) 
	 : Texture2D(stream, manager) {
		m_filename = stream->readString();
		Log(EInfo, "Unserializing texture \"%s\"", m_filename.leaf().c_str());
		m_gamma = stream->readFloat();
		m_format = static_cast<Bitmap::EFileFormat>(stream->readInt());
		m_filterType = (MIPMap::EFilterType) stream->readInt();
		m_wrapMode = (MIPMap::EWrapMode) stream->readUInt();
		m_maxAnisotropy = stream->readFloat();
		uint32_t size = stream->readUInt();
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
			m_stream->setPos(0);
		}
	}

	inline Float fromSRGBComponent(Float value) {
		if (value <= (Float) 0.04045)
			return value / (Float) 12.92;
		return std::pow((value + (Float) 0.055)
			/ (Float) (1.0 + 0.055), (Float) 2.4);
	}

	void initializeFrom(Bitmap *bitmap) {
		ref<Bitmap> corrected;
		m_bpp = bitmap->getBitsPerPixel();
		if (bitmap->getBitsPerPixel() == 128) {
			/* Nothing needs to be done */
			corrected = bitmap;
		} else {
			corrected = new Bitmap(bitmap->getWidth(), bitmap->getHeight(), 128);

			float tbl[256];
			if (m_gamma == -1) {
				for (int i=0; i<256; ++i) 
					tbl[i] = fromSRGBComponent((Float) i / (Float) 255);
			} else {
				for (int i=0; i<256; ++i)
					tbl[i] = std::pow((Float) i / (Float) 255, m_gamma);
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
			} else if (bitmap->getBitsPerPixel() == 1) {
				int pos = 0;
				for (int y=0; y<bitmap->getHeight(); ++y) {
					for (int x=0; x<bitmap->getWidth(); ++x) {
						int entry = pos / 8;
						int bit   = pos % 8;
						int value = (data[entry] & (1 << bit)) ? 255 : 0;
						float col = tbl[value];
						*flData++ = col;
						*flData++ = col;
						*flData++ = col;
						*flData++ = 1.0f;
						pos++;
					}
				}
			} else {
				Log(EError, "%i bpp images are currently not supported!", bitmap->getBitsPerPixel());
			}
		}

		m_mipmap = MIPMap::fromBitmap(corrected, m_filterType,
				m_wrapMode, m_maxAnisotropy);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture2D::serialize(stream, manager);
		stream->writeString(m_filename.file_string());
		stream->writeFloat(m_gamma);
		stream->writeInt(m_format);
		stream->writeInt(m_filterType);
		stream->writeUInt(m_wrapMode);
		stream->writeFloat(m_maxAnisotropy);

		if (m_stream.get()) {
			stream->writeUInt((uint32_t) m_stream->getSize());
			stream->write(m_stream->getData(), m_stream->getSize());
		} else {
			ref<Stream> mStream = new MemoryStream();
			ref<Stream> is = new FileStream(m_filename, FileStream::EReadOnly);
			stream->writeUInt((uint32_t) is->getSize());
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
		return m_mipmap->getAverage();
	}

	Spectrum getMaximum() const {
		return m_mipmap->getMaximum();
	}

	Spectrum getMinimum() const {
		return m_mipmap->getMinimum();
	}

	bool isConstant() const {
		return false;
	}

	bool usesRayDifferentials() const {
		return true;
	}

	Vector3i getResolution() const {
		return Vector3i(
			m_mipmap->getWidth(),
			m_mipmap->getHeight(),
			1
		);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "BitmapTexture[" << endl
			<< "  filename = \"" << m_filename << "\"," << endl
			<< "  bpp = " << m_bpp;
		if (m_bpp < 128) {
			oss << "," << endl
				<< "  gamma = " << m_gamma << endl;
		} else {
			oss << endl;
		}
		oss << "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
protected:
	ref<MIPMap> m_mipmap;
	ref<MemoryStream> m_stream;
	fs::path m_filename;
	Bitmap::EFileFormat m_format;
	MIPMap::EFilterType m_filterType;
	MIPMap::EWrapMode m_wrapMode;
	Float m_maxAnisotropy;
	Float m_gamma;
	int m_bpp;
};

// ================ Hardware shader implementation ================ 

class BitmapTextureShader : public Shader {
public:
	BitmapTextureShader(Renderer *renderer, std::string filename, ref<Bitmap> bitmap,
			const Point2 &uvOffset, const Vector2 &uvScale, MIPMap::EWrapMode wrapMode, 
			Float maxAnisotropy) 
		: Shader(renderer, ETextureShader), m_uvOffset(uvOffset), m_uvScale(uvScale) {
		m_gpuTexture = renderer->createGPUTexture(filename, bitmap);
		if (wrapMode == MIPMap::ERepeat)
			m_gpuTexture->setWrapType(GPUTexture::ERepeat);
		else
			m_gpuTexture->setWrapType(GPUTexture::EClampToEdge);
		m_gpuTexture->setMaxAnisotropy(maxAnisotropy);
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

Shader *BitmapTexture::createShader(Renderer *renderer) const {
	return new BitmapTextureShader(renderer, m_filename.leaf(),
			m_mipmap->getLDRBitmap(), m_uvOffset, m_uvScale,
			m_wrapMode, (m_filterType == MIPMap::EEWA)
			? m_maxAnisotropy : 1.0f);
}

MTS_IMPLEMENT_CLASS_S(BitmapTexture, false, Texture2D)
MTS_IMPLEMENT_CLASS(BitmapTextureShader, false, Shader)
MTS_EXPORT_PLUGIN(BitmapTexture, "Bitmap texture (EXR/JPG/PNG/TGA/BMP)");
MTS_NAMESPACE_END

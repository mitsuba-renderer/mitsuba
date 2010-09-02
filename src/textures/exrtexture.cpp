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
#include <mitsuba/core/mstream.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/mipmap.h>

MTS_NAMESPACE_BEGIN

/**
 * Simple linear (i.e. not gamma corrected) bitmap texture
 * using the EXR file format
 */
class EXRTexture : public Texture {
public:
	EXRTexture(const Properties &props) : Texture(props) {
		m_filename = props.getString("filename");
		m_filename = FileResolver::getInstance()->resolve(m_filename);
		Log(EInfo, "Loading texture \"%s\"", m_filename.c_str());

		ref<FileStream> fs = new FileStream(m_filename, FileStream::EReadOnly);
		ref<Bitmap> bitmap = new Bitmap(Bitmap::EEXR, fs);
		m_mipmap = MIPMap::fromBitmap(bitmap);
		m_average = m_mipmap->triangle(m_mipmap->getLevels()-1, 0, 0);
		m_maximum = m_mipmap->getMaximum();
	}
	
	EXRTexture(Stream *stream, InstanceManager *manager) 
	 : Texture(stream, manager) {
		m_filename = stream->readString();
		Log(EInfo, "Unserializing texture \"%s\"", m_filename.c_str());
		int size = stream->readInt();
		ref<MemoryStream> mStream = new MemoryStream(size);
		stream->copyTo(mStream, size);
		mStream->setPos(0);
		ref<Bitmap> bitmap = new Bitmap(Bitmap::EEXR, mStream);
		m_mipmap = MIPMap::fromBitmap(bitmap);
		m_average = m_mipmap->triangle(m_mipmap->getLevels()-1, 0, 0);
		m_maximum = m_mipmap->getMaximum();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture::serialize(stream, manager);
		stream->writeString(m_filename);
		ref<Stream> is = new FileStream(m_filename, FileStream::EReadOnly);
		stream->writeInt(is->getSize());
		is->copyTo(stream);
	}

	Spectrum getValue(const Intersection &its) const {
		return m_mipmap->getValue(its);
	}

	Spectrum getMaximum() const {
		return m_maximum;
	}

	Spectrum getAverage() const {
		return m_average;
	}

	bool usesRayDifferentials() const {
		return true;
	}
	
	std::string toString() const {
		std::ostringstream oss;
		oss << "EXRTexture[filename=\"" << m_filename << "\"]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
	ref<MIPMap> m_mipmap;
	std::string m_filename;
	Spectrum m_average, m_maximum;
};

MTS_IMPLEMENT_CLASS_S(EXRTexture, false, Texture)
MTS_EXPORT_PLUGIN(EXRTexture, "HDR texture (EXR)");
MTS_NAMESPACE_END

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
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/mipmap.h>

MTS_NAMESPACE_BEGIN

/**
 * Simple linear (i.e. not gamma corrected) bitmap texture
 * using the EXR file format
 */
class EXRTexture : public Texture2D {
public:
	EXRTexture(const Properties &props) : Texture2D(props) {
		m_filename = Thread::getThread()->getFileResolver()->resolve(
			props.getString("filename"));
		Log(EInfo, "Loading texture \"%s\"", m_filename.leaf().c_str());

		ref<FileStream> fs = new FileStream(m_filename, FileStream::EReadOnly);
		ref<Bitmap> bitmap = new Bitmap(Bitmap::EEXR, fs);
		m_mipmap = MIPMap::fromBitmap(bitmap);
		m_average = m_mipmap->triangle(m_mipmap->getLevels()-1, 0, 0);
		m_maximum = m_mipmap->getMaximum();
	}

	EXRTexture(Stream *stream, InstanceManager *manager) 
	 : Texture2D(stream, manager) {
		m_filename = stream->readString();
		Log(EInfo, "Unserializing texture \"%s\"", m_filename.leaf().c_str());
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
		Texture2D::serialize(stream, manager);
		stream->writeString(m_filename.file_string());
		ref<Stream> is = new FileStream(m_filename, FileStream::EReadOnly);
		stream->writeInt(is->getSize());
		is->copyTo(stream);
	}

	Spectrum getValue(const Point2 &uv) const {
		return m_mipmap->triangle(0, uv.x, uv.y);
	}

	Spectrum getValue(const Point2 &uv, Float dudx, 
			Float dudy, Float dvdx, Float dvdy) const {
		return m_mipmap->getValue(uv.x, uv.y, dudx, dudy, dvdx, dvdy);
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
		oss << "EXRTexture[filename=\"" << m_filename.file_string() << "\"]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
	ref<MIPMap> m_mipmap;
	fs::path m_filename;
	Spectrum m_average, m_maximum;
};

MTS_IMPLEMENT_CLASS_S(EXRTexture, false, Texture2D)
MTS_EXPORT_PLUGIN(EXRTexture, "HDR texture (EXR)");
MTS_NAMESPACE_END

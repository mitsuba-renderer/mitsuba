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

#if !defined(__BITMAP_H)
#define __BITMAP_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/stream.h>

MTS_NAMESPACE_BEGIN

/** \brief 1/8/24/32/128-Bit Raster ("Bitmap") data structure with
 * support for PNG storage and retrieval. When set to 128 bits per pixel,
 * the implementation switches to HDR and the EXR file format.
 *
 * This class can efficiently handle 1-bit masks
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE Bitmap : public Object {
public:
	enum EFileFormat {
		EPNG = 0,
		EEXR,
		ETGA,
		EBMP,
		EJPEG
	};

	/// Create a new bitmap
	Bitmap(int width = 512, int height = 512, int bpp = 24);

	/// Load a bitmap of the given file format
	Bitmap(EFileFormat format, Stream *stream);

	/// Create a copy of this bitmap
	Bitmap *clone() const;

	/// Bitmap equality operator for unit-tests etc.
	bool operator==(const Bitmap &bitmap) const;

	/// Clear the bitmap to zero
	void clear();


	/**
	 * Save the bitmap using the give file format. Where
	 * applicable, the \a compression parameter can be used
	 * to control amount of compression (1: lowest, 9: highest)
	 * */
	void save(EFileFormat format, Stream *stream, int compression = 5) const;

	/// Return the image's title identifier
	inline const std::string &getTitle() const { return m_title; }
	
	/// Return the image's author identifier
	inline const std::string &getAuthor() const { return m_author; }

	/// Return the image's comment identifier
	inline const std::string &getComment() const { return m_comment; }

	/// Return the image's gamma identifier (-1: sRGB)
	inline Float getGamma() const { return m_gamma; }

	/// Set the image's title identifier
	inline void setTitle(const std::string &title) { m_title = title; }

	/// Set the image's author identifier
	inline void setAuthor(const std::string &author) { m_author = author; }

	/// Set the image's comment identifier
	inline void setComment(const std::string &comment) { m_comment = comment; }

	/// Set the image's gamma identifier (-1: sRGB)
	inline void setGamma(Float gamma) { m_gamma = gamma; }

	/**
	 * \brief Access the underlying raster
	 * \remark This function is not exposed in the Python bindings
	 */
	inline unsigned char *getData() { return m_data; }
	
	/**
	 * \brief Access the underlying raster
	 * \remark This function is not exposed in the Python bindings
	 */
	inline const unsigned char *getData() const { return m_data; }

	/**
	 * \brief Access the underlying raster (128bpp images)
	 * \remark This function is not exposed in the Python bindings
	 */
	inline float *getFloatData() { return (float *) m_data; }

	/**
	 * \brief Access the underlying raster (128bpp images)
	 * \remark This function is not exposed in the Python bindings
	 */
	inline const float *getFloatData() const { return (const float *) m_data; }

	/// Return the bitmap width
	inline int getWidth() const { return m_width; }

	/// Return the bitmap height 
	inline int getHeight() const { return m_height; }

	/// Return the bitmap's bits per pixel
	inline int getBitsPerPixel() const { return m_bpp; }

	/// Return the bitmap size in bytes
	inline size_t getSize() const { return m_size; }

	/// Return some human-readable information about this bitmap
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Constructor for subclasses
	inline Bitmap(bool unused) {
	}

	/// Virtual destructor
	virtual ~Bitmap();

	/// Load a file stored using the PNG file format
	void loadPNG(Stream *stream);
	
	/// Load a file stored using the TGA file format
	void loadTGA(Stream *stream);
	
	/// Load a file stored using the BMP file format
	void loadBMP(Stream *stream);

	/// Load a file stored using the JPEG file format
	void loadJPEG(Stream *stream);

	/// Load a file stored using the EXR file format
	void loadEXR(Stream *stream);

	/// Save a file using the PNG file format
	void savePNG(Stream *stream, int compression) const;

	/// Save a file using the JPEG file format
	void saveJPEG(Stream *stream, int quality = 100) const;

	/// Save a file using the EXR file format
	void saveEXR(Stream *stream) const;
protected:
	int m_width;
	int m_height;
	int m_bpp;
	size_t m_size;
	unsigned char *m_data;
	std::string m_title;
	std::string m_author;
	std::string m_comment;
	Float m_gamma;
};

MTS_NAMESPACE_END

#endif /* __BITMAP_H */

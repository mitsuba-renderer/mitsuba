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

	/// Save the bitmap using the give file format
	void save(EFileFormat format, Stream *stream) const;

	/// Return the image's title identifier
	inline const std::string &getTile() const { return m_title; }
	
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

	/// Access the underlying raster
	inline unsigned char *getData() { return m_data; }
	
	/// Access the underlying bit raster
	inline const unsigned char *getData() const { return m_data; }

	/// Access the underlying raster (only meant for 128bpp images)
	inline float *getFloatData() { return (float *) m_data; }

	/// Access the underlying raster (only meant for 128bpp images)
	inline const float *getFloatData() const { return (const float *) m_data; }

	/// Return the bitmap width
	inline const int getWidth() const { return m_width; }

	/// Return the bitmap height 
	inline const int getHeight() const { return m_height; }

	/// Return the bitmap's bits per pixel
	inline const int getBitsPerPixel() const { return m_bpp; }

	/// Return the bitmap size in bytes
	inline const int getSize() const { return m_size; }

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
	void savePNG(Stream *stream) const;

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

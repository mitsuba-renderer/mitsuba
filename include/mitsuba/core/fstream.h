#if !defined(__FSTREAM_H)
#define __FSTREAM_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/stream.h>

MTS_NAMESPACE_BEGIN

/** \brief Simple Stream implementation which can be used to access files.
 *
 * This class uses streams on posix platforms and the native
 * WIN32 API when used on windows.
 */
class MTS_EXPORT_CORE FileStream : public Stream {
public:
	/// File opening modes
	enum EFileMode {
		/// rb
		EReadOnly = 0,
		/// rb+
		EReadWrite,
		/// wb
		ETruncWrite,
		/// wb+
		ETruncReadWrite,
		/// ab
		EAppendWrite,
		/// ab+
		EAppendReadWrite
	};

	/// Create a file stream class with no file open
	FileStream();

	/// Create a file stream class and open a file with a given EFileMode
	FileStream(const std::string &filename, EFileMode mode = EReadOnly);

	/// Return the file name
	inline const std::string &getFileName() const { return m_filename; }

	/// Check whether a file exists
	static bool exists(const std::string &filename);

	/// Open a file with a given open mode
	void open(const std::string &filename, EFileMode mode = EReadOnly);

	/// Close the current file
	void close();

	/// Remove the current file
	void remove();

	/// Return a string representation
	std::string toString() const;

	/* Stream implementation */
	void read(void *ptr, size_t size);
	void write(const void *ptr, size_t size);
	void setPos(size_t pos);
	size_t getPos() const;
	size_t getSize() const;
	void truncate(size_t size);
	void flush();
	bool canWrite() const;
	bool canRead() const;

	MTS_DECLARE_CLASS()
protected:
	/** \brief Virtual destructor
	 *
	 * The destructor frees all resources and closes
	 * the file if it is still open
	 */
	virtual ~FileStream();
protected:
#ifdef WIN32
	HANDLE m_file;
#else
	FILE* m_file;
#endif
	bool m_write;
	bool m_read;
	EFileMode m_mode;
	std::string m_filename;
};

MTS_NAMESPACE_END

#endif /* __FSTREAM_H */

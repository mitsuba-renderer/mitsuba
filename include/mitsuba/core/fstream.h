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

#if !defined(__FSTREAM_H)
#define __FSTREAM_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/stream.h>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

MTS_NAMESPACE_BEGIN

/** \brief Simple \ref Stream implementation for accessing files.
 *
 * This class uses POSIX streams on Linux and OSX and the native
 * WIN32 API when used on Windows.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE FileStream : public Stream {
public:
	/// Supported file opening modes
	enum EFileMode {
		EReadOnly = 0,   ///< rb
		EReadWrite,      ///< rb+
		ETruncWrite,     ///< wb
		ETruncReadWrite, ///< wb+
		EAppendWrite,    ///< ab
		EAppendReadWrite ///< ab+
	};

	// =============================================================
	//! @{ \name Constructors
	// =============================================================

	/// Create a file stream class with no file open
	FileStream();

	/// Create a file stream class and open a file with a given EFileMode
	explicit FileStream(const fs::path &path, EFileMode mode = EReadOnly);

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name File-specific features
	// =============================================================

	/// Return the file path
	inline const fs::path &getPath() const { return m_path; }

	/// Open a file with a given open mode
	void open(const fs::path &filename, EFileMode mode = EReadOnly);

	/// Close the current file
	void close();

	/// Remove the current file
	void remove();

	/// Return a string representation
	std::string toString() const;

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Stream interface
	// =============================================================

	void read(void *ptr, size_t size);
	void write(const void *ptr, size_t size);
	void setPos(size_t pos);
	size_t getPos() const;
	size_t getSize() const;
	void truncate(size_t size);
	void flush();
	bool canWrite() const;
	bool canRead() const;

	//! @}
	// =============================================================

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
	fs::path m_path;
};

MTS_NAMESPACE_END

#endif /* __FSTREAM_H */

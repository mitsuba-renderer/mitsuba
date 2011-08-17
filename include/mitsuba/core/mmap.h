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

#if !defined(__MMAP_H)
#define __MMAP_H

#include <mitsuba/core/fstream.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Basic cross-platform abstraction for memory mapped files
 * \ingroup libcore
 */
class MTS_EXPORT_CORE MemoryMappedFile : public Object {
public:
	/// Map the specified file into memory
	MemoryMappedFile(const fs::path &filename);

	/// Return a pointer to the file contents in memory
	inline void *getData() { return m_data; }

	/// Return a pointer to the file contents in memory (const version)
	inline const void *getData() const { return m_data; }

	/// Return the size of the mapped region
	inline size_t getSize() const { return m_size; }

	/// Release all resources
	virtual ~MemoryMappedFile();

	MTS_DECLARE_CLASS()
protected:
	fs::path m_filename;

#if defined(WIN32)
	HANDLE m_file;
	HANDLE m_fileMapping;
#endif
	size_t m_size;
	void *m_data;
};

MTS_NAMESPACE_END

#endif /* __MMAP_H */

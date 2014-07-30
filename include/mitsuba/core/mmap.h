/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

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

#pragma once
#if !defined(__MITSUBA_CORE_MMAP_H_)
#define __MITSUBA_CORE_MMAP_H_

#include <mitsuba/core/fstream.h>
#include <boost/scoped_ptr.hpp>

MTS_NAMESPACE_BEGIN

/**
 * \brief Basic cross-platform abstraction for memory mapped files
 * \ingroup libcore
 */
class MTS_EXPORT_CORE MemoryMappedFile : public Object {
public:
	/// Create a new memory-mapped file of the specified size
	MemoryMappedFile(const fs::path &filename, size_t size);

	/// Map the specified file into memory
	MemoryMappedFile(const fs::path &filename, bool readOnly = true);

	/// Return a pointer to the file contents in memory
	void *getData();

	/// Return a pointer to the file contents in memory (const version)
	const void *getData() const;

	/// Return the size of the mapped region
	size_t getSize() const;

	/**
	 * \brief Resize the memory-mapped file
	 *
	 * This involves remapping the file, which will
	 * generally change the pointer obtained via getData()
	 */
	void resize(size_t size);

	/// Return the associated filename
	const fs::path &getFilename() const;

	/// Return whether the mapped memory region is read-only
	bool isReadOnly() const;

	/// Return a string representation
	std::string toString() const;

	/**
	 * \brief Create a temporary memory-mapped file
	 *
	 * \remark When closing the mapping, the file is
	 * automatically deleted.
	 */
	static ref<MemoryMappedFile> createTemporary(size_t size);

	MTS_DECLARE_CLASS()
protected:
	/// Internal constructor
	MemoryMappedFile();

	/// Release all resources
	virtual ~MemoryMappedFile();
private:
	struct MemoryMappedFilePrivate;
	boost::scoped_ptr<MemoryMappedFilePrivate> d;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_MMAP_H_ */

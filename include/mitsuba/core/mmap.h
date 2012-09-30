/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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
	/// Map the specified file into memory
	MemoryMappedFile(const fs::path &filename);

	/// Create a new memory-mapped file of the specified size
	MemoryMappedFile(const fs::path &filename, size_t size);

	/// Return a pointer to the file contents in memory
	void *getData();

	/// Return a pointer to the file contents in memory (const version)
	const void *getData() const;

	/// Return the size of the mapped region
	size_t getSize() const;

	/// Release all resources
	virtual ~MemoryMappedFile();

	MTS_DECLARE_CLASS()
private:
	struct MemoryMappedFilePrivate;
	boost::scoped_ptr<MemoryMappedFilePrivate> d;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_MMAP_H_ */

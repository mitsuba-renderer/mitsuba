#include <mitsuba/core/mmap.h>

#if defined(__LINUX__) || defined(__OSX__)
#include <sys/mman.h>
#include <fcntl.h>
#endif

MTS_NAMESPACE_BEGIN

MemoryMappedFile::MemoryMappedFile(const fs::path &filename) : m_filename(filename), m_data(NULL) {
	if (!fs::exists(filename))
		Log(EError, "The file \"%s\" does not exist!", filename.file_string().c_str());
	m_size = (size_t) fs::file_size(filename);
	Log(ETrace, "Mapping \"%s\" into memory (%s)..", 
			filename.filename().c_str(), memString(m_size).c_str());
#if defined(__LINUX__) || defined(__OSX__)
	int fd = open(filename.file_string().c_str(), O_RDONLY);
	if (fd == -1)
		Log(EError, "Could not open \"%s\"!", m_filename.file_string().c_str());
	m_data = mmap(NULL, m_size, PROT_READ, MAP_SHARED, fd, 0);
	if (m_data == NULL)
		Log(EError, "Could not map \"%s\" to memory!", m_filename.file_string().c_str());
	if (close(fd) != 0)
		Log(EError, "close(): unable to close file!");
#elif defined(WIN32)
	m_file = CreateFile(filename.file_string().c_str(), GENERIC_READ, 
		FILE_SHARE_READ, NULL, OPEN_EXISTING, 
		FILE_ATTRIBUTE_NORMAL, NULL);
	if (m_file == INVALID_HANDLE_VALUE)
		Log(EError, "Could not open \"%s\": %s", m_filename.file_string().c_str(),
			lastErrorText().c_str());
	m_fileMapping = CreateFileMapping(m_file, NULL, PAGE_READONLY, 0, 0, NULL);
	if (m_fileMapping == NULL)
		Log(EError, "CreateFileMapping: Could not map \"%s\" to memory: %s", 
			m_filename.file_string().c_str(), lastErrorText().c_str());
	m_data = (void *) MapViewOfFile(m_fileMapping, FILE_MAP_READ, 0, 0, 0);
	if (m_data == NULL)
		Log(EError, "MapViewOfFile: Could not map \"%s\" to memory: %s", 
			m_filename.file_string().c_str(), lastErrorText().c_str());
#endif
}

MemoryMappedFile::~MemoryMappedFile() {
	if (m_data) {
		Log(ETrace, "Unmapping \"%s\" from memory", 
			m_filename.file_string().c_str()); 

		#if defined(__LINUX__) || defined(__OSX__)
			int retval = munmap(m_data, m_size);
			if (retval != 0)
				Log(EError, "munmap(): unable to unmap memory!");
		#elif defined(WIN32)
			if (!UnmapViewOfFile(m_data))
				Log(EError, "UnmapViewOfFile(): unable to unmap memory: %s", lastErrorText().c_str());
			if (!CloseHandle(m_fileMapping))
				Log(EError, "CloseHandle(): unable to close file mapping: %s", lastErrorText().c_str());
			if (!CloseHandle(m_file))
				Log(EError, "CloseHandle(): unable to close file: %s", lastErrorText().c_str());
		#endif
	}
}

MTS_IMPLEMENT_CLASS(MemoryMappedFile, false, Object)
MTS_NAMESPACE_END

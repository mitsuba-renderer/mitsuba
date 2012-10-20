#include <mitsuba/core/mmap.h>

#if defined(__LINUX__) || defined(__OSX__)
# include <sys/mman.h>
# include <fcntl.h>
#elif defined(WIN32)
# include <windows.h>
#endif

MTS_NAMESPACE_BEGIN

struct MemoryMappedFile::MemoryMappedFilePrivate
{
	fs::path filename;
#if defined(WIN32)
	HANDLE file;
	HANDLE fileMapping;
#endif
	size_t size;
	void *data;

	MemoryMappedFilePrivate(const fs::path & f, size_t s = 0) :
	filename(f), size(s), data(NULL) {}
};

MemoryMappedFile::MemoryMappedFile(const fs::path &filename, size_t size)
	: d(new MemoryMappedFilePrivate(filename, size)) {
	Log(ETrace, "Creating memory-mapped file \"%s\" (%s)..",
			filename.filename().string().c_str(), memString(d->size).c_str());
#if defined(__LINUX__) || defined(__OSX__)
	int fd = open(filename.string().c_str(), O_RDWR | O_CREAT | O_TRUNC, 0664);
	if (fd == -1)
		Log(EError, "Could not open \"%s\"!", d->filename.string().c_str());
	int result = lseek(fd, size-1, SEEK_SET);
	if (result == -1)
		Log(EError, "Could not set file size of \"%s\"!", d->filename.string().c_str());
	result = write(fd, "", 1);
	if (result != 1)
		Log(EError, "Could not write to \"%s\"!", d->filename.string().c_str());
	d->data = mmap(NULL, d->size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (d->data == NULL)
		Log(EError, "Could not map \"%s\" to memory!", d->filename.string().c_str());
	if (close(fd) != 0)
		Log(EError, "close(): unable to close file!");
#elif defined(_WIN32)
	d->file = CreateFile(filename.string().c_str(), GENERIC_WRITE | GENERIC_READ,
		FILE_SHARE_WRITE, NULL, CREATE_ALWAYS,
		FILE_ATTRIBUTE_NORMAL, NULL);
	if (d->file == INVALID_HANDLE_VALUE)
		Log(EError, "Could not open \"%s\": %s", d->filename.string().c_str(),
			lastErrorText().c_str());
	d->fileMapping = CreateFileMapping(d->file, NULL, PAGE_READWRITE, 0,
		static_cast<DWORD>(size), NULL);
	if (d->fileMapping == NULL)
		Log(EError, "CreateFileMapping: Could not map \"%s\" to memory: %s",
			d->filename.string().c_str(), lastErrorText().c_str());
	d->data = (void *) MapViewOfFile(d->fileMapping, FILE_MAP_WRITE, 0, 0, 0);
	if (d->data == NULL)
		Log(EError, "MapViewOfFile: Could not map \"%s\" to memory: %s",
			d->filename.string().c_str(), lastErrorText().c_str());
#endif
}


MemoryMappedFile::MemoryMappedFile(const fs::path &filename)
	: d(new MemoryMappedFilePrivate(filename)) {
	if (!fs::exists(filename))
		Log(EError, "The file \"%s\" does not exist!", filename.string().c_str());
	d->size = (size_t) fs::file_size(filename);
	Log(ETrace, "Mapping \"%s\" into memory (%s)..",
			filename.filename().string().c_str(), memString(d->size).c_str());
#if defined(__LINUX__) || defined(__OSX__)
	int fd = open(filename.string().c_str(), O_RDONLY);
	if (fd == -1)
		Log(EError, "Could not open \"%s\"!", d->filename.string().c_str());
	d->data = mmap(NULL, d->size, PROT_READ, MAP_SHARED, fd, 0);
	if (d->data == NULL)
		Log(EError, "Could not map \"%s\" to memory!", d->filename.string().c_str());
	if (close(fd) != 0)
		Log(EError, "close(): unable to close file!");
#elif defined(WIN32)
	d->file = CreateFile(filename.string().c_str(), GENERIC_READ,
		FILE_SHARE_READ, NULL, OPEN_EXISTING,
		FILE_ATTRIBUTE_NORMAL, NULL);
	if (d->file == INVALID_HANDLE_VALUE)
		Log(EError, "Could not open \"%s\": %s", d->filename.string().c_str(),
			lastErrorText().c_str());
	d->fileMapping = CreateFileMapping(d->file, NULL, PAGE_READONLY, 0, 0, NULL);
	if (d->fileMapping == NULL)
		Log(EError, "CreateFileMapping: Could not map \"%s\" to memory: %s",
			d->filename.string().c_str(), lastErrorText().c_str());
	d->data = (void *) MapViewOfFile(d->fileMapping, FILE_MAP_READ, 0, 0, 0);
	if (d->data == NULL)
		Log(EError, "MapViewOfFile: Could not map \"%s\" to memory: %s",
			d->filename.string().c_str(), lastErrorText().c_str());
#endif
}

MemoryMappedFile::~MemoryMappedFile() {
	if (d->data != NULL) {
		Log(ETrace, "Unmapping \"%s\" from memory",
			d->filename.string().c_str());

		#if defined(__LINUX__) || defined(__OSX__)
			int retval = munmap(d->data, d->size);
			if (retval != 0)
				Log(EError, "munmap(): unable to unmap memory!");
		#elif defined(WIN32)
			if (!UnmapViewOfFile(d->data))
				Log(EError, "UnmapViewOfFile(): unable to unmap memory: %s", lastErrorText().c_str());
			if (!CloseHandle(d->fileMapping))
				Log(EError, "CloseHandle(): unable to close file mapping: %s", lastErrorText().c_str());
			if (!CloseHandle(d->file))
				Log(EError, "CloseHandle(): unable to close file: %s", lastErrorText().c_str());
		#endif
	}
}

void * MemoryMappedFile::getData() {
	return d->data;
}

/// Return a pointer to the file contents in memory (const version)
const void * MemoryMappedFile::getData() const {
	return d->data;
}

size_t MemoryMappedFile::getSize() const {
	return d->size;
}

MTS_IMPLEMENT_CLASS(MemoryMappedFile, false, Object)
MTS_NAMESPACE_END

#include <mitsuba/core/mmap.h>

#if defined(__LINUX__) || defined(__OSX__)
# include <sys/mman.h>
# include <fcntl.h>
#elif defined(__WINDOWS__)
# include <windows.h>
#endif

MTS_NAMESPACE_BEGIN

struct MemoryMappedFile::MemoryMappedFilePrivate {
	fs::path filename;
#if defined(__WINDOWS__)
	HANDLE file;
	HANDLE fileMapping;
#endif
	size_t size;
	void *data;
	bool readOnly;

	MemoryMappedFilePrivate(const fs::path & f, size_t s = 0)
		: filename(f), size(s), data(NULL), readOnly(false) {}

	void create() {
		#if defined(__LINUX__) || defined(__OSX__)
			int fd = open(filename.string().c_str(), O_RDWR | O_CREAT | O_TRUNC, 0664);
			if (fd == -1)
				Log(EError, "Could not open \"%s\"!", filename.string().c_str());
			int result = lseek(fd, size-1, SEEK_SET);
			if (result == -1)
				Log(EError, "Could not set file size of \"%s\"!", filename.string().c_str());
			result = write(fd, "", 1);
			if (result != 1)
				Log(EError, "Could not write to \"%s\"!", filename.string().c_str());
			data = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
			if (data == NULL)
				Log(EError, "Could not map \"%s\" to memory!", filename.string().c_str());
			if (close(fd) != 0)
				Log(EError, "close(): unable to close file!");
		#elif defined(__WINDOWS__)
			file = CreateFile(filename.string().c_str(), GENERIC_WRITE | GENERIC_READ,
				FILE_SHARE_READ, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
			if (file == INVALID_HANDLE_VALUE)
				Log(EError, "Could not open \"%s\": %s", filename.string().c_str(),
					lastErrorText().c_str());
			fileMapping = CreateFileMapping(file, NULL, PAGE_READWRITE, 0,
				static_cast<DWORD>(size), NULL);
			if (fileMapping == NULL)
				Log(EError, "CreateFileMapping: Could not map \"%s\" to memory: %s",
					filename.string().c_str(), lastErrorText().c_str());
			data = (void *) MapViewOfFile(fileMapping, FILE_MAP_WRITE, 0, 0, 0);
			if (data == NULL)
				Log(EError, "MapViewOfFile: Could not map \"%s\" to memory: %s",
					filename.string().c_str(), lastErrorText().c_str());
		#endif
		readOnly = false;
	}

	void map() {
		if (!fs::exists(filename))
			Log(EError, "The file \"%s\" does not exist!", filename.string().c_str());
		size = (size_t) fs::file_size(filename);
		#if defined(__LINUX__) || defined(__OSX__)
			int fd = open(filename.string().c_str(), readOnly ? O_RDONLY : O_RDWR);
			if (fd == -1)
				Log(EError, "Could not open \"%s\"!", filename.string().c_str());
			data = mmap(NULL, size, PROT_READ | (readOnly ? 0 : PROT_WRITE), MAP_SHARED, fd, 0);
			if (data == NULL)
				Log(EError, "Could not map \"%s\" to memory!", filename.string().c_str());
			if (close(fd) != 0)
				Log(EError, "close(): unable to close file!");
		#elif defined(__WINDOWS__)
			file = CreateFile(filename.string().c_str(), GENERIC_READ | (readOnly ? 0 : GENERIC_WRITE),
				FILE_SHARE_WRITE|FILE_SHARE_READ, NULL, OPEN_EXISTING,
				FILE_ATTRIBUTE_NORMAL, NULL);
			if (file == INVALID_HANDLE_VALUE)
				Log(EError, "Could not open \"%s\": %s", filename.string().c_str(),
					lastErrorText().c_str());
			fileMapping = CreateFileMapping(file, NULL, readOnly ? PAGE_READONLY : PAGE_READWRITE, 0, 0, NULL);
			if (fileMapping == NULL)
				Log(EError, "CreateFileMapping: Could not map \"%s\" to memory: %s",
					filename.string().c_str(), lastErrorText().c_str());
			data = (void *) MapViewOfFile(fileMapping, readOnly ? FILE_MAP_READ : FILE_MAP_WRITE, 0, 0, 0);
			if (data == NULL)
				Log(EError, "MapViewOfFile: Could not map \"%s\" to memory: %s",
					filename.string().c_str(), lastErrorText().c_str());
		#endif
	}

	void unmap() {
		SLog(ETrace, "Unmapping \"%s\" from memory",
			filename.string().c_str());
		#if defined(__LINUX__) || defined(__OSX__)
			int retval = munmap(data, size);
			if (retval != 0)
				Log(EError, "munmap(): unable to unmap memory!");
		#elif defined(__WINDOWS__)
			if (!UnmapViewOfFile(data))
				Log(EError, "UnmapViewOfFile(): unable to unmap memory: %s", lastErrorText().c_str());
			if (!CloseHandle(fileMapping))
				Log(EError, "CloseHandle(): unable to close file mapping: %s", lastErrorText().c_str());
			if (!CloseHandle(file))
				Log(EError, "CloseHandle(): unable to close file: %s", lastErrorText().c_str());
		#endif
		data = NULL;
		size = 0;
	}
};

MemoryMappedFile::MemoryMappedFile(const fs::path &filename, size_t size)
	: d(new MemoryMappedFilePrivate(filename, size)) {
	SLog(ETrace, "Creating memory-mapped file \"%s\" (%s)..",
		filename.filename().string().c_str(), memString(d->size).c_str());
	d->create();
}


MemoryMappedFile::MemoryMappedFile(const fs::path &filename, bool readOnly)
	: d(new MemoryMappedFilePrivate(filename)) {
	d->readOnly = readOnly;
	d->map();
	Log(ETrace, "Mapped \"%s\" into memory (%s)..",
		filename.filename().string().c_str(), memString(d->size).c_str());
}

MemoryMappedFile::~MemoryMappedFile() {
	if (d->data) {
		try {
			d->unmap();
		} catch (std::exception &e) {
			/* Don't throw exceptions from a constructor */
			Log(EWarn, "%s", e.what());
		}
	}
}

void MemoryMappedFile::resize(size_t size) {
	if (!d->data)
		Log(EError, "Internal error in MemoryMappedFile::resize()!");
	d->unmap();
	fs::resize_file(d->filename, size);
	d->size = size;
	d->map();
}

void *MemoryMappedFile::getData() {
	return d->data;
}

/// Return a pointer to the file contents in memory (const version)
const void *MemoryMappedFile::getData() const {
	return d->data;
}

size_t MemoryMappedFile::getSize() const {
	return d->size;
}

bool MemoryMappedFile::isReadOnly() const {
	return d->readOnly;
}

std::string MemoryMappedFile::toString() const {
	std::ostringstream oss;
	oss << "MemoryMappedFile[filename=\""
		<< d->filename.string() << "\", size="
		<< memString(d->size) << "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(MemoryMappedFile, false, Object)
MTS_NAMESPACE_END

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

#include <mitsuba/core/fstream.h>
#include <cerrno>

#if !defined(WIN32)
#include <unistd.h>
#endif

MTS_NAMESPACE_BEGIN

FileStream::FileStream()
 : m_file(0) {
}

FileStream::FileStream(const fs::path &path, EFileMode mode)
 : m_file(0) {
	open(path, mode);
}


FileStream::~FileStream() {
	if (m_file != 0)
		close();
}

std::string FileStream::toString() const {
	std::ostringstream oss;
	oss << "FileStream[" << Stream::toString()
		<< ", path=\"" << m_path.file_string()
		<< "\", mode=" << m_mode << "]";
	return oss.str();
}

void FileStream::open(const fs::path &path, EFileMode mode) {
	AssertEx(m_file == 0, "A file has already been opened using this stream");

	Log(ETrace, "Opening \"%s\"", path.file_string().c_str());

	m_path = path;
	m_mode = mode;
	m_write = true;
	m_read = true;

#ifdef WIN32
	DWORD dwDesiredAccess = GENERIC_READ;
	DWORD dwCreationDisposition = OPEN_EXISTING;

	switch (m_mode) {
	case EReadOnly:
		m_write = false;
		break;
	case EReadWrite:
		dwDesiredAccess |= GENERIC_WRITE;
		break;
	case ETruncWrite:
		m_read = false;
		dwDesiredAccess = GENERIC_WRITE;
		dwCreationDisposition = CREATE_ALWAYS;
		break;
	case ETruncReadWrite:
		dwDesiredAccess |= GENERIC_WRITE;
		dwCreationDisposition = CREATE_ALWAYS;
		break;
	case EAppendWrite:
		m_read = false;
		dwDesiredAccess = GENERIC_WRITE;
		break;
	case EAppendReadWrite:
		dwDesiredAccess |= GENERIC_WRITE;
		break;
	default: 
		Log(EError, "Unknown file mode");
		break;
	}

	m_file = CreateFile(path.file_string().c_str(), dwDesiredAccess, 
		FILE_SHARE_WRITE | FILE_SHARE_READ, 0, 
		dwCreationDisposition, FILE_ATTRIBUTE_NORMAL, 0);

	if (m_file == INVALID_HANDLE_VALUE)
		Log(EError, "Error while trying to open file \"%s\": %s", 
			m_path.file_string().c_str(), lastErrorText().c_str());
	
	if (m_mode == EAppendWrite || m_mode == EAppendReadWrite)
		setPos(getSize());
#else
	const char *modeString = NULL;

	switch (m_mode) {
	case EReadOnly:
		modeString = "rb";
		m_write = false;
		break;
	case EReadWrite:
		modeString = "rb+";
		break;
	case ETruncWrite:
		modeString = "wb";
		m_read = false;
		break;
	case ETruncReadWrite:
		modeString = "wb+";
		break;
	case EAppendWrite:
		modeString = "ab";
		m_read = false;
		break;
	case EAppendReadWrite:
		modeString = "ab+";
		break;
	default:
		Log(EError, "Unknown file mode");
		break;
	};

	m_file = fopen(m_path.file_string().c_str(), modeString);

	if (m_file == NULL) {
		Log(EError, "Error while trying to open file \"%s\": %s", 
			m_path.file_string().c_str(), strerror(errno));
	}
#endif
}

void FileStream::close() {
	AssertEx(m_file != 0, "No file is currently open");
	Log(ETrace, "Closing \"%s\"", m_path.file_string().c_str());

#ifdef WIN32
	if (!CloseHandle(m_file)) {
		Log(EError, "Error while trying to close file \"%s\": %s", 
			m_path.file_string().c_str(), lastErrorText().c_str());
	}
#else
	if (fclose(m_file)) {
		Log(EError, "Error while trying to close file \"%s\": %s", 
			m_path.file_string().c_str(), strerror(errno));
	}
#endif
	m_file = 0;
}


void FileStream::remove() {
	close();
	Log(EDebug, "Removing \"%s\"", 	m_path.file_string().c_str());

	fs::remove(m_path);
}

void FileStream::setPos(size_t pos) {
	AssertEx(m_file != 0, "No file is currently open");
	
#ifdef WIN32
	LARGE_INTEGER fpos;
	fpos.QuadPart = pos;
	if (SetFilePointerEx(m_file, fpos, 0, FILE_BEGIN) == INVALID_SET_FILE_POINTER) {
		Log(EError, "Error while trying to seek to position %i in file \"%s\": %s", 
			pos, m_path.file_string().c_str(), lastErrorText().c_str());
	}
#else
	if (fseek(m_file, pos, SEEK_SET)) {
		Log(EError, "Error while trying to seek to position %i in file \"%s\": %s", 
			pos, m_path.file_string().c_str(), strerror(errno));
	}
#endif
}

size_t FileStream::getPos() const {
	AssertEx(m_file != 0, "No file is currently open");
#ifdef WIN32
	DWORD pos = SetFilePointer(m_file, 0, 0, FILE_CURRENT);
	if (pos == INVALID_SET_FILE_POINTER) {
		Log(EError, "Error while looking up the position in file \"%s\": %s",
			m_path.file_string().c_str(), lastErrorText().c_str());
	}
	return (size_t) pos;
#else
	long pos;
	pos = ftell(m_file);
	if (pos == -1) {
		Log(EError, "Error while looking up the position in file \"%s\": %s", 
			m_path.file_string().c_str(), strerror(errno));
	}
	return (size_t) pos;
#endif
}

size_t FileStream::getSize() const {
	AssertEx(m_file != 0, "No file is currently open");

#ifdef WIN32
	LARGE_INTEGER result;
	if (GetFileSizeEx(m_file, &result) == 0) {
		Log(EError, "Error while getting the file size of \"%s\": %s",
			m_path.file_string().c_str(), lastErrorText().c_str());
	}
	return (size_t) result.QuadPart;
#else
	size_t size, tmp;
	
	tmp = getPos();
	if (fseek(m_file, 0, SEEK_END)) {
		Log(EError, "Error while seeking within \"%s\": %s",
			m_path.file_string().c_str(), strerror(errno));
	}
	size = getPos();
	if (fseek(m_file, tmp, SEEK_SET)) {
		Log(EError, "Error while seeking within \"%s\": %s",
			m_path.file_string().c_str(), strerror(errno));
	}
	return size;
#endif
}

void FileStream::truncate(size_t size) {
	AssertEx(m_file != 0, "No file is currently open");
	AssertEx(m_write, "File is not open with write access");

	size_t pos = getPos();
	if (pos > size) 
		pos = size;

#ifdef WIN32
	setPos(size);
	if (!SetEndOfFile(m_file)) {
		Log(EError, "Error while truncating file \"%s\": %s",
			m_path.file_string().c_str(), lastErrorText().c_str());
	}
#else
	/* File truncation support blows on posix.. */
	setPos(pos);
	flush();

	if (ftruncate(fileno(m_file), size)) {
		Log(EError, "Error while truncating file \"%s\": %s",
			m_path.file_string().c_str(), strerror(errno));
	}
#endif
	setPos(pos);
}

void FileStream::flush() {
	AssertEx(m_file != 0, "No file is currently open");
	AssertEx(m_write, "File is not open with write access");
#ifdef WIN32
	if (!FlushFileBuffers(m_file)) {
		Log(EError, "Error while flusing the buffers of \"%s\": %s",
			m_path.file_string().c_str(), lastErrorText().c_str());
	}
#else
	if (fflush(m_file) != 0) {
		Log(EError, "Error while flusing the buffers of \"%s\": %s",
			m_path.file_string().c_str(), strerror(errno));
	}
#endif
}

void FileStream::read(void *pPtr, size_t size) {
	AssertEx(m_file != 0, "No file is currently open");
	AssertEx(m_read, "File is not open with read access");
	
	if (size == 0)
		return;
#ifdef WIN32
	DWORD lpNumberOfBytesRead;
	if (!ReadFile(m_file, pPtr, (DWORD) size, &lpNumberOfBytesRead, 0)) {
		Log(EError, "Error while reading from file \"%s\": %s",
			m_path.file_string().c_str(), lastErrorText().c_str());
	}
	if (lpNumberOfBytesRead != (DWORD) size) 
		throw EOFException(formatString("Read less data than expected (%i bytes required) "
			"from file \"%s\"", size, m_path.file_string().c_str()), (size_t) lpNumberOfBytesRead);
#else
	size_t bytesRead;
	if ((bytesRead = fread(pPtr, 1, size, m_file)) != size) {
		if (ferror(m_file) != 0) {
			Log(EError, "Error while reading from file \"%s\": %s",
				m_path.file_string().c_str(), strerror(errno));
		}
		throw EOFException(formatString("Read less data than expected (%i bytes required) "
			"from file \"%s\"", size, m_path.file_string().c_str()), bytesRead);
	}
#endif
}

void FileStream::write(const void *pPtr, size_t size) {
	AssertEx(m_file != 0, "No file is currently open");
	AssertEx(m_write, "File is not open with write access");

	if (size == 0)
		return;

#ifdef WIN32
	DWORD lpNumberOfBytesWritten;
	if (!WriteFile(m_file, pPtr, (DWORD) size, &lpNumberOfBytesWritten, 0)) {
		Log(EError, "Error while writing to file \"%s\": %s",
			m_path.file_string().c_str(), lastErrorText().c_str());
	}
	if (lpNumberOfBytesWritten != (DWORD) size) 
		throw EOFException(formatString("Wrote less data than expected (%i bytes required) "
			"to file \"%s\"", size, m_path.file_string().c_str()), (size_t) lpNumberOfBytesWritten);
#else
	size_t bytesWritten;
	if ((bytesWritten = fwrite(pPtr, 1, size, m_file)) != size) {
		if (ferror(m_file))
			Log(EError, "Error while writing to file \"%s\": %s",
				m_path.file_string().c_str(), strerror(errno));
		throw EOFException(formatString("Wrote less data than expected (%i bytes required) "
			"to file \"%s\"", size, m_path.file_string().c_str()), bytesWritten);
	}
#endif
}

bool FileStream::canRead() const {
	AssertEx(m_file != 0, "No file is currently open");
	return m_read;
}

bool FileStream::canWrite() const {
	AssertEx(m_file != 0, "No file is currently open");
	return m_write;
}

MTS_IMPLEMENT_CLASS(FileStream, false, Stream)
MTS_NAMESPACE_END

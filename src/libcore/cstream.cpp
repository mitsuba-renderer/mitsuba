#include <mitsuba/core/cstream.h>
#include <errno.h>

MTS_NAMESPACE_BEGIN

ConsoleStream::ConsoleStream() {
	setByteOrder(ENetworkByteOrder);
}

ConsoleStream::~ConsoleStream() {
}

std::string ConsoleStream::toString() const {
	return "ConsoleStream[]";
}

void ConsoleStream::setPos(size_t pos) {
	Log(EError, "Cannot seek within a console stream!");
}

size_t ConsoleStream::getPos() const {
	Log(EError, "Cannot determine the position within a console stream!");
	return 0;
}

size_t ConsoleStream::getSize() const {
	Log(EError, "Cannot determine the size of a console stream!");
	return 0;
}

void ConsoleStream::truncate(size_t) {
	Log(EError, "Cannot truncate a console stream!");
}

void ConsoleStream::flush() {
	if (fflush(stdout) == EOF)
		Log(EError, "Error in fflush(): %s!", strerror(errno));
}

void ConsoleStream::read(void *ptr, size_t size) {
	if (fread(ptr, size, 1, stdin) != 1) {
		if (feof(stdin))
			Log(EError, "Error in fread(): end of file!");
		else if (ferror(stdin))
			Log(EError, "Error in fread(): stream error!");
		/* Otherwise, ignore (strange, but seems to be required) */
	}

}

void ConsoleStream::write(const void *ptr, size_t size) {
	if (fwrite(ptr, size, 1, stdout) != 1) {
		if (feof(stdout))
			Log(EError, "Error in fwrite(): end of file!");
		else if (ferror(stdout))
			Log(EError, "Error in fwrite(): stream error!");
		/* Otherwise, ignore (strange, but seems to be required) */
	}
}

bool ConsoleStream::canRead() const {
	return true;
}

bool ConsoleStream::canWrite() const {
	return true;
}

MTS_IMPLEMENT_CLASS(ConsoleStream, false, Stream)
MTS_NAMESPACE_END

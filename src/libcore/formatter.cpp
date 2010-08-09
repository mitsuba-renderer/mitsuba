#include <mitsuba/mitsuba.h>
#include <mitsuba/core/thread.h>
#include <ctime>

MTS_NAMESPACE_BEGIN

DefaultFormatter::DefaultFormatter() 
 : m_haveDate(true) {
}

std::string DefaultFormatter::format(ELogLevel logLevel, const Class *theClass,
		const Thread *thread, const std::string &text, const char *file, int pLine) {
	std::ostringstream oss;
	char buffer[128];

	/* Date/Time */
	if (m_haveDate) {
		time_t theTime = std::time(NULL);
		strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S ", std::localtime(&theTime));
		oss << buffer;
	}

	/* Level */
	switch (logLevel) {
		case ETrace: oss << "TRACE "; break;
		case EDebug: oss << "DEBUG "; break;
		case EInfo:  oss << "INFO  "; break;
		case EWarn:  oss << "WARN  "; break;
		case EError: oss << "ERROR "; break;
		default:     oss << "CUSTM "; break;
	}

	/* Thread */
	if (thread) {
		oss << thread->getName();

		for (int i=0; i<(5 - (int) thread->getName().size()); i++)
			oss << ' ';
	}

	/* Class */
	if (theClass) {
		oss << "[" << theClass->getName() << "] ";
	} else if (pLine != -1 && file) {
		oss << "[" << file << ":" << pLine << "] ";
	}

	/* Text */
	oss << text;

	return oss.str();
}

MTS_IMPLEMENT_CLASS(Formatter, true, Object)
MTS_IMPLEMENT_CLASS(DefaultFormatter, false, Formatter)
MTS_NAMESPACE_END

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

Timer::Timer() {
	reset();
}

Timer::~Timer() {
}

void Timer::reset() {
#ifdef WIN32
	QueryPerformanceFrequency((PLARGE_INTEGER) &m_frequency);
	QueryPerformanceCounter((PLARGE_INTEGER) &m_start);
#else
	gettimeofday(&m_start, NULL);
#endif
}

unsigned int Timer::getMilliseconds() const {
#ifdef WIN32
	LARGE_INTEGER current;
	QueryPerformanceCounter(&current);

	return (int) ((current.QuadPart - m_start) * 1000 / m_frequency);
#else
	struct timeval current;

	gettimeofday(&current, NULL);

	return (current.tv_sec - m_start.tv_sec) * 1000 +
		   (current.tv_usec - m_start.tv_usec) / 1000;
#endif
}

unsigned int Timer::getMicroseconds() const {
#ifdef WIN32
	Log(EError, "Microsecond timer resolution is not available on WIN32");
	return 0;
#else
	struct timeval current;

	gettimeofday(&current, NULL);

	return (current.tv_sec - m_start.tv_sec) * 1000000 +
		   (current.tv_usec - m_start.tv_usec);
#endif
}

std::string Timer::toString() const {
	std::ostringstream oss;
	oss << "Timer[ms=" << getMilliseconds() << "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(Timer, false, Object)
MTS_NAMESPACE_END

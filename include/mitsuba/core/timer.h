#if !defined(__TIMER_H)
#define __TIMER_H

#ifndef WIN32
#include <sys/time.h>
#endif

MTS_NAMESPACE_BEGIN

/** \brief Platform independent milli/microsecond timer
 */
class MTS_EXPORT_CORE Timer : public Object {
public:
	/// Create a new timer and reset it
	Timer();

	/// Reset the timer
	void reset();

	/// Return the milliseconds which have passed since the last reset
	unsigned int getMilliseconds() const;

	/// Return the microseconds which have passed since the last reset
	unsigned int getMicroseconds() const;

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Timer();
private:
#ifdef WIN32
	__int64 m_start;
	__int64 m_frequency;
#else
	struct timeval m_start;
#endif
};

MTS_NAMESPACE_END

#endif /* __TIMER_H */

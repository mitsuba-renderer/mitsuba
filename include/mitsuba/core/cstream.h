#if !defined(__CSTREAM_H)
#define __CSTREAM_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/** \brief Console stream
 */
class MTS_EXPORT_CORE ConsoleStream : public Stream {
public:
	/// Create a new console stream
	ConsoleStream();

	/* Stream implementation */
	void read(void *ptr, size_t size);
	void write(const void *ptr, size_t size);
	void setPos(size_t pos);
	size_t getPos() const;
	size_t getSize() const;
	void truncate(size_t size);
	void flush();
	bool canWrite() const;
	bool canRead() const;

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/** \brief Virtual destructor
	 *
	 * The destructor frees all resources and closes
	 * the socket if it is still open
	 */
	virtual ~ConsoleStream();
};

MTS_NAMESPACE_END

#endif /* __CSTREAM_H */

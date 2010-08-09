#if !defined(__SESSION_H)
#define __SESSION_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

class Device;

/** \brief Abstract windowing environment session
 */
class MTS_EXPORT_HW Session : public Object {
	friend class Device;
public:
	/// Create a new session using the appropriate implementation
	static Session *create();

	/// Initialize the session
	virtual void init();

	/// Shut the session down
	virtual void shutdown();

	/// Process all events and call event callbacks
	virtual void processEvents() = 0;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Session() { }

	/// Create a new session
	Session();
protected:
	bool m_initialized;
	std::vector<Device *> m_devices;
};

MTS_NAMESPACE_END

#endif /* __SESSION_H */

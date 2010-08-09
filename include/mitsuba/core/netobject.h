#if !defined(__NETOBJECT_H)
#define __NETOBJECT_H

#include <mitsuba/core/sched.h>

MTS_NAMESPACE_BEGIN

/** \brief Interface of an abstract object referencing
 * globally shared resources. When it is serialized for use in a
 * parallel process executed on several machines, the object
 * is first given the opportunity to bind named resources to 
 * the process, which will then be distributed to all participating 
 * compute servers. Once unserialized on the remote side, 
 * <tt>wakeup</tt> is called to let the object re-associate
 * with the shared resources.
 */
class MTS_EXPORT_CORE NetworkedObject : public ConfigurableObject {
public:
	/// Constructor
	inline NetworkedObject(const Properties &props) : ConfigurableObject(props) { }

	/// Unserialize a configurable object
	inline NetworkedObject(Stream *stream, InstanceManager *manager) 
	 : ConfigurableObject(stream, manager) {
	}

	virtual void bindUsedResources(ParallelProcess *proc) const;
	virtual void wakeup(std::map<std::string, SerializableObject *> &params);

	/// Serialize this object to a stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~NetworkedObject() { }
};

MTS_NAMESPACE_END

#endif /* __NETOBJECT_H */

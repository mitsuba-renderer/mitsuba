#if !defined(__SERIALIZATION_H)
#define __SERIALIZATION_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/** \brief Serializable interface
 */
class MTS_EXPORT_CORE Serializable {
public:
	/// Serialize this object to a stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const = 0;
protected:
	/// Virtual deconstructor
	virtual ~Serializable() { }
};

/** \brief Base class of all reference-counted objects with serialization support.
 * To support unserialization from a stream, the implementation should use one of the
 * RTTI macros "MTS_IMPLEMENT_CLASS_*S".
 */
class MTS_EXPORT_CORE SerializableObject : public Object, public Serializable {
public:
	inline SerializableObject() { }

	/// Unserialize a serializable object
	SerializableObject(Stream *stream, InstanceManager *manager);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual deconstructor
	virtual ~SerializableObject() { }
};

/** \brief The instance manager coordinates the serialization/
 * unserialization process when link structures are encountered
 * in which an object is referenced by multiple other objects
 * within the same data stream. Cyclic dependencies can also
 * be handled.
 */
class MTS_EXPORT_CORE InstanceManager : public Object {
	friend class SerializableObject;
public:
	/** \brief Construct a new instance manager */
	InstanceManager();

	/// Retrieve/load an instance by ID
	SerializableObject *getInstance(Stream *stream);

	/// Store/skip an instance returning its ID
	void serialize(Stream *stream, const SerializableObject *inst);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~InstanceManager();

	/// Called from the unserialization constructor of SerializableObject
	void registerInstance(SerializableObject *object);
private:
	unsigned int m_counter, m_lastID;
	std::vector<SerializableObject *> m_fullyAllocated;
	std::map<unsigned int, SerializableObject *> m_idToObj;
	std::map<const SerializableObject *, unsigned int> m_objToId;
};

MTS_NAMESPACE_END

#endif /* __SERIALIZATION_H */

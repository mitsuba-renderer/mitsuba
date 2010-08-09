#if !defined(__GPUSYNC_H)
#define __GPUSYNC_H

#include <mitsuba/hw/renderer.h>

MTS_NAMESPACE_BEGIN

/** \brief Abstract GPU synchronization object implementing
 * a memory fence operation.
 */
class MTS_EXPORT_HW GPUSync : public Object {
public:
	/// Allocate memory for a new synchronization object
	GPUSync();

	/// Create the synchronization object on the GPU
	virtual void init() = 0;

	/// Wait on the fence (blocking)
	virtual void wait() = 0;

	/// Enqueue a wait command, but do not block
	virtual void enqueueWait() = 0;

	/// Remove the synchronization object
	virtual void cleanup() = 0;

	/// Return a string representation of this class
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GPUSync();
};

MTS_NAMESPACE_END

#endif /* __GPUSYNC_H */

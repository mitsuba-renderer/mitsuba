#if !defined(__GLSYNC_H)
#define __GLSYNC_H

#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/hw/gpusync.h>

MTS_NAMESPACE_BEGIN

/** \brief OpenGL-based GPUSync implementation
 */
class MTS_EXPORT_HW GLSync : public GPUSync {
public:
	/// Allocate memory for a new synchronization object
	GLSync();

	/// Create the synchronization object on the GL
	void init();

	/// Wait on the fence (blocking)
	void wait();

	/// Enqueue a wait command, but do not block
	void enqueueWait();

	/// Remove the synchronization object
	void cleanup();

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	~GLSync();
protected:
	GLsync m_sync;
};

MTS_NAMESPACE_END

#endif /* __GLSYNC_H */

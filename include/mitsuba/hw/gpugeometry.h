#if !defined(__GPUGEOMETRY_H)
#define __GPUGEOMETRY_H

#include <mitsuba/hw/renderer.h>

MTS_NAMESPACE_BEGIN

/** \brief Abstract geometry storage on a graphics card
 */
class MTS_EXPORT_HW GPUGeometry : public Object {
public:
	/// Create an empty program
	GPUGeometry(const TriMesh *mesh);

	/// Return the name of this geometry object
	inline const std::string &getName() { return m_mesh->getName(); }

	/// Upload the geometry object
	virtual void init() = 0;

	/// Refresh (re-upload) the geometry object
	virtual void refresh() = 0;

	/// Bind the geometry object
	virtual void bind() = 0;

	/// Unbind the geometry object
	virtual void unbind() = 0;
	
	/// Free the geometry object from GPU memory
	virtual void cleanup() = 0;

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GPUGeometry();
protected:
	ref<const TriMesh> m_mesh;
};

MTS_NAMESPACE_END

#endif /* __GPUGEOMETRY_H */

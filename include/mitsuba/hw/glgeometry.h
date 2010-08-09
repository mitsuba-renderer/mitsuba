#if !defined(__GLGEOMETRY_H)
#define __GLGEOMETRY_H

#include <mitsuba/hw/glrenderer.h>
#include <mitsuba/hw/gpugeometry.h>

MTS_NAMESPACE_BEGIN

/** \brief OpenGL-based GPUGeometry implementation
 */
class MTS_EXPORT_HW GLGeometry : public GPUGeometry {
	friend class GLRenderer;
public:
	/// Create a new GLGeometry instance for the given mesh
	GLGeometry(const TriMesh *mesh);

	/// Upload the geometry object
	void init();

	/// Refresh (re-upload) the geometry object
	void refresh();

	/// Bind the geometry object
	void bind();

	/// Unbind the geometry object
	void unbind();

	/// Free the geometry object from GPU memory
	void cleanup();
	
	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GLGeometry();
protected:
	GLuint m_vertexID;
	GLuint m_indexID;
	GLuint64 m_vertexAddr;
	GLuint64 m_indexAddr;
	GLuint m_vertexSize;
	GLuint m_indexSize;
};

MTS_NAMESPACE_END

#endif /* __GLGEOMETRY_H */

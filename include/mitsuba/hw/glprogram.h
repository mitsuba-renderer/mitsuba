#if !defined(__GLPROGRAM_H)
#define __GLPROGRAM_H

#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/hw/glrenderer.h>

MTS_NAMESPACE_BEGIN

/** \brief OpenGL shader class -- responsible from compiling
 * and linking GLSL fragments
 */
class MTS_EXPORT_HW GLProgram : public GPUProgram {
public:
	/// Construct a new (empty) shader
	GLProgram(const std::string &name = "default");

	/// Initialize the shader
	void init();

	/// Bind the shader
	void bind();

	/// Set the default shader program
	void unbind();

	/// Remove all related OpenGL handles
	void cleanup();

	/// Determine the ID number of a named parameter
	int getParameterID(const std::string &name, bool failIfMissing = true) const;
	
	/// Set a boolean parameter
	void setParameter(int id, bool value);

	/// Set a float parameter
	void setParameter(int id, Float value);

	/// Set a Vector parameter
	void setParameter(int id, const Vector &value);

	/// Set a Vector3i parameter
	void setParameter(int id, const Vector3i &value);

	/// Set a Vector2 parameter
	void setParameter(int id, const Vector2 &value);

	/// Set a Vector2i parameter
	void setParameter(int id, const Vector2i &value);

	/// Set a Vector4 parameter
	void setParameter(int id, const Vector4 &value);

	/// Set a Transform parameter
	void setParameter(int id, const Transform &value);

	/// Set a Spectrum parameter (will be converted to linear RGB)
	void setParameter(int id, const Spectrum &value);

	/** Set a GPUTexture parameter. Must be executed after
	    binding the texture to a texture unit */
	void setParameter(int id, const GPUTexture *value);

	MTS_DECLARE_CLASS()
protected:
	/// Create a shader from source code
	int createShader(int type, const std::string &source);

	/// Return the info log (shader)
	std::string getInfoLogShader(int id);

	/// Return the info log (program)
	std::string getInfoLogProgram();

	/// Virtual destructor
	virtual ~GLProgram();
private:
	GLuint m_id[3];
	GLuint m_program;
};

MTS_NAMESPACE_END

#endif /* __GLPROGRAM_H */

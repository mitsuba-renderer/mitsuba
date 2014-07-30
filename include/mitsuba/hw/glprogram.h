/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#if !defined(__MITSUBA_HW_GLPROGRAM_H_)
#define __MITSUBA_HW_GLPROGRAM_H_

#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/hw/glrenderer.h>

MTS_NAMESPACE_BEGIN

/** \brief OpenGL shader class -- responsible from compiling
 * and linking GLSL fragments
 * \ingroup libhw
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

	/// Set a integer parameter
	void setParameter(int id, int value);

	/// Set a unsigned integer parameter
	void setParameter(int id, uint32_t value);

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

	/// Set a Point parameter
	void setParameter(int id, const Point &value);

	/// Set a Point3i parameter
	void setParameter(int id, const Point3i &value);

	/// Set a Point2 parameter
	void setParameter(int id, const Point2 &value);

	/// Set a Point2i parameter
	void setParameter(int id, const Point2i &value);

	/// Set a Point4 parameter
	void setParameter(int id, const Point4 &value);

	/// Set a Matrix2x2 parameter
	void setParameter(int id, const Matrix2x2 &value);

	/// Set a Matrix3x3 parameter
	void setParameter(int id, const Matrix3x3 &value);

	/// Set a Matrix4x4 parameter
	void setParameter(int id, const Matrix4x4 &value);

	/// Set a Spectrum parameter (will be converted to linear RGB)
	void setParameter(int id, const Spectrum &value);

	/// Set a Color3 parameter
	void setParameter(int id, const Color3 &value);

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

#endif /* __MITSUBA_HW_GLPROGRAM_H_ */

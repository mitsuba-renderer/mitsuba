/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#if !defined(__GPUPROGRAM_H)
#define __GPUPROGRAM_H

#include <mitsuba/hw/renderer.h>

MTS_NAMESPACE_BEGIN

/** \brief Abstract shader program (for fragment/vertex shading)
 * \ingroup libhw
 */
class MTS_EXPORT_HW GPUProgram : public Object {
public:
	enum EType {
		EVertexProgram = 0,
		EFragmentProgram,
		EGeometryProgram
	};

	/// Create an empty program
	GPUProgram(const std::string &name = "default");

	/// Set the name of this program
	inline void setName(const std::string &name) { m_name = name; }

	/// Return the name of this program
	inline const std::string &getName() const { return m_name; }

	/// Set the source code of this program
	inline void setSource(EType type, const std::string &source) { m_source[type] = source; }
	
	/// Set the source code of this program by filename
	void setSourceFile(EType type, const fs::path &path);

	/// Get the source code of this program
	inline const std::string &getSource(EType type) const { return m_source[type]; }

	/// Upload to the GPU
	virtual void init() = 0;

	/// Bind the shader
	virtual void bind() = 0;

	/// Set the default shader program
	virtual void unbind() = 0;

	/// Remove all allocated handles
	virtual void cleanup() = 0;

	/// Determine the ID number of a named parameter
	virtual int getParameterID(const std::string &name, bool failIfMissing = true) const = 0;
	
	/// Set a boolean parameter by name
	inline void setParameter(const std::string &name, bool value, 
		bool failIfMissing = true) {
		setParameter(getParameterID(name, failIfMissing), value);
	}

	/// Set a integer parameter by name
	inline void setParameter(const std::string &name, int value, 
		bool failIfMissing = true) {
		setParameter(getParameterID(name, failIfMissing), value);
	}

	/// Set a float parameter by name
	inline void setParameter(const std::string &name, Float value, 
		bool failIfMissing = true) {
		setParameter(getParameterID(name, failIfMissing), value);
	}

	/// Set a Vector parameter by name
	inline void setParameter(const std::string &name, const Vector &value,
		bool failIfMissing = true) { 
		setParameter(getParameterID(name, failIfMissing), value);
	}
	
	/// Set a Vector3i parameter by name
	inline void setParameter(const std::string &name, const Vector3i &value,
		bool failIfMissing = true) { 
		setParameter(getParameterID(name, failIfMissing), value);
	}
	
	/// Set a Vector2 parameter by name
	inline void setParameter(const std::string &name, const Vector2 &value,
		bool failIfMissing = true) { 
		setParameter(getParameterID(name, failIfMissing), value);
	}

	/// Set a Vector2i parameter by name
	inline void setParameter(const std::string &name, const Vector2i &value,
		bool failIfMissing = true) { 
		setParameter(getParameterID(name, failIfMissing), value);
	}

	/// Set a Vector4 parameter by name
	inline void setParameter(const std::string &name, const Vector4 &value,
		bool failIfMissing = true) { 
		setParameter(getParameterID(name, failIfMissing), value);
	}

	/// Set a Point parameter by name
	inline void setParameter(const std::string &name, const Point &value,
		bool failIfMissing = true) { 
		setParameter(getParameterID(name, failIfMissing), value);
	}
	
	/// Set a Point3i parameter by name
	inline void setParameter(const std::string &name, const Point3i &value,
		bool failIfMissing = true) { 
		setParameter(getParameterID(name, failIfMissing), value);
	}
	
	/// Set a Point2 parameter by name
	inline void setParameter(const std::string &name, const Point2 &value,
		bool failIfMissing = true) { 
		setParameter(getParameterID(name, failIfMissing), value);
	}

	/// Set a Point2i parameter by name
	inline void setParameter(const std::string &name, const Point2i &value,
		bool failIfMissing = true) { 
		setParameter(getParameterID(name, failIfMissing), value);
	}

	/// Set a Point4 parameter by name
	inline void setParameter(const std::string &name, const Point4 &value,
		bool failIfMissing = true) { 
		setParameter(getParameterID(name, failIfMissing), value);
	}

	/// Set a Matrix2x2 parameter by name
	inline void setParameter(const std::string &name, const Matrix2x2 &value,
		bool failIfMissing = true) { 
		setParameter(getParameterID(name, failIfMissing), value);
	}

	/// Set a Matrix3x3 parameter by name
	inline void setParameter(const std::string &name, const Matrix3x3 &value,
		bool failIfMissing = true) { 
		setParameter(getParameterID(name, failIfMissing), value);
	}

	/// Set a Matrix4x4 parameter by name
	inline void setParameter(const std::string &name, const Matrix4x4 &value,
		bool failIfMissing = true) { 
		setParameter(getParameterID(name, failIfMissing), value);
	}

	/// Set a Transform parameter by name
	inline void setParameter(const std::string &name, const Transform &value,
		bool failIfMissing = true) {
		setParameter(getParameterID(name, failIfMissing), value.getMatrix());
	}

	/// Set a Spectrum parameter (will be converted to linear RGB) by name
	inline void setParameter(const std::string &name, const Spectrum &value,
		bool failIfMissing = true) { 
		setParameter(getParameterID(name, failIfMissing), value);
	}

	/** Set a GPUTexture parameter by name. Must be executed after
	    binding the texture to a texture unit */
	inline void setParameter(const std::string &name, const GPUTexture *value,
		bool failIfMissing = true) { 
		setParameter(getParameterID(name, failIfMissing), value);
	}

	/// Set a boolean parameter
	virtual void setParameter(int id, bool value) = 0;

	/// Set a float parameter
	virtual void setParameter(int id, Float value) = 0;

	/// Set a int parameter
	virtual void setParameter(int id, int value) = 0;

	/// Set a Vector parameter
	virtual void setParameter(int id, const Vector &value) = 0;
	
	/// Set a Vector3i parameter
	virtual void setParameter(int id, const Vector3i &value) = 0;

	/// Set a Vector2 parameter
	virtual void setParameter(int id, const Vector2 &value) = 0;
	
	/// Set a Vector2i parameter
	virtual void setParameter(int id, const Vector2i &value) = 0;

	/// Set a Vector4 parameter
	virtual void setParameter(int id, const Vector4 &value) = 0;

	/// Set a Point parameter
	virtual void setParameter(int id, const Point &value) = 0;
	
	/// Set a Point3i parameter
	virtual void setParameter(int id, const Point3i &value) = 0;

	/// Set a Point2 parameter
	virtual void setParameter(int id, const Point2 &value) = 0;
	
	/// Set a Point2i parameter
	virtual void setParameter(int id, const Point2i &value) = 0;

	/// Set a Point4 parameter
	virtual void setParameter(int id, const Point4 &value) = 0;

	/// Set a Matrix2x2 parameter
	virtual void setParameter(int id, const Matrix2x2 &value) = 0;

	/// Set a Matrix3x3 parameter
	virtual void setParameter(int id, const Matrix3x3 &value) = 0;

	/// Set a Matrix4x4 parameter
	virtual void setParameter(int id, const Matrix4x4 &value) = 0;

	/// Set a Transform parameter
	inline void setParameter(int id, const Transform &value) {
		setParameter(id, value.getMatrix());
	}

	/// Set a Spectrum parameter (will be converted to linear RGB)
	virtual void setParameter(int id, const Spectrum &value) = 0;

	/** Set a GPUTexture parameter. Must be executed after
	    binding the texture to a texture unit */
	virtual void setParameter(int id, const GPUTexture *value) = 0;

	/// Return a string representation of this class
	std::string toString() const;

	/// Returns the max. number of vertices generated by the geometry shader
	inline int getMaxVertices() const { return m_maxVertices; }

	/// Set the max. number of vertices generated by the geometry shader
	inline void setMaxVertices(int maxVertices) { m_maxVertices = maxVertices; }

	/** Return whether or not this program is currently bound */
	inline bool isBound() const { return m_bound; }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GPUProgram();
protected:
	std::string m_name;
	std::string m_source[3];
	int m_maxVertices;
	bool m_bound;
};

MTS_NAMESPACE_END

#endif /* __GPUPROGRAM_H */

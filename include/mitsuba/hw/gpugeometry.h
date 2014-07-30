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
#if !defined(__MITSUBA_HW_GPUGEOMETRY_H_)
#define __MITSUBA_HW_GPUGEOMETRY_H_

#include <mitsuba/hw/renderer.h>

MTS_NAMESPACE_BEGIN

class Shader;

/** \brief Abstract geometry storage on a graphics card
 * \ingroup libhw
 */
class MTS_EXPORT_HW GPUGeometry : public Object {
public:
	/// Create an empty geometry object
	GPUGeometry(const TriMesh *mesh);

	/// Return the name of this geometry object
	inline std::string getName() const { return m_mesh->getName(); }

	/// Return the associated triangle mesh
	inline const TriMesh *getTriMesh() const { return m_mesh.get(); }

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

	/// Return an (auxiliary) shader instance associated with the geometry
	inline Shader *getShader() { return m_shader; }

	/// Return an (auxiliary) shader instance associated with the geometry
	inline const Shader *getShader() const { return m_shader; }

	/// Set an (auxiliary) shader instance associated with the geometry
	inline void setShader(Shader *shader) { m_shader = shader; }

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GPUGeometry();
protected:
	ref<const TriMesh> m_mesh;
	Shader *m_shader;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_GPUGEOMETRY_H_ */

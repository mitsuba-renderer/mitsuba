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

#if !defined(__GPUGEOMETRY_H)
#define __GPUGEOMETRY_H

#include <mitsuba/hw/renderer.h>

MTS_NAMESPACE_BEGIN

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

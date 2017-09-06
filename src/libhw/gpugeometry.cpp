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

#include <mitsuba/hw/gpugeometry.h>

MTS_NAMESPACE_BEGIN

GPUGeometry::GPUGeometry(const TriMesh *mesh)
 : m_mesh(mesh), m_shader(NULL) {
}

GPUGeometry::~GPUGeometry() {
}

std::string GPUGeometry::toString() const {
    std::ostringstream oss;
    oss << "GPUGeometry[" << endl
        << "  mesh = " << indent(m_mesh->toString()) << endl
        << "]";
    return oss.str();
}

MTS_IMPLEMENT_CLASS(GPUGeometry, true, Object)
MTS_NAMESPACE_END

/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

GPUProgram::GPUProgram(const std::string &name)
 : m_name(name), m_maxVertices(0), m_bound(false) {
}

GPUProgram::~GPUProgram() {
}

std::string GPUProgram::toString() const {
	std::ostringstream oss;
	oss << "GPUProgram[name = '" << m_name<< "'";
	if (m_maxVertices != 0)
		oss << ", maxVertices=" << m_maxVertices;
	oss << "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(GPUProgram, true, Object)
MTS_NAMESPACE_END

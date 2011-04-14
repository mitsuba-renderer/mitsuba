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

#include <mitsuba/render/rfilter.h>

MTS_NAMESPACE_BEGIN

/**
 * Box filter -- fastest, but prone to aliasing.
 */
class BoxFilter : public ReconstructionFilter {
public:
	BoxFilter(const Properties &props) 
		: ReconstructionFilter(props) {
		m_size = Vector2(0.5f, 0.5f);
	}

	BoxFilter(Stream *stream, InstanceManager *manager) 
		: ReconstructionFilter(stream, manager) {
		m_size = Vector2(0.5f, 0.5f);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		ReconstructionFilter::serialize(stream, manager);
	}

	Float evaluate(Float x, Float y) const {
		return 1.0f;
	}
	
	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(BoxFilter, false, ReconstructionFilter);
MTS_EXPORT_PLUGIN(BoxFilter, "Box filter");
MTS_NAMESPACE_END

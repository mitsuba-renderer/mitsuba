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

#include <mitsuba/core/rfilter.h>

MTS_NAMESPACE_BEGIN

/**
 * Box filter: this is the fastest, but also about the worst possible
 * reconstruction filter, since it is extremely prone to aliasing.
 *
 * It is included mainly for completeness, though some rare situations
 * may warrant its use.
 */
class BoxFilter : public ReconstructionFilter {
public:
	BoxFilter(const Properties &props)
		: ReconstructionFilter(props) {
		/* Filter radius in pixels. A tiny epsilon is added, since some
		   samplers (Hammersley and Halton in particular) place samples
		   at positions like (0, 0). Without such an epsilon and rounding
		   errors, samples may end up not contributing to any pixel. */
		m_radius = props.getFloat("radius", 0.5f) + 1e-5f;
	}

	BoxFilter(Stream *stream, InstanceManager *manager)
		: ReconstructionFilter(stream, manager) {
		configure();
	}

	Float eval(Float x) const {
		return std::abs(x) <= m_radius ? 1.0f : 0.0f;
	}

	std::string toString() const {
		return formatString("BoxFilter[radius=%f]", m_radius);
	}

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(BoxFilter, false, ReconstructionFilter);
MTS_EXPORT_PLUGIN(BoxFilter, "Box filter");
MTS_NAMESPACE_END

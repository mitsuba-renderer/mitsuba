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
 * Simple tent (triangular) filter. This reconstruction filter never
 * suffers from ringing and usually causes less aliasing than a naive
 * box filter. When rendering scenes with sharp brightness discontinuities,
 * this may be useful; otherwise, negative-lobed filters may be preferable
 * (e.g. Mitchell-Netravali or Lanczos Sinc)
 */
class TentFilter : public ReconstructionFilter {
public:
    TentFilter(const Properties &props)
        : ReconstructionFilter(props) {
        m_radius = 1.0f;
    }

    TentFilter(Stream *stream, InstanceManager *manager)
        : ReconstructionFilter(stream, manager) {
        configure();
    }

    Float eval(Float x) const {
        return std::max((Float) 0.0f, 1.0f - std::abs(x / m_radius));
    }

    std::string toString() const {
        return formatString("TentFilter[radius=%f]", m_radius);
    }

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(TentFilter, false, ReconstructionFilter);
MTS_EXPORT_PLUGIN(TentFilter, "Tent filter");
MTS_NAMESPACE_END

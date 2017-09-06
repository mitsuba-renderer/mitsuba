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
 * This is a windowed version of the theoretically optimal low-pass filter.
 * It is generally one of the best available filters in terms of producing sharp
 * high-quality output. Its main disadvantage is that it produces ringing around
 * discontinuities, which can become a serious problem when rendering bright objects
 * with sharp edges (a directly visible light source will for instance have black
 * fringing artifacts around it).
 */
class LanczosSincFilter : public ReconstructionFilter {
public:
    LanczosSincFilter(const Properties &props)
        : ReconstructionFilter(props) {
        m_radius = (Float) props.getInteger("lobes", 3);
    }

    LanczosSincFilter(Stream *stream, InstanceManager *manager)
        : ReconstructionFilter(stream, manager) {
        configure();
    }

    Float eval(Float x) const {
        x = std::abs(x);

        if (x < Epsilon)
            return 1.0f;
        else if (x > m_radius)
            return 0.0f;

        Float x1 = M_PI * x;
        Float x2 = x1 / m_radius;

        return (std::sin(x1) * std::sin(x2)) / (x1 * x2);
    }

    std::string toString() const {
        return formatString("LanczosSincFilter[lobes=%f]", m_radius);
    }

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(LanczosSincFilter, false, ReconstructionFilter);
MTS_EXPORT_PLUGIN(LanczosSincFilter, "Lanczos Sinc filter");
MTS_NAMESPACE_END

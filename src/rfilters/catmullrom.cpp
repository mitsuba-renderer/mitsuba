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
 * Special version of the Mitchell-Netravali filter with constants B and C configured
 * to match the Catmull-Rom spline. It usually does a better job at at preserving sharp
 * features at the cost of more ringing.
 */
class CatmullRomFilter : public ReconstructionFilter {
public:
    CatmullRomFilter(const Properties &props)
        : ReconstructionFilter(props) {
        m_radius = 2.0f;
    }

    CatmullRomFilter(Stream *stream, InstanceManager *manager)
        : ReconstructionFilter(stream, manager) {
        configure();
    }

    Float eval(Float x) const {
        x = std::abs(x);

        Float x2 = x*x, x3 = x2*x;
        Float B = 0.0f, C = 0.5f;

        if (x < 1) {
            return 1.0f/6.0f * ((12-9*B-6*C)*x3
                    + (-18+12*B+6*C) * x2 + (6-2*B));
        } else if (x < 2) {
            return 1.0f/6.0f * ((-B-6*C)*x3 + (6*B+30*C) * x2
                    + (-12*B-48*C) * x + (8*B + 24*C));
        } else {
            return 0.0f;
        }
    }

    std::string toString() const {
        return formatString("CatmullRomFilter[radius=%f]", m_radius);
    }

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(CatmullRomFilter, false, ReconstructionFilter);
MTS_EXPORT_PLUGIN(CatmullRomFilter, "Catmull-Rom filter");
MTS_NAMESPACE_END

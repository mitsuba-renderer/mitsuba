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
 * Separable cubic spline reconstruction filter by Mitchell and Netravali.
 * This is often a good compromise between sharpness and ringing.
 *
 * D. Mitchell, A. Netravali, Reconstruction filters for computer graphics,
 * Proceedings of SIGGRAPH 88, Computer Graphics 22(4), pp. 221-228, 1988.
 */
class MitchellNetravaliFilter : public ReconstructionFilter {
public:
    MitchellNetravaliFilter(const Properties &props)
        : ReconstructionFilter(props) {
        /* Filter radius */
        m_radius = 2.0f;
        /* B parameter from the paper */
        m_B = props.getFloat("B", 1.0f / 3.0f);
        /* C parameter from the paper */
        m_C = props.getFloat("C", 1.0f / 3.0f);
    }

    MitchellNetravaliFilter(Stream *stream, InstanceManager *manager)
        : ReconstructionFilter(stream, manager) {
        m_B = stream->readFloat();
        m_C = stream->readFloat();
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        ReconstructionFilter::serialize(stream, manager);
        stream->writeFloat(m_B);
        stream->writeFloat(m_C);
    }

    Float eval(Float x) const {
        x = std::abs(x);

        Float x2 = x*x, x3 = x2*x;

        if (x < 1) {
            return 1.0f/6.0f * ((12-9*m_B-6*m_C)*x3
                    + (-18+12*m_B+6*m_C) * x2 + (6-2*m_B));
        } else if (x < 2) {
            return 1.0f/6.0f * ((-m_B-6*m_C)*x3 + (6*m_B+30*m_C) * x2
                    + (-12*m_B-48*m_C) * x + (8*m_B + 24*m_C));
        } else {
            return 0.0f;
        }
    }

    std::string toString() const {
        return formatString("MitchellNetravaliFilter[radius=%f, B=%f, C=%f]", m_radius, m_B, m_C);
    }

    MTS_DECLARE_CLASS()
protected:
    Float m_B, m_C;
};

MTS_IMPLEMENT_CLASS_S(MitchellNetravaliFilter, false, ReconstructionFilter);
MTS_EXPORT_PLUGIN(MitchellNetravaliFilter, "Mitchell-Netravali filter");
MTS_NAMESPACE_END

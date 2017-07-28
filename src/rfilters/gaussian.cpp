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
 * This is a windowed Gaussian filter with configurable standard deviation.
 * It often produces pleasing results, but may introduce too much blurring.
 *
 * When no reconstruction filter is explicitly requested, this is the default
 * choice in Mitsuba.
 */
class GaussianFilter : public ReconstructionFilter {
public:
    GaussianFilter(const Properties &props)
        : ReconstructionFilter(props) {
        /* Standard deviation */
        m_stddev = props.getFloat("stddev", 0.5f);

        /* Cut off after 4 standard deviations */
        m_radius = 4 * m_stddev;
    }

    GaussianFilter(Stream *stream, InstanceManager *manager)
        : ReconstructionFilter(stream, manager) {
        m_stddev = stream->readFloat();
        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        ReconstructionFilter::serialize(stream, manager);
        stream->writeFloat(m_stddev);
    }

    Float eval(Float x) const {
        Float alpha = -1.0f / (2.0f * m_stddev*m_stddev);
        return std::max((Float) 0.0f,
            math::fastexp(alpha * x * x) -
            math::fastexp(alpha * m_radius * m_radius));
    }

    std::string toString() const {
        return formatString("GaussianFilter[stddev=%f, radius=%f]", m_stddev, m_radius);
    }

    MTS_DECLARE_CLASS()
protected:
    Float m_stddev;
};

MTS_IMPLEMENT_CLASS_S(GaussianFilter, false, ReconstructionFilter);
MTS_EXPORT_PLUGIN(GaussianFilter, "Gaussian reconstruction filter");
MTS_NAMESPACE_END

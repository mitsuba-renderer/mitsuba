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
 * Windowed Gaussian filter with configurable extent
 * and standard deviation. Often produces pleasing 
 * results, but may introduce too much blurring.
 */
class GaussianFilter : public ReconstructionFilter {
public:
	GaussianFilter(const Properties &props) 
		: ReconstructionFilter(props) {
		/* Half filter size */
		Float halfSize = props.getFloat("halfSize", 2.0f);
		/* Standard deviation of the Gaussian */
		Float stddev = props.getFloat("stddev", 0.5f);

		/* Exponent multiplicator */
		m_alpha = 1 / (2*stddev*stddev);
		m_size = Vector2(halfSize, halfSize);

		/* Negative offset pre-computation */
		m_const = std::fastexp(-m_alpha * m_size.x * m_size.x);
	}

	GaussianFilter(Stream *stream, InstanceManager *manager) 
		: ReconstructionFilter(stream, manager) {
		Float halfSize = stream->readFloat();
		m_size = Vector2(halfSize, halfSize);
		m_alpha = stream->readFloat();
		m_const = stream->readFloat();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		ReconstructionFilter::serialize(stream, manager);
		stream->writeFloat(m_size.x);
		stream->writeFloat(m_alpha);
		stream->writeFloat(m_const);
	}

	virtual ~GaussianFilter() {
	}

	Float evaluate(Float x, Float y) const {
		return std::max((Float) 0.0f, std::fastexp(-m_alpha * x * x) - m_const)
			 * std::max((Float) 0.0f, std::fastexp(-m_alpha * y * y) - m_const);
	}

	MTS_DECLARE_CLASS()
protected:
	Float m_alpha, m_const;
};

MTS_IMPLEMENT_CLASS_S(GaussianFilter, false, ReconstructionFilter);
MTS_EXPORT_PLUGIN(GaussianFilter, "Gaussian reconstruction filter");
MTS_NAMESPACE_END

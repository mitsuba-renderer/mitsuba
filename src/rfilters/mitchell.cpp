#include <mitsuba/render/rfilter.h>

MTS_NAMESPACE_BEGIN

/**
 * Separable reconstruction filter by Mitchell and Netravali
 * 
 * D. Mitchell, A. Netravali, Reconstruction filters for computer graphics, 
 * Proceedings of SIGGRAPH 88, Computer Graphics 22(4), pp. 221-228, 1988.
 */
class MitchellNetravaliFilter : public ReconstructionFilter {
public:
	MitchellNetravaliFilter(const Properties &props) 
		: ReconstructionFilter(props) {
		/* Filter size in pixels */
		Float halfSize = props.getFloat("halfSize", 2.0f);
		m_size = Vector2(halfSize, halfSize);
		/* B parameter from the paper */
		m_B = props.getFloat("B", 1.0f / 3.0f);
		/* C parameter from the paper */
		m_C = props.getFloat("C", 1.0f / 3.0f);
	}

	MitchellNetravaliFilter(Stream *stream, InstanceManager *manager) 
		: ReconstructionFilter(stream, manager) {
		Float halfSize = stream->readFloat();
		m_size = Vector2(halfSize, halfSize);
		m_B = stream->readFloat();
		m_C = stream->readFloat();
	}
 
	void serialize(Stream *stream, InstanceManager *manager) const {
		ReconstructionFilter::serialize(stream, manager);
		stream->writeFloat(m_size.x);
		stream->writeFloat(m_B);
		stream->writeFloat(m_C);
	}

	Float evaluate(Float x, Float y) const {
		return mitchellNetravali(2.0f * x / m_size.x)
			 * mitchellNetravali(2.0f * y / m_size.y);
	}

	Float mitchellNetravali(Float x) const {
		x = std::abs(x);
		Float xSquared = x*x, xCubed = xSquared*x;

		if (x < 1) {
			return 1.0f/6.0f * ((12-9*m_B-6*m_C)*xCubed 
					+ (-18+12*m_B+6*m_C) * xSquared + (6-2*m_B));
		} else if (x < 2) {
			return 1.0f/6.0f * ((-m_B-6*m_C)*xCubed + (6*m_B+30*m_C) * xSquared
					+ (-12*m_B-48*m_C) * x + (8*m_B + 24*m_C));
		} else {
			return 0.0f;
		}
	}

	MTS_DECLARE_CLASS()
protected:
	Float m_B, m_C;
};

MTS_IMPLEMENT_CLASS_S(MitchellNetravaliFilter, false, ReconstructionFilter);
MTS_EXPORT_PLUGIN(MitchellNetravaliFilter, "Mitchell-Netravali filter");
MTS_NAMESPACE_END

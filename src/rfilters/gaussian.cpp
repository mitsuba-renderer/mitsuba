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
		m_const = std::exp(-m_alpha * m_size.x * m_size.x);
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
		return std::max((Float) 0.0f, std::exp(-m_alpha * x * x) - m_const)
			 * std::max((Float) 0.0f, std::exp(-m_alpha * y * y) - m_const);
	}

	MTS_DECLARE_CLASS()
protected:
	Float m_alpha, m_const;
};

MTS_IMPLEMENT_CLASS_S(GaussianFilter, false, ReconstructionFilter);
MTS_EXPORT_PLUGIN(GaussianFilter, "Gaussian reconstruction filter");
MTS_NAMESPACE_END

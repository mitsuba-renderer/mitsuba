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

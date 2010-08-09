#include <mitsuba/render/volume.h>

MTS_NAMESPACE_BEGIN

class ConstantDataSource : public VolumeDataSource {
public:
	ConstantDataSource(const Properties &props) 
		: VolumeDataSource(props) {
		m_type = props.getType("value");

		if (m_type == Properties::EFloat)
			m_float = props.getFloat("value");
		else if (m_type == Properties::EPoint)
			m_vector = props.getVector("value");
		else if (m_type == Properties::ESpectrum)
			m_spectrum = props.getSpectrum("value");
		else
			Log(EError, "The value of a 'constvolume' must have "
				"one of the following types: float, vector, spectrum");
	}

	ConstantDataSource(Stream *stream, InstanceManager *manager)
		: VolumeDataSource(stream, manager) {
		m_type = stream->readInt();
		if (m_type == Properties::EFloat)
			m_float = stream->readFloat();
		else if (m_type == Properties::EPoint)
			m_vector = Vector(stream);
		else if (m_type == Properties::ESpectrum)
			m_spectrum = Spectrum(stream);
		else
			Log(EError, "Internal error - unknown data type");
	}

	virtual ~ConstantDataSource() {
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		VolumeDataSource::serialize(stream, manager);
		stream->writeInt(m_type);
		if (m_type == Properties::ESpectrum)
			m_spectrum.serialize(stream);
		else if (m_type == Properties::EFloat)
			stream->writeFloat(m_float);
		else if (m_type == Properties::EPoint)
			m_vector.serialize(stream);
		else
			Log(EError, "Internal error - unknown data type");
	}

	Float lookupFloat(const Point &p) const {
		return m_float;
	}

	Spectrum lookupSpectrum(const Point &p) const {
		return m_spectrum;
	}
	
	Vector lookupVector(const Point &p) const {
		return m_vector;
	}

	bool supportsFloatLookups() const {
		return m_type == Properties::EFloat;
	}

	bool supportsSpectrumLookups() const {
		return m_type == Properties::ESpectrum;
	}

	bool supportsVectorLookups() const {
		return m_type == Properties::EPoint;
	}

	Float getStepSize() const {
		return std::numeric_limits<Float>::infinity();
	}

	MTS_DECLARE_CLASS()
protected:
	int m_type;
	Float m_float;
	Vector m_vector;
	Spectrum m_spectrum;
};

MTS_IMPLEMENT_CLASS_S(ConstantDataSource, false, VolumeDataSource);
MTS_EXPORT_PLUGIN(ConstantDataSource, "Constant data source");
MTS_NAMESPACE_END

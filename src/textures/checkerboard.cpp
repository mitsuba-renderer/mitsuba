#include <mitsuba/render/texture.h>

MTS_NAMESPACE_BEGIN

/**
 * Checkerboard texture
 */
class Checkerboard : public Texture {
public:
	Checkerboard(const Properties &props) : Texture(props) {
		m_reflectance = props.getSpectrum("reflectance", Spectrum(1.0f));
	}

	Checkerboard(Stream *stream, InstanceManager *manager) 
	 : Texture(stream, manager) {
		m_reflectance = Spectrum(stream);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture::serialize(stream, manager);
		m_reflectance.serialize(stream);
	}

	Spectrum getValue(const Intersection &its) const {
		int x = 2*(((int) its.uv.x) % 2) - 1, y = 2*(((int) its.uv.y) % 2) - 1;
	
		if (x*y == 1)
			return m_reflectance;
		else
			return Spectrum(0.0f);
	}
	
	bool usesRayDifferentials() const {
		return false;
	}
	
	Spectrum getAverage() const {
		return m_reflectance * .5f;
	}
	
	Spectrum getMaximum() const {
		return m_reflectance;
	}

	std::string toString() const {
		return "Checkerboard[]";
	}

	MTS_DECLARE_CLASS()
protected:
	Spectrum m_reflectance;
};

MTS_IMPLEMENT_CLASS_S(Checkerboard, false, Texture)
MTS_EXPORT_PLUGIN(Checkerboard, "Checkerboard texture");
MTS_NAMESPACE_END

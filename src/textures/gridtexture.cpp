#include <mitsuba/render/texture.h>

MTS_NAMESPACE_BEGIN

/**
 * Grid texture
 */
class GridTexture : public Texture {
public:
	GridTexture(const Properties &props) : Texture(props) {
		m_reflectance = props.getSpectrum("reflectance", Spectrum(.5f));
		m_width = props.getFloat("width", .01f);
	}

	GridTexture(Stream *stream, InstanceManager *manager) 
	 : Texture(stream, manager) {
		m_reflectance = Spectrum(stream);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture::serialize(stream, manager);
		m_reflectance.serialize(stream);
	}

	Spectrum getValue(const Intersection &its) const {
		Float x = its.uv.x - (int) its.uv.x;
		Float y = its.uv.y - (int) its.uv.y;

		if (x > .5)
			x-=1;
		if (y > .5)
			y-=1;

		if (std::abs(x) < m_width || std::abs(y) < m_width)
			return Spectrum(0.0f);
		else
			return m_reflectance;
	}

	bool usesRayDifferentials() const {
		return false;
	}

	Spectrum getMaximum() const {
		return m_reflectance;
	}

	Spectrum getAverage() const {
		return m_reflectance; // that's not quite right
	}

	std::string toString() const {
		return "GridTexture[]";
	}

	MTS_DECLARE_CLASS()
protected:
	Spectrum m_reflectance;
	Float m_width;
};

MTS_IMPLEMENT_CLASS_S(GridTexture, false, Texture)
MTS_EXPORT_PLUGIN(GridTexture, "Grid texture");
MTS_NAMESPACE_END

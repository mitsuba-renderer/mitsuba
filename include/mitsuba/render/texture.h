#if !defined(__TEXTURE_H)
#define __TEXTURE_H

#include <mitsuba/render/records.h>
#include <mitsuba/render/shader.h>

MTS_NAMESPACE_BEGIN

class MTS_EXPORT_RENDER Texture : public ConfigurableObject, public HWResource {
public:
	virtual Spectrum getValue(const Intersection &its) const = 0;
	virtual Spectrum getAverage() const = 0;
	virtual Spectrum getMaximum() const = 0; /* Component-wise maximum */
	virtual bool usesRayDifferentials() const = 0;
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	MTS_DECLARE_CLASS()
protected:
	Texture(const Properties &props);
	Texture(Stream *stream, InstanceManager *manager);

	virtual ~Texture();
};

class MTS_EXPORT_RENDER ConstantTexture : public Texture {
public:
	inline ConstantTexture(const Spectrum &value)
		: Texture(Properties()), m_value(value) {
	}

	ConstantTexture(Stream *stream, InstanceManager *manager);

	inline Spectrum getValue(const Intersection &its) const {
		return m_value;
	}
	
	inline Spectrum getAverage() const {
		return m_value;
	}
	
	inline Spectrum getMaximum() const {
		return m_value;
	}

	inline std::string toString() const {
		std::ostringstream oss;
		oss << "ConstantTexture[value=" << m_value.toString() << "]";
		return oss.str();
	}

	inline bool usesRayDifferentials() const {
		return false;
	}

	Shader *createShader(Renderer *renderer) const;

	void serialize(Stream *stream, InstanceManager *manager) const;

	MTS_DECLARE_CLASS()
protected:
	Spectrum m_value;
};

MTS_NAMESPACE_END

#endif /* __TEXTURE_H */

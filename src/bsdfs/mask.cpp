#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

MTS_NAMESPACE_BEGIN

/**
 * Applies a transparency mask to a nested BSDF
 */
class Mask : public BSDF {
public:
	Mask(const Properties &props) 
		: BSDF(props) {
		m_opacity = new ConstantTexture(props.getSpectrum("opacity", Spectrum(.5f)));
	}

	Mask(Stream *stream, InstanceManager *manager) 
		: BSDF(stream, manager) {
		m_opacity = static_cast<Texture *>(manager->getInstance(stream));
		m_nestedBSDF = static_cast<BSDF *>(manager->getInstance(stream));
		configure();
	}

	virtual ~Mask() {
		if (m_type)
			delete m_type;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_opacity.get());
		manager->serialize(stream, m_nestedBSDF.get());
	}

	void configure() {
		if (!m_nestedBSDF)
			Log(EError, "A child BSDF is required");
		m_combinedType = m_nestedBSDF->getType() | EDeltaTransmission;
		m_usesRayDifferentials = m_nestedBSDF->usesRayDifferentials();
		m_componentCount = m_nestedBSDF->getComponentCount() + 1;
		m_type = new unsigned int[m_componentCount];
		for (int i=0; i<m_nestedBSDF->getComponentCount(); ++i)
			m_type[i] = m_nestedBSDF->getType(i);
		m_type[m_nestedBSDF->getComponentCount()] = EDeltaTransmission;
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return m_nestedBSDF->getDiffuseReflectance(its) * (m_opacity->getValue(its).getLuminance());
	}

	Spectrum f(const BSDFQueryRecord &bRec) const {
		return m_nestedBSDF->f(bRec) * (m_opacity->getValue(bRec.its).getLuminance());
	}

	Spectrum fDelta(const BSDFQueryRecord &bRec) const {
		return Spectrum(1 - m_opacity->getValue(bRec.its).getLuminance());
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		return m_nestedBSDF->pdf(bRec) * (m_opacity->getValue(bRec.its).getLuminance());
	}

	Float pdfDelta(const BSDFQueryRecord &bRec) const {
		return (1 - m_opacity->getValue(bRec.its).getLuminance())
				* std::abs(Frame::cosTheta(bRec.wo));
	}

	inline void transmit(const Vector &wi, Vector &wo) const {
		wo = Vector(-wi.x, -wi.y, -wi.z);
	}

	Spectrum sample(BSDFQueryRecord &bRec) const {
		Float probBSDF = m_opacity->getValue(bRec.its).getLuminance();
		Spectrum result(0.0f);

		bool sampleTransmission = bRec.typeMask & EDeltaTransmission
			&& (bRec.sampledComponent == -1 || bRec.sampledComponent == 
				m_nestedBSDF->getComponentCount());
		bool sampleNested = bRec.sampledComponent == -1 || 
			bRec.sampledComponent < m_nestedBSDF->getComponentCount();

		if (sampleTransmission && sampleNested) {
			if (bRec.sample.x <= probBSDF) {
				bRec.sample.x /= probBSDF;
				result = m_nestedBSDF->sample(bRec);
			} else {
				transmit(bRec.wi, bRec.wo);
				bRec.sampledComponent = m_nestedBSDF->getComponentCount();
				bRec.sampledType = EDeltaTransmission;
				result = Spectrum(1/std::abs(Frame::cosTheta(bRec.wo)));
			}
		} else if (sampleTransmission) {
			transmit(bRec.wi, bRec.wo);
			bRec.sampledComponent = m_nestedBSDF->getComponentCount();
			bRec.sampledType = EDeltaTransmission;
			result = Spectrum(1 - probBSDF) / std::abs(Frame::cosTheta(bRec.wo));
		} else if (sampleNested) {
			result = m_nestedBSDF->sample(bRec);
		}

		return result;
	}

	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf) const {
		Float probBSDF = m_opacity->getValue(bRec.its).getLuminance();
		Spectrum result(0.0f);

		bool sampleTransmission = bRec.typeMask & EDeltaTransmission
			&& (bRec.sampledComponent == -1 || bRec.sampledComponent == 
				m_nestedBSDF->getComponentCount());
		bool sampleNested = bRec.sampledComponent == -1 || 
			bRec.sampledComponent < m_nestedBSDF->getComponentCount();

		if (sampleTransmission && sampleNested) {
			if (bRec.sample.x <= probBSDF) {
				bRec.sample.x /= probBSDF;
				result = m_nestedBSDF->sample(bRec, pdf) * probBSDF;
				pdf *= probBSDF;
			} else {
				transmit(bRec.wi, bRec.wo);
				bRec.sampledComponent = m_nestedBSDF->getComponentCount();
				bRec.sampledType = EDeltaTransmission;
				pdf = (1 - probBSDF) * std::abs(Frame::cosTheta(bRec.wo));
				result = Spectrum(1 - probBSDF);
			}
		} else if (sampleTransmission) {
			transmit(bRec.wi, bRec.wo);
			bRec.sampledComponent = m_nestedBSDF->getComponentCount();
			bRec.sampledType = EDeltaTransmission;
			pdf = std::abs(Frame::cosTheta(bRec.wo));
			result = Spectrum(1 - probBSDF);
		} else if (sampleNested) {
			result = m_nestedBSDF->sample(bRec, pdf);
		}

		return result;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(Texture::m_theClass) && name == "opacity") {
			m_opacity = static_cast<Texture *>(child);
		} else if (child->getClass()->derivesFrom(BSDF::m_theClass)) {
			m_nestedBSDF = static_cast<BSDF *>(child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	MTS_DECLARE_CLASS()
protected:
	ref<Texture> m_opacity;
	ref<BSDF> m_nestedBSDF;
};


MTS_IMPLEMENT_CLASS_S(Mask, false, BSDF)
MTS_EXPORT_PLUGIN(Mask, "Mask BSDF");
MTS_NAMESPACE_END

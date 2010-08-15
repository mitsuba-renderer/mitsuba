#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/**
 * Composite material, represents a linear combination of 
 * one or more BRDFs.
 */
class Composite : public BSDF {
public:
	Composite(const Properties &props) 
		: BSDF(props), m_bsdfWeight(NULL) {

		/* Parse the weight parameter */
		std::vector<std::string> weights = 
			tokenize(props.getString("weights", ""), " ,;");
		m_bsdfCount = weights.size();
		m_bsdfWeight = new Float[m_bsdfCount];
		m_bsdfOffset = new int[m_bsdfCount];

		Float totalWeight = 0;
		char *end_ptr = NULL;
		for (size_t i=0; i<weights.size(); ++i) {
			Float weight = strtod(weights[i].c_str(), &end_ptr);
			if (*end_ptr != '\0')
				SLog(EError, "Could not parse the BRDF weights!");
			if (weight < 0)
				SLog(EError, "Invalid BRDF weight!");
			m_bsdfWeight[i] = weight;
			totalWeight += weight;
		}

		if (totalWeight > 1)
			Log(EWarn, "Energy conservation is violated!");
	}

	Composite(Stream *stream, InstanceManager *manager) 
	 : BSDF(stream, manager), m_bsdfWeight(NULL) {
		m_bsdfCount = stream->readInt();
		m_bsdfWeight = new Float[m_bsdfCount];
		m_bsdfOffset = new int[m_bsdfCount];
		for (int i=0; i<m_bsdfCount; ++i) {
			m_bsdfWeight[i] = stream->readFloat();
			BSDF *bsdf = static_cast<BSDF *>(manager->getInstance(stream));
			bsdf->incRef();
			m_bsdfs.push_back(bsdf);
		}
		configure();
	}

	virtual ~Composite() {
		for (int i=0; i<m_bsdfCount; ++i)
			m_bsdfs[i]->decRef();
		if (m_type)
			delete[] m_type;
		if (m_bsdfWeight) {
			delete[] m_bsdfWeight;
			delete[] m_bsdfOffset;
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);
		
		stream->writeInt(m_bsdfCount);
		for (int i=0; i<m_bsdfCount; ++i) {
			stream->writeFloat(m_bsdfWeight[i]);
			manager->serialize(stream, m_bsdfs[i]);
		}
	}

	void configure() {
		m_combinedType = 0;
		m_usesRayDifferentials = false;
		m_componentCount = 0;

		if ((int) m_bsdfs.size() != m_bsdfCount)
			Log(EError, "BSDF count mismatch: %i bsdfs, but specified %i weights",
				(int) m_bsdfs.size(), m_bsdfCount);

		int offset = 0;
		for (int i=0; i<m_bsdfCount; ++i)
			m_componentCount = m_bsdfs[i]->getComponentCount();

		m_pdf = DiscretePDF(m_bsdfs.size());

		m_type = new unsigned int[m_componentCount];
		for (int i=0; i<m_bsdfCount; ++i) {
			BSDF *bsdf = m_bsdfs[i];
			m_bsdfOffset[i] = offset;
			for (int j=0; j<bsdf->getComponentCount(); ++j) {
				int componentType = bsdf->getType(j);
				m_type[offset+j] = componentType;
				m_combinedType |= componentType;
			}
			offset += bsdf->getComponentCount();
			m_usesRayDifferentials |= bsdf->usesRayDifferentials();
			m_pdf[i] = m_bsdfWeight[i];
		}
		m_pdf.build();
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		Spectrum result(0.0f);
		for (int i=0; i<m_bsdfCount; ++i)
			result+= m_bsdfs[i]->getDiffuseReflectance(its) * m_bsdfWeight[i];
		return result;
	}

	inline Spectrum f(const BSDFQueryRecord &bRec) const {
		Spectrum result(0.0f);

		if (bRec.component == -1) {
			for (int i=0; i<m_bsdfCount; ++i)
				result += m_bsdfs[i]->f(bRec) * m_bsdfWeight[i];
		} else {
			/* Pick out an individual component */
			for (int i=0; i<m_bsdfCount; ++i) {
				int component = bRec.component - m_bsdfOffset[i];
				if (component < 0 || component >= m_bsdfs[i]->getComponentCount())
					continue;

				BSDFQueryRecord bRec2(bRec);
				bRec2.component = component;
				return m_bsdfs[i]->f(bRec2) * m_bsdfWeight[i];
			}
		}

		return result;
	}

	Float pdf(const BSDFQueryRecord &bRec) const {
		Float result = 0.0f;

		if (bRec.component == -1) {
			for (int i=0; i<m_bsdfCount; ++i)
				result += m_bsdfs[i]->pdf(bRec) * m_pdf[i];
		} else {
			/* Pick out an individual component */
			for (int i=0; i<m_bsdfCount; ++i) {
				int component = bRec.component - m_bsdfOffset[i];
				if (component < 0 || component >= m_bsdfs[i]->getComponentCount())
					continue;

				BSDFQueryRecord bRec2(bRec);
				bRec2.component = component;
				return m_bsdfs[i]->pdf(bRec2);
			}
		}

		return result;
	}
	
	Spectrum sample(BSDFQueryRecord &bRec, Float &pdf) const {
		if (bRec.component == -1) {
			Float componentPDF;
			int entry = m_pdf.sampleReuse(bRec.sample.x, componentPDF);
			Spectrum result = m_bsdfs[entry]->sample(bRec, pdf);
			pdf *= componentPDF;
			return result;
		} else {
			/* Pick out an individual component */
			for (int i=0; i<m_bsdfCount; ++i) {
				int component = bRec.component - m_bsdfOffset[i];
				if (component < 0 || component >= m_bsdfs[i]->getComponentCount())
					continue;

				BSDFQueryRecord bRec2(bRec);
				bRec2.component = component;
				return m_bsdfs[i]->sample(bRec2, pdf);
			}
		}
		Log(EError, "Internal error!");
		return Spectrum(0.0f);
	}

	Spectrum sample(BSDFQueryRecord &bRec) const {
		if (bRec.component == -1) {
			Float componentPDF;
			int entry = m_pdf.sampleReuse(bRec.sample.x, componentPDF);
			Spectrum result = m_bsdfs[entry]->sample(bRec);
			result /= componentPDF;
			return result;
		} else {
			/* Pick out an individual component */
			for (int i=0; i<m_bsdfCount; ++i) {
				int component = bRec.component - m_bsdfOffset[i];
				if (component < 0 || component >= m_bsdfs[i]->getComponentCount())
					continue;

				BSDFQueryRecord bRec2(bRec);
				bRec2.component = component;
				return m_bsdfs[i]->sample(bRec2);
			}
		}
		Log(EError, "Internal error!");
		return Spectrum(0.0f);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(BSDF::m_theClass)) {
			BSDF *bsdf = static_cast<BSDF *>(child);
			m_bsdfs.push_back(bsdf);
			bsdf->incRef();
		} else {
			BSDF::addChild(name, child);
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Composite[" << endl
			<< "  bsdfs = {" << endl;
		for (size_t i=0; i<m_bsdfs.size(); ++i) 
			oss << indent(m_bsdfs[i]->toString()) << "," << endl;
		oss << "  }" << endl
			<< "  pdf = " << m_pdf.toString() << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	int m_bsdfCount;
	Float *m_bsdfWeight;
	int *m_bsdfOffset;
	std::vector<BSDF *> m_bsdfs;
	DiscretePDF m_pdf;
};

MTS_IMPLEMENT_CLASS_S(Composite, false, BSDF)
MTS_EXPORT_PLUGIN(Composite, "Composite BRDF");
MTS_NAMESPACE_END

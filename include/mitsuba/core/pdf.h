#if !defined(__PDF_H)
#define __PDF_H

MTS_NAMESPACE_BEGIN

/**
 * Utility class, which represents a discrete probability
 * distribution function along with its cumulative distribution
 * function. Can transform uniformly distributed samples so
 * that they match the stored distribution.
 */
class DiscretePDF {
public:
	/// Allocate a PDF with the given number of entries
	explicit inline DiscretePDF(int nEntries = 0) : m_ready(false) {
		m_pdf.resize(nEntries);
		m_cdf.resize(nEntries+1);
	}

	/// Append a PDF entry
	inline void put(Float pdfValue) {
		m_pdf.push_back(pdfValue);
		m_cdf.push_back(0.0f);
	}

	/// Return the amount of entries
	inline size_t size() const {
		return m_pdf.size();
	}

	/// Access a PDF entry
	inline Float &operator[](unsigned int entry) {
		return m_pdf[entry];
	}

	/// Access a PDF entry
	inline const Float &operator[](unsigned int entry) const {
		return m_pdf[entry];
	}

	/// Has the CDF been built?
	inline bool isReady() const {
		return m_ready;
	}

	/// Return the original (unnormalized) sum of all PDF entries
	inline Float getOriginalSum() const {
		return m_originalSum;
	}

	/**
	 * Normalize the PDF and build the associated cumulative 
	 * distribution function. Returns the sum of all unnormalized 
	 * PDF values.
	 */
	inline Float build() {
		SAssert(m_pdf.size() > 0 && !m_ready);
		m_cdf[0] = 0.0f;
		for (unsigned int i=1; i<m_cdf.size(); ++i)
			m_cdf[i] = m_cdf[i-1] + m_pdf[i-1];
		m_originalSum = m_cdf[m_cdf.size()-1];
		for (unsigned int i=0; i<m_pdf.size(); ++i) {
			m_cdf[i] /= m_originalSum;
			m_pdf[i] /= m_originalSum;
		}
		m_cdf[m_cdf.size()-1] = 1.0f;
		m_ready = true;
		return m_originalSum;
	}

	/**
	 * Transform a uniformly distributed sample. Returns the
	 * PDF entry index.
	 */
	inline int sample(Float sampleValue) const {
		std::vector<Float>::const_iterator entry = 
				std::lower_bound(m_cdf.begin(), m_cdf.end(), sampleValue);
		int index = std::max(0, (int) (entry - m_cdf.begin()) - 1);
		return std::min(index, (int) m_pdf.size()-1); // should sampleValue be > 1
	}

	/**
	 * Transform a uniformly distributed sample. Returns the
	 * PDF entry index and the probability value at that index
	 */
	inline int sample(Float sampleValue, Float &pdf) const {
		int index = sample(sampleValue);
		pdf = m_pdf[index];
		return index;
	}

	/**
	 * Transform a uniformly distributed sample. Returns the
	 * PDF entry index. The original sample is transformed so 
	 * that it can be re-used.
	 */
	inline int sampleReuse(Float &sampleValue) const {
		int index = sample(sampleValue);
		sampleValue = (sampleValue - m_cdf[index])
			/ (m_cdf[index + 1] - m_cdf[index]);
		return index;
	}

	/**
	 * Transform a uniformly distributed sample. Returns the
	 * PDF entry index and the probability value at that index.
	 * The original sample is transformed so that it can be re-used.
	 */
	inline int sampleReuse(Float &sampleValue, Float &pdf) const {
		int index = sample(sampleValue, pdf);
		sampleValue = (sampleValue - m_cdf[index])
			/ (m_cdf[index + 1] - m_cdf[index]);
		return index;
	}

	/**
	 * Print the underlying CDF
	 */
	std::string toString() const {
		std::ostringstream oss;
		oss << "DiscretePDF[originalSum=" << m_originalSum << ", ready=" 
			<< (int) m_ready << ", cdf={";
		for (size_t i=0; i<m_cdf.size(); ++i) {
			oss << m_cdf[i];
			if (i != m_cdf.size()-1)
				oss << ", ";
		}
		oss << "}]";
		return oss.str();
	}
private:
	std::vector<Float> m_pdf, m_cdf;
	Float m_originalSum;
	bool m_ready;
};

MTS_NAMESPACE_END

#endif /* __PDF_H */

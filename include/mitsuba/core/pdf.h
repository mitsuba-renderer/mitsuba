/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__PDF_H)
#define __PDF_H

MTS_NAMESPACE_BEGIN

/**
 * \brief Stores a discrete probability distribution
 * 
 * This class can be used to transform uniformly distributed samples
 * so that they match the stored distribution.
 * 
 * \ingroup libcore
 */
struct DiscretePDF {
public:
	/// Allocate a PDF with the given number of entries
	explicit inline DiscretePDF(size_t nEntries = 0) : m_ready(false) {
		m_pdf.resize(nEntries);
		m_cdf.resize(nEntries+1);
	}

	/// Reserve memory for a certain number of entries
	inline void reserve(size_t nEntries) {
		m_pdf.reserve(nEntries);
		m_cdf.reserve(nEntries+1);
	}

	/// Append a PDF entry
	inline void put(Float pdfValue) {
		m_pdf.push_back(pdfValue);
		m_cdf.push_back(0.0f);
	}

	/// Return the number of entries
	inline size_t size() const {
		return m_pdf.size();
	}

	/// Access a PDF entry
	inline Float &operator[](size_t entry) {
		return m_pdf[entry];
	}

	/// Access a PDF entry
	inline const Float &operator[](size_t entry) const {
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
	 * \brief Normalize the PDF and build the associated cumulative 
	 * distribution function.
	 * \return Sum of all unnormalized PDF values
	 */
	inline Float build() {
		SAssert(m_pdf.size() > 0);
		m_cdf[0] = 0.0f;
		for (size_t i=1; i<m_cdf.size(); ++i)
			m_cdf[i] = m_cdf[i-1] + m_pdf[i-1];
		m_originalSum = m_cdf[m_cdf.size()-1];
		for (size_t i=0; i<m_pdf.size(); ++i) {
			m_cdf[i] /= m_originalSum;
			m_pdf[i] /= m_originalSum;
		}
		m_cdf[m_cdf.size()-1] = 1.0f;
		m_ready = true;
		return m_originalSum;
	}

	/**
	 * \brief %Transform a uniformly distributed sample
	 * \param[in] sampleValue Uniform sample
	 * \return Sample index 
	 */
	inline size_t sample(Float sampleValue) const {
		std::vector<Float>::const_iterator entry = 
				std::lower_bound(m_cdf.begin(), m_cdf.end(), sampleValue);
		size_t index = (size_t) std::max((ptrdiff_t) 0, entry - m_cdf.begin() - 1);
		return std::min(index, m_pdf.size()-1); // should sampleValue be > 1
	}

	/**
	 * \brief %Transform a uniformly distributed sample. 
	 * \param[in] sampleValue Uniform sample
	 * \param[out] pdf Probability value of the sample
	 * \return Sample index 
	 */
	inline size_t sample(Float sampleValue, Float &pdf) const {
		size_t index = sample(sampleValue);
		pdf = m_pdf[index];
		return index;
	}

	/**
	 * \brief %Transform a uniformly distributed sample. 
	 * 
	 * The original sample is adjusted so that it can be reused.
	 * \param[in,out] sampleValue Uniform sample
	 * \return Sample index 
	 */
	inline size_t sampleReuse(Float &sampleValue) const {
		size_t index = sample(sampleValue);
		sampleValue = (sampleValue - m_cdf[index])
			/ (m_cdf[index + 1] - m_cdf[index]);
		return index;
	}

	/**
	 * \brief %Transform a uniformly distributed sample. 
	 * 
	 * The original sample is adjusted so that it can be reused.
	 * \param[in,out] sampleValue Uniform sample
	 * \param[out] pdf Probability value of the sample
	 * \return Sample index 
	 */
	inline size_t sampleReuse(Float &sampleValue, Float &pdf) const {
		size_t index = sample(sampleValue, pdf);
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

/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__SHEXP4D_H)
#define __SHEXP4D_H

#include <mitsuba/core/shvector.h>

MTS_NAMESPACE_BEGIN

/**
 * Stores a 4D function f(wi, wo) (such as a BRDF or phase function)
 * using a 2D table of spherical harmonics expansions. Discretizaiton
 * occurs in the 'wi' space. Later lookups interpolate amongst 
 * the 4 adjacent samples
 */
struct SHVector4D {
public:
	/// Construct an invalid expandion
	inline SHVector4D() : m_resTheta(0), m_resPhi(0), m_bands(0) {
	}

	/// Allocate memory for the specified resolution and number of bands
	inline SHVector4D(int resTheta, int resPhi, int bands) 
		: m_resTheta(resTheta), m_resPhi(resPhi), m_bands(bands), m_data((resTheta+1)*resPhi) {
		for (size_t i=0; i<m_data.size(); ++i)
			m_data[i] = SHVector(bands);
	}

	/// Unserialize from a stream
	inline SHVector4D(Stream *stream) {
		m_resTheta = stream->readInt();
		m_resPhi = stream->readInt();
		m_bands = stream->readInt();
		m_data.resize((m_resTheta+1)*m_resPhi);
		for (size_t i=0; i<m_data.size(); ++i)
			m_data[i] = SHVector(stream);
		SLog(EInfo, "Unserialized a %ix%i discretization (%i SH bands)",
			m_resTheta, m_resPhi, m_bands);
	}

	/// Project the given function f(wi, wo)
	template<typename Functor> void project(const Functor &f, int res = 32) {
		SAssert(res % 2 == 0);
		/* Nested composite Simpson's rule */
		Float hExt = M_PI / res,
			hInt = (2*M_PI)/(res*2);

		Vector *wi = new Vector[(m_resTheta+1) * m_resPhi];
		Float *values = new Float[(m_resTheta+1) * m_resPhi];
		for (int o=0; o<=m_resTheta; ++o)
			for (int p=0; p<m_resPhi; ++p)
				wi[o*m_resPhi+p] = sphericalDirection(o * M_PI / (Float) m_resTheta, 
					p * 2 * M_PI / (Float) m_resPhi);
		
		Float *sinPhi = (Float *) alloca(sizeof(Float)*m_bands),
			  *cosPhi = (Float *) alloca(sizeof(Float)*m_bands);

		for (int i=0; i<=res; ++i) {
			Float theta = hExt*i, cosTheta = std::cos(theta);
			Float weightExt = (i & 1) ? 4 : 2;
			if (i == 0 || i == res)
				weightExt = 1;
			SLog(EInfo, "4D projection: %i%% done", (int) (100*i/(Float) res));

			for (int j=0; j<=res*2; ++j) {
				Float phi = hInt*j;
				Float weightInt = (j & 1) ? 4 : 2;
				if (j == 0 || j == 2*res)
					weightInt = 1;

				for (int m=0; m<m_bands; ++m) {
					sinPhi[m] = std::sin((m+1)*phi);
					cosPhi[m] = std::cos((m+1)*phi);
				}

				Float weight = std::sin(theta)*weightInt*weightExt;
				Vector wo = sphericalDirection(theta, phi);

				for (int o=0; o<=m_resTheta; ++o)
					for (int p=0; p<m_resPhi; ++p)
						values[o*m_resPhi+p] = f(wi[o*m_resPhi+p], wo) * weight;

				for (int l=0; l<m_bands; ++l) {
					for (int m=1; m<=l; ++m) {
						Float L = SHVector::legendre(l, m, cosTheta) * SHVector::normalization(l, m);
						for (int o=0; o<=m_resTheta; ++o)
							for (int p=0; p<m_resPhi; ++p)
								m_data[o*m_resPhi+p](l, -m) += values[o*m_resPhi+p] * SQRT_TWO * sinPhi[m-1] * L;
						for (int o=0; o<=m_resTheta; ++o)
							for (int p=0; p<m_resPhi; ++p)
								m_data[o*m_resPhi+p](l, m) += values[o*m_resPhi+p] * SQRT_TWO * cosPhi[m-1] * L;
					}

					for (int o=0; o<=m_resTheta; ++o)
						for (int p=0; p<m_resPhi; ++p)
							m_data[o*m_resPhi+p](l, 0) += values[o*m_resPhi+p] 
								* SHVector::legendre(l, 0, cosTheta) * SHVector::normalization(l, 0);
				}
			}
		}

		for (int o=0; o<=m_resTheta; ++o)
			for (int p=0; p<m_resPhi; ++p)
				for (int l=0; l<m_bands; ++l)
					for (int m=-l; m<=l; ++m)
						m_data[o*m_resPhi+p](l,m) *= hExt*hInt/9;

		delete[] wi;
		delete[] values;
	}

	/// Offset all stored projections to mask negative areas
	inline void offsetNegativeRegions(int res) {
		for (int o=0; o<=m_resTheta; ++o) {
			for (int p=0; p<m_resPhi; ++p) {
				SHVector &vec = m_data[o*m_resPhi+p];
				Float minimum = vec.findMinimum(res);
				vec.offset(std::max((Float) 0, -minimum));
				cout << "Minimum now: " << minimum << endl;
				minimum = vec.findMinimum(res);
				cout << "Minimum was: " << minimum << endl;
			}
		}
	}

	/// Normalize all stored projections to ensure they are proper distributions
	inline void normalize() {
		for (int o=0; o<=m_resTheta; ++o)
			for (int p=0; p<m_resPhi; ++p)
				m_data[o*m_resPhi+p].normalize();
	}
	
	/// Offset all stored projections by a constant value
	inline void offset(Float value) {
		for (int o=0; o<=m_resTheta; ++o)
			for (int p=0; p<m_resPhi; ++p)
				m_data[o*m_resPhi+p].offset(value);
	}

	/// Perform an interpolated lookup and save to 'vec' (any content will be overwritten)
	inline void lookup(const Vector &wi, SHVector &vec) const {
		SAssert(vec.getBands() == m_bands);

		Point2 sph = toSphericalCoordinates(wi);
		if (sph.y<0)
			sph.y += 2*M_PI;
		Float theta = sph.x * m_resTheta/M_PI;
		Float phi = sph.y * m_resPhi/(2*M_PI);

		int i0 = (int) std::floor(theta), i1 = (int) std::ceil(theta);
		int j0 = (int) std::floor(phi), j1 = (int) std::ceil(phi);

		Float alpha = theta - (Float) i0;
		Float beta = phi - (Float) j0;

		i0 = std::min(std::max(i0, 0), m_resTheta);
		i1 = std::min(std::max(i1, 0), m_resTheta);
		j0 = std::min(std::max(j0, 0), m_resPhi);
		j1 = std::min(std::max(j1, 0), m_resPhi);

		if (j0 == m_resPhi) j0 = 0;
		if (j1 == m_resPhi) j1 = 0;

		vec.clear();
		vec.madd((1-alpha)*(1-beta), m_data[i0*m_resPhi+j0]);
		vec.madd(alpha*(1-beta), m_data[i1*m_resPhi+j0]);
		vec.madd((1-alpha)*beta, m_data[i0*m_resPhi+j1]);
		vec.madd(alpha*beta, m_data[i1*m_resPhi+j1]);
	}

	/// Serialize to a binary data stream
	inline void serialize(Stream *stream) const {
		stream->writeInt(m_resTheta);
		stream->writeInt(m_resPhi);
		stream->writeInt(m_bands);
		for (size_t i=0; i<m_data.size(); ++i)
			m_data[i].serialize(stream);
	}

	/// Return the number of SH coefficient bands
	inline int getBands() const {
		return m_bands;
	}

	/// Get the energy per band
	inline Float energy(int band) const {
		Float result = 0;
		for (int o=0; o<=m_resTheta; ++o)
			for (int p=0; p<m_resPhi; ++p)
				result += m_data[o*m_resPhi+p].energy(band);
		return result;
	}

	inline std::string toString() const {
		std::ostringstream oss;
		oss << "SHVector4D[resTheta=" << m_resTheta 
			<< ", resPhi=" << m_resPhi
			<< ", bands=" << m_bands
			<< ", size=" << m_bands*m_bands*m_data.size()*sizeof(Float) / 1024 
			<< " KiB]";
		return oss.str();
	}
private:
	int m_resTheta, m_resPhi, m_bands;
	std::vector<SHVector> m_data;
};

MTS_NAMESPACE_END

#endif /* __SHEXP4D_H */

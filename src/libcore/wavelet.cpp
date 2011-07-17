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

#include <mitsuba/core/wavelet.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN
	
/* ==================================================================== */
/*                        2D Wavelet Transform                          */
/* ==================================================================== */

Wavelet2D::Wavelet2D(const Bitmap *bitmap, int selectedChannel) {
	Assert(bitmap->getWidth() == bitmap->getHeight());
	Assert(isPowerOfTwo(bitmap->getWidth()));

	m_size = bitmap->getWidth();
	m_data = new float[m_size * m_size];
	m_temp = new float[m_size];

	int offset = selectedChannel, 
		stride = bitmap->getBitsPerPixel()/8;
	if (bitmap->getBitsPerPixel() != 128) {
		/* Byte image */
		const uint8_t *data = bitmap->getData() + offset;
		for (size_t i=0; i<m_size; i++) {
			for (size_t j=0; j<m_size; j++) {
				float value = *data / 255.0f;
				m_data[j+i*m_size] = value / m_size;
				data += stride;
			}
		}
	} else {
		/* Floating point image */
		const float *data = ((float *) bitmap->getData()) + offset;
		stride /= 4;
		for (size_t i=0; i<m_size; i++) {
			for (size_t j=0; j<m_size; j++) {
				m_data[j+i*m_size] = *data / m_size;
				data += stride;
			}
		}
	}
	nonstandardDecomposition();
}

Wavelet2D::Wavelet2D(const SparseWavelet2D *sw) {
	Assert(isPowerOfTwo(sw->getSize()));
	m_size = sw->getSize();
	m_data = new float[m_size * m_size];
	m_temp = new float[m_size];
	memset(m_data, 0, sizeof(float)*m_size*m_size);

	/* Extract the scaling function */
	m_data[0] = sw->getScalingFunction();

	/* Loop over the resolution hierarchy and extract the coefficients */
	uint8_t level = 0;
	uint16_t size = 1, offset = 1;
	while (size < m_size) {
		/* Loop over the different types of differentiation */
		for (uint8_t type=0; type<3; type++) {
			size_t offsetX = 0, offsetY = 0;
			switch (type) {
				case 0: offsetX = offset; break;
				case 1: offsetY = offset; break;
				case 2: offsetX = offset; offsetY = offset; break;
			}
			for (uint16_t j=0; j<size; j++) {
				for (uint16_t i=0; i<size; i++) {
					size_t pos = (j+offsetY)*m_size + i + offsetX;
					Assert(pos < m_size*m_size);
					m_data[pos] = sw->get(SparseWavelet2D::Key::create(level, type, i, j));
				}
			}
		}
		offset += size;
		size <<= 1;
		level += 1;
	}
}

SparseWavelet2D *Wavelet2D::toSparseWavelet() const {
	SparseWavelet2D *sparse = new SparseWavelet2D(m_size);
	sparse->setScalingFunction(m_data[0]);

	/* Loop over the resolution hierarchy and extract the coefficients */
	uint8_t level = 0;
	uint16_t size = 1, offset = 1;

	while (size < m_size) {
		/* Loop over the different types of differentiation */
		for (int type=0; type<3; type++) {
			size_t offsetX = 0, offsetY = 0;
			switch (type) {
				case 0: offsetX = offset; break;
				case 1: offsetY = offset; break;
				case 2: offsetX = offset; offsetY = offset; break;
			}

			for (uint16_t j=0; j<size; j++) {
				for (uint16_t i=0; i<size; i++) {
					size_t pos = (j+offsetY)*m_size + i + offsetX;
					Assert(pos < m_size*m_size);
					if (m_data[pos] == 0.0f)
						continue;
					sparse->put(SparseWavelet2D::Key::create(level, type, i, j),
						m_data[pos]);
				}
			}
		}
		offset += size;
		size <<= 1;
		level += 1;
	}
	return sparse;
}

void Wavelet2D::nonstandardDecomposition() {
	size_t size = m_size;
	while (size > 1) {
		for (size_t i=0; i<size; i++)
			forwardStepX(i, size);
		for (size_t i=0; i<size; i++)
			forwardStepY(i, size);
		size >>= 1;
	}
}

void Wavelet2D::decode(Bitmap *bitmap, float offset, float scale) {
	Assert(bitmap->getWidth() == bitmap->getHeight());
	Assert(bitmap->getBitsPerPixel() == 8
		|| bitmap->getBitsPerPixel() == 128);
	Assert((size_t) bitmap->getWidth() == m_size);

	/* Back up the current coefficients */
	float *backup = m_data;
	m_data = new float[m_size * m_size];
	for (size_t i=0; i<m_size*m_size; i++)
		m_data[i] = backup[i];

	/* Do an inverse non-standard wavelet transform */
	size_t size = 2;
	while (size <= m_size) {
		for (size_t i=0; i<size; i++)
			backwardStepX(i, size);
		for (size_t i=0; i<size; i++)
			backwardStepY(i, size);
		size <<= 1;
	}
	scale *= m_size;

	if (bitmap->getBitsPerPixel() == 8) {
		scale *= 255;
		/* Copy the image (byte image)*/
		for (size_t i=0; i<m_size; i++) {
			for (size_t j=0; j<m_size; j++) {
				bitmap->getData()[i+j*m_size] = std::min(255, 
					std::max(0, (int) ((m_data[i + j*m_size]+offset) * scale)));
			}
		}
	} else {
		float *data = bitmap->getFloatData();
		for (size_t i=0; i<m_size; i++) {
			for (size_t j=0; j<m_size; j++) {
				float value = std::max(0.0f, (m_data[j + i*m_size] + offset)*scale);
				*data++ = value; *data++ = value;
				*data++ = value; *data++ = 1.0f;
			}
		}
	}

	/* Restore the coefficients */
	delete[] m_data;
	m_data = backup;
}

void Wavelet2D::discard(Float fraction) {
	Assert(fraction >= 0 && fraction <= 1.0f);

	size_t nEntries = m_size*m_size;
	size_t toDiscard = (int) (fraction * nEntries);

	/* Find the minimum value cutoff */
	std::vector<coeff_t> coeffs(nEntries);
	for (size_t i=0; i<nEntries; i++)
		coeffs[i] = coeff_t(i, m_data[i]);
	std::nth_element(coeffs.begin(), coeffs.begin() + toDiscard,
		coeffs.end(), std::less<coeff_t>());

	/* Set to zero */
	float minValue = coeffs[toDiscard].value;
	for (size_t i=0; i<nEntries; i++) {
		if (std::abs(m_data[i]) <= minValue) {
			m_data[i] = 0;
		}
	}
}

Float Wavelet2D::compress(Float maxError) {
	size_t nEntries = m_size*m_size;

	std::vector<coeff_t> coeffs(nEntries);
	for (size_t i=0; i<nEntries; i++) 
		coeffs[i] = coeff_t(i, m_data[i]);

	std::sort(coeffs.begin(), coeffs.end(), std::less<coeff_t>());
	float t, tMin = std::abs(coeffs[0].value),
			 tMax = std::abs(coeffs[coeffs.size()-1].value);
	int discarded = 0;

	maxError = std::pow(m_data[0]*maxError, 2);

	do {
		t = .5f * (tMin+tMax);
		Float s = 0;

		for (size_t i=0; i<nEntries; ++i) {
			float value = coeffs[i].value;
			if (std::abs(value) > t)
				break;
			s = s + value*value;
		}

		if (s < maxError)
			tMin = t;
		else
			tMax = t;

	} while (tMax - tMin > Epsilon);

	for (size_t i=0; i<nEntries; ++i) {
		if (std::abs(coeffs[i].value) > t)
			break;
		m_data[coeffs[i].index] = 0;
		discarded++;
	}

	return discarded/(Float) nEntries;
}

void Wavelet2D::forwardStepX(size_t y, size_t size) {
	size_t halfSize = size/2;
	size_t offset = y*m_size;
	for (size_t i=0; i<halfSize; i++) {
		/* Averaging step */
		m_temp[i] =
			(m_data[offset + 2*i]
		  +  m_data[offset + 2*i+1]) * INV_SQRT_TWO;
		/* Feature extraction step */
		m_temp[halfSize + i] =
			(m_data[offset + 2*i]
		  -  m_data[offset + 2*i+1]) * INV_SQRT_TWO;	
	}
	for (size_t i=0; i<size; i++)
		m_data[offset+i] = m_temp[i];
}

void Wavelet2D::forwardStepY(size_t x, size_t size) {
	size_t halfSize = size/2;
	for (size_t i=0; i<halfSize; i++) {
		/* Averaging step */
		m_temp[i] =
			(m_data[x + (2*i)*m_size]
		  +  m_data[x + (2*i+1)*m_size]) * INV_SQRT_TWO;
		/* Feature extraction step */
		m_temp[halfSize + i] =
			(m_data[x + (2*i)*m_size]
		  -  m_data[x + (2*i+1)*m_size]) * INV_SQRT_TWO;	
	}
	for (size_t i=0; i<size; i++)
		m_data[x+i*m_size] = m_temp[i];
}

void Wavelet2D::backwardStepX(size_t y, size_t size) {
	size_t halfSize = size/2;
	size_t offset = y*m_size;
	for (size_t i=0; i<halfSize; i++) {
		m_temp[2*i]   = (m_data[offset + i] 
			+ m_data[halfSize + offset + i]) * INV_SQRT_TWO;
		m_temp[2*i+1] = (m_data[offset + i] 
			- m_data[halfSize + offset + i]) * INV_SQRT_TWO;
	}
	for (size_t i=0; i<size; i++)
		m_data[offset+i] = m_temp[i];
}

void Wavelet2D::backwardStepY(size_t x, size_t size) {
	size_t halfSize = size/2;
	for (size_t i=0; i<halfSize; i++) {
		m_temp[2*i]   = (m_data[x + i*m_size] 
			+ m_data[x + (halfSize + i)*m_size]) * INV_SQRT_TWO;
		m_temp[2*i+1] = (m_data[x + i*m_size] 
			- m_data[x + (halfSize + i)*m_size]) * INV_SQRT_TWO;
	}
	for (size_t i=0; i<size; i++)
		m_data[x+i*m_size] = m_temp[i];
}

Wavelet2D::~Wavelet2D() {
	delete[] m_data;
	delete[] m_temp;
}

/* ==================================================================== */
/*                         Sparse 2D Wavelet                            */
/* ==================================================================== */

SparseWavelet2D::SparseWavelet2D(size_t size) 
 : m_scalingFunction(0), m_size(size) {
	Assert(sizeof(Key) == 4);
#if defined(USE_GOOGLE_DENSE_HASHMAP)
	m_data.set_empty_key(0xFFFFFFFFFFFFFFFFULL);
#endif
	m_maxLevel = log2i(m_size)-1;
}

SparseWavelet2D::SparseWavelet2D(const SparseWavelet2D *sw) 
 : m_data(sw->m_data), m_scalingFunction(sw->m_scalingFunction), m_size(sw->m_size),
   m_maxLevel(sw->m_maxLevel) {
}

SparseWavelet2D::SparseWavelet2D(Stream *stream, InstanceManager *Manager) {
	m_size = stream->readSize();
	m_scalingFunction = stream->readSingle();
	size_t coefficientCount = stream->readSize();
#if defined(USE_GOOGLE_DENSE_HASHMAP)
	m_data.set_empty_key(0xFFFFFFFFFFFFFFFFULL);
#endif
	for (size_t i=0; i<coefficientCount; i++) {
		uint64_t key = stream->readULong();
		m_data[key] = stream->readSingle();
	}
	m_maxLevel = log2i(m_size)-1;
}

void SparseWavelet2D::clear() {
	m_data.clear();
#if defined(USE_GOOGLE_DENSE_HASHMAP)
	m_data.resize(0);
#endif
	m_scalingFunction = 0.0f;
}


Float SparseWavelet2D::getPixel(const Point2i &pt) const {
	Key key = Key::create(0, 0, pt.x, pt.y);
	Float value = 0;

	for (int i=m_maxLevel; i>=0; i--) {
		int offsetX = key.i & 1, offsetY = key.j & 1;
		key.i >>= 1; key.j >>= 1; key.level = i;
		for (int j=0; j<3; ++j) {
			key.type = j;
			float coeff = get(key) * key.quadrantSign(offsetX, offsetY);
			value += coeff * std::pow(2.0f, (float) i);
		}
	}
	value += m_scalingFunction;

	return value;
}

void SparseWavelet2D::serialize(Stream *stream, InstanceManager *Manager) const {
	stream->writeSize(m_size);
	stream->writeSingle(m_scalingFunction);
	stream->writeSize(m_data.size());

	for (CoefficientIterator it = m_data.begin(); it != m_data.end(); ++it) {
		stream->writeULong((*it).first);
		stream->writeSingle((*it).second);
	}
}

Float SparseWavelet2D::lineIntegral(Point2 start, Point2 end) const {
	Vector2 d = end-start;
	Float accum = 0, maxt = d.length();
	size_t res = m_size;
	Key key;

	d/= maxt;
	key.empty = 0;

	for (int level=m_maxLevel; level>=0; --level) {
		key.level = level;

		Vector2i dpos(
			std::max(0, std::min((int) start.x, (int) res-1)),
			std::max(0, std::min((int) start.y, (int) res-1))
		);

		Vector2 delta, next;
		Vector2i step;
		for (int i=0; i<2; ++i) {
			if (std::abs(d[i]) < Epsilon) {
				delta[i] = 0; step[i] = 0; next[i] =
					std::numeric_limits<Float>::infinity();
			} else if (d[i] > 0) {
				delta[i] = 1/d[i];
				next[i] = (dpos[i] + 1 - start[i]) / d[i];
				step[i] = 1;
			} else {
				delta[i] = -1/d[i];
				next[i] = (dpos[i] - start[i]) / d[i];
				step[i] = -1;
			}
		}

		Float t = 0, nextT;

		while (true) {
			nextT = std::min(std::min(next.x, next.y), maxt);
			key.i = dpos.x >> 1; key.j = dpos.y >> 1;
			float temp = 0;

			for (int type=0; type<3; ++type) {
				key.type = type;
				const float coeff = get(key),
					sign = key.quadrantSign(dpos.x - (key.i << 1), dpos.y - (key.j << 1));
				temp += coeff*sign;
			}

			accum += .5f * temp * (nextT-t);
			t = nextT;

			if (nextT == maxt) {
				break;
			} else if (next.x <= next.y) {
				dpos.x += step.x;
				next.x += delta.x;
			} else {
				dpos.y += step.y;
				next.y += delta.y;
			}
		}
		res /= 2; maxt /= 2; start /= 2; end /= 2;
	}
	accum += maxt * m_scalingFunction;

	return accum;
}

std::string SparseWavelet2D::toString() const {
	std::ostringstream oss;
	oss << "SparseWavelet2D[size=" << m_size 
		<< ", numCoefficents=" << m_data.size() << "]";
	return oss.str();
}

/* ==================================================================== */
/*                        3D Wavelet Transform                          */
/* ==================================================================== */

Wavelet3D::Wavelet3D(const float *data, size_t resolution) {
	Assert(isPowerOfTwo(resolution));

	m_size = resolution;
	size_t nEntries = m_size * m_size * m_size;
	m_data = new float[nEntries];
	m_temp = new float[m_size];
	m_slab = resolution*resolution;

	float normFactor = 1/(resolution * std::sqrt((float) resolution));

	for (size_t i=0; i<nEntries; ++i)
		m_data[i] = data[i] * normFactor;

	nonstandardDecomposition();
}


SparseWaveletOctree *Wavelet3D::toOctree() const {
	SparseWaveletOctree *octree = new SparseWaveletOctree(m_size, m_data[0]);

	/* Loop over the resolution hierarchy and extract the coefficients */
	size_t level = log2i(m_size)-1;
	size_t size = m_size/2;

	while (size > 0) {
		float factor = std::pow(2.0f, 3*level/2.0f);

		/* Loop over the different types of wavelet basis functions */
		for (size_t k=0; k<size; k++) {
			for (size_t j=0; j<size; j++) {
				for (size_t i=0; i<size; i++) {
					float coeff[7];
					bool active = false;

					for (int type=0; type<7; type++) {
						size_t offsetX = 0, offsetY = 0, offsetZ = 0;
						switch (type) {
							case 0: offsetX = size; break;
							case 1: offsetY = size; break;
							case 2: offsetX = size; offsetY = size; break;
							case 3: offsetZ = size; break;
							case 4: offsetZ = size; offsetX = size; break;
							case 5: offsetZ = size; offsetY = size; break;
							case 6: offsetZ = size; offsetX = size; offsetY = size; break;
						}

						coeff[type] = factor * m_data[(k + offsetZ) * m_slab 
							+ (j+offsetY)*m_size + i + offsetX];
						active |= (coeff[type] != 0);
					}

					if (active)
						octree->put(level, i, j, k, coeff);
				}
			}
		}
		size >>= 1;
		level -= 1;
	}

	return octree;
}


void Wavelet3D::nonstandardDecomposition() {
	size_t size = m_size;

	while (size > 1) {
		for (size_t z=0; z<size; z++) 
			for (size_t y=0; y<size; y++)
				forwardStepX(y, z, size);
	
		for (size_t z=0; z<size; z++) 
			for (size_t x=0; x<size; x++)
				forwardStepY(x, z, size);

		for (size_t y=0; y<size; y++)
			for (size_t x=0; x<size; x++)
				forwardStepZ(x, y, size);
		size >>= 1;
	}
}

void Wavelet3D::decode(float *target) {
	float *original = m_data;
	size_t nEntries = m_size * m_size * m_size;

	memcpy(target, m_data, nEntries * sizeof(float));
	m_data = target;

	size_t size = 2;
	while (size <= m_size) {
		for (size_t z=0; z<size; z++)
			for (size_t y=0; y<size; y++)
				backwardStepX(y, z, size);
		for (size_t z=0; z<size; z++)
			for (size_t x=0; x<size; x++)
				backwardStepY(x, z, size);
		for (size_t y=0; y<size; y++)
			for (size_t x=0; x<size; x++)
				backwardStepZ(x, y, size);
		size <<= 1;
	}

	float normFactor = m_size * std::sqrt((float) m_size);

	for (size_t i=0; i<nEntries; ++i)
		m_data[i] *= normFactor; 

	/* Restore the coefficients */
	m_data = original;
}

void Wavelet3D::discard(Float fraction) {
	Assert(fraction >= 0 && fraction <= 1.0f);

	size_t nEntries = m_size*m_size*m_size;
	size_t toDiscard = (int) (fraction * nEntries);

	/* Find the minimum value cutoff */
	std::vector<coeff_t> coeffs(nEntries);
	for (size_t i=0; i<nEntries; i++)
		coeffs[i] = coeff_t(i, m_data[i]);
	std::nth_element(coeffs.begin(), coeffs.begin() + toDiscard,
		coeffs.end(), std::less<coeff_t>());

	/* Set to zero */
	float minValue = coeffs[toDiscard].value;
	for (size_t i=0; i<nEntries; i++) {
		if (std::abs(m_data[i]) <= minValue)
			m_data[i] = 0;
	}
}

Float Wavelet3D::compress(Float maxError) {
	/* From: "Wavelets for computer graphics: A primer, part 1" by
	   Eric J. Stollnitz, Tony D. DeRose, and David H. Salesin
	   (IEEE Computer Graphics and Applications, May 1995) */

	size_t nEntries = m_size*m_size*m_size;

	std::vector<coeff_t> coeffs(nEntries);
	for (size_t i=0; i<nEntries; i++) 
		coeffs[i] = coeff_t(i, m_data[i]);

	std::sort(coeffs.begin(), coeffs.end(), std::less<coeff_t>());
	float t, tMin = std::abs(coeffs[0].value),
			 tMax = std::abs(coeffs[coeffs.size()-1].value);
	int discarded = 0;

	maxError = std::pow(m_data[0]*maxError, 2);

	do {
		t = .5f * (tMin+tMax);
		Float s = 0;

		for (size_t i=0; i<nEntries; ++i) {
			float value = coeffs[i].value;
			if (std::abs(value) > t)
				break;
			s = s + value*value;
		}

		if (s < maxError)
			tMin = t;
		else
			tMax = t;
	} while (tMax - tMin > Epsilon);

	for (size_t i=0; i<nEntries; ++i) {
		if (std::abs(coeffs[i].value) > t)
			break;
		m_data[coeffs[i].index] = 0;
		discarded++;
	}

	return discarded/(Float) nEntries;
}

void Wavelet3D::forwardStepX(size_t y, size_t z, size_t size) {
	size_t halfSize = size/2;
	size_t offset = y*m_size + z * m_slab;
	for (size_t i=0; i<halfSize; i++) {
		/* Averaging step */
		m_temp[i] =
			(m_data[offset + 2*i]
		  +  m_data[offset + 2*i+1]) * INV_SQRT_TWO;

		/* Detail extraction step */
		m_temp[halfSize + i] =
			(m_data[offset + 2*i]
		  -  m_data[offset + 2*i+1]) * INV_SQRT_TWO;	
	}
	for (size_t i=0; i<size; i++)
		m_data[offset+i] = m_temp[i];
}

void Wavelet3D::forwardStepY(size_t x, size_t z, size_t size) {
	size_t halfSize = size/2;
	size_t offset = z * m_slab + x;
	for (size_t i=0; i<halfSize; i++) {
		/* Averaging step */
		m_temp[i] =
			(m_data[offset + (2*i)*m_size]
		  +  m_data[offset + (2*i+1)*m_size]) * INV_SQRT_TWO;

		/* Detail extraction step */
		m_temp[halfSize + i] =
			(m_data[offset + (2*i)*m_size]
		  -  m_data[offset + (2*i+1)*m_size]) * INV_SQRT_TWO;	
	}
	for (size_t i=0; i<size; i++)
		m_data[offset+i*m_size] = m_temp[i];
}

void Wavelet3D::forwardStepZ(size_t x, size_t y, size_t size) {
	size_t halfSize = size/2;
	size_t offset = y * m_size + x;
	for (size_t i=0; i<halfSize; i++) {
		/* Averaging step */
		m_temp[i] =
			(m_data[offset + (2*i)*m_slab]
		  +  m_data[offset + (2*i+1)*m_slab]) * INV_SQRT_TWO;

		/* Detail extraction step */
		m_temp[halfSize + i] =
			(m_data[offset + (2*i)*m_slab]
		  -  m_data[offset + (2*i+1)*m_slab]) * INV_SQRT_TWO;	
	}
	for (size_t i=0; i<size; i++)
		m_data[offset+i*m_slab] = m_temp[i];
}

void Wavelet3D::backwardStepX(size_t y, size_t z, size_t size) {
	size_t halfSize = size/2;
	size_t offset = y*m_size + z*m_slab;
	for (size_t i=0; i<halfSize; i++) {
		m_temp[2*i]   = (m_data[offset + i] 
			+ m_data[halfSize + offset + i]) * INV_SQRT_TWO;
		m_temp[2*i+1] = (m_data[offset + i] 
			- m_data[halfSize + offset + i]) * INV_SQRT_TWO;
	}
	for (size_t i=0; i<size; i++)
		m_data[offset+i] = m_temp[i];
}

void Wavelet3D::backwardStepY(size_t x, size_t z, size_t size) {
	size_t halfSize = size/2;
	size_t offset = x + z*m_slab;
	for (size_t i=0; i<halfSize; i++) {
		m_temp[2*i]   = (m_data[offset + i*m_size] 
			+ m_data[offset + (halfSize + i)*m_size]) * INV_SQRT_TWO;
		m_temp[2*i+1] = (m_data[offset + i*m_size] 
			- m_data[offset + (halfSize + i)*m_size]) * INV_SQRT_TWO;
	}
	for (size_t i=0; i<size; i++)
		m_data[offset+i*m_size] = m_temp[i];
}

void Wavelet3D::backwardStepZ(size_t x, size_t y, size_t size) {
	size_t halfSize = size/2;
	size_t offset = x + y*m_size;
	for (size_t i=0; i<halfSize; i++) {
		m_temp[2*i]   = (m_data[offset + i*m_slab] 
			+ m_data[offset + (halfSize + i)*m_slab]) * INV_SQRT_TWO;
		m_temp[2*i+1] = (m_data[offset + i*m_slab] 
			- m_data[offset + (halfSize + i)*m_slab]) * INV_SQRT_TWO;
	}
	for (size_t i=0; i<size; i++)
		m_data[offset+i*m_slab] = m_temp[i];
}

Wavelet3D::~Wavelet3D() {
	delete[] m_data;
	delete[] m_temp;
}

/* ==================================================================== */
/*                Sparse 3D Wavelet - Octree representation             */
/* ==================================================================== */

SparseWaveletOctree::SparseWaveletOctree(size_t size, float scalingFunction) : m_size(size) {
	m_maxLevel = log2i(m_size)-1;
	m_nodes.push_back(Node(0));
	m_nodes.push_back(Node(scalingFunction));
}
	
void SparseWaveletOctree::put(int level, int i, int j, int k, float coeff[7]) {
	const float quadrantSigns[7][8] = {
		{1, 1, 1, 1, -1, -1, -1, -1}, /*   X differencing */
		{1, 1, -1, -1, 1, 1, -1, -1}, /*   Y differencing */
		{1, 1, -1, -1, -1, -1, 1, 1}, /*  XY differencing */
		{1, -1, 1, -1, 1, -1, 1, -1}, /*   Z differencing */
		{1, -1, 1, -1, -1, 1, -1, 1}, /*  XZ differencing */
		{1, -1, -1, 1, 1, -1, -1, 1}, /*  YZ differencing */
		{1, -1, -1, 1, -1, 1, 1, -1}  /* XYZ differencing */
	};

	int currentLevel = 0;
	int node = 1;

	while (currentLevel != level) {
		int shift = level - currentLevel - 1;
		int index = (k >> shift) + ((j >> shift) << 1) + ((i >> shift) << 2);

		if (m_nodes[node].child[index] == 0) {
			m_nodes[node].child[index] = m_nodes.size();
			node = m_nodes.size();
			m_nodes.push_back(Node(0));
		} else {
			node = m_nodes[node].child[index];
			if (node < 0)
				Log(EError, "Internal error while building wavelet octree -"
				" expected coefficients in order of decreasing support");
		}

		i -= ((i >> shift) << shift);
		j -= ((j >> shift) << shift);
		k -= ((k >> shift) << shift);
		++currentLevel;
	}

	for (int i=0; i<8; ++i) {
		float value = 0;
		for (int type=0; type<7; ++type)
			value += quadrantSigns[type][i]*coeff[type];

		if (m_nodes[node].child[i] == 0) {
			/* Sacrifice one bit of accuracy */
			union {
				float floatValue;
				int intValue;
			} a;
			a.floatValue = value;
			m_nodes[node].child[i] = (a.intValue >> 1) | 0x80000000;
		} else {
			if (m_nodes[node].child[i] < 0)
				Log(EError, "Internal error while building wavelet octree -"
				" expected coefficients in order of decreasing support");
			m_nodes[m_nodes[node].child[i]].value += value;
		}
	}
}

/*
  Octree traversal algorithm presented in
  "An Efficient Parametric Algorithm for Octree Traversal"
  by J. Revelles, C. Urena, M. Lastra
*/

/* Tables 1 and 2 */
inline int first_node(Float tx0, Float ty0, Float tz0, Float txm, Float tym, Float tzm) {
	uint8_t result = 0;

	if (tx0 > ty0 && tx0 > tz0) {
		if (tym < tx0) result |= 2;
		if (tzm < tx0) result |= 1;
	} else if (ty0 > tz0) {
		if (txm < ty0) result |= 4;
		if (tzm < ty0) result |= 1;
	} else {
		if (txm < tz0) result |= 4;
		if (tym < tz0) result |= 2;
	}

	return result;
}

/* Table 3 */
inline int new_node(Float t1, int a, Float t2, int b, Float t3, int c) {
	return ((t1 < t2 && t1 < t3) ? a : (t2 < t3 ? b : c));
}

Float SparseWaveletOctree::lineIntegral(Point start, Point end) const {
	start /= (Float) m_size;
	end /= (Float) m_size;

	Ray ray(start, normalize(end-start), 0.0f);
	
	uint8_t a = 0;
	if (ray.d.x < 0) {
		ray.d.x = -ray.d.x;
		ray.dRcp.x = -ray.dRcp.x;
		ray.o.x = 1-ray.o.x;
		a |= 4;
	}
	if (ray.d.y < 0) {
		ray.d.y = -ray.d.y;
		ray.dRcp.y = -ray.dRcp.y;
		ray.o.y = 1-ray.o.y;
		a |= 2;
	}
	if (ray.d.z < 0) {
		ray.d.z = -ray.d.z;
		ray.dRcp.z = -ray.dRcp.z;
		ray.o.z = 1-ray.o.z;
		a |= 1;
	}

	Float tx0 = -ray.o.x*ray.dRcp.x,
		  ty0 = -ray.o.y*ray.dRcp.y,
		  tz0 = -ray.o.z*ray.dRcp.z,
		  tx1 = (1-ray.o.x)*ray.dRcp.x,
		  ty1 = (1-ray.o.y)*ray.dRcp.y,
		  tz1 = (1-ray.o.z)*ray.dRcp.z,
		  mint = std::max(std::max(tx0, ty0), tz0),
		  maxt = std::min(std::min(tx1, ty1), tz1);

	if (mint >= maxt)
		return 0.0f;

	return lineIntegral(1, tx0, ty0, tz0, tx1, ty1, tz1, a);
}

Float SparseWaveletOctree::lineIntegral(int32_t idx,
	Float tx0, Float ty0, Float tz0,
	Float tx1, Float ty1, Float tz1, uint8_t a) const {
	if (idx == 0)
		return 0;

	Float overlap = std::min(std::min(tx1, ty1), tz1) 
		- std::max(std::max(tx0, ty0), tz0);

	if (idx < 0) {
		/* Leaf node */
		union {
			float floatValue;
			int intValue;
		} a;
		a.intValue = idx << 1;
		return a.floatValue * overlap;
	}

	const Node &node = m_nodes[idx];
	Float result = node.value * overlap;
	const Float txm = .5f * (tx0+tx1);
	const Float tym = .5f * (ty0+ty1);
	const Float tzm = .5f * (tz0+tz1);
	unsigned int currNode = first_node(tx0, ty0, tz0, txm, tym, tzm);

	do {
		switch (currNode) {
			case 0:
				result += lineIntegral(node.child[a], tx0, ty0, tz0, txm, tym, tzm, a);
				currNode = new_node(txm, 4, tym, 2, tzm, 1);
				break;

			case 1:
				result += lineIntegral(node.child[1^a], tx0, ty0, tzm, txm, tym, tz1, a);
				currNode = new_node(txm, 5, tym, 3, tz1, 8);
				break;

			case 2:
				result += lineIntegral(node.child[2^a], tx0, tym, tz0, txm, ty1, tzm, a);
				currNode = new_node(txm, 6, ty1, 8, tzm, 3);
				break;

			case 3:
				result += lineIntegral(node.child[3^a], tx0, tym, tzm, txm, ty1, tz1, a);
				currNode = new_node(txm, 7, ty1, 8, tz1, 8);
				break;

			case 4:
				result += lineIntegral(node.child[4^a], txm, ty0, tz0, tx1, tym, tzm, a);
				currNode = new_node(tx1, 8, tym, 6, tzm, 5);
				break;

			case 5:
				result += lineIntegral(node.child[5^a], txm, ty0, tzm, tx1, tym, tz1, a);
				currNode = new_node(tx1, 8, tym, 7, tz1, 8);
				break;

			case 6:
				result += lineIntegral(node.child[6^a], txm, tym, tz0, tx1, ty1, tzm, a);
				currNode = new_node(tx1, 8, ty1, 8, tzm, 7);
				break;

			case 7:
				result += lineIntegral(node.child[7^a], txm, tym, tzm, tx1, ty1, tz1, a);
				currNode = 8;
				break;
		}
	} while (currNode < 8);

	return result;
}

std::string SparseWaveletOctree::toString() const {
	std::ostringstream oss;
	oss << "SparseWaveletTree[size=" 
		<< m_nodes.size() * sizeof(Node) / 1024 << " KiB]" << endl;
	return oss.str();
}

MTS_IMPLEMENT_CLASS(Wavelet2D, false, Object)
MTS_IMPLEMENT_CLASS(Wavelet3D, false, Object)
MTS_IMPLEMENT_CLASS_S(SparseWavelet2D, false, SerializableObject)
MTS_IMPLEMENT_CLASS(SparseWaveletOctree, false, Object)
MTS_NAMESPACE_END

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

#include <mitsuba/core/statistics.h>
#include <mitsuba/render/mipmap.h>

MTS_NAMESPACE_BEGIN

static StatsCounter mipmapLookups("Texture", "Mip-map texture lookups");
static StatsCounter ewaLookups("Texture", "EWA texture lookups");

/* Isotropic/anisotropic EWA mip-map texture map class based on PBRT */
MIPMap::MIPMap(int width, int height, Spectrum *pixels, 
	EFilterType filterType, EWrapMode wrapMode, Float maxAnisotropy) 
		: m_width(width), m_height(height), m_filterType(filterType), 
		  m_wrapMode(wrapMode), m_maxAnisotropy(maxAnisotropy) {
	Spectrum *texture = pixels;

	if (filterType != ENone && (!isPowerOfTwo(width) || !isPowerOfTwo(height))) {
		m_width = (int) roundToPowerOfTwo((uint32_t) width);
		m_height = (int) roundToPowerOfTwo((uint32_t) height);

		/* The texture needs to be up-sampled */
		Spectrum *texture1 = new Spectrum[m_width*height];

		/* Re-sample into the X direction */
		ResampleWeight *weights = resampleWeights(width, m_width);
		for (int y=0; y<height; y++) {
			for (int x=0; x<m_width; x++) {
				texture1[x+m_width*y] = Spectrum(0.0f);
				for (int j=0; j<4; j++) {
					int pos = weights[x].firstTexel + j;
					if (pos < 0 || pos >= height) {
						if (wrapMode == ERepeat) 
							pos = modulo(pos, width);
						else if (wrapMode == EClamp)
							pos = clamp(pos, 0, width-1);
					}
					if (pos >= 0 && pos < width)
						texture1[x+m_width*y] += pixels[pos+y*width] 
							* weights[x].weight[j];
				}
			}
		}
		delete[] weights;
		delete[] pixels;

		/* Re-sample into the Y direction */
		texture = new Spectrum[m_width*m_height];
		weights = resampleWeights(height, m_height);
		memset(texture, 0, sizeof(Spectrum)*m_width*m_height);
		for (int x=0; x<m_width; x++) {
			for (int y=0; y<m_height; y++) {
				for (int j=0; j<4; j++) {
					int pos = weights[y].firstTexel + j;
					if (pos < 0 || pos >= height) {
						if (wrapMode == ERepeat) 
							pos = modulo(pos, height);
						else if (wrapMode == EClamp)
							pos = clamp(pos, 0, height-1);
					}
					if (pos >= 0 && pos < height)
						texture[x+m_width*y] += texture1[x+pos*m_width]
							* weights[y].weight[j];
				}
			}
		}
		for (int y=0; y<m_height; y++)
			for (int x=0; x<m_width; x++)
			texture[x+m_width*y].clampNegative();
		delete[] weights;
		delete[] texture1;
	}

	if (m_filterType != ENone)
		m_levels = 1 + log2i((uint32_t) std::max(width, height));
	else
		m_levels = 1;

	m_pyramid = new Spectrum*[m_levels];
	m_pyramid[0] = texture;
	m_levelWidth = new int[m_levels];
	m_levelHeight= new int[m_levels];
	m_levelWidth[0] = m_width;
	m_levelHeight[0] = m_height;
	m_maximum = Spectrum(-std::numeric_limits<Float>::infinity());
	m_minimum = Spectrum(std::numeric_limits<Float>::infinity());
	m_average = Spectrum(0.0f);
	
	if (m_levels > 1) {
		/* Generate the mip-map hierarchy */
		for (int i=1; i<m_levels; i++) {
			m_levelWidth[i]  = std::max(1, m_levelWidth[i-1]/2);
			m_levelHeight[i] = std::max(1, m_levelHeight[i-1]/2);
			m_pyramid[i] = new Spectrum[m_levelWidth[i] * m_levelHeight[i]];
	
			if (i == 1) {
				for (int y = 0; y < m_levelHeight[i]; y++) {
					for (int x = 0; x < m_levelWidth[i]; x++) {
						Spectrum t00 = getTexel(i-1, 2*x, 2*y),
								 t10 = getTexel(i-1, 2*x+1, 2*y),
								 t01 = getTexel(i-1, 2*x, 2*y+1),
								 t11 = getTexel(i-1, 2*x+1, 2*y+1);
	
						/* Compute minima and maxima while processing level 0 */
						for (int k=0; k<SPECTRUM_SAMPLES; ++k) {
							m_maximum[k] = std::max(m_maximum[k],
								std::max(std::max(t00[k], t10[k]), 
										 std::max(t01[k], t11[k])));
							m_minimum[k] = std::min(m_minimum[k],
								std::min(std::min(t00[k], t10[k]), 
										 std::min(t01[k], t11[k])));
						}
	
						m_pyramid[i][x+y*m_levelWidth[i]] = 
							(t00 + t10 + t01 + t11) * 0.25f;
					}
				}
			} else {
				for (int y = 0; y < m_levelHeight[i]; y++) {
					for (int x = 0; x < m_levelWidth[i]; x++) {
						m_pyramid[i][x+y*m_levelWidth[i]] = (
							getTexel(i-1, 2*x, 2*y) + 
							getTexel(i-1, 2*x+1, 2*y) + 
							getTexel(i-1, 2*x, 2*y+1) + 
							getTexel(i-1, 2*x+1, 2*y+1)) * 0.25f;
					}
				}
			}
		}
		m_average = m_pyramid[m_levels-1][0];
	} else {
		/* Nearest filtering, no hierarchy needed -- 
		   still compute average/min/max values */
		int width = m_levelWidth[0], height = m_levelHeight[0];
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				Spectrum value = m_pyramid[0][x + width * y];

				/* Compute minima and maxima while processing level 0 */
				for (int k=0; k<SPECTRUM_SAMPLES; ++k) {
					m_maximum[k] = std::max(m_maximum[k], value[k]);
					m_minimum[k] = std::min(m_minimum[k], value[k]);
				}

				m_average += value;
			}
		}

		m_average /= width*height;
	}
			  
	if (m_wrapMode == EBlack || m_wrapMode == EWhite) {
		Float value = (m_wrapMode == EBlack) ? 0.0f : 1.0f;
		for (int k=0; k<3; ++k) {
			m_maximum[k] = std::max(m_maximum[k], value);
			m_minimum[k] = std::min(m_minimum[k], value);
		}
	}

	if (m_filterType == EEWA) {
		m_weightLut = static_cast<Float *>(allocAligned(sizeof(Float)*MIPMAP_LUTSIZE));
		for (int i=0; i<MIPMAP_LUTSIZE; ++i) {
			Float pos = (Float) i / (Float) (MIPMAP_LUTSIZE-1);
			m_weightLut[i] = std::fastexp(-2.0f * pos) - std::fastexp(-2.0f);
		}
	}
}

MIPMap::~MIPMap() {
	if (m_filterType == EEWA) 
		freeAligned(m_weightLut);
	for (int i=0; i<m_levels; i++)
		delete[] m_pyramid[i];
	delete[] m_levelHeight;
	delete[] m_levelWidth;
	delete[] m_pyramid;
}

ref<MIPMap> MIPMap::fromBitmap(Bitmap *bitmap, EFilterType filterType,
		EWrapMode wrapMode, Float maxAnisotropy,
		Spectrum::EConversionIntent intent) {
	int width = bitmap->getWidth();
	int height = bitmap->getHeight();
	float *data = bitmap->getFloatData();
	Spectrum s, *pixels = new Spectrum[width*height];

	for (int y=0; y<height; y++) {
		for (int x=0; x<width; x++) {
			float r = data[(y*width+x)*4+0];
			float g = data[(y*width+x)*4+1];
			float b = data[(y*width+x)*4+2];
			/* Convert to a spectral representation */
			s.fromLinearRGB(r, g, b, intent);
			s.clampNegative();
			pixels[y*width+x] = s;
		}
	}

	return new MIPMap(width, height, pixels,
		filterType, wrapMode, maxAnisotropy);
}

MIPMap::ResampleWeight *MIPMap::resampleWeights(int oldRes, int newRes) const {
	/* Resample using a Lanczos windowed sinc reconstruction filter */
	Assert(newRes >= oldRes);
	Float filterWidth = 2.0f;
	ResampleWeight *weights = new ResampleWeight[newRes];
	for (int i=0; i<newRes; i++) {
		Float center = (i + .5f) * oldRes / newRes;
		weights[i].firstTexel = floorToInt(center - filterWidth + (Float) 0.5f);
		Float weightSum = 0;
		for (int j=0; j<4; j++) {
			Float pos = weights[i].firstTexel + j + .5f;
			Float weight = lanczosSinc((pos - center) / filterWidth);
			weightSum += weight;
			weights[i].weight[j] = weight;
		}

		/* Normalize */
		Float invWeights = 1.0f / weightSum;
		for (int j=0; j<4; j++)
			weights[i].weight[j] *= invWeights;
	}
	return weights;
}
	
Spectrum MIPMap::getTexel(int level, int x, int y) const {
	int levelWidth = m_levelWidth[level];
	int levelHeight = m_levelHeight[level];

	if (x < 0 || y < 0 || x >= levelWidth || y >= levelHeight) {
		switch (m_wrapMode) {
			case ERepeat:
				x = modulo(x, levelWidth);
				y = modulo(y, levelHeight);
				break;
			case EClamp:
				x = clamp(x, 0, levelWidth - 1);
				y = clamp(y, 0, levelHeight - 1);
				break;
			case EBlack:
				return Spectrum(0.0f);
			case EWhite:
				return Spectrum(1.0f);
		}
	}

	return m_pyramid[level][x + levelWidth*y];
}
	
Spectrum MIPMap::triangle(int level, Float x, Float y) const {
	if (m_filterType == ENone) {
		int xPos = floorToInt(x*m_levelWidth[0]),
			yPos = floorToInt(y*m_levelHeight[0]);
		return getTexel(0, xPos, yPos);
	} else {
		level = clamp(level, 0, m_levels - 1);
		x = x * m_levelWidth[level] - 0.5f;
		y = y * m_levelHeight[level] - 0.5f;
		int xPos = floorToInt(x), yPos = floorToInt(y);
		Float dx = x - xPos, dy = y - yPos;
		return getTexel(level, xPos, yPos) * (1.0f - dx) * (1.0f - dy)
			+ getTexel(level, xPos, yPos + 1) * (1.0f - dx) * dy
			+ getTexel(level, xPos + 1, yPos) * dx * (1.0f - dy)
			+ getTexel(level, xPos + 1, yPos + 1) * dx * dy;
	}	
}
		
Spectrum MIPMap::getValue(Float u, Float v, 
		Float dudx, Float dudy, Float dvdx, Float dvdy) const {
	if (m_filterType == ETrilinear) {
		++mipmapLookups;
		/* Conservatively estimate a square lookup region */
		Float width = 2.0f * std::max(
			std::max(std::abs(dudx), std::abs(dudy)),
			std::max(std::abs(dvdx), std::abs(dvdy)));
		Float mipmapLevel = m_levels - 1 + 
			log2(std::max(width, (Float) 1e-8f));

		if (mipmapLevel < 0) {
			/* The lookup is smaller than one pixel */
			return triangle(0, u, v);
		} else if (mipmapLevel >= m_levels - 1) {
			/* The lookup is larger than the whole texture */
			return getTexel(m_levels - 1, 0, 0);
		} else {
			/* Tri-linear interpolation */
			int level = (int) mipmapLevel;
			Float delta = mipmapLevel - level;
			return triangle(level, u, v) * (1.0f - delta)
				+ triangle(level, u, v) * delta;
		}
	} else if (m_filterType == EEWA) {
		if (dudx*dudx + dudy*dudy < dvdx*dvdx + dvdy*dvdy) {
			std::swap(dudx, dvdx);
			std::swap(dudy, dvdy);
		}

		Float majorLength = std::sqrt(dudx * dudx + dudy * dudy);
		Float minorLength = std::sqrt(dvdx * dvdx + dvdy * dvdy);

		if (minorLength * m_maxAnisotropy < majorLength && minorLength > 0.0f) {
			Float scale = majorLength / (minorLength * m_maxAnisotropy);
			dvdx *= scale; dvdy *= scale;
			minorLength *= scale;
		}

		if (minorLength == 0)
			return triangle(0, u, v);
	
		// The min() below avoids overflow in the int conversion when lod=inf
		Float lod = 
			std::min(std::max((Float) 0, m_levels - 1 + log2(minorLength)), 
				(Float) (m_levels-1));
		int ilod = floorToInt(lod);
		Float d = lod - ilod;

		return EWA(u, v, dudx, dudy, dvdx, dvdy, ilod)   * (1-d) +
			   EWA(u, v, dudx, dudy, dvdx, dvdy, ilod+1) * d;
	} else {
		int xPos = floorToInt(u*m_levelWidth[0]),
			yPos = floorToInt(v*m_levelHeight[0]);
		return getTexel(0, xPos, yPos);
	}
}
	
Spectrum MIPMap::EWA(Float u, Float v, Float dudx, Float dudy, Float dvdx, 
	Float dvdy, int level) const {
	++ewaLookups;
	if (level >= m_levels)
		return getTexel(m_levels-1, 0, 0);

	Spectrum result(0.0f);
	Float denominator = 0.0f;
	u = u * m_levelWidth[level]; v = v * m_levelHeight[level];
	dudx = dudx * m_levelWidth[level]; dudy = dudy * m_levelHeight[level];
	dvdx = dvdx * m_levelWidth[level]; dvdy = dvdy * m_levelHeight[level];

	Float A = dudy * dudy + dvdy * dvdy + 1.0f;
	Float B = -2.0f * (dudx * dudy + dvdx * dvdy);
	Float C = dudx * dudx + dvdx * dvdx + 1.0f;
	Float F = A * C - B * B * 0.25f;
	Float du = std::sqrt(C), dv = std::sqrt(A);
	int u0 = (int) std::ceil(u - du);
	int u1 = (int) std::floor(u + du);
	int v0 = (int) std::ceil(v - dv);
	int v1 = (int) std::floor(v + dv);
	Float invF = 1.0f / F;
	A *= invF; B *= invF; C *= invF;

	for (int ut = u0; ut <= u1; ++ut) {
		const Float uu = ut - u;
		for (int vt = v0; vt <= v1; ++vt) {
			const Float vv = vt - v;
			const Float r2 = A*uu*uu + B*uu*vv + C*vv*vv;
			if (r2 < 1) {
				const Float weight = m_weightLut[
					std::max(0, std::min((int) (r2 * MIPMAP_LUTSIZE), MIPMAP_LUTSIZE - 1))];
				result += getTexel(level, ut, vt) * weight;
				denominator += weight;
			}
		}
	}
	return result / denominator;
}

Bitmap *MIPMap::getBitmap() const {
	Bitmap *bitmap = new Bitmap(m_width, m_height, 128);
	float *floatData = bitmap->getFloatData();
	Spectrum *specData = m_pyramid[0];

	for (int y=0; y<m_height; ++y) {
		for (int x=0; x<m_width; ++x) {
			Float r, g, b;
			(specData++)->toLinearRGB(r, g, b);
			*floatData++ = r;
			*floatData++ = g;
			*floatData++ = b;
			*floatData++ = 1.0f;
		}
	}
	return bitmap;
}

Bitmap *MIPMap::getLDRBitmap() const {
	Bitmap *bitmap = new Bitmap(m_width, m_height, 24);
	uint8_t *data = bitmap->getData();
	Spectrum *specData = m_pyramid[0];

	for (int y=0; y<m_height; ++y) {
		for (int x=0; x<m_width; ++x) {
			Float r, g, b;
			(specData++)->toLinearRGB(r, g, b);
			*data++ = (uint8_t) std::min(255, std::max(0, (int) (r*255)));
			*data++ = (uint8_t) std::min(255, std::max(0, (int) (g*255)));
			*data++ = (uint8_t) std::min(255, std::max(0, (int) (b*255)));
		}
	}
	return bitmap;
}

MTS_IMPLEMENT_CLASS(MIPMap, false, Object)
MTS_NAMESPACE_END

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

#include <mitsuba/core/statistics.h>
#include <mitsuba/render/mipmap.h>

MTS_NAMESPACE_BEGIN

static StatsCounter mipmapLookups("Texture", "Mip-map texture lookups");
static StatsCounter ewaLookups("Texture", "EWA texture lookups");

/* Isotropic/anisotropic EWA mip-map texture map class based on PBRT */
MIPMap::MIPMap(int width, int height, Spectrum *pixels, 
	bool isotropic, EWrapMode wrapMode, Float maxAnisotropy) 
		: m_width(width), m_height(height), m_isotropic(isotropic), 
		  m_wrapMode(wrapMode), m_maxAnisotropy(maxAnisotropy) {
	Spectrum *texture = pixels;

	if (!isPowerOfTwo(width) || !isPowerOfTwo(height)) {
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

	m_levels = 1 + log2i((uint32_t) std::max(width, height));
	m_pyramid = new Spectrum*[m_levels];
	m_pyramid[0] = texture;
	m_levelWidth = new int[m_levels];
	m_levelHeight= new int[m_levels];
	m_levelWidth[0] = m_width;
	m_levelHeight[0] = m_height;

	/* Generate the mip-map hierarchy */
	for (int i=1; i<m_levels; i++) {
		m_levelWidth[i]  = std::max(1, m_levelWidth[i-1]/2);
		m_levelHeight[i] = std::max(1, m_levelHeight[i-1]/2);
		m_pyramid[i] = new Spectrum[m_levelWidth[i] * m_levelHeight[i]];
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

	if (!isotropic) {
		m_weightLut = static_cast<Float *>(allocAligned(sizeof(Float)*MIPMAP_LUTSIZE));
		for (int i=0; i<MIPMAP_LUTSIZE; ++i) {
			Float pos = (Float) i / (Float) (MIPMAP_LUTSIZE-1);
			m_weightLut[i] = std::exp(-2.0f * pos) - std::exp(-2.0f);
		}
	}
}

MIPMap::~MIPMap() {
	if (!m_isotropic)
		freeAligned(m_weightLut);
	for (int i=0; i<m_levels; i++)
		delete[] m_pyramid[i];
	delete[] m_levelHeight;
	delete[] m_levelWidth;
	delete[] m_pyramid;
}

Spectrum MIPMap::getMaximum() const {
	Spectrum max(0.0f);
	int height = m_levelHeight[0];
	int width = m_levelWidth[0];
	Spectrum *pixels = m_pyramid[0];
	for (int y=0; y<height; ++y) {
		for (int x=0; x<width; ++x) {
			Spectrum value = *pixels++;
			for (int j=0; j<SPECTRUM_SAMPLES; ++j)
				max[j] = std::max(max[j], value[j]);
		}
	}
	if (m_wrapMode == EWhite) {
		for (int i=0; i<SPECTRUM_SAMPLES; ++i)
			max[i] = std::max(max[i], (Float) 1.0f);
	}
	return max;
}
	
ref<MIPMap> MIPMap::fromBitmap(Bitmap *bitmap, bool isotropic, 
		EWrapMode wrapMode, Float maxAnisotropy) {
	int width = bitmap->getWidth();
	int height = bitmap->getHeight();
	float *data = bitmap->getFloatData();
	Spectrum s, *pixels = new Spectrum[width*height];

	for (int y=0; y<height; y++) {
		for (int x=0; x<width; x++) {
			float r = data[(y*width+x)*4+0];
			float g = data[(y*width+x)*4+1];
			float b = data[(y*width+x)*4+2];
			s.fromLinearRGB(r, g, b);
			s.clampNegative();
			/* Convert to a spectral representation */
			pixels[y*width+x] = s;
		}
	}

	return new MIPMap(width, height, pixels,
		isotropic, wrapMode, maxAnisotropy);
}

MIPMap::ResampleWeight *MIPMap::resampleWeights(int oldRes, int newRes) const {
	/* Resample using a Lanczos windowed sinc reconstruction filter */
	Assert(newRes >= oldRes);
	Float filterWidth = 2.0f;
	ResampleWeight *weights = new ResampleWeight[newRes];
	for (int i=0; i<newRes; i++) {
		Float center = (i + .5f) * oldRes / newRes;
		weights[i].firstTexel = (int) std::floor(center - filterWidth + (Float) 0.5f);
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

	if (x <= 0 || y < 0 || x >= levelWidth || y >= levelHeight) {
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
	level = clamp(level, 0, m_levels - 1);
	x = x * m_levelWidth[level] - 0.5f;
	y = y * m_levelHeight[level] - 0.5f;
	int xPos = (int) std::floor(x), yPos = (int) std::floor(y);
	Float dx = x - xPos, dy = y - yPos;
	return getTexel(level, xPos, yPos) * (1.0f - dx) * (1.0f - dy)
		+ getTexel(level, xPos, yPos + 1) * (1.0f - dx) * dy
		+ getTexel(level, xPos + 1, yPos) * dx * (1.0f - dy)
		+ getTexel(level, xPos + 1, yPos + 1) * dx * dy;
}
		
Spectrum MIPMap::getValue(Float u, Float v, 
		Float dudx, Float dudy, Float dvdx, Float dvdy) const {
	if (m_isotropic) {
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
	} else {
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
		int ilod = (int) std::floor(lod);
		Float d = lod - ilod;

		return EWA(u, v, dudx, dudy, dvdx, dvdy, ilod)   * (1-d) +
			   EWA(u, v, dudx, dudy, dvdx, dvdy, ilod+1) * d;
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

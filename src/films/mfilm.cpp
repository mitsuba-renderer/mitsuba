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

#include <mitsuba/render/film.h>
#include <mitsuba/core/bitmap.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

/*!\plugin{mfilm}{MATLAB M-file film}
 * \parameters{
 *     \parameter{width, height}{\Integer}{
 *       Width and height of the camera sensor in pixels
 *       \default{768, 576}
 *     }
 *     \parameter{cropOffsetX, cropOffsetY, cropWidth, cropHeight}{\Integer}{
 *       These parameter can optionally be provided to render a sub-rectangle
 *       of the output \default{Unused}
 *     }
 *     \parameter{alpha}{\Boolean}{Include an alpha channel in the output
 *        image? \default{\code{true}}
 *     }
 *     \parameter{banner}{\Boolean}{Include a small Mitsuba banner in the 
 *         output image? \default{\code{true}}
 *     }
 * }
 * 
 * This plugin provides a camera film that exports luminance
 * values as a matrix using the MATLAB M-file format. This is
 * useful when running Mitsuba as simulation step as part of a 
 * larger virtual experiment. It can also come in handy when
 * verifying parts of the renderer using a test suite.
 *
 * When Mitsuba is started with the ``test case mode'' parameter
 * (\code{-t}), this class will write triples consisting of
 * the luminance, variance, and sample count for every pixel
 * (instead of just the luminance).
 */
class MFilm : public Film {
protected:
	/* Pixel data structure */
	struct Pixel {
		Spectrum spec;
		Spectrum variance;
		Float weight;
		int nSamples;

		Pixel() : spec(0.0f), variance(0.0f), weight(0.0f) {
		}
	};

	Pixel *m_pixels;
	bool m_hasVariances;
	bool m_exportSpectra;
public:
	MFilm(const Properties &props) : Film(props) {
		m_pixels = new Pixel[m_cropSize.x * m_cropSize.y];
		/* Export luminance by default */
		m_exportSpectra = props.getBoolean("spectra", false);
	}

	MFilm(Stream *stream, InstanceManager *manager) 
		: Film(stream, manager) {
		m_pixels = new Pixel[m_cropSize.x * m_cropSize.y];
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Film::serialize(stream, manager);
	}

	virtual ~MFilm() {
		if (m_pixels)
			delete[] m_pixels;
	}
	
	void clear() {
		memset(m_pixels, 0, sizeof(Pixel) * m_cropSize.x * m_cropSize.y);
		m_hasVariances = false;
	}

	void fromBitmap(const Bitmap *bitmap) {
		Assert(bitmap->getWidth() == m_cropSize.x 
			&& bitmap->getHeight() == m_cropSize.y);
		Assert(bitmap->getBitsPerPixel() == 128);
		unsigned int lastIndex = m_cropSize.x * m_cropSize.y;

		for (unsigned int index=0; index<lastIndex; ++index) {
			Pixel &pixel = m_pixels[index];
			const float 
				r = bitmap->getFloatData()[index*4+0],
				g = bitmap->getFloatData()[index*4+1],
				b = bitmap->getFloatData()[index*4+2];
			pixel.spec.fromLinearRGB(r, g, b);
			pixel.weight = 1.0f;
		}
		m_hasVariances = false;
	}

	void toBitmap(Bitmap *bitmap) const {
		Assert(bitmap->getWidth() == m_cropSize.x 
			&& bitmap->getHeight() == m_cropSize.y);
		Assert(bitmap->getBitsPerPixel() == 128);
		unsigned int lastIndex = m_cropSize.x * m_cropSize.y;
		Float r, g, b;

		for (unsigned int index=0; index<lastIndex; ++index) {
			Pixel &pixel = m_pixels[index];
			Float invWeight = pixel.weight > 0 ? 1/pixel.weight : 0;
			pixel.spec.toLinearRGB(r, g, b);
			bitmap->getFloatData()[index*4+0] = r*invWeight;
			bitmap->getFloatData()[index*4+1] = g*invWeight;
			bitmap->getFloatData()[index*4+2] = b*invWeight;
			bitmap->getFloatData()[index*4+3] = 1.0f;
		}
	}

	Spectrum getValue(int xPixel, int yPixel) {
		xPixel -= m_cropOffset.x; yPixel -= m_cropOffset.y;
		if (!(xPixel >= 0 && xPixel < m_cropSize.x && yPixel >= 0 && yPixel < m_cropSize.y )) {
			Log(EWarn, "Pixel out of range : %i,%i", xPixel, yPixel); 
			return Spectrum(0.0f);
		}
		Pixel &pixel = m_pixels[xPixel + yPixel * m_cropSize.x];
		return pixel.spec / pixel.weight;
	}

	void putImageBlock(const ImageBlock *block) {
		int entry=0, imageY = block->getOffset().y - 
			block->getBorder() - m_cropOffset.y - 1;

		if (!block->collectStatistics()) {
			for (int y=0; y<block->getFullSize().y; ++y) {
				if (++imageY < 0 || imageY >= m_cropSize.y) {
					/// Skip a row if it is outside of the crop region
					entry += block->getFullSize().x;
					continue;
				}

				int imageX = block->getOffset().x - block->getBorder()
					- m_cropOffset.x - 1;
				for (int x=0; x<block->getFullSize().x; ++x) {
					if (++imageX < 0 || imageX >= m_cropSize.x) {
						++entry;
						continue;
					}

					Pixel &pixel = m_pixels[imageY * m_cropSize.x + imageX];

					pixel.spec += block->getPixel(entry);
					pixel.weight += block->getWeight(entry++);
				}
			}
		} else {
			Assert(block->getBorder() == 0);
			m_hasVariances = true;
			for (int y=0; y<block->getFullSize().y; ++y) {
				if (++imageY < 0 || imageY >= m_cropSize.y) {
					/// Skip a row if it is outside of the crop region
					entry += block->getFullSize().x;
					continue;
				}

				int imageX = block->getOffset().x - block->getBorder()
					- m_cropOffset.x - 1;
				for (int x=0; x<block->getFullSize().x; ++x) {
					if (++imageX < 0 || imageX >= m_cropSize.x) {
						++entry;
						continue;
					}

					Pixel &pixel = m_pixels[imageY * m_cropSize.x + imageX];

					pixel.spec += block->getPixel(entry);
					pixel.weight = block->getWeight(entry);
					pixel.nSamples = block->getSampleCount(entry);
					pixel.variance = block->getVariance(entry++);
				}
			}
		}
	}

	void develop(const fs::path &destFile) {
		fs::path filename = destFile;
		std::string extension = boost::to_lower_copy(fs::extension(filename));
		if (extension != ".m")
			filename.replace_extension(".m");

		Log(EInfo, "Writing image to \"%s\" ..", filename.leaf().c_str());
	
		FILE *f = fopen(filename.file_string().c_str(), "w");
		if (!f)
			Log(EError, "Output file cannot be created!");

		fprintf(f, "[");
		int pos = 0;
		for (int y=0; y<m_cropSize.y; y++) {
			for (int x=0; x<m_cropSize.x; x++) {
				Pixel &pixel = m_pixels[pos];
				if (m_exportSpectra) {
					Float invWeight = pixel.weight > 0 ? 1/pixel.weight : 1;
					Spectrum spec = invWeight*pixel.spec;
					for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
						if (!m_hasVariances)
							fprintf(f, "%f", spec[i]);
						else
							fprintf(f, "%f %f %i", spec[i], pixel.variance[i], pixel.nSamples);
						if (i+1 < SPECTRUM_SAMPLES)
							fprintf(f, " ");
					}
					if (x + 1 < m_cropSize.x)
						fprintf(f, ", ");
				} else {
					Float invWeight = pixel.weight > 0 ? 1/pixel.weight : 1;
					Float luminance = invWeight*pixel.spec.getLuminance();
					if (!m_hasVariances)
						fprintf(f, "%f", luminance);
					else
						fprintf(f, "%f %f %i", luminance, pixel.variance.getLuminance(), pixel.nSamples);
					if (x + 1 < m_cropSize.x)
						fprintf(f, ", ");
				}
				++pos;
			}
			if (y + 1 < m_cropSize.y)
				fprintf(f, ";\n ");
		}

		fprintf(f, "]\n");
		fclose(f);
	}
	
	bool destinationExists(const fs::path &baseName) const {
		fs::path filename = baseName;
		if (boost::to_lower_copy(filename.extension()) != ".m")
			filename.replace_extension(".m");
		return fs::exists(filename);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MFilm[" << std::endl
			<< "  size = " << m_size.toString() << "," << std::endl
			<< "  cropOffset = " << m_cropOffset.toString() << "," << std::endl
			<< "  cropSize = " << m_cropSize.toString() << std::endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(MFilm, false, Film)
MTS_EXPORT_PLUGIN(MFilm, "MATLAB film");
MTS_NAMESPACE_END

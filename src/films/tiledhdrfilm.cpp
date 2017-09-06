/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

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
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/version.h>
#include <boost/algorithm/string.hpp>

#if defined(_MSC_VER)
#pragma warning(disable : 4231) // nonstandard extension used : 'extern' before template explicit instantiation
#endif

#include <ImfTiledOutputFile.h>
#include <ImfChannelList.h>
#include <ImfStringAttribute.h>
#include <ImfFrameBuffer.h>
#include <ImfStandardAttributes.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{tiledhdrfilm}{Tiled high dynamic range film}
 * \order{2}
 * \parameters{
 *     \parameter{width, height}{\Integer}{
 *       Width and height of the camera sensor in pixels
 *       \default{768, 576}
 *     }
 *     \parameter{cropOffsetX, cropOffsetY, cropWidth, cropHeight}{\Integer}{
 *       These parameters can optionally be provided to select a sub-rectangle
 *       of the output. In this case, Mitsuba will only render the requested
 *       regions. \default{Unused}
 *     }
 *     \parameter{pixelFormat}{\String}{Specifies the desired pixel format
 *         for OpenEXR output images. The options are \code{luminance},
 *         \code{luminanceAlpha}, \code{rgb}, \code{rgba}, \code{xyz},
 *         \code{xyza}, \code{spectrum}, and \code{spectrumAlpha}. In the latter two cases,
 *         the number of written channels depends on the value assigned to
 *         \code{SPECTRUM\_SAMPLES} during compilation (see Section~\ref{sec:compiling}
 *         for details)
 *         \default{\code{rgb}}
 *     }
 *     \parameter{componentFormat}{\String}{Specifies the desired floating
 *         point component format used for the output. The options are
 *         \code{float16}, \code{float32}, or \code{uint32}
 *         \default{\code{float16}}
 *     }
 *
 *     \parameter{\Unnamed}{\RFilter}{Reconstruction filter that should
 *     be used by the film. \default{\code{gaussian}, a windowed Gaussian filter}}
 * }
 *
 * This plugin implements a camera film that stores the captured image
 * as a \emph{tiled} high dynamic-range OpenEXR file. It is very similar to
 * \pluginref{hdrfilm}, the main difference being that it does not keep
 * the rendered image in memory. Instead, image tiles are directly written
 * to disk as they are being rendered, which enables renderings of extremely
 * large output images that would otherwise not fit into memory (e.g.
 * 100K$\times$100K).
 *
 * When the image can fit into memory, usage of this plugin is discouraged:
 * due to the extra overhead of tracking image tiles, the rendering process
 * will be slower, and the output files also generally do not compress as
 * well as those produced by \pluginref{hdrfilm}.
 *
 * Based on the provided parameter values, the film will either write a luminance,
 * luminance/alpha, RGB(A), XYZ(A) tristimulus, or spectrum/spectrum-alpha-based
 * bitmap having a \code{float16}, \code{float32}, or \code{uint32}-based
 * internal representation. The default is RGB and \code{float16}.
 * Note that the spectral output options only make sense when using a
 * custom compiled Mitsuba distribution that has spectral rendering
 * enabled. This is not the case for the downloadable release builds.
 *
 * When RGB output is selected, the measured spectral power distributions are
 * converted to linear RGB based on the CIE 1931 XYZ color matching curves and
 * the ITU-R Rec. BT.709 primaries with a D65 white point.
 *
 * \remarks{
 *    \item This film is only meant for command line-based rendering. When
 *    used with \texttt{mtsgui}, the preview image will be black.
 *    \item This plugin is slower than \pluginref{hdrfilm}, and therefore should
 *    only be used when the output image is too large to fit into system memory.
 * }
 */

class TiledHDRFilm : public Film {
public:
    TiledHDRFilm(const Properties &props) : Film(props), m_output(NULL), m_frameBuffer(NULL) {
        std::vector<std::string> pixelFormats = tokenize(boost::to_lower_copy(
            props.getString("pixelFormat", "rgb")), " ,");
        std::vector<std::string> channelNames = tokenize(
            props.getString("channelNames", ""), ", ");
        std::string componentFormat = boost::to_lower_copy(
            props.getString("componentFormat", "float16"));

        if (pixelFormats.empty())
            Log(EError, "At least one pixel format must be specified!");

        if ((pixelFormats.size() != 1 && channelNames.size() != pixelFormats.size()) ||
            (pixelFormats.size() == 1 && channelNames.size() > 1))
            Log(EError, "Number of channel names must match the number of specified pixel formats!");

        for (size_t i=0; i<pixelFormats.size(); ++i) {
            std::string pixelFormat = pixelFormats[i];
            std::string name = i < channelNames.size() ? (channelNames[i] + std::string(".")) : "";

            if (pixelFormat == "luminance") {
                m_pixelFormats.push_back(Bitmap::ELuminance);
                m_channelNames.push_back(name + "Y");
            } else if (pixelFormat == "luminancealpha") {
                m_pixelFormats.push_back(Bitmap::ELuminanceAlpha);
                m_channelNames.push_back(name + "Y");
                m_channelNames.push_back(name + "A");
            } else if (pixelFormat == "rgb") {
                m_pixelFormats.push_back(Bitmap::ERGB);
                m_channelNames.push_back(name + "R");
                m_channelNames.push_back(name + "G");
                m_channelNames.push_back(name + "B");
            } else if (pixelFormat == "rgba") {
                m_pixelFormats.push_back(Bitmap::ERGBA);
                m_channelNames.push_back(name + "R");
                m_channelNames.push_back(name + "G");
                m_channelNames.push_back(name + "B");
                m_channelNames.push_back(name + "A");
            } else if (pixelFormat == "xyz") {
                m_pixelFormats.push_back(Bitmap::EXYZ);
                m_channelNames.push_back(name + "X");
                m_channelNames.push_back(name + "Y");
                m_channelNames.push_back(name + "Z");
            } else if (pixelFormat == "xyza") {
                m_pixelFormats.push_back(Bitmap::EXYZA);
                m_channelNames.push_back(name + "X");
                m_channelNames.push_back(name + "Y");
                m_channelNames.push_back(name + "Z");
                m_channelNames.push_back(name + "A");
            } else if (pixelFormat == "spectrum") {
                m_pixelFormats.push_back(Bitmap::ESpectrum);
                for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
                    std::pair<Float, Float> coverage = Spectrum::getBinCoverage(i);
                    m_channelNames.push_back(name + formatString("%.2f-%.2fnm", coverage.first, coverage.second));
                }
            } else if (pixelFormat == "spectrumalpha") {
                m_pixelFormats.push_back(Bitmap::ESpectrumAlpha);
                for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
                    std::pair<Float, Float> coverage = Spectrum::getBinCoverage(i);
                    m_channelNames.push_back(name + formatString("%.2f-%.2fnm", coverage.first, coverage.second));
                }
                m_channelNames.push_back(name + "A");
            } else {
                Log(EError, "The \"pixelFormat\" parameter must either be equal to "
                    "\"luminance\", \"luminanceAlpha\", \"rgb\", \"rgba\", \"xyz\", \"xyza\", "
                    "\"spectrum\", or \"spectrumAlpha\"!");
            }
        }

        for (size_t i=0; i<m_pixelFormats.size(); ++i) {
            if (SPECTRUM_SAMPLES == 3 && (m_pixelFormats[i] == Bitmap::ESpectrum || m_pixelFormats[i] == Bitmap::ESpectrumAlpha))
                Log(EError, "You requested to render a spectral image, but Mitsuba is currently "
                    "configured for a RGB flow (i.e. SPECTRUM_SAMPLES = 3). You will need to recompile "
                    "it with a different configuration. Please see the documentation for details.");
        }

        if (componentFormat == "float16") {
            m_componentFormat = Bitmap::EFloat16;
        } else if (componentFormat == "float32") {
            m_componentFormat = Bitmap::EFloat32;
        } else if (componentFormat == "uint32") {
            m_componentFormat = Bitmap::EUInt32;
        } else {
            Log(EError, "The \"componentFormat\" parameter must either be "
                "equal to \"float16\", \"float32\", or \"uint32\"!");
        }

        if (m_highQualityEdges)
            Log(EError, "The 'highQualityEdges' parameter is incompatible with the "
                "tiled EXR film. Please disable it.");
    }

    TiledHDRFilm(Stream *stream, InstanceManager *manager)
        : Film(stream, manager), m_output(NULL), m_frameBuffer(NULL) {
        m_pixelFormats.resize((size_t) stream->readUInt());
        for (size_t i=0; i<m_pixelFormats.size(); ++i)
            m_pixelFormats[i] = (Bitmap::EPixelFormat) stream->readUInt();
        m_channelNames.resize((size_t) stream->readUInt());
        for (size_t i=0; i<m_channelNames.size(); ++i)
            m_channelNames[i] = stream->readString();
        m_componentFormat = (Bitmap::EComponentFormat) stream->readUInt();
    }

    virtual ~TiledHDRFilm() {
        develop(NULL, 0);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Film::serialize(stream, manager);
        for (size_t i=0; i<m_pixelFormats.size(); ++i)
            stream->writeUInt(m_pixelFormats[i]);
        stream->writeUInt((uint32_t) m_channelNames.size());
        for (size_t i=0; i<m_channelNames.size(); ++i)
            stream->writeString(m_channelNames[i]);
        stream->writeUInt(m_componentFormat);
    }

    void setDestinationFile(const fs::path &destFile, uint32_t blockSize) {
        if (m_output)
            develop(NULL, 0);

        fs::path filename = destFile;
        std::string extension = boost::to_lower_copy(filename.extension().string());
        if (extension != ".exr")
            filename.replace_extension(".exr");

        Log(EInfo, "Commencing creation of a tiled EXR image at \"%s\" ..", filename.string().c_str());

        Imf::Header header(m_size.x, m_size.y);
        header.setTileDescription(Imf::TileDescription(blockSize, blockSize, Imf::ONE_LEVEL));
        header.insert("generated-by", Imf::StringAttribute("Mitsuba version " MTS_VERSION));

        if (m_pixelFormats.size() == 1) {
            /* Write a chromaticity tag when this is possible */
            Bitmap::EPixelFormat pixelFormat = m_pixelFormats[0];
            if (pixelFormat == Bitmap::EXYZ || pixelFormat == Bitmap::EXYZA) {
                Imf::addChromaticities(header, Imf::Chromaticities(
                    Imath::V2f(1.0f, 0.0f),
                    Imath::V2f(0.0f, 1.0f),
                    Imath::V2f(0.0f, 0.0f),
                    Imath::V2f(1.0f/3.0f, 1.0f/3.0f)));
            } else if (pixelFormat == Bitmap::ERGB || pixelFormat == Bitmap::ERGBA) {
                Imf::addChromaticities(header, Imf::Chromaticities());
            }
        }

        Imf::PixelType compType;
        size_t compStride;

        if (m_componentFormat == Bitmap::EFloat16) {
            compType = Imf::HALF;
            compStride = 2;
        } else if (m_componentFormat == Bitmap::EFloat32) {
            compType = Imf::FLOAT;
            compStride = 4;
        } else if (m_componentFormat == Bitmap::EUInt32) {
            compType = Imf::UINT;
            compStride = 4;
        } else {
            Log(EError, "Invalid component type (must be "
                "float16, float32, or uint32)");
            return;
        }

        Imf::ChannelList &channels = header.channels();
        for (size_t i=0; i<m_channelNames.size(); ++i)
            channels.insert(m_channelNames[i].c_str(), Imf::Channel(compType));

        m_output = new Imf::TiledOutputFile(filename.string().c_str(), header);
        m_frameBuffer = new Imf::FrameBuffer();
        m_blockSize = (int) blockSize;
        m_blocksH = (m_size.x + blockSize - 1) / blockSize;
        m_blocksV = (m_size.y + blockSize - 1) / blockSize;

        m_pixelStride = m_channelNames.size() * compStride;
        m_rowStride = m_pixelStride * m_blockSize;

        if (m_pixelFormats.size() == 1) {
            m_tile = new Bitmap(m_pixelFormats[0], m_componentFormat,
                    Vector2i(m_blockSize, m_blockSize));
        } else {
            m_tile = new Bitmap(Bitmap::EMultiChannel, m_componentFormat,
                    Vector2i(m_blockSize, m_blockSize), (uint8_t) m_channelNames.size());
            m_tile->setChannelNames(m_channelNames);
        }

        char *ptr = (char *) m_tile->getUInt8Data();

        for (size_t i=0; i<m_channelNames.size(); ++i) {
            m_frameBuffer->insert(m_channelNames[i].c_str(),
                Imf::Slice(compType, ptr, m_pixelStride, m_rowStride));
            ptr += compStride;
        }

        m_output->setFrameBuffer(*m_frameBuffer);
        m_peakUsage = 0;
    }

    void put(const ImageBlock *block) {
        Assert(m_output != NULL);

        if ((block->getOffset().x % m_blockSize) != 0 ||
            (block->getOffset().y % m_blockSize) != 0)
            Log(EError, "Encountered an unaligned block!");

        if (block->getSize().x > m_blockSize ||
            block->getSize().y > m_blockSize)
            Log(EError, "Encountered an oversized block!");

        int x = block->getOffset().x / (int) m_blockSize;
        int y = block->getOffset().y / (int) m_blockSize;

        /* Create two copies: a clean one, and one that is used for accumulation */
        ref<ImageBlock> copy1, copy2;
        if (m_freeBlocks.size() > 0) {
            copy1 = m_freeBlocks.back();
            block->copyTo(copy1);
            m_freeBlocks.pop_back();
        } else {
            copy1 = block->clone();
            copy1->incRef();
            ++m_peakUsage;
        }

        if (m_freeBlocks.size() > 0) {
            copy2 = m_freeBlocks.back();
            block->copyTo(copy2);
            m_freeBlocks.pop_back();
        } else {
            copy2 = block->clone();
            copy2->incRef();
            ++m_peakUsage;
        }

        uint32_t idx = (uint32_t) x + (uint32_t) y * m_blocksH;
        m_origBlocks[idx]   = copy1;
        m_mergedBlocks[idx] = copy2;

        for (int yo = -1; yo <= 1; ++yo)
            for (int xo = -1; xo <= 1; ++xo)
                potentiallyWrite(x + xo, y + yo);
    }

    void setBitmap(const Bitmap *bitmap, Float multiplier) {
        Log(EError, "setBitmap(): Global image updates are permitted by this film, "
            "which operates strictly on tiles! Please either switch to a compatible "
            "rendering technique or use a non-tiled film. (e.g. 'hdrfilm')");
    }

    void addBitmap(const Bitmap *bitmap, Float multiplier) {
        Log(EError, "addBitmap(): Global image updates are permitted by this film, "
            "which operates strictly on tiles! Please either switch to a compatible "
            "rendering technique or use a non-tiled film. (e.g. 'hdrfilm')");
    }

    void potentiallyWrite(int x, int y) {
        if (x < 0 || y < 0 || x >= m_blocksH || y >= m_blocksV)
            return;

        uint32_t idx = (uint32_t) x + (uint32_t) y * m_blocksH;
        std::map<uint32_t, ImageBlock *>::iterator it = m_origBlocks.find(idx);
        if (it == m_origBlocks.end())
            return;

        ImageBlock *origBlock = it->second;
        if (origBlock == NULL)
            return;

        /* This could be accelerated using some counters */
        for (int yo = -1; yo <= 1; ++yo) {
            for (int xo = -1; xo <= 1; ++xo) {
                int xp = x + xo, yp = y + yo;
                if (xp < 0 || yp < 0 || xp >= m_blocksH || yp >= m_blocksV
                   || (xp == x && yp == y))
                    continue;

                uint32_t idx2 = (uint32_t) xp + (uint32_t) yp * m_blocksH;
                if (m_origBlocks.find(idx2) == m_origBlocks.end())
                    return; /* Not all neighboring blocks are there yet */
            }
        }

        ImageBlock *mergedBlock = m_mergedBlocks[idx];
        if (mergedBlock == NULL)
            return;

        /* All neighboring blocks are there -- join overlapping regions */
        for (int yo = -1; yo <= 1; ++yo) {
            for (int xo = -1; xo <= 1; ++xo) {
                int xp = x + xo, yp = y + yo;
                if (xp < 0 || yp < 0 || xp >= m_blocksH || yp >= m_blocksV
                   || (xp == x && yp == y))
                    continue;
                uint32_t idx2 = (uint32_t) xp + (uint32_t) yp * m_blocksH;
                ImageBlock *origBlock2   = m_origBlocks[idx2];
                ImageBlock *mergedBlock2 = m_mergedBlocks[idx2];
                if (!origBlock2 || !mergedBlock2)
                    continue;

                mergedBlock->put(origBlock2);
                mergedBlock2->put(origBlock);
            }
        }

        const Bitmap *source = mergedBlock->getBitmap();

        size_t sourceBpp = source->getBytesPerPixel();
        size_t targetBpp = m_tile->getBytesPerPixel();

        const uint8_t *sourceData = source->getUInt8Data()
            + mergedBlock->getBorderSize() * sourceBpp * (1 + source->getWidth());
        uint8_t *targetData = m_tile->getUInt8Data();

        const FormatConverter *cvt = FormatConverter::getInstance(
            std::make_pair(Bitmap::EFloat, m_tile->getComponentFormat())
        );

        for (int i=0; i<m_blockSize; ++i) {
            if (m_pixelFormats.size() == 1)
                cvt->convert(source->getPixelFormat(), 1.0f, sourceData,
                    m_tile->getPixelFormat(), m_tile->getGamma(), targetData,
                    m_tile->getWidth());
            else
                Bitmap::convertMultiSpectrumAlphaWeight(source, sourceData,
                    m_tile, targetData, m_pixelFormats, m_componentFormat, m_tile->getWidth());

            sourceData += source->getWidth() * sourceBpp;
            targetData += m_tile->getWidth() * targetBpp;
        }

        /* Commit to disk */
        size_t ptrOffset = mergedBlock->getOffset().x * m_pixelStride +
            mergedBlock->getOffset().y * m_rowStride;

        for (Imf::FrameBuffer::Iterator it = m_frameBuffer->begin();
            it != m_frameBuffer->end(); ++it)
            it.slice().base -= ptrOffset;

        m_output->setFrameBuffer(*m_frameBuffer);
        m_output->writeTile(x, y);

        for (Imf::FrameBuffer::Iterator it = m_frameBuffer->begin();
            it != m_frameBuffer->end(); ++it)
            it.slice().base += ptrOffset;

        /* Release the block */
        m_freeBlocks.push_back(origBlock);
        m_freeBlocks.push_back(mergedBlock);
        m_origBlocks[idx] = NULL;
        m_mergedBlocks[idx] = NULL;
    }

    bool develop(const Point2i &sourceOffset, const Vector2i &size,
            const Point2i &targetOffset, Bitmap *target) const {
        target->fillRect(targetOffset, size, Spectrum(0.0f));
        return false; /* Not supported by the tiled EXR film! */
    }

    void develop(const Scene *scene, Float renderTime) {
        if (m_output) {
            Log(EInfo, "Closing EXR file (%u tiles in total, peak memory usage: %u tiles)..",
                m_blocksH * m_blocksV, m_peakUsage);
            delete m_output;
            delete m_frameBuffer;
            m_output = NULL;
            m_frameBuffer = NULL;
            m_tile = NULL;

            for (std::vector<ImageBlock *>::iterator it = m_freeBlocks.begin();
                it != m_freeBlocks.end(); ++it)
                (*it)->decRef();
            m_freeBlocks.clear();

            for (std::map<uint32_t, ImageBlock *>::iterator it = m_origBlocks.begin();
                it != m_origBlocks.end(); ++it) {
                if ((*it).second)
                    (*it).second->decRef();
            }
            m_origBlocks.clear();

            for (std::map<uint32_t, ImageBlock *>::iterator it = m_mergedBlocks.begin();
                it != m_mergedBlocks.end(); ++it) {
                if ((*it).second)
                    (*it).second->decRef();
            }
            m_mergedBlocks.clear();
        }
    }

    void clear() { /* Do nothing */ }

    bool hasAlpha() const {
        for (size_t i=0; i<m_pixelFormats.size(); ++i) {
            if (m_pixelFormats[i] == Bitmap::ELuminanceAlpha ||
                m_pixelFormats[i] == Bitmap::ERGBA ||
                m_pixelFormats[i] == Bitmap::EXYZA ||
                m_pixelFormats[i] == Bitmap::ESpectrumAlpha)
                return true;
        }
        return false;
    }

    bool destinationExists(const fs::path &baseName) const {
        fs::path filename = baseName;
        if (boost::to_lower_copy(filename.extension().string()) != ".exr")
            filename.replace_extension(".exr");
        return fs::exists(filename);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "TiledHDRFilm[" << endl
            << "  size = " << m_size.toString() << "," << endl
            << "  pixelFormat = ";
        for (size_t i=0; i<m_pixelFormats.size(); ++i)
            oss << m_pixelFormats[i] << ", ";
        oss << endl
            << "  channelNames = ";
        for (size_t i=0; i<m_channelNames.size(); ++i)
            oss << "\"" << m_channelNames[i] << "\"" << ", ";
        oss << endl
            << "  componentFormat = " << m_componentFormat << "," << endl
            << "  cropOffset = " << m_cropOffset.toString() << "," << endl
            << "  cropSize = " << m_cropSize.toString() << "," << endl
            << "  filter = " << indent(m_filter->toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
protected:
    std::vector<Bitmap::EPixelFormat> m_pixelFormats;
    std::vector<std::string> m_channelNames;
    Bitmap::EComponentFormat m_componentFormat;
    std::vector<ImageBlock *> m_freeBlocks;
    std::map<uint32_t, ImageBlock *> m_origBlocks, m_mergedBlocks;
    Imf::TiledOutputFile *m_output;
    Imf::FrameBuffer *m_frameBuffer;
    ref<Bitmap> m_tile;
    size_t m_pixelStride, m_rowStride;
    int m_blocksH, m_blocksV, m_peakUsage;
    int m_blockSize;
};

MTS_IMPLEMENT_CLASS_S(TiledHDRFilm, false, Film)
MTS_EXPORT_PLUGIN(TiledHDRFilm, "Tiled high dynamic range film");
MTS_NAMESPACE_END

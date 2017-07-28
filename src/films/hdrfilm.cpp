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
#include <mitsuba/core/statistics.h>
#include <boost/algorithm/string.hpp>
#include "banner.h"
#include "annotations.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{hdrfilm}{High dynamic range film}
 * \order{1}
 * \parameters{
 *     \parameter{width, height}{\Integer}{
 *       Width and height of the camera sensor in pixels
 *       \default{768, 576}
 *     }
 *     \parameter{fileFormat}{\String}{
 *       Denotes the desired output file format. The options
 *       are \code{openexr} (for ILM's OpenEXR format),
 *       \code{rgbe} (for Greg Ward's RGBE format),
 *       or \code{pfm} (for the Portable Float Map format)
 *       \default{\code{openexr}}
 *     }
 *     \parameter{pixelFormat}{\String}{Specifies the desired pixel format
 *         of output images. The options are \code{luminance},
 *         \code{luminanceAlpha}, \code{rgb}, \code{rgba}, \code{xyz},
 *         \code{xyza}, \code{spectrum}, and \code{spectrumAlpha}.
 *         For the \code{spectrum*} options, the number of written channels depends on
 *         the value assigned to \code{SPECTRUM\_SAMPLES} during compilation
 *         (see \secref{compiling} for details)
 *         \default{\code{rgb}}
 *     }
 *     \parameter{componentFormat}{\String}{Specifies the desired floating
 *         point component format of output images. The options are
 *         \code{float16}, \code{float32}, or \code{uint32}.
 *         \default{\code{float16}}
 *     }
 *     \parameter{cropOffsetX, cropOffsetY, cropWidth, cropHeight}{\Integer}{
 *       These parameters can optionally be provided to select a sub-rectangle
 *       of the output. In this case, Mitsuba will only render the requested
 *       regions. \default{Unused}
 *     }
 *     \parameter{attachLog}{\Boolean}{Mitsuba can optionally attach
 *         the entire rendering log file as a metadata field so that this
 *         information is permanently saved.
 *         \default{\code{true}, i.e. attach it}
 *     }
 *     \parameter{banner}{\Boolean}{Include a small Mitsuba banner in the
 *         output image? \default{\code{true}}
 *     }
 *     \parameter{highQualityEdges}{\Boolean}{
 *        If set to \code{true}, regions slightly outside of the film
 *        plane will also be sampled. This may improve the image
 *        quality at the edges, especially when using very large
 *        reconstruction filters. In general, this is not needed though.
 *        \default{\code{false}, i.e. disabled}
 *     }
 *     \parameter{\Unnamed}{\RFilter}{Reconstruction filter that should
 *     be used by the film. \default{\code{gaussian}, a windowed Gaussian filter}}
 * }
 *
 * This is the default film plugin that is used when none is explicitly
 * specified. It stores the captured image as a high dynamic range OpenEXR file
 * and tries to preserve the rendering as much as possible by not performing any
 * kind of post processing, such as gamma correction---the output file
 * will record linear radiance values.
 *
 * When writing OpenEXR files, the film will either produce a luminance, luminance/alpha,
 * RGB(A), XYZ(A) tristimulus, or spectrum/spectrum-alpha-based bitmap having a
 * \code{float16}, \code{float32}, or \code{uint32}-based internal representation
 * based on the chosen parameters.
 * The default configuration is RGB with a \code{float16} component format,
 * which is appropriate for most purposes. Note that the spectral output options
 * only make sense when using a custom build of Mitsuba that has spectral
 * rendering enabled (this is not the case for the downloadable release builds).
 * For OpenEXR files, Mitsuba also supports fully general multi-channel output;
 * refer to the \pluginref{multichannel} plugin for details on how this works.
 *
 * The plugin can also write RLE-compressed files in the Radiance RGBE format
 * pioneered by Greg Ward (set \code{fileFormat=rgbe}), as well as the
 * Portable Float Map format (set \code{fileFormat=pfm}).
 * In the former case,
 * the \code{componentFormat} and \code{pixelFormat} parameters are ignored,
 * and the output is ``\code{float8}''-compressed RGB data.
 * PFM output is restricted to \code{float32}-valued images using the
 * \code{rgb} or \code{luminance} pixel formats.
 * Due to the superior accuracy and adoption of OpenEXR, the use of these
 * two alternative formats is discouraged however.
 *
 * When RGB(A) output is selected, the measured spectral power distributions are
 * converted to linear RGB based on the CIE 1931 XYZ color matching curves and
 * the ITU-R Rec. BT.709-3 primaries with a D65 white point.
 *
 * \begin{xml}[caption=Instantiation of a film that writes a full-HD RGBA OpenEXR file without the Mitsuba banner]
 * <film type="hdrfilm">
 *     <string name="pixelFormat" value="rgba"/>
 *     <integer name="width" value="1920"/>
 *     <integer name="height" value="1080"/>
 *     <boolean name="banner" value="false"/>
 * </film>
 * \end{xml}
 *
 * \subsubsection*{Render-time annotations:}
 * \label{sec:film-annotations}
 * The \pluginref{ldrfilm} and \pluginref{hdrfilm} plugins support a
 * feature referred to as \emph{render-time annotations} to facilitate
 * record keeping.
 * Annotations are used to embed useful information inside a rendered image so
 * that this information is later available to anyone viewing the image.
 * Exemplary uses of this feature might be to store the frame or take number,
 * rendering time, memory usage, camera parameters, or other relevant scene
 * information.
 *
 * Currently, two different types are supported: a \code{metadata} annotation
 * creates an entry in the metadata table of the image, which is preferable
 * when the image contents should not be touched. Alternatively, a \code{label}
 * annotation creates a line of text that is overlaid on top of the image. Note
 * that this is only visible when opening the output file (i.e. the line is not
 * shown in the interactive viewer).
 * The syntax of this looks as follows:
 *
 * \begin{xml}
 * <film type="hdrfilm">
 *  <!-- Create a new metadata entry 'my_tag_name' and set it to the
 *       value 'my_tag_value' -->
 *  <string name="metadata['key_name']" value="Hello!"/>
 *
 *  <!-- Add the label 'Hello' at the image position X=50, Y=80 -->
 *  <string name="label[50, 80]" value="Hello!"/>
 * </film>
 * \end{xml}
 *
 * The \code{value="..."} argument may also include certain keywords that will be
 * evaluated and substituted when the rendered image is written to disk. A list all available
 * keywords is provided in Table~\ref{tbl:film-keywords}.
 *
 * Apart from querying the render time,
 * memory usage, and other scene-related information, it is also possible
 * to `paste' an existing parameter that was provided to another plugin---for instance,
 * the camera transform matrix would be obtained as \code{\$sensor['toWorld']}. The name of
 * the active integrator plugin is given by \code{\$integrator['type']}, and so on.
 * All of these can be mixed to build larger fragments, as following example demonstrates.
 * The result of this annotation is shown in Figure~\ref{fig:annotation-example}.
 * \begin{xml}[mathescape=false]
 * <string name="label[10, 10]" value="Integrator: $integrator['type'],
 *   $film['width']x$film['height'], $sampler['sampleCount'] spp,
 *   render time: $scene['renderTime'], memory: $scene['memUsage']"/>
 * \end{xml}
 * \vspace{1cm}
 * \renderings{
 * \fbox{\includegraphics[width=.8\textwidth]{images/annotation_example}}\hfill\,
 * \caption{\label{fig:annotation-example}A demonstration of the label annotation feature
 *  given the example string shown above.}
 * }
 * \vspace{2cm}
 * \begin{table}[htb]
 * \centering
 * \begin{savenotes}
 * \begin{tabular}{ll}
 * \toprule
 * \code{\$scene['renderTime']}& Image render time, use \code{renderTimePrecise} for more digits.\\
 * \code{\$scene['memUsage']}& Mitsuba memory usage\footnote{The definition of this quantity unfortunately
 * varies a bit from platform to platform. On Linux and Windows, it denotes the total
 * amount of allocated RAM and disk-based memory that is private to the process (i.e. not
 * shared or shareable), which most intuitively captures the amount of memory required for
 * rendering. On OSX, it denotes the working set size---roughly speaking, this is the
 * amount of RAM apportioned to the process (i.e. excluding disk-based memory).}.
 * Use \code{memUsagePrecise} for more digits.\\
 * \code{\$scene['coreCount']}& Number of local and remote cores working on the rendering job\\
 * \code{\$scene['blockSize']}& Block size used to parallelize up the rendering workload\\
 * \code{\$scene['sourceFile']}& Source file name\\
 * \code{\$scene['destFile']}& Destination file name\\
 * \code{\$integrator['..']}& Copy a named integrator parameter\\
 * \code{\$sensor['..']}& Copy a named sensor parameter\\
 * \code{\$sampler['..']}& Copy a named sampler parameter\\
 * \code{\$film['..']}& Copy a named film parameter\\
 * \bottomrule
 * \end{tabular}
 * \end{savenotes}
 * \caption{\label{tbl:film-keywords}A list of all special
 * keywords supported by the annotation feature}
 * \end{table}
 *
 */

class HDRFilm : public Film {
public:
    HDRFilm(const Properties &props) : Film(props) {
        /* Should an Mitsuba banner be added to the output image? */
        m_banner = props.getBoolean("banner", true);
        /* Attach the log file as the EXR comment attribute? */
        m_attachLog = props.getBoolean("attachLog", true);

        std::string fileFormat = boost::to_lower_copy(
            props.getString("fileFormat", "openexr"));
        std::vector<std::string> pixelFormats = tokenize(boost::to_lower_copy(
            props.getString("pixelFormat", "rgb")), " ,");
        std::vector<std::string> channelNames = tokenize(
            props.getString("channelNames", ""), ", ");
        std::string componentFormat = boost::to_lower_copy(
            props.getString("componentFormat", "float16"));

        if (fileFormat == "openexr") {
            m_fileFormat = Bitmap::EOpenEXR;
        } else if (fileFormat == "rgbe") {
            m_fileFormat = Bitmap::ERGBE;
        } else if (fileFormat == "pfm") {
            m_fileFormat = Bitmap::EPFM;
        } else {
            Log(EError, "The \"fileFormat\" parameter must either be "
                "equal to \"openexr\", \"pfm\", or \"rgbe\"!");
        }

        if (pixelFormats.empty())
            Log(EError, "At least one pixel format must be specified!");

        if ((pixelFormats.size() != 1 && channelNames.size() != pixelFormats.size()) ||
            (pixelFormats.size() == 1 && channelNames.size() > 1))
            Log(EError, "Number of channel names must match the number of specified pixel formats!");

        if (pixelFormats.size() != 1 && m_fileFormat != Bitmap::EOpenEXR)
            Log(EError, "General multi-channel output is only supported when writing OpenEXR files!");

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

        if (m_fileFormat == Bitmap::ERGBE) {
            /* RGBE output; override pixel & component format if necessary */
            if (m_pixelFormats.size() != 1)
                Log(EError, "The RGBE format does not support general multi-channel images!");
            if (m_pixelFormats[0] != Bitmap::ERGB) {
                Log(EWarn, "The RGBE format only supports pixelFormat=\"rgb\". Overriding..");
                m_pixelFormats[0] = Bitmap::ERGB;
            }
            if (m_componentFormat != Bitmap::EFloat32) {
                Log(EWarn, "The RGBE format only supports componentFormat=\"float32\". Overriding..");
                m_componentFormat = Bitmap::EFloat32;
            }
        } else if (m_fileFormat == Bitmap::EPFM) {
            /* PFM output; override pixel & component format if necessary */
            if (m_pixelFormats.size() != 1)
                Log(EError, "The PFM format does not support general multi-channel images!");
            if (m_pixelFormats[0] != Bitmap::ERGB && m_pixelFormats[0] != Bitmap::ELuminance) {
                Log(EWarn, "The PFM format only supports pixelFormat=\"rgb\" or \"luminance\"."
                    " Overriding (setting to \"rgb\")..");
                m_pixelFormats[0] = Bitmap::ERGB;
            }
            if (m_componentFormat != Bitmap::EFloat32) {
                Log(EWarn, "The PFM format only supports componentFormat=\"float32\". Overriding..");
                m_componentFormat = Bitmap::EFloat32;
            }
        }

        std::vector<std::string> keys = props.getPropertyNames();
        for (size_t i=0; i<keys.size(); ++i) {
            std::string key = boost::to_lower_copy(keys[i]);
            key.erase(std::remove_if(key.begin(), key.end(), ::isspace), key.end());

            if ((boost::starts_with(key, "metadata['") && boost::ends_with(key, "']")) ||
                (boost::starts_with(key, "label[") && boost::ends_with(key, "]")))
                props.markQueried(keys[i]);
        }

        if (m_pixelFormats.size() == 1) {
            m_storage = new ImageBlock(Bitmap::ESpectrumAlphaWeight, m_cropSize);
        } else {
            m_storage = new ImageBlock(Bitmap::EMultiSpectrumAlphaWeight, m_cropSize,
                NULL, (int) (SPECTRUM_SAMPLES * m_pixelFormats.size() + 2));
        }
    }

    HDRFilm(Stream *stream, InstanceManager *manager)
        : Film(stream, manager) {
        m_banner = stream->readBool();
        m_attachLog = stream->readBool();
        m_fileFormat = (Bitmap::EFileFormat) stream->readUInt();
        m_pixelFormats.resize((size_t) stream->readUInt());
        for (size_t i=0; i<m_pixelFormats.size(); ++i)
            m_pixelFormats[i] = (Bitmap::EPixelFormat) stream->readUInt();
        m_channelNames.resize((size_t) stream->readUInt());
        for (size_t i=0; i<m_channelNames.size(); ++i)
            m_channelNames[i] = stream->readString();
        m_componentFormat = (Bitmap::EComponentFormat) stream->readUInt();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Film::serialize(stream, manager);
        stream->writeBool(m_banner);
        stream->writeBool(m_attachLog);
        stream->writeUInt(m_fileFormat);
        stream->writeUInt((uint32_t) m_pixelFormats.size());
        for (size_t i=0; i<m_pixelFormats.size(); ++i)
            stream->writeUInt(m_pixelFormats[i]);
        stream->writeUInt((uint32_t) m_channelNames.size());
        for (size_t i=0; i<m_channelNames.size(); ++i)
            stream->writeString(m_channelNames[i]);
        stream->writeUInt(m_componentFormat);
    }

    void clear() {
        m_storage->clear();
    }

    void put(const ImageBlock *block) {
        m_storage->put(block);
    }

    void setBitmap(const Bitmap *bitmap, Float multiplier) {
        bitmap->convert(m_storage->getBitmap(), multiplier);
    }

    void addBitmap(const Bitmap *bitmap, Float multiplier) {
        /* Currently, only accumulating spectrum-valued floating point images
           is supported. This function basically just exists to support the
           somewhat peculiar film updates done by BDPT */

        Vector2i size = bitmap->getSize();
        if (bitmap->getPixelFormat() != Bitmap::ESpectrum ||
            bitmap->getComponentFormat() != Bitmap::EFloat ||
            bitmap->getGamma() != 1.0f ||
            size != m_storage->getSize() ||
            m_pixelFormats.size() != 1) {
            Log(EError, "addBitmap(): Unsupported bitmap format!");
        }

        size_t nPixels = (size_t) size.x * (size_t) size.y;
        const Float *source = bitmap->getFloatData();
        Float *target = m_storage->getBitmap()->getFloatData();
        for (size_t i=0; i<nPixels; ++i) {
            Float weight = target[SPECTRUM_SAMPLES + 1];
            if (weight == 0)
                weight = target[SPECTRUM_SAMPLES + 1] = 1;
            weight *= multiplier;
            for (size_t j=0; j<SPECTRUM_SAMPLES; ++j)
                *target++ += *source++ * weight;
            target += 2;
        }
    }

    bool develop(const Point2i &sourceOffset, const Vector2i &size,
            const Point2i &targetOffset, Bitmap *target) const {
        const Bitmap *source = m_storage->getBitmap();
        const FormatConverter *cvt = FormatConverter::getInstance(
            std::make_pair(Bitmap::EFloat, target->getComponentFormat())
        );

        size_t sourceBpp = source->getBytesPerPixel();
        size_t targetBpp = target->getBytesPerPixel();

        const uint8_t *sourceData = source->getUInt8Data()
            + (sourceOffset.x + sourceOffset.y * source->getWidth()) * sourceBpp;
        uint8_t *targetData = target->getUInt8Data()
            + (targetOffset.x + targetOffset.y * target->getWidth()) * targetBpp;

        if (EXPECT_NOT_TAKEN(m_pixelFormats.size() != 1)) {
            /* Special case for general multi-channel images -- just develop the first component(s) */
            for (int i=0; i<size.y; ++i) {
                for (int j=0; j<size.x; ++j) {
                    Float weight = *((Float *) (sourceData + (j+1)*sourceBpp - sizeof(Float)));
                    Float invWeight = weight != 0 ? ((Float) 1 / weight) : (Float) 0;
                    cvt->convert(Bitmap::ESpectrum, 1.0f, sourceData + j*sourceBpp,
                        target->getPixelFormat(), target->getGamma(), targetData + j * targetBpp,
                        1, invWeight);
                }

                sourceData += source->getWidth() * sourceBpp;
                targetData += target->getWidth() * targetBpp;
            }

        } else if (size.x == m_cropSize.x && target->getWidth() == m_storage->getWidth()) {
            /* Develop a connected part of the underlying buffer */
            cvt->convert(source->getPixelFormat(), 1.0f, sourceData,
                target->getPixelFormat(), target->getGamma(), targetData,
                size.x*size.y);
        } else {
            /* Develop a rectangular subregion */
            for (int i=0; i<size.y; ++i) {
                cvt->convert(source->getPixelFormat(), 1.0f, sourceData,
                    target->getPixelFormat(), target->getGamma(), targetData,
                    size.x);

                sourceData += source->getWidth() * sourceBpp;
                targetData += target->getWidth() * targetBpp;
            }
        }

        return true;
    }

    void setDestinationFile(const fs::path &destFile, uint32_t blockSize) {
        m_destFile = destFile;
    }

    void develop(const Scene *scene, Float renderTime) {
        if (m_destFile.empty())
            return;

        Log(EDebug, "Developing film ..");

        ref<Bitmap> bitmap;
        if (m_pixelFormats.size() == 1) {
            bitmap = m_storage->getBitmap()->convert(m_pixelFormats[0], m_componentFormat);
            bitmap->setChannelNames(m_channelNames);
        } else {
            bitmap = m_storage->getBitmap()->convertMultiSpectrumAlphaWeight(m_pixelFormats,
                    m_componentFormat, m_channelNames);
        }

        if (m_banner && m_cropSize.x > bannerWidth+5 && m_cropSize.y > bannerHeight + 5 && m_pixelFormats.size() == 1) {
            int xoffs = m_cropSize.x - bannerWidth - 5,
                yoffs = m_cropSize.y - bannerHeight - 5;
            for (int y=0; y<bannerHeight; y++) {
                for (int x=0; x<bannerWidth; x++) {
                    if (banner[x+y*bannerWidth])
                        continue;
                    bitmap->setPixel(Point2i(x+xoffs, y+yoffs), Spectrum(1024));
                }
            }
        }

        fs::path filename = m_destFile;
        std::string properExtension;
        if (m_fileFormat == Bitmap::EOpenEXR)
            properExtension = ".exr";
        else if (m_fileFormat == Bitmap::ERGBE)
            properExtension = ".rgbe";
        else
            properExtension = ".pfm";

        std::string extension = boost::to_lower_copy(filename.extension().string());
        if (extension != properExtension)
            filename.replace_extension(properExtension);

        Log(EInfo, "Writing image to \"%s\" ..", filename.string().c_str());
        ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);

        if (m_pixelFormats.size() == 1)
            annotate(scene, m_properties, bitmap, renderTime, 1.0f);

        /* Attach the log file to the image if this is requested */
        Logger *logger = Thread::getThread()->getLogger();
        std::string log;
        if (m_attachLog && logger->readLog(log)) {
            log += "\n\n";
            log += Statistics::getInstance()->getStats();
            bitmap->setMetadataString("log", log);
        }

        bitmap->write(m_fileFormat, stream);
    }

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
        std::string properExtension;
        if (m_fileFormat == Bitmap::EOpenEXR)
            properExtension = ".exr";
        else if (m_fileFormat == Bitmap::ERGBE)
            properExtension = ".rgbe";
        else
            properExtension = ".pfm";

        fs::path filename = baseName;
        if (boost::to_lower_copy(filename.extension().string()) != properExtension)
            filename.replace_extension(properExtension);
        return fs::exists(filename);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "HDRFilm[" << endl
            << "  size = " << m_size.toString() << "," << endl
            << "  fileFormat = " << m_fileFormat << "," << endl
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
            << "  banner = " << m_banner << "," << endl
            << "  filter = " << indent(m_filter->toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
protected:
    Bitmap::EFileFormat m_fileFormat;
    std::vector<Bitmap::EPixelFormat> m_pixelFormats;
    std::vector<std::string> m_channelNames;
    Bitmap::EComponentFormat m_componentFormat;
    bool m_banner;
    bool m_attachLog;
    fs::path m_destFile;
    ref<ImageBlock> m_storage;
};

MTS_IMPLEMENT_CLASS_S(HDRFilm, false, Film)
MTS_EXPORT_PLUGIN(HDRFilm, "High dynamic range film");
MTS_NAMESPACE_END

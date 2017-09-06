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

/*!\plugin{ldrfilm}{Low dynamic range film}
 * \order{3}
 * \parameters{
 *     \parameter{width, height}{\Integer}{
 *       Camera sensor resolution in pixels
 *       \default{768, 576}
 *     }
 *     \parameter{fileFormat}{\Integer}{
 *       The desired output file format:
 *       \code{png} or \code{jpeg}. \default{\code{png}}
 *     }
 *     \parameter{pixelFormat}{\String}{Specifies the pixel format
 *         of the generated image. The options are \code{luminance},
 *         \code{luminanceAlpha}, \code{rgb} or \code{rgba} for PNG output
 *         and \code{rgb} or \code{luminance} for JPEG output.
 *     }
 *     \parameter{tonemapMethod}{\String}{
 *       Method used to tonemap recorded radiance values
 *       \vspace{-2mm}
 *       \begin{enumerate}[(i)]
 *          \item \code{gamma}: Exposure and gamma correction (default)
 *          \vspace{-1mm}
 *          \item \code{reinhard}: Apply the
 *          tonemapping technique by Reinhard et al. \cite{Reinhard2002Photographic}
 *          followd by gamma correction.
 *       \vspace{-4mm}
 *       \end{enumerate}
 *     }
 *     \parameter{gamma}{\Float}{
 *       The gamma curve applied to correct the output image,
 *       where the special value -1 indicates sRGB. \default{-1}
 *     }
 *     \parameter{exposure}{\Float}{
 *       When \code{gamma} tonemapping is active, this parameter specifies
 *       an exposure factor in f-stops that is applied to the image before
 *       gamma correction (scaling the radiance values by $2^{\,\text{exposure}}$).
 *       \default{0, i.e. do not change the exposure}
 *     }
 *     \parameter{key}{\Float}{
 *       When \code{reinhard} tonemapping is active, this parameter in $(0,1]$ specifies
 *       whether a low-key or high-key image is desired.
 *       \default{0.18, corresponding to a middle-grey}
 *     }
 *     \parameter{burn}{\Float}{
 *       When \code{reinhard} tonemapping is active, this parameter in $[0,1]$ specifies how much
 *       highlights can burn out. \default{0, i.e. map all luminance values into the displayable range}
 *     }
 *     \parameter{banner}{\Boolean}{Include a banner in the
 *         output image?\default{\code{true}}
 *     }
 *     \parameter{cropOffsetX, cropOffsetY, cropWidth, cropHeight}{\Integer}{
 *       These parameters can optionally be provided to select a sub-rectangle
 *       of the output. In this case, Mitsuba will only render the requested
 *       regions. \default{Unused}
 *     }
 *     \parameter{highQualityEdges}{\Boolean}{
 *        If set to \code{true}, regions slightly outside of the film
 *        plane will also be sampled. This may improve image
 *        quality at the edges, but is not needed in general.
 *        \default{\code{false}}
 *     }
 *     \parameter{\Unnamed}{\RFilter}{Reconstruction filter that should
 *     be used by the film. \default{\code{gaussian}, a windowed Gaussian filter}}
 * }
 * This plugin implements a low dynamic range film that can write out 8-bit PNG
 * and JPEG images in various configurations. It provides basic tonemapping techniques
 * to map recorded radiance values into a reasonable displayable range. An alpha (opacity)
 * channel can be written if desired. By default, the plugin writes gamma-corrected
 * PNG files using the sRGB color space and no alpha channel.
 *
 * This film is a good choice when low dynamic range output is desired
 * and the rendering setup can be configured to capture the relevant portion
 * of the dynamic range reliably enough so that the original HDR data can safely
 * be discarded. When this is not the case, it may be easier to use \pluginref{hdrfilm}
 * along with the batch tonemapper (\secref{tonemapper}).
 *
 * By default, the plugin assumes that no special tonemapping needs to be done and simply
 * applies an exposure multiplier and sRGB gamma curve to the recorded radiance values
 * before converting them to 8 bit. When the dynamic range varies greatly, it may be
 * preferable to use the photographic tonemapping technique by Reinhard et al.
 * \cite{Reinhard2002Photographic}, which can be activated by setting \code{tonemapMethod=reinhard}.
 *
 * Note that the interactive tonemapper that is available in the graphical user interface \code{mtsgui}
 * interoperates with this plugin. In particular, when saving the scene
 * (\emph{File}$\to$\emph{Save}), the currently active tonemapper
 * settings are automatically exported into the updated scene file.
 *
 * The RGB values exported by this plugin correspond to the ITU-R Rec. BT. 709-3
 * primaries with a D65 white point. When $\texttt{gamma}$ is set to $\code{-1}$ (the default),
 * the output is in the sRGB color space and will display as intended on compatible devices.
 *
 * Note that this plugin supports render-time \emph{annotations}, which
 * are described on page~\pageref{sec:film-annotations}.
 */
class LDRFilm : public Film {
public:
    enum ETonemapMethod {
        EGamma,
        EReinhard
    };

    LDRFilm(const Properties &props) : Film(props) {
        /* Should an Mitsuba banner be added to the output image? */
        m_hasBanner = props.getBoolean("banner", true);

        std::string fileFormat = boost::to_lower_copy(
            props.getString("fileFormat", "png"));
        std::string pixelFormat = boost::to_lower_copy(
            props.getString("pixelFormat", "rgb"));
        std::string tonemapMethod = boost::to_lower_copy(
            props.getString("tonemapMethod", "gamma"));

        if (fileFormat == "png") {
            m_fileFormat = Bitmap::EPNG;
        } else if (fileFormat == "jpg" || fileFormat == "jpeg") {
            m_fileFormat = Bitmap::EJPEG;
        } else {
            Log(EError, "The \"fileFormat\" parameter must either be equal to "
                "\"png\" or \"jpeg\"!");
        }

        if (tonemapMethod == "gamma") {
            m_tonemapMethod = EGamma;
        } else if (tonemapMethod == "reinhard") {
            m_tonemapMethod = EReinhard;
        } else {
            Log(EError, "The \"method\" parameter must either be equal to "
                "\"gamma\" or \"reinhard\"!");
        }

        if (pixelFormat == "luminance") {
            m_pixelFormat = Bitmap::ELuminance;
        } else if (pixelFormat == "luminancealpha") {
            m_pixelFormat = Bitmap::ELuminanceAlpha;
        } else if (pixelFormat == "rgb") {
            m_pixelFormat = Bitmap::ERGB;
        } else if (pixelFormat == "rgba") {
            m_pixelFormat = Bitmap::ERGBA;
        } else {
            Log(EError, "The \"pixelFormat\" parameter must either be equal to "
                "\"luminance\", \"luminanceAlpha\", \"rgb\", or \"rgba\"!");
        }

        if (m_pixelFormat == Bitmap::ELuminanceAlpha && m_fileFormat == Bitmap::EJPEG) {
            Log(EWarn, "Ignoring request to write an alpha channel, since the JPEG format does not support it.");
            m_pixelFormat = Bitmap::ELuminance;
        } else if (m_pixelFormat == Bitmap::ERGBA && m_fileFormat == Bitmap::EJPEG) {
            Log(EWarn, "Ignoring request to write an alpha channel, since the JPEG format does not support it.");
            m_pixelFormat = Bitmap::ERGB;
        }

        m_gamma = props.getFloat("gamma", -1);
        m_exposure = props.getFloat("exposure", 0.0f);
        m_reinhardKey = props.getFloat("key", 0.18f);
        m_reinhardBurn = props.getFloat("burn", 0.0);

        std::vector<std::string> keys = props.getPropertyNames();
        for (size_t i=0; i<keys.size(); ++i) {
            std::string key = boost::to_lower_copy(keys[i]);
            key.erase(std::remove_if(key.begin(), key.end(), ::isspace), key.end());

            if ((boost::starts_with(key, "metadata['") && boost::ends_with(key, "']")) ||
                (boost::starts_with(key, "label[") && boost::ends_with(key, "]")))
                props.markQueried(keys[i]);
        }

        m_storage = new ImageBlock(Bitmap::ESpectrumAlphaWeight, m_cropSize);
    }

    LDRFilm(Stream *stream, InstanceManager *manager)
        : Film(stream, manager) {
        m_hasBanner = stream->readBool();
        m_pixelFormat = (Bitmap::EPixelFormat) stream->readUInt();
        m_fileFormat = (Bitmap::EFileFormat) stream->readUInt();
        m_gamma = stream->readFloat();
        m_tonemapMethod = (ETonemapMethod) stream->readUInt();
        m_exposure = stream->readFloat();
        m_reinhardKey = stream->readFloat();
        m_reinhardBurn = stream->readFloat();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Film::serialize(stream, manager);
        stream->writeBool(m_hasBanner);
        stream->writeUInt(m_pixelFormat);
        stream->writeUInt(m_fileFormat);
        stream->writeFloat(m_gamma);
        stream->writeUInt(m_tonemapMethod);
        stream->writeFloat(m_exposure);
        stream->writeFloat(m_reinhardKey);
        stream->writeFloat(m_reinhardBurn);
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
            size != m_storage->getSize()) {
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

        if (size.x == m_cropSize.x && target->getWidth() == m_storage->getWidth()) {
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

        ref<Bitmap> bitmap = m_storage->getBitmap();
        Float multiplier = 1.0f;

        if (m_tonemapMethod == EReinhard) {
            bitmap = bitmap->convert(m_pixelFormat, Bitmap::EFloat);

            Float logAvgLuminance = 0, maxLuminance = 0; /* Unused */
            bitmap->tonemapReinhard(logAvgLuminance, maxLuminance,
                m_reinhardKey, m_reinhardBurn);
            Log(EInfo, "Tonemapping finished (log-avg luminance=%f, max luminance=%f)",
                logAvgLuminance, maxLuminance);
        } else {
            multiplier = std::pow((Float) 2, (Float) m_exposure);
        }

        bitmap = bitmap->convert(m_pixelFormat, Bitmap::EUInt8, m_gamma, multiplier);

        if (m_hasBanner && m_cropSize.x > bannerWidth+5 && m_cropSize.y > bannerHeight + 5) {
            int xoffs = m_cropSize.x - bannerWidth - 5,
                yoffs = m_cropSize.y - bannerHeight - 5;
            for (int y=0; y<bannerHeight; y++) {
                for (int x=0; x<bannerWidth; x++) {
                    if (banner[x+y*bannerWidth])
                        continue;
                    bitmap->setPixel(Point2i(x+xoffs, y+yoffs), Spectrum(1));
                }
            }
        }

        fs::path filename = m_destFile;
        std::string extension = boost::to_lower_copy(filename.extension().string());
        std::string expectedExtension;
        switch (m_fileFormat) {
            case Bitmap::EPNG: expectedExtension = ".png"; break;
            case Bitmap::EJPEG: expectedExtension = ".jpg"; break;
            default:
                Log(EError, "Unknown file format!");
        }
        if (extension != expectedExtension)
            filename.replace_extension(expectedExtension);

        Log(EInfo, "Writing image to \"%s\" ..", filename.string().c_str());
        ref<FileStream> stream = new FileStream(filename, FileStream::ETruncWrite);

        annotate(scene, m_properties, bitmap, renderTime, m_gamma);

        bitmap->write(m_fileFormat, stream);
    }

    bool hasAlpha() const {
        return
            m_pixelFormat == Bitmap::ELuminanceAlpha ||
            m_pixelFormat == Bitmap::ERGBA;
    }

    bool destinationExists(const fs::path &baseName) const {
        fs::path filename = baseName;
        std::string extension;
        switch (m_fileFormat) {
            case Bitmap::EPNG: extension = ".png"; break;
            case Bitmap::EJPEG: extension = ".jpg"; break;
            default:
                Log(EError, "Unknown file format!");
                return false;
        }
        if (boost::to_lower_copy(filename.extension().string()) != extension)
            filename.replace_extension(extension);
        return fs::exists(filename);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "LDRFilm[" << endl
            << "  size = " << m_size.toString() << "," << endl
            << "  fileFormat = " << m_fileFormat << "," << endl
            << "  pixelFormat = " << m_pixelFormat << "," << endl
            << "  gamma = " << m_gamma << "," << endl
            << "  cropOffset = " << m_cropOffset.toString() << "," << endl
            << "  cropSize = " << m_cropSize.toString() << "," << endl
            << "  banner = " << m_hasBanner << "," << endl
            << "  method = " << ((m_tonemapMethod == EGamma) ? "gamma" : "reinhard") << "," << endl
            << "  exposure = " << m_exposure << "," << endl
            << "  reinhardKey = " << m_reinhardKey << "," << endl
            << "  reinhardBurn = " << m_reinhardBurn << "," << endl
            << "  filter = " << indent(m_filter->toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
protected:
    Bitmap::EFileFormat m_fileFormat;
    Bitmap::EPixelFormat m_pixelFormat;
    bool m_hasBanner;
    fs::path m_destFile;
    Float m_gamma;
    ref<ImageBlock> m_storage;
    ETonemapMethod m_tonemapMethod;
    Float m_exposure, m_reinhardKey, m_reinhardBurn;
};

MTS_IMPLEMENT_CLASS_S(LDRFilm, false, Film)
MTS_EXPORT_PLUGIN(LDRFilm, "Low dynamic range film");
MTS_NAMESPACE_END

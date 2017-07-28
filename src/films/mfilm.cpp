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
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/plugin.h>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/fstream.hpp>
#include <iomanip>
#include "cnpy.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{mfilm}{MATLAB / Mathematica / NumPy film}
 * \order{4}
 * \parameters{
 *     \parameter{width, height}{\Integer}{
 *       Width and height of the sensor in pixels
 *       \default{1, 1}
 *     }
 *     \parameter{cropOffsetX, cropOffsetY, cropWidth, cropHeight}{\Integer}{
 *       These parameters can optionally be provided to select a sub-rectangle
 *       of the output. In this case, Mitsuba will only render the requested
 *       regions. \default{Unused}
 *     }
 *     \parameter{fileFormat}{\String}{
 *       Specifies the desired output format; must be one of
 *       \code{matlab}, \code{mathematica}, or \code{numpy}. \default{\code{matlab}}
 *     }
 *     \parameter{digits}{\Integer}{
 *       Number of significant digits to be written \default{4}
 *     }
 *     \parameter{variable}{\String}{
 *       Name of the target variable \default{\code{"data"}}
 *     }
 *     \parameter{pixelFormat}{\String}{Specifies the desired pixel format
 *         of the generated image. The options are \code{luminance},
 *         \code{luminanceAlpha}, \code{rgb}, \code{rgba}, \code{spectrum},
 *         and \code{spectrumAlpha}. In the latter two cases,
 *         the number of written channels depends on the value assigned to
 *         \code{SPECTRUM\_SAMPLES} during compilation (see Section~\ref{sec:compiling}
 *         for details) \default{\code{luminance}}
 *     }
 *     \parameter{highQualityEdges}{\Boolean}{
 *        If set to \code{true}, regions slightly outside of the film
 *        plane will also be sampled. This may improve the image
 *        quality at the edges, especially when using very large
 *        reconstruction filters. In general (and particularly using the
 *        default box filter), this is not needed though.
 *        \default{\code{false}, i.e. disabled}
 *     }
 *     \parameter{\Unnamed}{\RFilter}{Reconstruction filter that should
 *     be used by the film. \default{\code{box}, a simple box filter}}
 * }
 *
 * \renderings{
 *     \rendering{Importing and tonemapping an image in Mathematica}{film_mfilm_mathematica.jpg}
 * }
 *
 * This plugin provides a camera film that exports spectrum, RGB, XYZ, or
 * luminance values as a matrix to a MATLAB or Mathematica ASCII file or a NumPy binary file.
 * This is useful when running Mitsuba as simulation step as part of a
 * larger virtual experiment. It can also come in handy when
 * verifying parts of the renderer using an automated test suite.
 */
class MFilm : public Film {
public:
    enum EMode {
        EMATLAB = 0,
        EMathematica,
        ENumPy
    };

    MFilm(const Properties &props) : Film(props) {
        std::string pixelFormat = boost::to_lower_copy(
            props.getString("pixelFormat", "luminance"));

        std::string fileFormat = boost::to_lower_copy(
            props.getString("fileFormat", "matlab"));

        if (pixelFormat == "luminance") {
            m_pixelFormat = Bitmap::ELuminance;
        } else if (pixelFormat == "luminancealpha") {
            m_pixelFormat = Bitmap::ELuminanceAlpha;
        } else if (pixelFormat == "rgb") {
            m_pixelFormat = Bitmap::ERGB;
        } else if (pixelFormat == "rgba") {
            m_pixelFormat = Bitmap::ERGBA;
        } else if (pixelFormat == "xyz") {
            m_pixelFormat = Bitmap::EXYZ;
        } else if (pixelFormat == "xyza") {
            m_pixelFormat = Bitmap::EXYZA;
        } else if (pixelFormat == "spectrum") {
            m_pixelFormat = Bitmap::ESpectrum;
        } else if (pixelFormat == "spectrumalpha") {
            m_pixelFormat = Bitmap::ESpectrumAlpha;
        } else {
            Log(EError, "The \"pixelFormat\" parameter must either be equal to "
                "\"luminance\", \"luminanceAlpha\", \"rgb\", \"rgba\", \"xyz\", \"xyza\", "
                "\"spectrum\", or \"spectrumAlpha\"!");
        }

        if (SPECTRUM_SAMPLES == 3 && (m_pixelFormat == Bitmap::ESpectrum || m_pixelFormat == Bitmap::ESpectrumAlpha))
            Log(EError, "You requested to render a spectral image, but Mitsuba is currently "
                "configured for a RGB flow (i.e. SPECTRUM_SAMPLES = 3). You will need to recompile "
                "it with a different configuration. Please see the documentation for details.");

        if (fileFormat == "matlab") {
            m_fileFormat = EMATLAB;
        } else if (fileFormat == "mathematica") {
            m_fileFormat = EMathematica;
        } else if (fileFormat == "numpy") {
            m_fileFormat = ENumPy;
        } else {
            Log(EError, "The \"fileFormat\" parameter must either be equal to "
                "\"matlab\" or \"mathematica\" or \"numpy\" or \"numpycompressed\"!");
        }

        m_digits = props.getInteger("digits", 4);
        m_variable = props.getString("variable", "data");

        m_storage = new ImageBlock(Bitmap::ESpectrumAlphaWeight, m_cropSize);
    }

    MFilm(Stream *stream, InstanceManager *manager)
        : Film(stream, manager) {
        m_pixelFormat = (Bitmap::EPixelFormat) stream->readUInt();
        m_fileFormat = (EMode) stream->readUInt();
        m_digits = stream->readInt();
        m_variable = stream->readString();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Film::serialize(stream, manager);
        stream->writeUInt(m_pixelFormat);
        stream->writeUInt(m_fileFormat);
        stream->writeInt(m_digits);
        stream->writeString(m_variable);
    }

    void configure() {
        if (m_filter == NULL) {
            /* No reconstruction filter has been selected. Load a box filter by default */
            m_filter = static_cast<ReconstructionFilter *> (PluginManager::getInstance()->
                    createObject(MTS_CLASS(ReconstructionFilter), Properties("box")));
            m_filter->configure();
        }

        Film::configure();
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

        fs::path filename = m_destFile;
        std::string extension = boost::to_lower_copy(filename.extension().string());
        std::string expectedExtension;
        if (m_fileFormat == EMathematica || m_fileFormat == EMATLAB) {
            expectedExtension = ".m";
        } else if (m_fileFormat == ENumPy) {
            expectedExtension = ".npy";
        } else {
            Log(EError, "Invalid file format!");
        }
        if (extension != expectedExtension)
            filename.replace_extension(expectedExtension);

        ref<Bitmap> bitmap = m_storage->getBitmap()->convert(
            m_pixelFormat, Bitmap::EFloat);

        Log(EInfo, "Writing image to \"%s\" ..", filename.filename().string().c_str());

        if (m_fileFormat == EMathematica || m_fileFormat == EMATLAB) {
            fs::ofstream os(filename);
            if (!os.good() || os.fail())
                Log(EError, "Output file cannot be created!");

            os << std::setprecision(m_digits);

            int rowSize = bitmap->getWidth();

            for (int ch=0; ch<bitmap->getChannelCount(); ++ch) {
                if (m_fileFormat == EMATLAB) {
                    if (ch == 0) {
                        os << m_variable << " = [";
                    } else {
                        os << endl << m_variable << "(:, :, " << ch + 1 << ") = [";
                    }
                } else {
                    if (ch == 0) {
                        if (bitmap->getChannelCount() == 1)
                            os << m_variable << " = {{";
                        else
                            os << m_variable << " = Transpose[{{{";
                    }
                }
                Float *ptr = bitmap->getFloatData();
                ptr += ch;

                for (int y=0; y < bitmap->getHeight(); y++) {
                    for (int x=0; x < rowSize; x++) {
                        if (m_fileFormat == EMATLAB) {
                            os << *ptr;
                        } else {
                            /* Mathematica uses the peculiar '*^' notation rather than the standard 'e' notation. */
                            std::ostringstream oss;
                            oss << std::setprecision(m_digits);
                            oss << *ptr;
                            std::string str = oss.str();
                            boost::replace_first(str, "e", "*^");
                            os << str;
                        }

                        ptr += bitmap->getChannelCount();
                        if (x + 1 < rowSize) {
                            os << ", ";
                        } else {
                            if (m_fileFormat == EMATLAB) {
                                if (y + 1 < bitmap->getHeight())
                                    os << ";" << endl << "\t";
                                else
                                    os << "];" << endl;
                            } else {
                                if (y + 1 < bitmap->getHeight()) {
                                    os << "}," << endl << "\t{";
                                } else if (ch + 1 == bitmap->getChannelCount()) {
                                    if (bitmap->getChannelCount() == 1)
                                        os << "}};" << endl;
                                    else
                                        os << "}}}, {3,1,2}];" << endl;
                                } else {
                                    os << "}}," << endl << endl << "\t{{";
                                }
                            }
                        }
                    }
                }
            }
        } else {
            unsigned int shape[] = {
                (unsigned int) bitmap->getHeight(),
                (unsigned int) bitmap->getWidth(),
                (unsigned int) bitmap->getChannelCount()
            };
            unsigned int N = 3, *shape_ptr = shape;

            if (bitmap->getChannelCount() == 1)
                N = 2;

            const Float *data = bitmap->getFloatData();
            cnpy::npy_save(filename.string(), data, shape_ptr, N, "w");
        }
    }

    bool destinationExists(const fs::path &baseName) const {
        fs::path filename = baseName;
        std::string expectedExtension;
        if (m_fileFormat == EMathematica || m_fileFormat == EMATLAB) {
            expectedExtension = ".m";
        } else if (m_fileFormat == ENumPy) {
            expectedExtension = ".npy";
        } else {
            Log(EError, "Invalid file format!");
        }
        if (boost::to_lower_copy(filename.extension().string()) != expectedExtension)
            filename.replace_extension(expectedExtension);
        return fs::exists(filename);
    }

    bool hasAlpha() const {
        return
            m_pixelFormat == Bitmap::ELuminanceAlpha ||
            m_pixelFormat == Bitmap::ERGBA ||
            m_pixelFormat == Bitmap::EXYZA ||
            m_pixelFormat == Bitmap::ESpectrumAlpha;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "MFilm[" << endl
            << "  size = " << m_size.toString() << "," << endl
            << "  pixelFormat = " << m_pixelFormat << "," << endl
            << "  digits = " << m_digits << "," << endl
            << "  variable = \"" << m_variable << "\"," << endl
            << "  cropOffset = " << m_cropOffset.toString() << "," << endl
            << "  cropSize = " << m_cropSize.toString() << "," << endl
            << "  filter = " << indent(m_filter->toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
protected:
    Bitmap::EPixelFormat m_pixelFormat;
    EMode m_fileFormat;
    fs::path m_destFile;
    ref<ImageBlock> m_storage;
    std::string m_variable;
    int m_digits;
};

MTS_IMPLEMENT_CLASS_S(MFilm, false, Film)
MTS_EXPORT_PLUGIN(MFilm, "MATLAB / Mathematica / NumPy film");
MTS_NAMESPACE_END

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

#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/mipmap.h>
#include <mitsuba/hw/renderer.h>
#include <mitsuba/hw/gputexture.h>
#include <mitsuba/hw/gpuprogram.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

/*!\plugin{bitmap}{Bitmap texture}
 * \order{1}
 * \parameters{
 *     \parameter{filename}{\String}{
 *       Filename of the bitmap to be loaded
 *     }
 *     \parameter{wrapMode, wrapModeU, wrapModeV}{\String}{
 *       Behavior of texture lookups outside of the $[0,1]$ $uv$ range.\vspace{-1mm}
 *       \begin{enumerate}[(i)]
 *           \item \code{repeat}: Repeat the texture indefinitely\vspace{-1mm}
 *           \item \code{mirror}: Mirror the texture along its boundaries\vspace{-1mm}
 *           \item \code{clamp}: Clamp $uv$ coordinates to $[0,1]$ before a lookup\vspace{-1mm}
 *           \item \code{zero}: Switch to a zero-valued texture \vspace{-1mm}
 *           \item \code{one}: Switch to a one-valued texture \vspace{-1mm}
 *       \end{enumerate}
 *       Default: \code{repeat}. The parameter \code{wrapMode} is a shortcut for
 *       setting both \code{wrapModeU} and \code{wrapModeV} at the same time.
 *     }
 *     \parameter{gamma}{\Float}{
 *       Optional parameter to override the gamma value of the source bitmap,
 *       where 1 indicates a linear color space and the special value -1
 *       corresponds to sRGB. \default{automatically detect based on the
 *       image type and metadata}
 *     }
 *     \parameter{filterType}{\String}{
 *       Specifies the texture filturing that should be used for lookups\vspace{-1mm}
 *       \begin{enumerate}[(i)]
 *           \item \code{ewa}: Elliptically weighted average (a.k.a.
 *           anisotropic filtering). This produces the best quality\vspace{-1mm}
 *           \item \code{trilinear}: Simple trilinear (isotropic) filtering.\vspace{-1mm}
 *           \item \code{nearest}: No filtering, do nearest neighbor lookups.\vspace{-1mm}
 *       \end{enumerate}
 *       Default: \code{ewa}.
 *     }
 *     \parameter{maxAnisotropy}{\Float}{
 *        Specific to \code{ewa} filtering, this parameter limits the
 *        anisotropy (and thus the computational cost) of filtured texture lookups. The
 *        default of 20 is a good compromise.
 *     }
 *     \parameter{cache}{\Boolean}{
 *        Preserve generated MIP map data in a cache file? This will cause a file named
 *        \emph{filename}\code{.mip} to be created.
 *        \default{automatic---use caching for textures larger than 1M pixels.}
 *     }
 *     \parameter{uoffset, voffset}{\Float}{
 *       Numerical offset that should be applied to UV lookups
 *     }
 *     \parameter{uscale, vscale}{\Float}{
 *       Multiplicative factors that should be applied to UV lookups
 *     }
 *     \parameter{channel}{\String}{
 *       Create a monochromatic texture based on one of the image channels
 *       (e.g. \texttt{r}, \texttt{g}, \texttt{b}, \texttt{a}, \texttt{x},
 *       \texttt{y}, \texttt{z} etc.). \default{use all channels}
 *     }
 * }
 * This plugin provides a bitmap-backed texture source that supports \emph{filtered}
 * texture lookups on\footnote{Some of these may not be available depending on how
 * Mitsuba was compiled.} JPEG, PNG, OpenEXR, RGBE, TGA, and BMP files. Filtered
 * lookups are useful to avoid aliasing when rendering textures that contain high
 * frequencies (see the next page for an example).
 *
 * The plugin operates as follows: when loading a bitmap file, it is first converted
 * into a linear color space. Following this, a MIP map is constructed that is necessary
 * to perform filtered lookups during rendering. A \emph{MIP map} is a hierarchy of
 * progressively lower resolution versions of the input image, where the resolution of
 * adjacent levels differs by a factor of two. Mitsuba creates this hierarchy using
 * Lanczos resampling to obtain very high quality results.
 * Note that textures may have an arbitrary resolution and are not limited to powers of two.
 * Three different filtering modes are supported:
 *
 * \begin{enumerate}[(i)]
 * \item Nearest neighbor lookups effectively disable filtering and always query the
 * highest-resolution version of the texture without any kind of interpolation. This is
 * fast and requires little memory (no MIP map is created), but results in visible aliasing.
 * Only a single pixel value is accessed.
 *
 * \item The trilinear filter performs bilinear interpolation on two adjacent MIP levels
 * and blends the results. Because it cannot do anisotropic (i.e. slanted) lookups in texture space,
 * it must compromise either on the side of blurring or aliasing. The implementation in Mitsuba
 * chooses blurring over aliasing (though note that (\textbf{b}) is an extreme case).
 * Only 8 pixel values are accessed.
 *
 * \item The EWA filter performs anisotropicically filtered lookups on two adjacent MIP map levels
 * and blends them. This produces the best quality, but at the expense of computation time.
 * Generally, 20-40 pixel values must be read for a single EWA texture lookup. To limit
 * the number of pixel accesses, the \code{maxAnisotropy} parameter can be used to bound
 * the amount of anisotropy that a texture lookup is allowed to have.
 * \end{enumerate}
 * \renderings{
 *     \rendering{Nearest-neighbor filter. Note the aliasing}{tex_bitmap_nearest}
 *     \rendering{Trilinear filter. Note the blurring}{tex_bitmap_trilinear}
 *     \vspace{-5mm}
 * }
 * \renderings{
 *     \setcounter{subfigure}{2}
 *     \rendering{EWA filter}{tex_bitmap_ewa}
 *     \rendering{Ground truth (512 samples per pixel)}{tex_bitmap_gt}
 *     \caption{A somewhat contrived comparison of the different filters when rendering a high-frequency
 *     checkerboard pattern using four samples per pixel. The EWA method (the default)
 *     pre-filters the texture anisotropically to limit blurring and aliasing, but has a
 *     higher computational cost than the other filters.}
 * }
 * \paragraph{Caching and memory requirements:}
 * When a texture is read, Mitsuba internally converts it into an uncompressed linear format
 * using a half precision (\code{float16})-based representation. This is convenient for
 * rendering but means that textures require copious amounts of memory (in particular, the
 * size of the occupied memory region might be orders of magnitude greater than that of the
 * original input file).
 *
 * For instance, a basic 10 megapixel image requires as much as 76 MiB of memory! Loading,
 * color space transformation, and MIP map construction require up to several seconds in this case.
 * To reduce these overheads, Mitsuba 0.4.0 introduced MIP map caches. When a large
 * texture is loaded for the first time, a MIP map cache file with the name \emph{filename}\code{.mip}
 * is generated. This is essentially a verbatim copy of the in-memory representation created
 * during normal rendering. Storing this information as a separate file has two advantages:
 *
 * \begin{enumerate}[(i)]
 *    \item MIP maps do not have to be regenerated in subsequent Mitsuba runs,
 *     which substantially reduces scene loading times.
 *    \item Because the texture storage is entirely disk-backed and can be \emph{memory-mapped},
 *    Mitsuba is able to work with truly massive textures that would otherwise exhaust the main system memory.
 * \end{enumerate}
 *
 * The texture caches are automatically regenerated when the input texture is modified.
 * Of course, the cache files can be cumbersome when they are not needed anymore. On Linux
 * or Mac OS, they can safely be deleted by executing the following command within a scene directory.
 *
 * \begin{shell}
 * $\code{\$}$ find . -name "*.mip" -delete
 * \end{shell}
 */

class BitmapTexture : public Texture2D {
public:
    /* Store texture data using half precision, but perform computations in
       single/double precision based on compilation flags. The following
       generates efficient implementations for both luminance and RGB data */
    typedef TSpectrum<Float, 1> Color1;
    typedef TSpectrum<Float, 3> Color3;
    typedef TSpectrum<half, 1>  Color1h;
    typedef TSpectrum<half, 3>  Color3h;
    typedef TMIPMap<Color1, Color1h> MIPMap1;
    typedef TMIPMap<Color3, Color3h> MIPMap3;

    BitmapTexture(const Properties &props) : Texture2D(props) {
        uint64_t timestamp = 0;
        bool tryReuseCache = false;
        fs::path cacheFile;
        ref<Bitmap> bitmap;

        m_channel = boost::to_lower_copy(props.getString("channel", ""));

        if (props.hasProperty("bitmap")) {
            /* Support initialization via raw data passed from another plugin */
            bitmap = reinterpret_cast<Bitmap *>(props.getData("bitmap").ptr);
        } else {
            m_filename = Thread::getThread()->getFileResolver()->resolve(
                props.getString("filename"));

            Log(EInfo, "Loading texture \"%s\"", m_filename.filename().string().c_str());
            if (!fs::exists(m_filename))
                Log(EError, "Texture file \"%s\" could not be found!", m_filename.string().c_str());

            boost::system::error_code ec;
            timestamp = (uint64_t) fs::last_write_time(m_filename, ec);
            if (ec.value())
                Log(EError, "Could not determine modification time of \"%s\"!", m_filename.string().c_str());

            cacheFile = m_filename;

            if (m_channel.empty())
                cacheFile.replace_extension(".mip");
            else
                cacheFile.replace_extension(formatString(".%s.mip", m_channel.c_str()));

            tryReuseCache = fs::exists(cacheFile) && props.getBoolean("cache", true);
        }

        std::string filterType = boost::to_lower_copy(props.getString("filterType", "ewa"));
        std::string wrapMode = props.getString("wrapMode", "repeat");
        m_wrapModeU = parseWrapMode(props.getString("wrapModeU", wrapMode));
        m_wrapModeV = parseWrapMode(props.getString("wrapModeV", wrapMode));

        m_gamma = props.getFloat("gamma", 0);

        if (filterType == "ewa")
            m_filterType = EEWA;
        else if (filterType == "bilinear")
            m_filterType = EBilinear;
        else if (filterType == "trilinear")
            m_filterType = ETrilinear;
        else if (filterType == "nearest")
            m_filterType = ENearest;
        else
            Log(EError, "Unknown filter type '%s' -- must be "
                "'ewa', 'trilinear', or 'nearest'!", filterType.c_str());

        m_maxAnisotropy = props.getFloat("maxAnisotropy", 20);

        if (m_filterType != EEWA)
            m_maxAnisotropy = 1.0f;

        if (tryReuseCache && MIPMap3::validateCacheFile(cacheFile, timestamp,
                Bitmap::ERGB, m_wrapModeU, m_wrapModeV, m_filterType, m_gamma)) {
            /* Reuse an existing MIP map cache file */
            m_mipmap3 = new MIPMap3(cacheFile, m_maxAnisotropy);
        } else if (tryReuseCache && MIPMap1::validateCacheFile(cacheFile, timestamp,
                Bitmap::ELuminance, m_wrapModeU, m_wrapModeV, m_filterType, m_gamma)) {
            /* Reuse an existing MIP map cache file */
            m_mipmap1 = new MIPMap1(cacheFile, m_maxAnisotropy);
        } else {
            if (bitmap == NULL) {
                /* Load the input image if necessary */
                ref<Timer> timer = new Timer();
                ref<FileStream> fs = new FileStream(m_filename, FileStream::EReadOnly);
                bitmap = new Bitmap(Bitmap::EAuto, fs);
                if (m_gamma != 0)
                    bitmap->setGamma(m_gamma);
                Log(EDebug, "Loaded \"%s\" in %i ms", m_filename.filename().string().c_str(),
                    timer->getMilliseconds());
            }

            Bitmap::EPixelFormat pixelFormat;
            if (!m_channel.empty()) {
                /* Create a texture from a certain channel of an image */
                pixelFormat = Bitmap::ELuminance;
                bitmap = bitmap->extractChannel(findChannel(bitmap, m_channel));
                if (m_channel == "a")
                    bitmap->setGamma(1.0f);
            } else {
                switch (bitmap->getPixelFormat()) {
                    case Bitmap::ELuminance:
                    case Bitmap::ELuminanceAlpha:
                        pixelFormat = Bitmap::ELuminance;
                        break;
                    case Bitmap::ERGB:
                    case Bitmap::ERGBA:
                        pixelFormat = Bitmap::ERGB;
                        break;
                    default:
                        Log(EError, "The input image has an unsupported pixel format!");
                        return;
                }
            }

            /* (Re)generate the MIP map hierarchy; downsample using a
                2-lobed Lanczos reconstruction filter */
            Properties rfilterProps("lanczos");
            rfilterProps.setInteger("lobes", 2);
            ref<ReconstructionFilter> rfilter = static_cast<ReconstructionFilter *> (
                PluginManager::getInstance()->createObject(
                MTS_CLASS(ReconstructionFilter), rfilterProps));
            rfilter->configure();

            /* Potentially create a new MIP map cache file */
            bool createCache = !cacheFile.empty() && props.getBoolean("cache",
                bitmap->getSize().x * bitmap->getSize().y > 1024*1024);

            if (pixelFormat == Bitmap::ELuminance)
                m_mipmap1 = new MIPMap1(bitmap, pixelFormat, Bitmap::EFloat,
                    rfilter, m_wrapModeU, m_wrapModeV, m_filterType, m_maxAnisotropy,
                    createCache ? cacheFile : fs::path(), timestamp);
            else
                m_mipmap3 = new MIPMap3(bitmap, pixelFormat, Bitmap::EFloat,
                    rfilter, m_wrapModeU, m_wrapModeV, m_filterType, m_maxAnisotropy,
                    createCache ? cacheFile : fs::path(), timestamp);
        }
    }

    static int findChannel(const Bitmap *bitmap, const std::string channel) {
        int found = -1;
        std::string channelNames;
        for (int i=0; i<bitmap->getChannelCount(); ++i) {
            std::string name = boost::to_lower_copy(bitmap->getChannelName(i));
            if (name == channel)
                found = i;
            channelNames += name;
            if (i + 1 < bitmap->getChannelCount())
                channelNames += std::string(", ");
        }

        if (found == -1) {
            Log(EError, "Channel \"%s\" not found! Must be one of: [%s]",
                channel.c_str(), channelNames.c_str());
        }

        return found;
    }

    inline ReconstructionFilter::EBoundaryCondition parseWrapMode(const std::string &wrapMode) {
        if (wrapMode == "repeat")
            return ReconstructionFilter::ERepeat;
        else if (wrapMode == "clamp")
            return ReconstructionFilter::EClamp;
        else if (wrapMode == "mirror")
            return ReconstructionFilter::EMirror;
        else if (wrapMode == "zero" || wrapMode == "black")
            return ReconstructionFilter::EZero;
        else if (wrapMode == "one" || wrapMode == "white")
            return ReconstructionFilter::EOne;
        else
            Log(EError, "Unknown wrap mode '%s' -- must be "
                "'repeat', 'clamp', 'black', or 'white'!", wrapMode.c_str());
        return ReconstructionFilter::EZero; // make gcc happy
    }

    BitmapTexture(Stream *stream, InstanceManager *manager)
     : Texture2D(stream, manager) {
        m_filename = stream->readString();
        Log(EDebug, "Unserializing texture \"%s\"", m_filename.filename().string().c_str());
        m_filterType = (EMIPFilterType) stream->readUInt();
        m_wrapModeU = (ReconstructionFilter::EBoundaryCondition) stream->readUInt();
        m_wrapModeV = (ReconstructionFilter::EBoundaryCondition) stream->readUInt();
        m_gamma = stream->readFloat();
        m_maxAnisotropy = stream->readFloat();
        m_channel = stream->readString();

        size_t size = stream->readSize();
        ref<MemoryStream> mStream = new MemoryStream(size);
        stream->copyTo(mStream, size);
        mStream->seek(0);
        ref<Bitmap> bitmap = new Bitmap(Bitmap::EAuto, mStream);
        if (m_gamma != 0)
            bitmap->setGamma(m_gamma);

        /* Downsample using a 2-lobed Lanczos reconstruction filter */
        Properties rfilterProps("lanczos");
        rfilterProps.setInteger("lobes", 2);
        ref<ReconstructionFilter> rfilter = static_cast<ReconstructionFilter *> (
            PluginManager::getInstance()->createObject(
            MTS_CLASS(ReconstructionFilter), rfilterProps));
        rfilter->configure();

        Bitmap::EPixelFormat pixelFormat;
        if (!m_channel.empty()) {
            /* Create a texture from a certain channel of an image */
            pixelFormat = Bitmap::ELuminance;
            bitmap = bitmap->extractChannel(findChannel(bitmap, m_channel));
            if (m_channel == "a")
                bitmap->setGamma(1.0f);
        } else {
            switch (bitmap->getPixelFormat()) {
                case Bitmap::ELuminance:
                case Bitmap::ELuminanceAlpha:
                    pixelFormat = Bitmap::ELuminance;
                    break;
                case Bitmap::ERGB:
                case Bitmap::ERGBA:
                    pixelFormat = Bitmap::ERGB;
                    break;
                default:
                    Log(EError, "The input image has an unsupported pixel format!");
                    return;
            }
        }

        if (pixelFormat == Bitmap::ELuminance)
            m_mipmap1 = new MIPMap1(bitmap, pixelFormat, Bitmap::EFloat,
                rfilter, m_wrapModeU, m_wrapModeV, m_filterType, m_maxAnisotropy,
                fs::path(), 0);
        else
            m_mipmap3 = new MIPMap3(bitmap, pixelFormat, Bitmap::EFloat,
                rfilter, m_wrapModeU, m_wrapModeV, m_filterType, m_maxAnisotropy,
                fs::path(), 0);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Texture2D::serialize(stream, manager);
        stream->writeString(m_filename.string());
        stream->writeUInt(m_filterType);
        stream->writeUInt(m_wrapModeU);
        stream->writeUInt(m_wrapModeV);
        stream->writeFloat(m_gamma);
        stream->writeFloat(m_maxAnisotropy);

        if (!m_filename.empty() && fs::exists(m_filename)) {
            /* We still have access to the original image -- use that, since
               it is probably much smaller than the in-memory representation */
            ref<Stream> is = new FileStream(m_filename, FileStream::EReadOnly);
            stream->writeString(m_channel);
            stream->writeSize(is->getSize());
            is->copyTo(stream);
        } else {
            /* No access to the original image anymore. Create an EXR image
               from the top MIP map level and serialize that */
            ref<MemoryStream> mStream = new MemoryStream();
            ref<Bitmap> bitmap = m_mipmap1.get() ?
                m_mipmap1->toBitmap() : m_mipmap3->toBitmap();
            bitmap->write(Bitmap::EOpenEXR, mStream);

            stream->writeString("");
            stream->writeSize(mStream->getSize());
            stream->write(mStream->getData(), mStream->getSize());
        }
    }

    Spectrum eval(const Point2 &uv) const {
        /* There are no ray differentials to do any kind of
           prefiltering. Evaluate the full-resolution texture */

        Spectrum result;
        if (m_mipmap3.get()) {
            Color3 value;
            if (m_mipmap3->getFilterType() != ENearest)
                value = m_mipmap3->evalBilinear(0, uv);
            else
                value = m_mipmap3->evalBox(0, uv);
            result.fromLinearRGB(value[0], value[1], value[2]);
        } else {
            Color1 value;
            if (m_mipmap1->getFilterType() != ENearest)
                value = m_mipmap1->evalBilinear(0, uv);
            else
                value = m_mipmap1->evalBox(0, uv);
            result = Spectrum(value[0]);
        }
        stats::filteredLookups.incrementBase();

        return result;
    }

    void evalGradient(const Point2 &uv, Spectrum *gradient) const {
        /* There are no ray differentials to do any kind of
           prefiltering. Evaluate the full-resolution texture */

        if (m_mipmap3.get()) {
            Color3 result[2];
            if (m_mipmap3->getFilterType() != ENearest) {
                m_mipmap3->evalGradientBilinear(0, uv, result);
                gradient[0].fromLinearRGB(result[0][0], result[0][1], result[0][2]);
                gradient[1].fromLinearRGB(result[1][0], result[1][1], result[1][2]);
            } else {
                gradient[0] = gradient[1] = Spectrum(0.0f);
            }
        } else {
            Color1 result[2];
            if (m_mipmap1->getFilterType() != ENearest) {
                m_mipmap1->evalGradientBilinear(0, uv, result);
                gradient[0] = Spectrum(result[0][0]);
                gradient[1] = Spectrum(result[1][0]);
            } else {
                gradient[0] = gradient[1] = Spectrum(0.0f);
            }
        }
        stats::filteredLookups.incrementBase();
    }

    ref<Bitmap> getBitmap(const Vector2i &/* unused */) const {
        return m_mipmap1.get() ? m_mipmap1->toBitmap() : m_mipmap3->toBitmap();
    }

    Spectrum eval(const Point2 &uv, const Vector2 &d0, const Vector2 &d1) const {
        stats::filteredLookups.incrementBase();
        ++stats::filteredLookups;

        Spectrum result;
        if (m_mipmap3.get()) {
            Color3 value = m_mipmap3->eval(uv, d0, d1);
            result.fromLinearRGB(value[0], value[1], value[2]);
        } else {
            Color1 value = m_mipmap1->eval(uv, d0, d1);
            result = Spectrum(value[0]);
        }
        return result;
    }

    Spectrum getAverage() const {
        Spectrum result;
        if (m_mipmap3.get()) {
            Color3 value = m_mipmap3->getAverage();
            result.fromLinearRGB(value[0], value[1], value[2]);
        } else {
            Color1 value = m_mipmap1->getAverage();
            result = Spectrum(value[0]);
        }
        return result;
    }

    Spectrum getMaximum() const {
        Spectrum result;
        if (m_mipmap3.get()) {
            Color3 value = m_mipmap3->getMaximum();
            result.fromLinearRGB(value[0], value[1], value[2]);
        } else {
            Color1 value = m_mipmap1->getMaximum();
            result = Spectrum(value[0]);
        }
        return result;
    }

    Spectrum getMinimum() const {
        Spectrum result;
        if (m_mipmap3.get()) {
            Color3 value = m_mipmap3->getMinimum();
            result.fromLinearRGB(value[0], value[1], value[2]);
        } else {
            Color1 value = m_mipmap1->getMinimum();
            result = Spectrum(value[0]);
        }
        return result;
    }

    bool isConstant() const {
        return false;
    }

    bool usesRayDifferentials() const {
        return true;
    }

    bool isMonochromatic() const {
        return m_mipmap1.get() != NULL;
    }

    Vector3i getResolution() const {
        if (m_mipmap3.get()) {
            return Vector3i(
                m_mipmap3->getWidth(),
                m_mipmap3->getHeight(),
                1
            );
        } else {
            return Vector3i(
                m_mipmap1->getWidth(),
                m_mipmap1->getHeight(),
                1
            );
        }
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "BitmapTexture[" << endl
            << "  filename = \"" << m_filename.string() << "\"," << endl;

        if (m_mipmap3.get())
            oss << "  mipmap = " << indent(m_mipmap3.toString()) << endl;
        else
            oss << "  mipmap = " << indent(m_mipmap1.toString()) << endl;

        oss << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
protected:
    ref<MIPMap1> m_mipmap1;
    ref<MIPMap3> m_mipmap3;
    EMIPFilterType m_filterType;
    ReconstructionFilter::EBoundaryCondition m_wrapModeU;
    ReconstructionFilter::EBoundaryCondition m_wrapModeV;
    Float m_gamma, m_maxAnisotropy;
    std::string m_channel;
    fs::path m_filename;
};

// ================ Hardware shader implementation ================
class BitmapTextureShader : public Shader {
public:
    BitmapTextureShader(Renderer *renderer, const std::string &filename,
            const BitmapTexture::MIPMap1* mipmap1,
            const BitmapTexture::MIPMap3* mipmap3,
            const Point2 &uvOffset, const Vector2 &uvScale,
            ReconstructionFilter::EBoundaryCondition wrapModeU,
            ReconstructionFilter::EBoundaryCondition wrapModeV,
            Float maxAnisotropy)
        : Shader(renderer, ETextureShader), m_uvOffset(uvOffset), m_uvScale(uvScale) {

        ref<Bitmap> bitmap = mipmap1 ? mipmap1->toBitmap() : mipmap3->toBitmap();
        m_gpuTexture = renderer->createGPUTexture(filename, bitmap);

        switch (wrapModeU) {
            case ReconstructionFilter::EClamp:
                m_gpuTexture->setWrapType(GPUTexture::EClampToEdge);
                break;
            case ReconstructionFilter::EMirror:
                m_gpuTexture->setWrapType(GPUTexture::EMirror);
                break;
            case ReconstructionFilter::ERepeat:
                m_gpuTexture->setWrapType(GPUTexture::ERepeat);
                break;
            case ReconstructionFilter::EZero:
                m_gpuTexture->setWrapType(GPUTexture::EClampToBorder);
                m_gpuTexture->setBorderColor(Color3(0.0f));
                break;
            case ReconstructionFilter::EOne:
                m_gpuTexture->setWrapType(GPUTexture::EClampToBorder);
                m_gpuTexture->setBorderColor(Color3(1.0f));
                break;
            default:
                Log(EError, "Unknown wrap mode!");
        }

        switch (mipmap1 ? mipmap1->getFilterType() : mipmap3->getFilterType()) {
            case ENearest:
                m_gpuTexture->setFilterType(GPUTexture::ENearest);
                break;
            case EBilinear:
                m_gpuTexture->setFilterType(GPUTexture::ELinear);
                m_gpuTexture->setMipMapped(false);
                break;
            default:
                m_gpuTexture->setFilterType(GPUTexture::EMipMapLinear);
                break;
        }

        m_gpuTexture->setMaxAnisotropy(maxAnisotropy);
        m_gpuTexture->setMaxAnisotropy(maxAnisotropy);
        m_gpuTexture->initAndRelease();
    }

    void cleanup(Renderer *renderer) {
        m_gpuTexture->cleanup();
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "uniform sampler2D " << evalName << "_texture;" << endl
            << "uniform vec2 " << evalName << "_uvOffset;" << endl
            << "uniform vec2 " << evalName << "_uvScale;" << endl
            << endl
            << "vec3 " << evalName << "(vec2 uv) {" << endl
            << "    return texture2D(" << evalName << "_texture, vec2(" << endl
            << "          uv.x * " << evalName << "_uvScale.x + " << evalName << "_uvOffset.x," << endl
            << "          uv.y * " << evalName << "_uvScale.y + " << evalName << "_uvOffset.y)).rgb;" << endl
            << "}" << endl;
    }

    void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
        parameterIDs.push_back(program->getParameterID(evalName + "_texture", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_uvOffset", false));
        parameterIDs.push_back(program->getParameterID(evalName + "_uvScale", false));
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs,
        int &textureUnitOffset) const {
        m_gpuTexture->bind(textureUnitOffset++);
        program->setParameter(parameterIDs[0], m_gpuTexture.get());
        program->setParameter(parameterIDs[1], m_uvOffset);
        program->setParameter(parameterIDs[2], m_uvScale);
    }

    void unbind() const {
        m_gpuTexture->unbind();
    }

    MTS_DECLARE_CLASS()
private:
    ref<GPUTexture> m_gpuTexture;
    Point2 m_uvOffset;
    Vector2 m_uvScale;
};

Shader *BitmapTexture::createShader(Renderer *renderer) const {
    return new BitmapTextureShader(renderer, m_filename.filename().string(),
            m_mipmap1.get(), m_mipmap3.get(), m_uvOffset, m_uvScale,
            m_wrapModeU, m_wrapModeV, m_maxAnisotropy);
}

MTS_IMPLEMENT_CLASS_S(BitmapTexture, false, Texture2D)
MTS_IMPLEMENT_CLASS(BitmapTextureShader, false, Shader)
MTS_EXPORT_PLUGIN(BitmapTexture, "Bitmap texture");
MTS_NAMESPACE_END

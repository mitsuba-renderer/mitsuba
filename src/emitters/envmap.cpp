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

#include <mitsuba/render/scene.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/render/mipmap.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/hw/gputexture.h>

MTS_NAMESPACE_BEGIN

#if SPECTRUM_SAMPLES == 3
# define ENVMAP_PIXELFORMAT Bitmap::ERGB
#else
# define ENVMAP_PIXELFORMAT Bitmap::ESpectrum
#endif

/*!\plugin{envmap}{Environment emitter}
 * \icon{emitter_envmap}
 * \order{9}
 * \parameters{
 *     \parameter{filename}{\String}{
 *       Filename of the radiance-valued input image to be loaded;
 *       must be in latitude-longitude format.
 *     }
 *     \parameter{scale}{\Float}{
 *         A scale factor that is applied to the
 *         radiance values stored in the input image. \default{1}
 *     }
 *     \parameter{toWorld}{\Transform}{
 *	      Specifies an optional linear emitter-to-world space rotation.
 *        \default{none (i.e. emitter space $=$ world space)}
 *     }
 *     \parameter{gamma}{\Float}{
 *       Optional parameter to override the gamma value of the source bitmap,
 *       where 1 indicates a linear color space and the special value -1
 *       corresponds to sRGB. \default{automatically detect based on the
 *       image type and metadata}
 *     }
 *     \parameter{cache}{\Boolean}{
 *        Preserve generated MIP map data in a cache file? This will cause a file named
 *        \emph{filename}\code{.mip} to be created.
 *        \default{automatic---use caching for images larger than 1M pixels.}
 *     }
 *     \parameter{samplingWeight}{\Float}{
 *         Specifies the relative amount of samples
 *         allocated to this emitter. \default{1}
 *     }
 * }
 * \renderings{
 *   \rendering{The museum environment map by Bernhard Vogl that is used
 *   in many example renderings in this document}{emitter_envmap_example}
 *   \unframedrendering{Coordinate conventions used when mapping the input
 *   image onto the sphere.}{emitter_envmap_axes}
 * }
 *
 * This plugin provides a HDRI (high dynamic range imaging) environment map,
 * which is a type of light source that is well-suited for representing ``natural''
 * illumination. Many images in this document are made using the environment map
 * shown in (a).
 *
 * The implementation loads a captured illumination environment from a image in
 * latitude-longitude format and turns it into an infinitely distant emitter.
 * The image could either be a processed photograph or a rendering made using the
 * \pluginref{spherical} sensor. The direction conventions of this transformation
 * are shown in (b).
 * The plugin can work with all types of images that are natively supported by Mitsuba
 * (i.e. JPEG, PNG, OpenEXR, RGBE, TGA, and BMP). In practice, a good environment
 * map will contain high-dynamic range data that can only be represented using the
 * OpenEXR or RGBE file formats.
 * High quality free light probes are available on Paul Debevec's website
 * (\url{http://gl.ict.usc.edu/Data/HighResProbes}) and Bernhard Vogl's website
 * (\url{http://dativ.at/lightprobes/}).
 *
 * Like the \pluginref{bitmap} texture, this plugin generates a cache file
 * named \emph{filename}\code{.mip} when given a large input image. This
 * significantly accelerates the loading times of subsequent renderings. When this
 * is not desired, specify \code{cache=false} to the plugin.
 */
class EnvironmentMap : public Emitter {
public:
	/* Store the environment in a blocked MIP map using half precision */
	typedef TSpectrum<half, SPECTRUM_SAMPLES> SpectrumHalf;
	typedef TMIPMap<Spectrum, SpectrumHalf> MIPMap;

	EnvironmentMap(const Properties &props) : Emitter(props),
			m_mipmap(NULL), m_cdfRows(NULL), m_cdfCols(NULL), m_rowWeights(NULL) {
		m_type |= EOnSurface | EEnvironmentEmitter;
		uint64_t timestamp = 0;
		bool tryReuseCache = false;
		fs::path cacheFile;
		ref<Bitmap> bitmap;

		if (props.hasProperty("bitmap")) {
			/* Support initialization via raw data passed from another plugin */
			bitmap = reinterpret_cast<Bitmap *>(props.getData("bitmap").ptr);
		} else {
			m_filename = Thread::getThread()->getFileResolver()->resolve(
				props.getString("filename"));

			Log(EInfo, "Loading environment map \"%s\"", m_filename.filename().string().c_str());
			if (!fs::exists(m_filename))
				Log(EError, "Environment map file \"%s\" could not be found!", m_filename.string().c_str());

			boost::system::error_code ec;
			timestamp = (uint64_t) fs::last_write_time(m_filename, ec);
			if (ec.value())
				Log(EError, "Could not determine modification time of \"%s\"!", m_filename.string().c_str());

			/* Create MIP map a cache when the environment map is large, and
			   reuse cache files that have been created previously */
			cacheFile = m_filename;
			cacheFile.replace_extension(".mip");
			tryReuseCache = fs::exists(cacheFile) && props.getBoolean("cache", true);
		}

		/* Gamma override */
		m_gamma = props.getFloat("gamma", 0);

		/* These are reasonable MIP map defaults for environment maps, I don't
		   think there is a need to expose them through plugin parameters */
		EMIPFilterType filterType = EEWA;
		Float maxAnisotropy = 10.0f;

		if (tryReuseCache && MIPMap::validateCacheFile(cacheFile, timestamp,
				ENVMAP_PIXELFORMAT, ReconstructionFilter::ERepeat,
				ReconstructionFilter::EClamp, filterType, m_gamma)) {
			/* Reuse an existing MIP map cache file */
			m_mipmap = new MIPMap(cacheFile, maxAnisotropy);
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

			if (std::max(bitmap->getWidth(), bitmap->getHeight()) > 0xFFFF)
				Log(EError, "Environment maps images must be smaller than 65536 "
					" pixels in width and height");

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

			m_mipmap = new MIPMap(bitmap, ENVMAP_PIXELFORMAT, Bitmap::EFloat,
				rfilter, ReconstructionFilter::ERepeat, ReconstructionFilter::EClamp,
				filterType, maxAnisotropy, createCache ? cacheFile : fs::path(), timestamp,
				std::numeric_limits<Float>::infinity(), Spectrum::EIlluminant);
		}

		if (props.hasProperty("intensityScale"))
			Log(EError, "The 'intensityScale' parameter has been deprecated and is now called scale.");

		/* Scale factor */
		m_scale = props.getFloat("scale", 1.0f);
	}

	EnvironmentMap(Stream *stream, InstanceManager *manager) : Emitter(stream, manager),
			m_mipmap(NULL), m_cdfRows(NULL), m_cdfCols(NULL), m_rowWeights(NULL) {
		m_filename = stream->readString();
		Log(EDebug, "Unserializing texture \"%s\"", m_filename.filename().string().c_str());
		m_gamma = stream->readFloat();
		m_scale = stream->readFloat();
		m_sceneBSphere = BSphere(stream);
		m_geoBSphere = BSphere(stream);

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

		m_mipmap = new MIPMap(bitmap, ENVMAP_PIXELFORMAT, Bitmap::EFloat, rfilter,
			ReconstructionFilter::ERepeat, ReconstructionFilter::EClamp, EEWA, 10.0f,
			fs::path(), 0, std::numeric_limits<Float>::infinity(), Spectrum::EIlluminant);

		configure();
	}

	virtual ~EnvironmentMap() {
		if (m_mipmap)
			delete m_mipmap;
		if (m_cdfRows)
			delete[] m_cdfRows;
		if (m_cdfCols)
			delete[] m_cdfCols;
		if (m_rowWeights)
			delete[] m_rowWeights;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Emitter::serialize(stream, manager);
		stream->writeString(m_filename.string());
		stream->writeFloat(m_gamma);
		stream->writeFloat(m_scale);
		m_sceneBSphere.serialize(stream);
		m_geoBSphere.serialize(stream);

		if (!m_filename.empty() && fs::exists(m_filename)) {
			/* We still have access to the original image -- use that, since
			   it is probably much smaller than the in-memory representation */
			ref<Stream> is = new FileStream(m_filename, FileStream::EReadOnly);
			stream->writeSize(is->getSize());
			is->copyTo(stream);
		} else {
			/* No access to the original image anymore. Create an EXR image
			   from the top MIP map level and serialize that */
			ref<MemoryStream> mStream = new MemoryStream();
			ref<Bitmap> bitmap = m_mipmap->toBitmap();
			bitmap->write(Bitmap::EOpenEXR, mStream);

			stream->writeSize(mStream->getSize());
			stream->write(mStream->getData(), mStream->getSize());
		}
	}

	void configure() {
		Emitter::configure();

		if (!m_rowWeights) {
			/// Build CDF tables to sample the environment map
			const MIPMap::Array2DType &array = m_mipmap->getArray();
			m_size = array.getSize();

			size_t nEntries = (size_t) (m_size.x + 1) * (size_t) m_size.y,
				totalStorage = sizeof(float) * (m_size.x + 1 + nEntries);

			Log(EInfo, "Precomputing data structures for environment map sampling (%s)",
				memString(totalStorage).c_str());

			ref<Timer> timer = new Timer();
			m_cdfCols = new float[nEntries];
			m_cdfRows = new float[m_size.y + 1];
			m_rowWeights = new Float[m_size.y];

			size_t colPos = 0, rowPos = 0;
			Float rowSum = 0.0f;

			/* Build a marginal & conditional cumulative distribution
			   function over luminances weighted by sin(theta) */
			m_cdfRows[rowPos++] = 0;
			for (int y=0; y<m_size.y; ++y) {
				Float colSum = 0;

				m_cdfCols[colPos++] = 0;
				for (int x=0; x<m_size.x; ++x) {
					Spectrum value(array(x, y));

					colSum += value.getLuminance();
					m_cdfCols[colPos++] = (float) colSum;
				}

				float normalization = 1.0f / (float) colSum;
				for (int x=1; x<m_size.x; ++x)
					m_cdfCols[colPos-x-1] *= normalization;
				m_cdfCols[colPos-1] = 1.0f;

				Float weight = std::sin((y + 0.5f) * M_PI / m_size.y);
				m_rowWeights[y] = weight;
				rowSum += colSum * weight;
				m_cdfRows[rowPos++] = (float) rowSum;
			}

			float normalization = 1.0f / (float) rowSum;
			for (int y=1; y<m_size.y; ++y)
				m_cdfRows[rowPos-y-1] *= normalization;
			m_cdfRows[rowPos-1] = 1.0f;

			if (rowSum == 0)
				Log(EError, "The environment map is completely black -- this is not allowed.");
			else if (!std::isfinite(rowSum))
				Log(EError, "The environment map contains an invalid floating"
					" point value (nan/inf) -- giving up.");

			m_normalization = 1.0f / (rowSum *
				(2 * M_PI / m_size.x) * (M_PI / m_size.y));

			/* Size of a pixel in spherical coordinates */
			m_pixelSize = Vector2(2 * M_PI / m_size.x, M_PI / m_size.y);

			Log(EInfo, "Done (took %i ms)", timer->getMilliseconds());
		}
		Float surfaceArea = 4 * M_PI * m_sceneBSphere.radius * m_sceneBSphere.radius;
		m_invSurfaceArea = 1 / surfaceArea;
		m_power = surfaceArea * m_scale / m_normalization;
	}

	ref<Shape> createShape(const Scene *scene) {
		/* Create a bounding sphere that surrounds the scene */
		BSphere sceneBSphere(scene->getAABB().getBSphere());
		sceneBSphere.radius = std::max(Epsilon, sceneBSphere.radius * 1.5f);
		BSphere geoBSphere(scene->getKDTree()->getAABB().getBSphere());

		if (sceneBSphere != m_sceneBSphere || geoBSphere != m_geoBSphere) {
			m_sceneBSphere = sceneBSphere;
			m_geoBSphere = geoBSphere;
			configure();
		}

		Transform trafo =
			Transform::translate(Vector(m_sceneBSphere.center)) *
			Transform::scale(Vector(m_sceneBSphere.radius));

		Properties props("sphere");
		props.setTransform("toWorld", trafo);
		props.setBoolean("flipNormals", true);
		Shape *shape = static_cast<Shape *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Shape), props));
		shape->addChild(this);
		shape->configure();

		return shape;
	}

	bool fillDirectSamplingRecord(DirectSamplingRecord &dRec, const Ray &ray) const {
		Float nearT, farT;

		if (!m_sceneBSphere.rayIntersect(ray, nearT, farT) || nearT > 0 || farT < 0) {
			Log(EWarn, "fillDirectSamplingRecord(): internal error!");
			return false;
		}

		dRec.p = ray(farT);
		dRec.n = normalize(m_sceneBSphere.center - dRec.p);
		dRec.measure = ESolidAngle;
		dRec.object = this;
		dRec.d = ray.d;
		dRec.dist = farT;

		return true;
	}

	Spectrum eval(const Intersection &its, const Vector &d) const {
		return evalEnvironment(RayDifferential(its.p, -d, its.time));
	}

	Spectrum evalEnvironment(const RayDifferential &ray) const {
		const Transform &trafo = m_worldTransform->eval(ray.time);
		Vector v = trafo.inverse()(ray.d);

		/* Convert to latitude-longitude texture coordinates */
		Point2 uv(
			std::atan2(v.x, -v.z) * INV_TWOPI,
			math::safe_acos(v.y) * INV_PI
		);

		Spectrum value;
		if (!ray.hasDifferentials) {
			/* No differentials - perform a lookup at the highest level */
			value = m_mipmap->evalBilinear(0, uv);
		} else {
			/* Compute texture-space partials and perform a filtered lookup */
			Vector dvdx = trafo.inverse()(ray.rxDirection) - v,
			       dvdy = trafo.inverse()(ray.ryDirection) - v;

			Float t1 = INV_TWOPI / (v.x*v.x+v.z*v.z),
			      t2 = -INV_PI / std::max(math::safe_sqrt(1.0f-v.y*v.y), Epsilon);

			Vector2 dudx(t1 * (dvdx.z*v.x - dvdx.x*v.z), t2 * dvdx.y),
			        dudy(t1 * (dvdy.z*v.x - dvdy.x*v.z), t2 * dvdy.y);

			++stats::filteredLookups;
			value = m_mipmap->eval(uv, dudx, dudy);
		}
		stats::filteredLookups.incrementBase();
		return value * m_scale;
	}

	Spectrum samplePosition(PositionSamplingRecord &pRec, const Point2 &sample,
			const Point2 *extra) const {
		Vector d = warp::squareToUniformSphere(sample);

		pRec.p = m_sceneBSphere.center + d * m_sceneBSphere.radius;
		pRec.n = -d;
		pRec.measure = EArea;
		pRec.pdf = m_invSurfaceArea;

		return Spectrum(m_power);
	}

	Spectrum evalPosition(const PositionSamplingRecord &pRec) const {
		return Spectrum(m_power * m_invSurfaceArea);
	}

	Float pdfPosition(const PositionSamplingRecord &pRec) const {
		return m_invSurfaceArea;
	}

	/**
	 * A note regarding sampleArea()/sampleDirection():
	 *
	 * Mitsuba's emitter/sensor API permits sampling the spatial and directional
	 * components separately (in that order!). As nice as this is for bidirectional
	 * techniques, it is unfortunately very difficult to reconcile with the default
	 * environment map strategy of importance sampling a direction *first* from the
	 * environment, followed by (*second*) uniformly picking a position on a hyperplane
	 * perpendicular to that direction.
	 *
	 * Therefore, the following compromise is implemented:
	 * 1. sampleArea(): uniformly sample a position on the scene's bounding sphere
	 * 2. sampleDirection(): importance sample a direction from the environment map
	 *
	 * This is not that great, as many useless samples will be generated. Note however
	 * that rendering scenes by tracing from the environment side is usually a
	 * terrible strategy in any case, so at least not much is lost by making it worse.
	 *
	 * The direct illumination sampling code (\ref sampleDirect()), and the combined
	 * position + direction code (\ref sampleRay) both does the right thing and are
	 * not affected by any of this. Many of the bidirectional rendering
	 * algorithms in Mitsuba can make use of direct illumination sampling strategies.
	 */
	Spectrum sampleDirection(DirectionSamplingRecord &dRec,
			PositionSamplingRecord &pRec,
			const Point2 &sample,
			const Point2 *extra) const {
		const Transform &trafo = m_worldTransform->eval(pRec.time);

		/* Sample a direction from the environment map */
		Spectrum value; Vector d; Float pdf;
		internalSampleDirection(sample, d, value, pdf);

		dRec.measure = ESolidAngle;
		dRec.pdf = pdf;
		dRec.d = trafo(-d);

		/* Be wary of roundoff errors */
		if (value.isZero() || pdf == 0)
			return Spectrum(0.0f);
		else
			return (value * m_normalization) / (pdf * m_scale);
	}

	Float pdfDirection(const DirectionSamplingRecord &dRec,
			const PositionSamplingRecord &pRec) const {
		const Transform &trafo = m_worldTransform->eval(pRec.time);
		return internalPdfDirection(-trafo.inverse()(dRec.d));
	}

	Spectrum evalDirection(const DirectionSamplingRecord &dRec,
			const PositionSamplingRecord &pRec) const {
		const Transform &trafo = m_worldTransform->eval(pRec.time);
		Vector v = -trafo.inverse()(dRec.d);

		/* Convert to latitude-longitude texture coordinates */
		Point2 uv(
			std::atan2(v.x, -v.z) * INV_TWOPI,
			math::safe_acos(v.y) * INV_PI
		);

		stats::filteredLookups.incrementBase();

		return m_mipmap->evalBilinear(0, uv) * m_normalization;
	}

	Spectrum sampleRay(Ray &ray,
			const Point2 &spatialSample,
			const Point2 &directionalSample,
			Float time) const {
		const Transform &trafo = m_worldTransform->eval(time);
		Vector d; Spectrum value; Float pdf;
		internalSampleDirection(directionalSample, d, value, pdf);
		d = -trafo(d);
		Point2 offset = warp::squareToUniformDiskConcentric(spatialSample);
		Vector perpOffset = Frame(d).toWorld(Vector(offset.x, offset.y, 0));

		ray.setOrigin(m_geoBSphere.center + (perpOffset - d) * m_geoBSphere.radius);
		ray.setDirection(d);
		ray.setTime(time);

		return value * M_PI * m_geoBSphere.radius * m_geoBSphere.radius / pdf;
	}

	Spectrum sampleDirect(DirectSamplingRecord &dRec, const Point2 &sample) const {
		const Transform &trafo = m_worldTransform->eval(dRec.time);

		/* Sample a direction from the environment map */
		Spectrum value; Vector d; Float pdf;
		internalSampleDirection(sample, d, value, pdf);

		/* Intersect against the scene's bounding sphere. This may
		   seem somewhat excessive, but it's needed by the bidirectional
		   integrators that expect all the different sampling methods in
		   this class to be consistent with respect to each other. */
		Ray ray(dRec.ref, trafo(d), 0);
		Float nearT, farT;
		if (value.isZero() || pdf == 0 || !m_sceneBSphere.rayIntersect(ray, nearT, farT)
			|| nearT >= 0 || farT <= 0) {
			dRec.pdf = 0.0f;
			return Spectrum(0.0f);
		}

		dRec.pdf = pdf;
		dRec.p = ray(farT);
		dRec.n = normalize(m_sceneBSphere.center - dRec.p);
		dRec.dist = farT;
		dRec.d = ray.d;
		dRec.measure = ESolidAngle;

		return value / pdf;
	}

	Float pdfDirect(const DirectSamplingRecord &dRec) const {
		const Transform &trafo = m_worldTransform->eval(dRec.time);
		Float pdfSA = internalPdfDirection(trafo.inverse()(dRec.d));

		if (dRec.measure == ESolidAngle)
			return pdfSA;
		else if (dRec.measure == EArea)
			return pdfSA * absDot(dRec.d, dRec.n)
				/ (dRec.dist * dRec.dist);
		else
			return 0.0f;
	}

	AABB getAABB() const {
		/* The scene sets its bounding box so that it contains all shapes and
		   emitters, but this particular emitter always wants to be *a little*
		   bigger than the scene. To avoid a silly recursion, just return a
		   point here. */
		return AABB(m_sceneBSphere.center);
	}

	/// Helper function that samples a direction from the environment map
	void internalSampleDirection(Point2 sample, Vector &d, Spectrum &value, Float &pdf) const {
		/* Sample a discrete pixel position */
		uint32_t row = sampleReuse(m_cdfRows, m_size.y, sample.y),
		         col = sampleReuse(m_cdfCols + row * (m_size.x+1), m_size.x, sample.x);

		/* Using the remaining bits of precision to shift the sample by an offset
		   drawn from a tent function. This effectively creates a sampling strategy
		   for a linearly interpolated environment map */
		Point2 pos = Point2((Float) col, (Float) row) + warp::squareToTent(sample);

		/* Bilinearly interpolate colors from the adjacent four neighbors */
		int xPos = math::floorToInt(pos.x), yPos = math::floorToInt(pos.y);
		Float dx1 = pos.x - xPos, dx2 = 1.0f - dx1,
		      dy1 = pos.y - yPos, dy2 = 1.0f - dy1;

		Spectrum value1 = m_mipmap->evalTexel(0, xPos, yPos) * dx2 * dy2
		                + m_mipmap->evalTexel(0, xPos + 1, yPos) * dx1 * dy2;
		Spectrum value2 = m_mipmap->evalTexel(0, xPos, yPos + 1) * dx2 * dy1
		                + m_mipmap->evalTexel(0, xPos + 1, yPos + 1) * dx1 * dy1;
		stats::filteredLookups.incrementBase();

		/* Compute the final color and probability density of the sample */
		value = (value1 + value2) * m_scale;
		pdf = (value1.getLuminance() * m_rowWeights[math::clamp(yPos,   0, m_size.y-1)] +
		       value2.getLuminance() * m_rowWeights[math::clamp(yPos+1, 0, m_size.y-1)]) * m_normalization;

		/* Turn into a proper direction on the sphere */
		Float sinPhi, cosPhi, sinTheta, cosTheta;
		math::sincos(m_pixelSize.x * (pos.x + 0.5f), &sinPhi, &cosPhi);
		math::sincos(m_pixelSize.y * (pos.y + 0.5f), &sinTheta, &cosTheta);

		d = Vector(sinPhi*sinTheta, cosTheta, -cosPhi*sinTheta);
		pdf /= std::max(std::abs(sinTheta), Epsilon);
	}

	/// Helper function that computes the solid angle density of \ref internalSampleDirection()
	Float internalPdfDirection(const Vector &d) const {
		/* Convert to latitude-longitude texture coordinates */
		Point2 uv(
			std::atan2(d.x, -d.z) * INV_TWOPI,
			math::safe_acos(d.y) * INV_PI
		);

		if (EXPECT_NOT_TAKEN(!std::isfinite(uv.x) || !std::isfinite(uv.y))) {
			Log(EWarn, "pdfDirect(): encountered a NaN!");
			return 0.0f;
		}

		/* Convert to fractional pixel coordinates on the specified level */
		Float u = uv.x * m_size.x - 0.5f, v = uv.y * m_size.y - 0.5f;

		/* Bilinearly interpolate colors from the adjacent four neighbors */
		int xPos = math::floorToInt(u), yPos = math::floorToInt(v);
		Float dx1 = u - xPos, dx2 = 1.0f - dx1,
		      dy1 = v - yPos, dy2 = 1.0f - dy1;

		Spectrum value1 = m_mipmap->evalTexel(0, xPos, yPos) * dx2 * dy2
		                + m_mipmap->evalTexel(0, xPos + 1, yPos) * dx1 * dy2;
		Spectrum value2 = m_mipmap->evalTexel(0, xPos, yPos + 1) * dx2 * dy1
		                + m_mipmap->evalTexel(0, xPos + 1, yPos + 1) * dx1 * dy1;
		stats::filteredLookups.incrementBase();

		Float sinTheta = math::safe_sqrt(1-d.y*d.y);
		return (value1.getLuminance() * m_rowWeights[math::clamp(yPos,   0, m_size.y-1)] +
		        value2.getLuminance() * m_rowWeights[math::clamp(yPos+1, 0, m_size.y-1)])
			* m_normalization / std::max(std::abs(sinTheta), Epsilon);
	}

	ref<Bitmap> getBitmap(const Vector2i &/* unused */) const {
		return m_mipmap->toBitmap();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "EnvironmentMap[" << endl
			<< "  filename = \"" << m_filename.string() << "\"," << endl
			<< "  samplingWeight = " << m_samplingWeight << "," << endl
			<< "  bsphere = " << m_sceneBSphere.toString() << "," << endl
			<< "  worldTransform = " << indent(m_worldTransform.toString()) << "," << endl
			<< "  mipmap = " << indent(m_mipmap->toString()) << "," << endl
			<< "  medium = " << indent(m_medium.toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	/// Sample from an array using the inversion method
	inline uint32_t sampleReuse(float *cdf, uint32_t size, Float &sample) const {
		float *entry = std::lower_bound(cdf, cdf+size+1, (float) sample);
		uint32_t index = std::min((uint32_t) std::max((ptrdiff_t) 0, entry - cdf - 1), size-1);
		sample = (sample - (Float) cdf[index]) / (Float) (cdf[index+1] - cdf[index]);
		return index;
	}
private:
	MIPMap *m_mipmap;
	float *m_cdfRows, *m_cdfCols;
	Float *m_rowWeights;
	fs::path m_filename;
	Float m_gamma, m_scale;
	Float m_normalization;
	Float m_power;
	Float m_invSurfaceArea;
	BSphere m_geoBSphere;
	BSphere m_sceneBSphere;
	Vector2i m_size;
	Vector2 m_pixelSize;
};

// ================ Hardware shader implementation ================

class EnvironmentMapShader : public Shader {
public:
	EnvironmentMapShader(Renderer *renderer, const fs::path &filename, ref<Bitmap> bitmap,
			const Transform &worldToEmitter, Float scale) : Shader(renderer, EEmitterShader) {
		std::string name = filename.filename().string();
		if (name.empty())
			name = "Environment map";
		m_gpuTexture = renderer->createGPUTexture(name, bitmap);
		m_gpuTexture->setWrapTypeU(GPUTexture::ERepeat);
		m_gpuTexture->setWrapTypeV(GPUTexture::EClamp);
		m_gpuTexture->setMaxAnisotropy(8);
		m_gpuTexture->initAndRelease();
		m_worldToEmitter = worldToEmitter;
		m_scale = scale;
		m_useCustomTextureFiltering =
			renderer->getCapabilities()->isSupported(RendererCapabilities::ECustomTextureFiltering);
	}

	void cleanup(Renderer *renderer) {
		m_gpuTexture->cleanup();
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_texture", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_worldToEmitter", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_scale", false));
	}

	void generateCode(std::ostringstream &oss, const std::string &evalName,
			const std::vector<std::string> &depNames) const {

		oss << "uniform sampler2D " << evalName << "_texture;" << endl
			<< "uniform mat4 " << evalName << "_worldToEmitter;" << endl
			<< "uniform float " << evalName << "_scale;" << endl
			<< endl
			<< "vec3 " << evalName << "_dir(vec3 wo) {" << endl
			<< "   return vec3(1.0);" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_background(vec3 v_) {" << endl
			<< "   const float inv_pi = 0.318309886183791, inv_twopi = 0.159154943091895;" << endl
			<< "   vec3 v = normalize((" << evalName << "_worldToEmitter * vec4(v_, 0.0)).xyz);" << endl
			<< "   vec2 coords = vec2(atan(v.x, -v.z) * inv_twopi, acos(v.y) * inv_pi);" << endl;

		if (m_useCustomTextureFiltering) {
			oss << "   /* Manually compute derivatives to work around discontinuities */" << endl
				<< "   vec3 dvdx = dFdx(v), dvdy = dFdy(v);" << endl
				<< "   float t1 = inv_twopi / (v.x*v.x+v.z*v.z)," << endl
				<< "         t2 = -inv_pi / sqrt(1.0-v.y*v.y);" << endl
				<< "   vec2 dudx = vec2(t1 * (dvdx.z*v.x - dvdx.x*v.z), t2 * dvdx.y)," << endl
				<< "        dudy = vec2(t1 * (dvdy.z*v.x - dvdy.x*v.z), t2 * dvdy.y);" << endl
				<< "   if (abs(v.y) > 0.99) { dudx = dudy = vec2(0.0); /* Don't blur over vertical edges */ }" << endl
				<< "   return texture2DGrad(" << evalName << "_texture, coords, dudx, dudy).rgb * " << evalName << "_scale;" << endl;
		} else {
			oss << "   return texture2D(" << evalName << "_texture, coords).rgb * " << evalName << "_scale;" << endl;
		}

		oss << "}" << endl;
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs,
			int &textureUnitOffset) const {
		m_gpuTexture->bind(textureUnitOffset++);
		program->setParameter(parameterIDs[0], m_gpuTexture.get());
		program->setParameter(parameterIDs[1], m_worldToEmitter);
		program->setParameter(parameterIDs[2], m_scale);
	}

	void unbind() const {
		m_gpuTexture->unbind();
	}

	MTS_DECLARE_CLASS()
private:
	ref<GPUTexture> m_gpuTexture;
	bool m_useCustomTextureFiltering;
	Transform m_worldToEmitter;
	Float m_scale;
};

Shader *EnvironmentMap::createShader(Renderer *renderer) const {
	Transform trafo = m_worldTransform->eval(0).inverse();
	ref<Bitmap> bitmap = m_mipmap->toBitmap();
#if SPECTRUM_SAMPLES != 3
	bitmap = bitmap->convert(Bitmap::ERGB, bitmap->getComponentFormat());
#endif
	return new EnvironmentMapShader(renderer, m_filename, bitmap, trafo, m_scale);
}

MTS_IMPLEMENT_CLASS(EnvironmentMapShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(EnvironmentMap, false, Emitter)
MTS_EXPORT_PLUGIN(EnvironmentMap, "Environment map");
MTS_NAMESPACE_END

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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/volume.h>
#include <fstream>

MTS_NAMESPACE_BEGIN

/*
 * When inverting an integral of the form f(t):=\int_0^t [...] dt
 * using composite Simpson's rule quadrature, the following constant
 * specified how many times the step size will be reduced when we
 * have almost reached the desired value.
 */
#define NUM_REFINES 8

/**
 * Allow to stop integrating densities when the resulting segment
 * has a throughput of less than 'Epsilon'
 */
#define HET_EARLY_EXIT 1

/**
 * Flexible heterogeneous medium implementation, which acquires its data from 
 * nested <tt>Volume</tt> instances. These can be constant, use a procedural 
 * function, or fetch data from disk, e.g. using a memory-mapped density grid.
 *
 * Instead of allowing separate volumes to be provided for the scattering
 * parameters sigma_s and sigma_t, this class instead takes the approach of 
 * enforcing a spectrally uniform sigma_t, which must be provided using a 
 * nested scalar-valued volume named 'density'.
 *
 * A nested spectrum-valued 'albedo' volume must also be provided, which is 
 * used to compute the parameter sigma_s using the expression
 * "sigma_s = density * albedo" (i.e. 'albedo' contains the single-scattering
 * albedo of the medium).
 *
 * Optionally, one can also provide an vector-valued 'orientation' volume,
 * which contains local particle orientations that will be passed to
 * scattering models such as Kajiya-Kay phase function.
 */
class HeterogeneousMedium : public Medium {
public:
	HeterogeneousMedium(const Properties &props) 
		: Medium(props) {
		m_kdTree = new ShapeKDTree();
		m_stepSize = props.getFloat("stepSize", 0);
	}

	/* Unserialize from a binary data stream */
	HeterogeneousMedium(Stream *stream, InstanceManager *manager) 
		: Medium(stream, manager) {
		m_kdTree = new ShapeKDTree();
		size_t shapeCount = stream->readUInt();
		for (size_t i=0; i<shapeCount; ++i) 
			addChild("", static_cast<Shape *>(manager->getInstance(stream)));
		m_densities = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_albedo = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_orientations = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_stepSize = stream->readFloat();
		configure();
	}

	virtual ~HeterogeneousMedium() {
		for (size_t i=0; i<m_shapes.size(); ++i)
			m_shapes[i]->decRef();
	}

	/* Serialize the volume to a binary data stream */
	void serialize(Stream *stream, InstanceManager *manager) const {
		Medium::serialize(stream, manager);
		stream->writeUInt((uint32_t) m_shapes.size());
		for (size_t i=0; i<m_shapes.size(); ++i)
			manager->serialize(stream, m_shapes[i]);
		manager->serialize(stream, m_densities.get());
		manager->serialize(stream, m_albedo.get());
		manager->serialize(stream, m_orientations.get());
		stream->writeFloat(m_stepSize);
	}

	void configure() {
		Medium::configure();
		if (m_densities.get() == NULL)
			Log(EError, "No densities specified!");
		if (m_albedo.get() == NULL)
			Log(EError, "No albedo specified!");

		if (m_shapes.size() != 0) {
			m_kdTree->build();
			m_aabb = m_kdTree->getAABB();
		} else {
			m_aabb = m_densities->getAABB();
		}

		if (m_stepSize == 0) {
			m_stepSize = std::min(
				m_densities->getStepSize(), m_albedo->getStepSize());
			if (m_orientations != NULL)
				m_stepSize = std::min(m_stepSize,
					m_orientations->getStepSize());

			if (m_stepSize == std::numeric_limits<Float>::infinity()) 
				Log(EError, "Unable to infer step size, please specify!");
		}
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Shape))) {
			Shape *shape = static_cast<Shape *>(child);
			if (shape->isCompound()) {
				int ctr = 0;
				while (true) {
					ref<Shape> childShape = shape->getElement(ctr++);
					if (!childShape)
						break;
					addChild("", childShape);
				}
			} else {
				m_kdTree->addShape(shape);
				shape->incRef();
				m_shapes.push_back(shape);
			}
		} else if (child->getClass()->derivesFrom(MTS_CLASS(VolumeDataSource))) {
			VolumeDataSource *volume = static_cast<VolumeDataSource *>(child);

			if (name == "albedo") {
				Assert(volume->supportsSpectrumLookups());
				m_albedo = volume;
			} else if (name == "density") {
				Assert(volume->supportsFloatLookups());
				m_densities = volume;
			} else if (name == "orientations") {
				Assert(volume->supportsVectorLookups());
				m_orientations = volume;
			}
		} else {
			Medium::addChild(name, child);
		}
	}
	
	/**
	 * \brief Intersect a ray with the stencil and 
	 * determine whether it started on the inside
	 */
	inline bool rayIntersect(const Ray &ray, Float &t, bool &inside) const {
		if (m_shapes.size() != 0) {
			Intersection its;
			if (!m_kdTree->rayIntersect(ray, its)) 
				return false;
			inside = dot(its.geoFrame.n, ray.d) > 0;
			t = its.t;
		} else {
			Float mint, maxt;
			if (!m_aabb.rayIntersect(ray, mint, maxt))
				return false;
			if (mint <= ray.mint && maxt <= ray.mint)
				return false;
			inside = mint <= 0 && maxt > 0;
			t = (mint <= 0) ? maxt : mint;
		}
		return true;
	}

	/// Integrate densities using the composite Simpson's rule
	inline Float integrateDensity(Ray ray, Float length) const {
		if (length == 0)
			return 0.0f;
		int nParts = (int) std::ceil(length/m_stepSize);
		nParts += nParts % 2;
		const Float stepSize = length/nParts;
		const Vector increment = ray.d * stepSize;

		Float m = 4;
		Float accumulatedDensity = m_densities->lookupFloat(ray.o) 
			+ m_densities->lookupFloat(ray(length));
#ifdef HET_EARLY_EXIT
		const Float stopAfterDensity = -std::log(Epsilon);
		const Float stopValue = stopAfterDensity*3/(stepSize
				* m_densityMultiplier);
#endif

		ray.o += increment;
		for (int i=1; i<nParts; ++i) {
			Float value = m_densities->lookupFloat(ray.o);
			accumulatedDensity += value*m;
			Point tmp = ray.o;
			ray.o += increment;

			if (ray.o == tmp) {
				Log(EWarn, "integrateDensity(): failed to make forward progress -- "
						"round-off error issues? The step size was %f", stepSize);
				break;
			}

#ifdef HET_EARLY_EXIT
			if (accumulatedDensity > stopValue) // Stop early
				return std::numeric_limits<Float>::infinity();
#endif
			m = 6 - m;
		}
		accumulatedDensity *= stepSize/3;

		return accumulatedDensity * m_densityMultiplier;
	}

	/**
	 * \brief Integrate the density covered by the specified ray
	 * segment, while clipping the volume to the underlying
	 * stencil volume.
	 */
	inline Float integrateDensity(const Ray &r) const {
		Ray ray(r(r.mint), r.d,
			0, std::numeric_limits<Float>::infinity(), 0.0f);
		Float integral = 0, remaining = r.maxt - r.mint;
		int iterations = 0;
		bool inside = false;
		Float t = 0;

		while (remaining > 0 && rayIntersect(ray, t, inside)) {
			if (inside) 
				integral += integrateDensity(ray, std::min(remaining, t));

			remaining -= t;
			ray.o = ray(t);
			ray.mint = Epsilon;

			if (++iterations > 100) {
				/// Just a precaution..
				Log(EWarn, "integrateDensity(): round-off error issues?");
				break;
			}
		}
		return integral;
	}

	/**
	 * Attempts to solve the following 1D integral equation for 't'
	 * \int_0^t density(ray.o + x * ray.d) * dx + accumulatedDensity == desiredDensity.
	 * When no solution can be found in [0, maxDist] the function returns
	 * false. For convenience, the function returns the current values of sigmaT 
	 * and the albedo, as well as the 3D position 'ray.o+t*ray.d' upon
	 * success.
	 */
	bool invertDensityIntegral(Ray ray, Float maxDist, 
			Float &accumulatedDensity, Float desiredDensity,
			Float &currentSigmaT, Spectrum &currentAlbedo, 
			Point &currentPoint) const {
		if (maxDist == 0)
			return std::numeric_limits<Float>::infinity();
		int nParts = (int) std::ceil(maxDist/m_stepSize);
		Float stepSize = maxDist/nParts;
		Vector fullIncrement = ray.d * stepSize,
			   halfIncrement = fullIncrement * .5f;

		Float node1 = m_densities->lookupFloat(ray.o);
		Float t = 0;
		int numRefines = 0;
		bool success = false;

		while (t <= maxDist) {
			Float node2 = m_densities->lookupFloat(ray.o + halfIncrement),
				  node3 = m_densities->lookupFloat(ray.o + fullIncrement);
			const Float newDensity = accumulatedDensity 
				+ (node1+node2*4+node3) * (stepSize/6) * m_densityMultiplier;

			if (newDensity > desiredDensity) {
				if (t+stepSize/2 <= maxDist) {
					/* Record the last "good" scattering event */
					success = true;
					if (EXPECT_TAKEN(node2 != 0)) {
						currentPoint = ray.o + halfIncrement;
						currentSigmaT = node2;
					} else if (node3 != 0 && t+stepSize <= maxDist) {
						currentPoint = ray.o + fullIncrement;
						currentSigmaT = node3;
					} else if (node1 != 0) {
						currentPoint = ray.o;
						currentSigmaT = node1;
					}
				}

				if (++numRefines > NUM_REFINES)
					break;

				stepSize *= .5f;
				fullIncrement = halfIncrement;
				halfIncrement *= .5f;
				continue;
			}
			
			if (ray.o+fullIncrement == ray.o) /* Unable to make forward progress - stop */
				break;

			accumulatedDensity = newDensity;
			node1 = node3;
			t += stepSize;
			ray.o += fullIncrement;
		}

		if (success) 
			currentAlbedo = m_albedo->lookupSpectrum(currentPoint);

		return success;
	}

	/// Evaluate the attenuation along a ray segment
	Spectrum tau(const Ray &ray) const {
		return Spectrum(std::exp(-integrateDensity(ray)));
	}

	bool sampleDistance(const Ray &r, MediumSamplingRecord &mRec,
			Sampler *sampler) const {
		Ray ray(r(r.mint), r.d, 0.0f);
		Float distSurf       = r.maxt - r.mint,
			  desiredDensity     = -std::log(1-sampler->next1D()),
			  accumulatedDensity = 0.0f,
			  currentSigmaT  = 0.0f;
		Spectrum currentAlbedo(0.0f);
		int iterations = 0;
		bool inside = false, success = false;
		Float t = 0;
		Float sigmaTOrigin = m_densities->lookupFloat(ray.o) * m_densityMultiplier;

		while (rayIntersect(ray, t, inside)) {
			if (inside) {
				Point currentPoint(0.0f);
				success = invertDensityIntegral(ray, std::min(t, distSurf), 
						accumulatedDensity, desiredDensity, currentSigmaT, currentAlbedo,
						currentPoint);
				if (success) {
					/* A medium interaction occurred */
					mRec.p = currentPoint;
					mRec.t = (mRec.p-r.o).length();
					success = true;
					break;
				}
			}

			distSurf -= t;

			if (distSurf < 0) {
				/* A surface interaction occurred */
				break;
			}
			ray.o = ray(t);
			ray.mint = Epsilon;

			if (++iterations > 100) {
				/// Just a precaution..
				Log(EWarn, "sampleDistance(): round-off error issues?");
				break;
			}
		}
		Float expVal = std::max(Epsilon, std::exp(-accumulatedDensity));

		mRec.pdfFailure = expVal;
		mRec.pdfSuccess = expVal * std::max((Float) Epsilon, currentSigmaT);
		mRec.pdfSuccessRev = expVal * std::max((Float) Epsilon, sigmaTOrigin);
		mRec.attenuation = Spectrum(expVal);
		
		if (!success)
			return false;

		mRec.sigmaS = currentAlbedo * currentSigmaT;
		mRec.sigmaA = Spectrum(currentSigmaT) - mRec.sigmaS;
		mRec.albedo = currentAlbedo.max();
		mRec.orientation = m_orientations != NULL 
			? m_orientations->lookupVector(mRec.p) : Vector(0.0f);
		mRec.medium = this;
		return true;
	}

	void pdfDistance(const Ray &ray, Float t, MediumSamplingRecord &mRec) const {
		Float expVal = std::exp(-integrateDensity(Ray(ray, 0, t)));

		mRec.attenuation = Spectrum(expVal);
		mRec.pdfFailure = expVal;
		mRec.pdfSuccess = expVal * std::max((Float) Epsilon, 
			m_densities->lookupFloat(ray(t)) * m_densityMultiplier);
		mRec.pdfSuccessRev = expVal * std::max((Float) Epsilon, 
			m_densities->lookupFloat(ray(ray.mint)) * m_densityMultiplier);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HeterogeneousMedium[" << endl
			<< "  albedo=" << indent(m_albedo.toString()) << "," << endl
			<< "  orientations=" << indent(m_orientations.toString()) << "," << endl
			<< "  densities=" << indent(m_densities.toString()) << "," << endl
			<< "  stepSize=" << m_stepSize << "," << endl
			<< "  densityMultiplier=" << m_densityMultiplier << endl
			<< "]";
		return oss.str();
	}
	MTS_DECLARE_CLASS()
private:
	ref<VolumeDataSource> m_densities;
	ref<VolumeDataSource> m_albedo;
	ref<VolumeDataSource> m_orientations;
	ref<ShapeKDTree> m_kdTree;
	std::vector<Shape *> m_shapes;
	Float m_stepSize;
};

MTS_IMPLEMENT_CLASS_S(HeterogeneousMedium, false, Medium)
MTS_EXPORT_PLUGIN(HeterogeneousMedium, "Heterogeneous medium");
MTS_NAMESPACE_END

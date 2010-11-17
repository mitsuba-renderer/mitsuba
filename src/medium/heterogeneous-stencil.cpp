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
#define NUM_REFINES 10

class HeterogeneousStencilMedium : public Medium {
public:
	HeterogeneousStencilMedium(const Properties &props) 
		: Medium(props) {
		m_kdTree = new KDTree();
		m_stepSize = props.getFloat("stepSize", 0);
	}

	/* Unserialize from a binary data stream */
	HeterogeneousStencilMedium(Stream *stream, InstanceManager *manager) 
		: Medium(stream, manager) {
		m_kdTree = new KDTree();
		size_t shapeCount = stream->readUInt();
		for (size_t i=0; i<shapeCount; ++i) 
			addChild("", static_cast<Shape *>(manager->getInstance(stream)));
		m_densities = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_albedo = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_orientations = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_stepSize = stream->readFloat();
		configure();
	}

	virtual ~HeterogeneousStencilMedium() {
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
		if (m_orientations.get() == NULL)
			Log(EError, "No orientations specified!");

		if (m_shapes.size() != 0) {
			m_kdTree->build();
			m_aabb = m_kdTree->getAABB();
		} else {
			m_aabb = m_densities->getAABB();
		}

		if (m_stepSize == 0) {
			m_stepSize = std::min(std::min(
				m_densities->getStepSize(),
				m_albedo->getStepSize()),
				m_orientations->getStepSize());

			if (m_stepSize == std::numeric_limits<Float>::infinity()) 
				Log(EError, "Unable to infer step size, please specify!");
		}
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(Shape::m_theClass)) {
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
		} else if (child->getClass()->derivesFrom(VolumeDataSource::m_theClass)) {
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

	Float distanceToMediumEntry(const Ray &ray) const {
		if (m_shapes.size() != 0) {
			Ray r(ray, Epsilon, std::numeric_limits<Float>::infinity());
			Intersection its;
			if (!m_kdTree->rayIntersect(r, its)) 
				return std::numeric_limits<Float>::infinity(); 
			return dot(ray.d, its.geoFrame.n) > 0 ? 0 : its.t;
		} else {
			Float mint, maxt;
			if (!m_aabb.rayIntersect(ray, mint, maxt))
				return std::numeric_limits<Float>::infinity(); 
			if (mint <= Epsilon && maxt <= Epsilon)
				return std::numeric_limits<Float>::infinity(); 
			else
				return (mint < 0) ? 0 : mint;
		}
	}

	Float distanceToMediumExit(const Ray &ray) const {
		if (m_shapes.size() != 0) {
			Ray r(ray, Epsilon, std::numeric_limits<Float>::infinity());
			Intersection its;
			if (!m_kdTree->rayIntersect(r, its)) 
				return 0;
			return dot(ray.d, its.geoFrame.n) < 0 ? 0 : its.t;
		} else {
			Float mint, maxt;
			if (!m_aabb.rayIntersect(ray, mint, maxt)) 
				return 0;
			if (mint < Epsilon && maxt > 0)
				return maxt;
			else
				return 0;
		}
	}

	Spectrum tau(const Ray &r) const {
		Float dLength = r.d.length();
		Ray ray(r(r.mint), r.d / dLength, 0.0f);
		Float remaining = (r.maxt - r.mint) * dLength;
		Float integral = 0.0f;
		int iterations = 0;

		Float entry = distanceToMediumEntry(ray);
		while (entry < remaining && entry != std::numeric_limits<Float>::infinity()) {
			ray.o = ray(entry);
			remaining -= entry;
			Float exit = distanceToMediumExit(ray);

			if (exit != std::numeric_limits<Float>::infinity())
				integral += integrateDensities(ray, std::min(remaining, exit));

			remaining -= exit;
			if (remaining <= 0) 
				break;

			ray.o = ray(exit);
			entry = distanceToMediumEntry(ray);
			if (++iterations > 10) {
				/// Just a precaution..
				Log(EWarn, "tau(): round-off error issues?");
				break;
			}
		}

		return Spectrum(integral * m_sizeMultiplier);
	}
	
	Float integrateDensities(Ray ray, Float length) const {
		int nParts = (int) std::ceil(length/m_stepSize);
		nParts += nParts % 2;
		const Float stepSize = length/nParts;
		const Vector increment = ray.d * stepSize;

		Float m=4;
		Float accumulatedTau = 
			m_densities->lookupFloat(ray.o) + m_densities->lookupFloat(ray(length));

		ray.o += increment;
		for (int i=1; i<nParts; ++i) {
			Float value = m_densities->lookupFloat(ray.o);
			accumulatedTau += value*m;
			ray.o += increment;
			m = 6-m;
		}
		accumulatedTau *= stepSize/3;

		return accumulatedTau;
	}

	Float invertDensityIntegral(Ray ray, Float maxDist, 
		Float &accumulatedTau, 
		Float &currentSigmaT, 
		Spectrum &currentAlbedo, 
		Float desiredTau) const {

		int nParts = (int) std::ceil(maxDist/m_stepSize);
		Float stepSize = maxDist/nParts;
		Vector fullIncrement = ray.d * stepSize,
			   halfIncrement = fullIncrement * .5f;
	
		Float lastNode = m_densities->lookupFloat(ray.o);
		Float t = 0;

		int numRefines = 0;
		while (t <= maxDist) {
			const Float node1 = lastNode,
						node2 = m_densities->lookupFloat(ray.o + halfIncrement),
						node3 = m_densities->lookupFloat(ray.o + fullIncrement);
			const Float newDensity = accumulatedTau 
				+ (node1+node2*4+node3) * (stepSize/6);

			if (newDensity > desiredTau) {
				if (++numRefines > NUM_REFINES) {
					currentAlbedo = m_albedo->lookupSpectrum(ray.o+fullIncrement);
					currentSigmaT = node3;
					return t;
				}

				stepSize *= .5f;
				fullIncrement *= .5f;
				halfIncrement *= .5f;
				continue;
			}

			accumulatedTau = newDensity;
			lastNode = node3;
			t += stepSize;
			ray.o += fullIncrement;
		}
		return std::numeric_limits<Float>::infinity();
	}

	bool sampleDistance(const Ray &r, Float maxDist, 
			MediumSamplingRecord &mRec,  Sampler *sampler) const {
		Float dLength = r.d.length();
		Ray ray(r(r.mint), r.d / dLength, 0.0f);
		Float remaining      = (maxDist - r.mint) * dLength,
			  desiredTau     = -std::log(1-sampler->next1D())/m_sizeMultiplier,
			  accumulatedTau = 0.0f,
			  currentSigmaT  = 0.0f;
		Spectrum currentAlbedo(0.0f), integral(0.0f);
		int iterations = 0;
		bool success = false;

		mRec.miWeight = 1;
		Float entry = distanceToMediumEntry(ray);
		while (entry < remaining && entry != std::numeric_limits<Float>::infinity()) {
			ray.o = ray(entry);
			remaining -= entry;
			Float exit = distanceToMediumExit(ray);

			if (exit != std::numeric_limits<Float>::infinity()) {
				Float t = invertDensityIntegral(ray, std::min(remaining, exit), 
						accumulatedTau, currentSigmaT, currentAlbedo, desiredTau);
				if (t != std::numeric_limits<Float>::infinity()) {
					mRec.p = ray(t);
					mRec.t = (mRec.p-r.o).length();
					success = true;
					break;
				}
			} else {
				break;
			}

			remaining -= exit;

			if (remaining < 0) {
				break;
			}

			ray.o = ray(exit);
			entry = distanceToMediumEntry(ray);
			if (++iterations > 10) {
				/// Just a precaution..
				Log(EWarn, "sampleDistance(): round-off error issues?");
				break;
			}
		}
		accumulatedTau *= m_sizeMultiplier;
		currentSigmaT *= m_sizeMultiplier;

		if (!success) {
			/* Could not achieve the desired density - hit the next surface instead */
			mRec.pdf = std::max((Float) Epsilon, std::exp(-accumulatedTau));
			mRec.attenuation = (-Spectrum(accumulatedTau)).exp();
			return false;
		}

		mRec.pdf = std::exp(-accumulatedTau) * std::max((Float) Epsilon, currentSigmaT);
		mRec.attenuation = (-Spectrum(accumulatedTau)).exp();

		mRec.sigmaS = currentAlbedo * currentSigmaT;
		mRec.sigmaA = Spectrum(currentSigmaT) - mRec.sigmaS;
		mRec.albedo = currentAlbedo.max();
		mRec.orientation = m_orientations->lookupVector(mRec.p);
		mRec.medium = this;
		return true;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HeterogeneousStencilMedium[" << endl
			<< "  albedo=" << indent(m_albedo.toString()) << endl
			<< "  orientations=" << indent(m_orientations.toString()) << endl
			<< "  densities=" << indent(m_densities.toString()) 
			<< "]";
		return oss.str();
	}
	MTS_DECLARE_CLASS()
private:
	ref<VolumeDataSource> m_densities;
	ref<VolumeDataSource> m_albedo;
	ref<VolumeDataSource> m_orientations;
	ref<KDTree> m_kdTree;
	std::vector<Shape *> m_shapes;
	Float m_stepSize;
};

MTS_IMPLEMENT_CLASS_S(HeterogeneousStencilMedium, false, Medium)
MTS_EXPORT_PLUGIN(HeterogeneousStencilMedium, "Stenciled heterogeneous medium");
MTS_NAMESPACE_END

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
#include <mitsuba/render/mipmap3d.h>
#include <fstream>

MTS_NAMESPACE_BEGIN

#define TEST_RESULT(i) if (i != 1) Log(EError, "File format error!");

class AdaptiveMedium : public Medium {
public:
	AdaptiveMedium(const Properties &props) : Medium(props) {
		std::string volData = props.getString("filename");
		volData = FileResolver::getInstance()->resolve(volData);

		m_mediumToWorld = props.getTransform("toWorld", Transform());
		m_worldToMedium = m_mediumToWorld.inverse();
		m_scaleFactor = m_mediumToWorld(Vector(1,0,0)).length();

		Float proposedSigma = std::numeric_limits<Float>::infinity();
		for (int i=0; i<SPECTRUM_SAMPLES; ++i)
			proposedSigma = std::min(proposedSigma, m_sigmaT[i]);
		
		/* Can be used to override the extinction coefficient used to sample distances 
			in the in-scatter line integral. By default, the smallest spectral sample of 
			<tt>sigmaA+sigmaT</tt> is used. */
		m_sigma = props.getFloat("sigma", proposedSigma);


		Log(EInfo, "Loading medium data from \"%s\" ..", volData.c_str());
		FILE *f = fopen(volData.c_str(), "r");
		if (f == NULL)
			Log(EError, "Invalid medium data file specified");

		Vector3i oldRes; AABB aabb;
		/* [Fedkiw et al.] volume density format */
		TEST_RESULT(fscanf(f, "%i", &oldRes.x));
		TEST_RESULT(fscanf(f, "%i", &oldRes.y));
		TEST_RESULT(fscanf(f, "%i", &oldRes.z));
		TEST_RESULT(fscanf(f, "%f", &aabb.min.x));
		TEST_RESULT(fscanf(f, "%f", &aabb.min.y));
		TEST_RESULT(fscanf(f, "%f", &aabb.min.z));
		TEST_RESULT(fscanf(f, "%f", &aabb.max.x));
		TEST_RESULT(fscanf(f, "%f", &aabb.max.y));
		TEST_RESULT(fscanf(f, "%f", &aabb.max.z));

		int res = roundToPowerOfTwo(std::max(std::max(oldRes.x, oldRes.y), oldRes.z));
		size_t nElems = res*res*res;
		float *data = new float[nElems];
		int slice = res*res;
		memset(data, 0, nElems*sizeof(float));

		for (int z=0; z<oldRes.x; ++z) {
			for (int y=0; y<oldRes.x; ++y) {
				for (int x=0; x<oldRes.x; ++x) {
					Float tmp;
					TEST_RESULT(fscanf(f, "%f", &tmp));
					data[x+y*res+z*slice] = std::max((Float) 0, tmp);
				}
			}
		}
		fclose(f);

		Float maxError = props.getFloat("maxError", .1f);
		Float fudgeFactor = props.getFloat("fudgeFactor", .01f); // Let's be honest..

		m_aabb.reset();
		for (int i=0; i<8; ++i)
			m_aabb.expandBy(m_mediumToWorld(aabb.getCorner(i)));

		m_mipmap = new SparseMipmap3D(aabb, res, data, maxError, fudgeFactor);

		delete[] data;
	}

	/* Unserialize from a binary data stream */
	AdaptiveMedium(Stream *stream, InstanceManager *manager) 
		: Medium(stream, manager) {
		m_mediumToWorld = Transform(stream);
		m_sigma = stream->readFloat();
		m_mipmap = static_cast<SparseMipmap3D *>(manager->getInstance(stream));
		m_worldToMedium = m_mediumToWorld.inverse();
		m_scaleFactor = m_mediumToWorld(Vector(1,0,0)).length();
		configure();
	}

	/* Serialize to a binary data stream */
	void serialize(Stream *stream, InstanceManager *manager) const {
		Medium::serialize(stream, manager);

		m_mediumToWorld.serialize(stream);
		stream->writeFloat(m_sigma);
		manager->serialize(stream, m_mipmap.get());
	}

	Spectrum tau(const Ray &r) const {
		Ray ray;
		m_worldToMedium(r, ray);
		return m_sigmaT * (m_mipmap->lineIntegral(ray) * m_scaleFactor);
	}

	bool sampleDistance(const Ray &r, Float maxDist, 
			MediumSamplingRecord &mRec,  Sampler *sampler) const {
		Ray ray;
		m_worldToMedium(r, ray);

		Float t, accumulatedDensity, sampleDensity,
			desiredDensity = -std::log(1-sampler->next1D())/(m_sigma * m_scaleFactor);

		mRec.miWeight = 1;
		if (!m_mipmap->invertLineIntegral(ray, desiredDensity, 
			accumulatedDensity, t, sampleDensity)) {
			accumulatedDensity *= m_scaleFactor;

			mRec.pdf = std::exp(-accumulatedDensity * m_sigma);
			mRec.attenuation = (m_sigmaT * (-accumulatedDensity)).exp();
			return false;
		} else {
			accumulatedDensity *= m_scaleFactor;

			mRec.t = t * m_scaleFactor;
			mRec.p = r(mRec.t);
			mRec.pdf = std::exp(-accumulatedDensity * m_sigma) * sampleDensity * m_sigma;
			mRec.attenuation = (m_sigmaT * (-accumulatedDensity)).exp();
			mRec.sigmaA = m_sigmaA * sampleDensity;
			mRec.sigmaS = m_sigmaS * sampleDensity;
			mRec.albedo = m_albedo;
			mRec.medium = this;
			return true;
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "AdaptiveMedium[" << endl
			<< "  sigmaA = " << m_sigmaA.toString() << "," << std::endl
			<< "  sigmaS = " << m_sigmaS.toString() << "," << std::endl
			<< "  sigmaT = " << m_sigmaT.toString() << "," << std::endl
			<< "  phase = " << indent(m_phaseFunction->toString()) << "," << std::endl
			<< "  mipmap = " << indent(m_mipmap->toString()) << std::endl
			<< "]";
		return oss.str();
	}
	MTS_DECLARE_CLASS()
private:
	Transform m_worldToMedium;
	Transform m_mediumToWorld;
	Float m_scaleFactor, m_sigma;
	ref<SparseMipmap3D> m_mipmap;
};

MTS_IMPLEMENT_CLASS_S(AdaptiveMedium, false, Medium)
MTS_EXPORT_PLUGIN(AdaptiveMedium, "Adaptive heterogeneous medium");
MTS_NAMESPACE_END

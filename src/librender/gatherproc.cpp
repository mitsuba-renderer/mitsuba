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

#include <mitsuba/render/gatherproc.h>

MTS_NAMESPACE_BEGIN

/**
 * This stores a number of photons, which can be sent over 
 * the wire as needed. Used to implement parallel photon
 * tracing passes.
 */
class PhotonVector : public WorkResult {
public:
	PhotonVector() { }

	inline void put(const PhotonMap::Photon &p) {
		m_photons.push_back(p);
	}

	inline size_t size() const {
		return m_photons.size();
	}

	inline void clear() {
		m_photons.clear();
	}

	inline const PhotonMap::Photon &operator[](size_t index) const {
		return m_photons[index];
	}

	void load(Stream *stream) {
		m_photons.clear();
		size_t count = stream->readUInt();
		m_photons.resize(count);
		for (size_t i=0; i<count; ++i)
			m_photons[i] = PhotonMap::Photon(stream);
	}

	void save(Stream *stream) const {
		stream->writeUInt((unsigned int) m_photons.size());
		for (size_t i=0; i<m_photons.size(); ++i)
			m_photons[i].serialize(stream);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "PhotonVector[size=" << m_photons.size() << "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
	// Virtual destructor
	virtual ~PhotonVector() { }
private:
	std::vector<PhotonMap::Photon> m_photons;
};

/**
 * This class does the actual photon tracing work
 */
class GatherPhotonWorker : public ParticleTracer {
public:
	GatherPhotonWorker(GatherPhotonProcess::EGatherType type, unsigned int granularity,
		int maxDepth, int rrDepth)
		: ParticleTracer(maxDepth, true, rrDepth),
		m_type(type), m_granularity(granularity) {
	}

	GatherPhotonWorker(Stream *stream, InstanceManager *manager) 
	 : ParticleTracer(stream, manager) {
		m_type = (GatherPhotonProcess::EGatherType) stream->readInt();
		m_granularity = stream->readUInt();
	}

	ref<WorkProcessor> clone() const {
		return new GatherPhotonWorker(m_type, m_granularity, m_maxDepth,
			m_rrDepth);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		ParticleTracer::serialize(stream, manager);
		stream->writeInt(m_type);
		stream->writeUInt(m_granularity);
	}

	ref<WorkResult> createWorkResult() const {
		return new PhotonVector();
	}

	void process(const WorkUnit *workUnit, WorkResult *workResult, 
		const bool &stop) {
		m_workResult = static_cast<PhotonVector *>(workResult);
		m_workResult->clear();
		ParticleTracer::process(workUnit, workResult, stop);
		m_workResult = NULL;
	}

	void handleSurfaceInteraction(int depth, bool caustic, const Intersection &its, const Spectrum &weight) {
		int bsdfType = its.shape->getBSDF()->getType();
		if (!(bsdfType & BSDF::EDiffuseReflection) && !(bsdfType & BSDF::EGlossyReflection))
			return;

		if ((m_type == GatherPhotonProcess::ECausticPhotons && depth > 1 && caustic)
		 || (m_type == GatherPhotonProcess::ESurfacePhotons && depth > 1 && !caustic)
		 || (m_type == GatherPhotonProcess::EAllSurfacePhotons)) 
			m_workResult->put(PhotonMap::Photon(its.p, its.geoFrame.n, -its.toWorld(its.wi), weight, depth));
	}

	void handleMediumInteraction(int depth, bool caustic, const MediumSamplingRecord &mRec, Float time, 
			const Vector &wi, const Spectrum &weight) {
		if (m_type == GatherPhotonProcess::EVolumePhotons && depth > 1)
			m_workResult->put(PhotonMap::Photon(mRec.p, Normal(), -wi, weight, depth));
	}

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GatherPhotonWorker() { }
protected:
	GatherPhotonProcess::EGatherType m_type;
	unsigned int m_granularity;
	ref<PhotonVector> m_workResult;
};

GatherPhotonProcess::GatherPhotonProcess(EGatherType type, size_t photonCount, 
	unsigned int granularity, int maxDepth, int rrDepth, const void *progressReporterPayload) 
	: ParticleProcess(ParticleProcess::EGather, photonCount, granularity, "Gathering photons", 
	  progressReporterPayload), m_type(type), m_maxDepth(maxDepth), m_rrDepth(rrDepth), 
	  m_excess(0), m_numShot(0) {
	m_photonMap = new PhotonMap(photonCount);
}

ref<WorkProcessor> GatherPhotonProcess::createWorkProcessor() const {
	return new GatherPhotonWorker(m_type, (unsigned int) m_granularity, m_maxDepth, 
		m_rrDepth);
}

void GatherPhotonProcess::processResult(const WorkResult *wr, bool cancelled) {
	if (cancelled)
		return;
	const PhotonVector &vec = *static_cast<const PhotonVector *>(wr);
	m_resultMutex->lock();
	increaseResultCount(vec.size());
	for (size_t i=0; i<vec.size(); ++i) {
		if (!m_photonMap->storePhoton(vec[i])) {
			m_excess += vec.size() - i;
			break;
		}
	}
	m_numShot += m_granularity;
	m_resultMutex->unlock();
}


MTS_IMPLEMENT_CLASS(GatherPhotonProcess, false, ParticleProcess) 
MTS_IMPLEMENT_CLASS_S(GatherPhotonWorker, false, ParticleTracer) 
MTS_IMPLEMENT_CLASS(PhotonVector, false, WorkResult) 
MTS_NAMESPACE_END

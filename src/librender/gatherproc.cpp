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

#include <mitsuba/render/gatherproc.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief This work result implementation stores a sequence of photons, which can be
 * sent over the wire as needed.
 *
 * It is used to implement parallel networked photon tracing passes.
 */
class PhotonVector : public WorkResult {
public:
	PhotonVector() { }

	inline void nextParticle() {
		m_particleIndices.push_back((uint32_t) m_photons.size());
	}

	inline void put(const Photon &p) {
		m_photons.push_back(p);
	}

	inline size_t size() const {
		return m_photons.size();
	}

	inline size_t getParticleCount() const {
		return m_particleIndices.size()-1;
	}

	inline size_t getParticleIndex(size_t idx) const {
		return m_particleIndices.at(idx);
	}

	inline void clear() {
		m_photons.clear();
		m_particleIndices.clear();
	}

	inline const Photon &operator[](size_t index) const {
		return m_photons[index];
	}

	void load(Stream *stream) {
		clear();
		size_t count = (size_t) stream->readUInt();
		m_particleIndices.resize(count);
		stream->readUIntArray(&m_particleIndices[0], count);
		count = (size_t) stream->readUInt();
		m_photons.resize(count);
		for (size_t i=0; i<count; ++i)
			m_photons[i] = Photon(stream);
	}

	void save(Stream *stream) const {
		stream->writeUInt((uint32_t) m_particleIndices.size());
		stream->writeUIntArray(&m_particleIndices[0], m_particleIndices.size());
		stream->writeUInt((uint32_t) m_photons.size());
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
	std::vector<Photon> m_photons;
	std::vector<uint32_t> m_particleIndices;
};

/**
 * This class does the actual photon tracing work
 */
class GatherPhotonWorker : public ParticleTracer {
public:
	GatherPhotonWorker(GatherPhotonProcess::EGatherType type, size_t granularity,
		int maxDepth, int rrDepth) : ParticleTracer(maxDepth, rrDepth, false),
		m_type(type), m_granularity(granularity) { }

	GatherPhotonWorker(Stream *stream, InstanceManager *manager)
	 : ParticleTracer(stream, manager) {
		m_type = (GatherPhotonProcess::EGatherType) stream->readInt();
		m_granularity = stream->readSize();
	}

	ref<WorkProcessor> clone() const {
		return new GatherPhotonWorker(m_type, m_granularity, m_maxDepth,
			m_rrDepth);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		ParticleTracer::serialize(stream, manager);
		stream->writeInt(m_type);
		stream->writeSize(m_granularity);
	}

	ref<WorkResult> createWorkResult() const {
		return new PhotonVector();
	}

	void process(const WorkUnit *workUnit, WorkResult *workResult,
		const bool &stop) {
		m_workResult = static_cast<PhotonVector *>(workResult);
		m_workResult->clear();
		ParticleTracer::process(workUnit, workResult, stop);
		m_workResult->nextParticle();
		m_workResult = NULL;
	}

	void handleNewParticle() {
		m_workResult->nextParticle();
	}

	void handleSurfaceInteraction(int depth_, int nullInteractions, bool delta,
			const Intersection &its, const Medium *medium,
			const Spectrum &weight) {
		int bsdfType = its.getBSDF()->getType(), depth = depth_ - nullInteractions;
		if (!(bsdfType & BSDF::EDiffuseReflection) && !(bsdfType & BSDF::EGlossyReflection))
			return;

		if ((m_type == GatherPhotonProcess::ECausticPhotons && depth > 1 && delta)
		 || (m_type == GatherPhotonProcess::ESurfacePhotons && depth > 1 && !delta)
		 || (m_type == GatherPhotonProcess::EAllSurfacePhotons))
			m_workResult->put(Photon(its.p, its.geoFrame.n, -its.toWorld(its.wi), weight, depth));
	}

	void handleMediumInteraction(int depth, int nullInteractions, bool delta,
			const MediumSamplingRecord &mRec, const Medium *medium,
			const Vector &wi, const Spectrum &weight) {
		if (m_type == GatherPhotonProcess::EVolumePhotons)
			m_workResult->put(Photon(mRec.p, Normal(0.0f, 0.0f, 0.0f),
				-wi, weight, depth-nullInteractions));
	}

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GatherPhotonWorker() { }
protected:
	GatherPhotonProcess::EGatherType m_type;
	size_t m_granularity;
	ref<PhotonVector> m_workResult;
};

GatherPhotonProcess::GatherPhotonProcess(EGatherType type, size_t photonCount,
	size_t granularity, int maxDepth, int rrDepth, bool isLocal, bool autoCancel,
	const void *progressReporterPayload)
	: ParticleProcess(ParticleProcess::EGather, photonCount, granularity, "Gathering photons",
	  progressReporterPayload), m_type(type), m_photonCount(photonCount), m_maxDepth(maxDepth),
	  m_rrDepth(rrDepth),  m_isLocal(isLocal), m_autoCancel(autoCancel), m_excess(0), m_numShot(0) {
	m_photonMap = new PhotonMap(photonCount);
}

bool GatherPhotonProcess::isLocal() const {
	return m_isLocal;
}

ref<WorkProcessor> GatherPhotonProcess::createWorkProcessor() const {
	return new GatherPhotonWorker(m_type, m_granularity, m_maxDepth, m_rrDepth);
}

void GatherPhotonProcess::processResult(const WorkResult *wr, bool cancelled) {
	if (cancelled)
		return;
	const PhotonVector &vec = *static_cast<const PhotonVector *>(wr);
	LockGuard lock(m_resultMutex);

	size_t nParticles = 0;
	for (size_t i=0; i<vec.getParticleCount(); ++i) {
		size_t start = vec.getParticleIndex(i),
			   end   = vec.getParticleIndex(i+1);
		++nParticles;
		bool full = false;
		for (size_t j=start; j<end; ++j) {
			if (!m_photonMap->tryAppend(vec[j])) {
				m_excess += vec.size() - j;
				full = true;
				break;
			}
		}
		if (full)
			break;
	}
	m_numShot += nParticles;
	increaseResultCount(vec.size());
}

ParallelProcess::EStatus GatherPhotonProcess::generateWork(WorkUnit *unit, int worker) {
	/* Use the same approach as PBRT for auto canceling */
	LockGuard lock(m_resultMutex);
	if (m_autoCancel && m_numShot > 100000
			&& unsuccessful(m_photonCount, m_photonMap->size(), m_numShot)) {
		Log(EInfo, "Not enough photons could be collected, giving up");
		return EFailure;
	}

	return ParticleProcess::generateWork(unit, worker);
}

MTS_IMPLEMENT_CLASS(GatherPhotonProcess, false, ParticleProcess)
MTS_IMPLEMENT_CLASS_S(GatherPhotonWorker, false, ParticleTracer)
MTS_IMPLEMENT_CLASS(PhotonVector, false, WorkResult)
MTS_NAMESPACE_END

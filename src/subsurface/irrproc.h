/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#if !defined(__DIPOLE_PROC_H)
#define __DIPOLE_PROC_H

#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/**
 * Abstract storage class for representing irradiance samples on
 * translucent surfaces
 */
class IrradianceSample {
public:
	/// Default (empty) constructor
	inline IrradianceSample() { }

	/// Unserialize an irradiance sample from a binary data stream
	inline IrradianceSample(Stream *stream) {
		p = Point(stream);
		E = Spectrum(stream);
		area = stream->readFloat();
	}

	/**
	 * \param p The sample point on the surface
	 * \param E The irradiance value at this point
	 * \param area The surface area associated with that sample
	 */
	inline IrradianceSample(const Point &p, const Spectrum &E, Float area) 
		: p(p), E(E), area(area) { }

	/// Serialize an irradiance sample to a binary data stream
	inline void serialize(Stream *stream) const {
		p.serialize(stream);
		E.serialize(stream);
		stream->writeFloat(area);
	}

	Point p;
	Spectrum E;
	Float area;
};


/**
 * This stores a number of irradiance samples, which can be sent
 * over the wire as needed. Used to implement parallel irradiance
 * sampling for the dipole BSSRDF.
 */
class IrradianceRecordVector : public WorkResult {
public:
	IrradianceRecordVector() { }

	inline void put(const IrradianceSample &rec) {
		m_samples.push_back(rec);
	}

	inline size_t size() const {
		return m_samples.size();
	}

	inline void clear() {
		m_samples.clear();
	}

	inline const IrradianceSample &operator[](size_t index) const {
		return m_samples[index];
	}

	/* WorkUnit interface */
	void load(Stream *stream);
	void save(Stream *stream) const;
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	// Virtual destructor
	virtual ~IrradianceRecordVector() { }
private:
	std::vector<IrradianceSample> m_samples;
};

/**
 * Parallel process for performing distributed irradiance sampling
 */
class IrradianceSamplingProcess : public ParallelProcess {
public:
	IrradianceSamplingProcess(size_t sampleCount, size_t granularity, 
		int ssIndex, int irrSamples, 
		bool irrIndirect, const void *progressReporterPayload);

	inline const IrradianceRecordVector *getSamples() const {
		return m_samples.get();
	}

	/* ParallelProcess implementation */
	ref<WorkProcessor> createWorkProcessor() const; 
	void processResult(const WorkResult *wr, bool cancelled);
	ParallelProcess::EStatus generateWork(WorkUnit *unit, int worker);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~IrradianceSamplingProcess();
private:
	size_t m_sampleCount, m_samplesRequested, m_resultCount, m_granularity;
	int m_ssIndex;
	int m_irrSamples;
	bool m_irrIndirect;
	ref<Mutex> m_resultMutex;
	ref<IrradianceRecordVector> m_samples;
	ProgressReporter *m_progress;
};

MTS_NAMESPACE_END

#endif /* __DIPOLE_PROC_H */

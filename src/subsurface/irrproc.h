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

#if !defined(__DIPOLE_PROC_H)
#define __DIPOLE_PROC_H

#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/**
 * Data structure representing position samples for
 * which irradiance values should be computed
 */
class PositionSample {
public:
    /// Default (empty) constructor
    inline PositionSample() { }

    /// Unserialize an Position sample from a binary data stream
    inline PositionSample(Stream *stream) {
        p = Point(stream);
        n = Normal(stream);
        shapeIndex = stream->readInt();
    }
    inline PositionSample(const Point &p, const Normal &n, int shapeIndex)
        : p(p), n(n), shapeIndex(shapeIndex) { }

    /// Serialize an position sample to a binary data stream
    inline void serialize(Stream *stream) const {
        p.serialize(stream);
        n.serialize(stream);
        stream->writeInt(shapeIndex);
    }

    Point p;
    Normal n;
    int shapeIndex;
};

/**
 * \brief This class stores a number of position samples, which can be sent
 * over the wire as needed.
 *
 * Used to implement parallel irradiance sampling for the dipole BSSRDF.
 */
class PositionSampleVector : public WorkUnit {
public:
    PositionSampleVector() { }

    inline void put(const PositionSample &rec) {
        m_samples.push_back(rec);
    }

    inline size_t size() const {
        return m_samples.size();
    }

    inline void clear() {
        m_samples.clear();
    }

    inline std::vector<PositionSample> &get() {
        return m_samples;
    }

    inline void reserve(size_t size) {
        m_samples.reserve(size);
    }

    inline const PositionSample &operator[](size_t index) const {
        return m_samples[index];
    }

    /* WorkUnit interface */
    void load(Stream *stream);
    void save(Stream *stream) const;
    void set(const WorkUnit *workUnit);
    std::string toString() const;

    MTS_DECLARE_CLASS()
protected:
    // Virtual destructor
    virtual ~PositionSampleVector() { }
private:
    std::vector<PositionSample> m_samples;
};

/**
 * Data structure for representing irradiance samples on
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
     * \param E The irradiance value at this position
     */
    inline IrradianceSample(const Point &p, const Spectrum &E)
        : p(p), E(E) { }

    /// Serialize an irradiance sample to a binary data stream
    inline void serialize(Stream *stream) const {
        p.serialize(stream);
        E.serialize(stream);
        stream->writeFloat(area);
    }

    /// Return the position (used by the octree code)
    inline const Point &getPosition() const {
        return p;
    }

    Point p;
    Spectrum E;
    Float area;    //!< total surface area represented by this sample
    uint8_t label; //!< used by the octree construction code
};

/**
 * \brief This class stores a number of irradiance samples, which can be sent
 * over the wire as needed.
 *
 * Used to implement parallel irradiance sampling for the dipole BSSRDF.
 */
class IrradianceSampleVector : public WorkResult {
public:
    IrradianceSampleVector() { }

    inline void put(const IrradianceSample &rec) {
        m_samples.push_back(rec);
    }

    inline size_t size() const {
        return m_samples.size();
    }

    inline void clear() {
        m_samples.clear();
    }

    inline std::vector<IrradianceSample> &get() {
        return m_samples;
    }

    inline void reserve(size_t size) {
        m_samples.reserve(size);
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
    virtual ~IrradianceSampleVector() { }
private:
    std::vector<IrradianceSample> m_samples;
};

/**
 * Parallel process for performing distributed irradiance sampling
 */
class IrradianceSamplingProcess : public ParallelProcess {
public:
    IrradianceSamplingProcess(PositionSampleVector *positions,
        size_t granularity, int irrSamples, bool irrIndirect,
        Float time, const void *data);

    inline IrradianceSampleVector *getIrradianceSampleVector() {
        return m_irradianceSamples.get();
    }

    inline PositionSampleVector *getPositionSampleVector() {
        return m_positionSamples.get();
    }

    inline const AABB &getAABB() const {
        return m_aabb;
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
    ref<PositionSampleVector> m_positionSamples;
    ref<IrradianceSampleVector> m_irradianceSamples;
    size_t m_samplesRequested, m_granularity;
    int m_irrSamples;
    bool m_irrIndirect;
    Float m_time;
    ref<Mutex> m_resultMutex;
    ProgressReporter *m_progress;
    AABB m_aabb;
};

MTS_NAMESPACE_END

#endif /* __DIPOLE_PROC_H */

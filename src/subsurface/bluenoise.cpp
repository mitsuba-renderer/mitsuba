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

#include "bluenoise.h"
#include <mitsuba/core/statistics.h>
#include <boost/unordered_map.hpp>

#if defined(MTS_OPENMP)
# include <omp.h>
#endif

#if defined(__GNUC__) && defined(MTS_OPENMP) && __GNUC__ >= 4 && __GNUC_MINOR__ >= 3
# define MTS_PARALLEL_SORT 1
# include <parallel/algorithm>
#else
# define MTS_PARALLEL_SORT 0
#endif

MTS_NAMESPACE_BEGIN

/// Used to store a white noise position sample along with its cell ID
struct UniformSample {
    Point p;
    Normal n;
    int64_t cellID;
    int shapeIndex;

    inline UniformSample() { }
    inline UniformSample(const Point &p, const Normal &n, int cellID, int shapeIndex)
        : p(p), n(n), cellID(cellID), shapeIndex(shapeIndex) { }
};

/// Functor for sorting 'UniformSample' instances with respect to their cell ID
struct CellIDOrdering {
    inline bool operator()(const UniformSample &s1, const UniformSample &s2) const {
        return s1.cellID < s2.cellID;
    }
};

/// Stores the first UniformSample that falls in this cell and the chosen one (if any)
struct Cell {
    int firstIndex;
    int sample;

    inline Cell() { }
    inline Cell(int firstIndex, int sample = -1)
        : firstIndex(firstIndex), sample(sample) { }
};

void blueNoisePointSet(const Scene *scene, const std::vector<Shape *> &shapes,
        Float radius, PositionSampleVector *target, Float &sa, AABB &aabb,
        const void *data) {
    int kmax = 8; /* Perform 8 trial runs */

    #if defined(MTS_OPENMP)
        int nproc = mts_omp_get_max_threads();
    #else
        int nproc = 1;
    #endif
    ProgressReporter rep("Generating sample positions", 27*kmax+5, data);

    /* Create a random number generator for each thread */
    ref<Random> rootRNG = new Random();
    ref_vector<Random> t_rng(nproc);
    std::vector<AABB> t_aabb(nproc);
    for (int i=0; i<nproc; ++i) {
        t_rng[i] = new Random(rootRNG);
        t_aabb[i].reset();
    }

    DiscreteDistribution areaDistr;
    std::vector<int> shapeMap(shapes.size());
    for (size_t i=0; i<shapes.size(); ++i) {
        shapeMap[i] = -1;
        for (size_t j=0; j<scene->getShapes().size(); ++j) {
            if (scene->getShapes()[j].get() == shapes[i]) {
                shapeMap[i] = (int) j;
                break;
            }
        }
        SAssert(shapeMap[i] != -1);
        areaDistr.append(shapes[i]->getSurfaceArea());
    }
    areaDistr.normalize();
    sa = areaDistr.getSum();

    /* Start with a fairly dense set of white noise points */
    int nsamples = math::ceilToInt(15 * sa / (M_PI * radius * radius));

    SLog(EInfo, "Creating a blue noise point set (radius=%f, "
        "surface area = %f)", radius, sa);
    SLog(EInfo, "  phase 1: creating dense white noise (%i samples)", nsamples);
    ref<Timer> timer = new Timer();
    std::vector<UniformSample> samples(nsamples);
    rep.update(0);

    #if defined(MTS_OPENMP)
        #pragma omp parallel for schedule(static)
    #endif
    for (int i=0; i<nsamples; ++i) {
        #if defined(MTS_OPENMP)
            int tid = mts_omp_get_thread_num();
        #else
            int tid = 0;
        #endif

        Random *random = t_rng[tid].get();
        Point2 sample(random->nextFloat(), random->nextFloat());
        int shapeIndex = (int) areaDistr.sampleReuse(sample.x);
        Shape *shape = shapes[shapeIndex];

        PositionSamplingRecord pRec(0);
        shape->samplePosition(pRec, sample);

        samples[i] = UniformSample(pRec.p, pRec.n, 0, shapeMap[shapeIndex]);
        t_aabb[tid].expandBy(pRec.p);
    }
    SLog(EInfo, "    done (took %i ms, %s)" , timer->getMilliseconds(),
        memString(sizeof(PositionSample) * samples.size()).c_str());
    rep.update(1);
    timer->reset();

    Float cellWidth = radius / std::sqrt(3.0f),
          invCellWidth = 1.0f / cellWidth;

    aabb.reset();
    for (int i=0; i<nproc; ++i)
        aabb.expandBy(t_aabb[i]);
    Vector extents = aabb.getExtents();

    Vector3i cellCount;
    for (int i=0; i<3; ++i)
        cellCount[i] = std::max(1, math::ceilToInt(extents[i] * invCellWidth));

    SLog(EInfo, "  phase 2: computing cell indices ..");
    #if defined(MTS_OPENMP)
        #pragma omp parallel for schedule(static)
    #endif
    for (int i=0; i<nsamples; ++i) {
        Vector rel = samples[i].p - aabb.min;
        Vector3i idx;
        for (int j=0; j<3; ++j)
            idx[j] = std::min((int) (rel[j] * invCellWidth), cellCount[j]-1);
        samples[i].cellID = idx[0] + (int64_t) cellCount[0] *
            (idx[1] + idx[2] * (int64_t) cellCount[1]);
    }
    SLog(EInfo, "    done (took %i ms)" , timer->getMilliseconds());
    rep.update(2);
    timer->reset();

    SLog(EInfo, "  phase 3: sorting ..");
    #if MTS_PARALLEL_SORT
        __gnu_parallel::sort(samples.begin(), samples.end(), CellIDOrdering());
    #else
        std::sort(samples.begin(), samples.end(), CellIDOrdering());
    #endif
    SLog(EInfo, "    done (took %i ms)" , timer->getMilliseconds());
    rep.update(3);
    timer->reset();

    SLog(EInfo, "  phase 4: establishing valid cells and phase groups ..");
    typedef boost::unordered_map<int64_t, Cell> CellMap;

    CellMap cells(samples.size());
    std::vector<std::vector<int64_t> > phaseGroups(27);
    for (int i=0; i<27; ++i)
        phaseGroups[i].reserve(cells.size() / 27);

    int64_t last = -1;
    for (int i=0; i<nsamples; ++i) {
        int64_t id = samples[i].cellID;
        if (id != last) {
            cells[id] = Cell(i);
            last = id;

            /* Schedule this cell wrt. the corresponding phase group */
            int64_t tmp = id;
            int64_t z = tmp / (cellCount[0] * cellCount[1]);
            tmp -= z * (cellCount[0] * cellCount[1]);
            int64_t y = tmp / cellCount[0];
            int64_t x = tmp - y * cellCount[0];
            int phaseID = x % 3 + (y % 3) * 3 + (z % 3) * 9;
            phaseGroups[phaseID].push_back(id);
        }
    }

    SLog(EInfo, "    done (took %i ms), got %i cells, avg. samples per cell: %f",
        timer->getMilliseconds(), (int) cells.size(),
        samples.size() / (Float) cells.size());
    rep.update(4);
    timer->reset();

    SLog(EInfo, "  phase 5: parallel sampling ..");
    for (int trial=0; trial<kmax; ++trial) {
        for (int phase=0; phase<27; ++phase) {
            const std::vector<int64_t> &phaseGroup = phaseGroups[phase];

            #if defined(MTS_OPENMP)
                #pragma omp parallel for
            #endif
            for (int i=0; i < (int) phaseGroup.size(); ++i) {
                int64_t cellID = phaseGroup[i];
                Cell &cell = cells[cellID];
                int arrayIndex = cell.firstIndex + trial;

                if (arrayIndex >= (int) samples.size() ||
                    samples[arrayIndex].cellID != cellID ||
                    cell.sample != -1)
                    continue;

                const UniformSample &sample = samples[arrayIndex];

                bool conflict = false;

                for (int z=-2; z<3; ++z) {
                    for (int y=-2; y<3; ++y) {
                        for (int x=-2; x<3; ++x) {
                            int64_t neighborCellID = cellID + x
                                + (int64_t) cellCount[0] * (y + z * (int64_t) cellCount[1]);

                            CellMap::iterator it = cells.find(neighborCellID);

                            if (it != cells.end()) {
                                const Cell &neighbor = it->second;
                                if (neighbor.sample != -1) {
                                    const UniformSample &sample2 = samples[neighbor.sample];

                                    if ((sample.p-sample2.p).lengthSquared() < radius*radius) {
                                        conflict = true;
                                        goto bailout;
                                    }
                                }
                            }
                        }
                    }
                }

            bailout:
                if (!conflict)
                    cell.sample = arrayIndex;
            }
            rep.update(5+trial*27+phase);
        }
    }
    SLog(EInfo, "    done (took %i ms)" , timer->getMilliseconds());
    timer->reset();

    for (CellMap::const_iterator it = cells.begin();
            it != cells.end(); ++it) {
        const Cell &cell = it->second;
        if (cell.sample == -1)
            continue;
        const UniformSample &sample = samples[cell.sample];
        target->put(PositionSample(sample.p, sample.n, sample.shapeIndex));
    }

    SLog(EInfo, "Sampling finished (obtained %i blue noise samples)", (int) target->size());
}

MTS_NAMESPACE_END

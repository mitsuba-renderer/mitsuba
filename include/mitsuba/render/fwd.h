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

#pragma once
#if !defined(__MITSUBA_RENDER_FWD_H_)
#define __MITSUBA_RENDER_FWD_H_

MTS_NAMESPACE_BEGIN

class BlockedImageProcess;
class BlockedRenderProcess;
class BlockListener;
class BSDF;
struct BSDFSamplingRecord;
struct DirectionSamplingRecord;
struct DirectSamplingRecord;
class Emitter;
class Film;
class GatherPhotonProcess;
class HemisphereSampler;
class HWResource;
class ImageBlock;
class Instanced;
class Integrator;
struct Intersection;
class IrradianceCache;
template <typename AABBType> class KDTreeBase;
template <typename AABBType, typename TreeConstructionHeuristic, typename Derived> class GenericKDTree;
template <typename Derived> class SAHKDTree3D;
class ShapeKDTree;
class LocalWorker;
struct LuminaireSamplingRecord;
class Medium;
struct MediumSamplingRecord;
class MIPMap;
class MonteCarloIntegrator;
class ParticleProcess;
class ParticleTracer;
struct PhaseFunctionSamplingRecord;
class PhaseFunction;
class PhotonMap;
struct PositionSamplingRecord;
class PreviewWorker;
class ProjectiveCamera;
struct RadianceQueryRecord;
struct PositionSamplingRecord;
class Random;
class RangeWorkUnit;
class ReconstructionFilter;
class RectangularWorkUnit;
class RenderJob;
class RenderListener;
class RenderQueue;
class SamplingIntegrator;
class Sampler;
class Sensor;
class Scene;
class SceneHandler;
class Shader;
class Shape;
class SparseMipmap3D;
class Spiral;
class Subsurface;
class Texture;
struct TriAccel;
struct TriAccel4;
class TriMesh;
class Utility;
class VolumeDataSource;
struct VPL;

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_FWD_H_ */

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

#if !defined(__RENDER_FWD_H)
#define __RENDER_FWD_H

MTS_NAMESPACE_BEGIN

class AbstractAnimationTrack;
template <typename T> class AnimationTrack;
typedef AnimationTrack<Float> FloatTrack;
typedef AnimationTrack<Quaternion> QuatTrack;
typedef AnimationTrack<Vector> VectorTrack;
typedef AnimationTrack<Point> PointTrack;
class AnimatedTransform;
class BlockedImageProcess;
class BlockedRenderProcess;
class BlockListener;
class BSDF;
struct BSDFQueryRecord;
class Camera;
struct EmissionRecord;
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
class Luminaire;
struct LuminaireSamplingRecord;
class Medium;
struct MediumSamplingRecord;
class MIPMap;
class MonteCarloIntegrator;
class ParticleProcess;
class ParticleTracer;
struct PhaseFunctionQueryRecord;
class PhaseFunction;
class PhotonMap;
class PreviewWorker;
class ProjectiveCamera;
struct RadianceQueryRecord;
class Random;
class RangeWorkUnit;
class ReconstructionFilter;
class RectangularWorkUnit;
class RenderJob;
class RenderListener;
class RenderQueue;
class SampleIntegrator;
class Sampler;
class Scene;
class SceneHandler;
class Shader;
class Shape;
struct ShapeSamplingRecord;
class SparseMipmap3D;
class Spiral;
class Subsurface;
class TabulatedFilter;
class Texture;
struct TriAccel;
struct TriAccel4;
class TriMesh;
class Utility;
class VolumeDataSource;
struct VPL;

MTS_NAMESPACE_END

#endif /* __RENDER_FWD_H */

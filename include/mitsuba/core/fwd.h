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

#if !defined(__CORE_FWD_H)
#define __CORE_FWD_H

MTS_NAMESPACE_BEGIN

struct AABB;
class Appender;
class Bitmap;
class BlackBodySpectrum;
struct BSphere;
class ConfigurableObject;
struct CacheLineCounter;
class Class;
class ConditionVariable;
class ConsoleStream;
class DefaultFormatter;
struct DiscretePDF;
class FileResolver;
class FileStream;
class Formatter;
struct Frame;
class InstanceManager;
class InterpolatedSpectrum;
class LocalWorker;
class Logger;
struct Matrix4x4;
class MemoryStream;
class Mutex;
class NetworkedObject;
struct Normal;
class Object;
class ParallelProcess;
class Plugin;
class PluginManager;
class ProgressReporter;
class Properties;
class Random;
struct Ray;
struct RayDifferential;
class RemoteProcess;
class RemoteWorker;
class RemoteWorkerReader;
class Scheduler;
class Serializable;
class SerializableObject;
struct SHRotation;
class SHSampler;
struct SHVector;
struct SHVector4D;
class SmoothSpectrum;
class SocketStream;
class SparseWavelet2D;
class SparseWaveletOctree;
struct Spectrum;
class SSHStream;
class Statistics;
class StatsCounter;
class Stream;
class StreamAppender;
class StreamBackend;
class Thread;
class Timer;
struct Transform;
struct Triangle;
class UnbufferedAppender;
template <typename T> struct TVector2;
template <typename T> struct TVector3;
template <typename T> struct TVector4;
template <typename T> struct TPoint2;
template <typename T> struct TPoint3;
template <typename T> struct TPoint4;
template <typename T> struct TQuaternion;
typedef TVector2<Float>    Vector2;
typedef TVector2<int>      Vector2i;
typedef TVector2<float>    Vector2f;
typedef TVector2<double>   Vector2d;
typedef TVector3<Float>    Vector;
typedef TVector3<Float>    Vector3;
typedef TVector3<int>      Vector3i;
typedef TVector3<float>    Vector3f;
typedef TVector3<double>   Vector3d;
typedef TVector4<Float>    Vector4;
typedef TVector4<int>      Vector4i;
typedef TVector4<float>    Vector4f;
typedef TVector4<double>   Vector4d;
typedef TPoint2<Float>     Point2;
typedef TPoint2<int>       Point2i;
typedef TPoint2<float>     Point2f;
typedef TPoint2<double>    Point2d;
typedef TPoint3<Float>     Point;
typedef TPoint3<Float>     Point3;
typedef TPoint3<int>       Point3i;
typedef TPoint3<float>     Point3f;
typedef TPoint3<double>    Point3d;
typedef TPoint4<Float>     Point4;
typedef TPoint4<int>       Point4i;
typedef TPoint4<float>     Point4f;
typedef TPoint4<double>    Point4d;
typedef TQuaternion<Float> Quaternion;
class WaitFlag;
class Wavelet2D;
class Wavelet3D;
class Worker;
class WorkProcessor;
class WorkResult;
class WorkUnit;
class ZStream;

MTS_NAMESPACE_END

#endif /* __CORE_FWD_H */

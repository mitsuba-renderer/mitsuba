#define BOOST_PYTHON_USE_GCC_SYMBOL_VISIBILITY 1
#include "base.h"
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/appender.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/aabb.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/sched_remote.h>
#include <mitsuba/core/netobject.h>
#include <mitsuba/core/sstream.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/cstream.h>
#include <mitsuba/core/qmc.h>
#include <mitsuba/core/vmf.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/sshstream.h>
#include <mitsuba/render/scenehandler.h>
#include <mitsuba/render/scene.h>
#include <boost/algorithm/string.hpp>
#include <boost/python/tuple.hpp>

using namespace mitsuba;

bool check_python_exception() {
    AcquireGIL gil;

	PyObject *exception_ = NULL, *value_ = NULL, *traceback_ = NULL;
	PyErr_Fetch(&exception_, &value_, &traceback_);

	if (!exception_)
		return false;

	PyErr_NormalizeException(&exception_, &value_, &traceback_);

	bp::object exception(bp::handle<>(bp::allow_null(exception_)));
	bp::object value(bp::handle<>(bp::allow_null(value_)));
	bp::object traceback(bp::handle<>(bp::allow_null(traceback_)));

	PyErr_Clear();

	bp::object traceback_package(bp::import("traceback"));
	bp::object format_exception = traceback_package.attr("format_exception");
	bp::object formatted_list = format_exception(exception, value, traceback);
	bp::object formatted = bp::str("\n").join(formatted_list);

	SLog(EWarn, "Caught a Python exception: %s",
		bp::extract<std::string>(formatted)().c_str());

	return true;
}

static void initializeFramework() {
	/* Initialize the core framework */
	Class::staticInitialization();
	Object::staticInitialization();
	PluginManager::staticInitialization();
	Statistics::staticInitialization();
	FileStream::staticInitialization();
	Thread::staticInitialization();
	Logger::staticInitialization();
	Spectrum::staticInitialization();
	Bitmap::staticInitialization();
	Scheduler::staticInitialization();
	SHVector::staticInitialization();
	SceneHandler::staticInitialization();
	Thread::registerCrashHandler(&check_python_exception);
}

static void shutdownFramework() {
	/* Shutdown the core framework */
	SceneHandler::staticShutdown();
	SHVector::staticShutdown();
	Scheduler::staticShutdown();
	Bitmap::staticShutdown();
	Spectrum::staticShutdown();
	Logger::staticShutdown();
	Thread::staticShutdown();
	FileStream::staticShutdown();
	Statistics::staticShutdown();
	PluginManager::staticShutdown();
	Object::staticShutdown();
	Class::staticShutdown();
}

StringVector *StringVector_fromList(bp::list list) {
	StringVector *result = new StringVector(bp::len(list));
	for (int i=0; i<bp::len(list); ++i)
		result->operator[](i) = bp::extract<std::string>(list[i]);
	return result;
}

template <typename SpectrumType> class SpectrumWrapper {
public:
	static Float get(const SpectrumType &spec, int i) {
		if (i < 0 || i >= SpectrumType::dim) {
			SLog(EError, "Index %i is out of range!", i);
			return 0.0f;
		}
		return spec[i];
	}

	static void set(SpectrumType &spec, int i, Float value) {
		if (i < 0 || i >= SpectrumType::dim)
			SLog(EError, "Index %i is out of range!", i);
		else
			spec[i] = value;
	}

	static int len(SpectrumType &) {
		return SpectrumType::dim;
	}
};

class properties_wrapper {
public:
	static bp::object get(const Properties &props, std::string name) {
		if (!props.hasProperty(name))
			SLog(EError, "Properties: keyword \"%s\" not found!", name.c_str());

		switch (props.getType(name)) {
			case Properties::EBoolean:
				return bp::object(props.getBoolean(name));
			case Properties::EString:
				return bp::object(props.getString(name));
			case Properties::EInteger:
				return bp::object(props.getInteger(name));
			case Properties::EFloat:
				return bp::object(props.getFloat(name));
			case Properties::EVector:
				return bp::object(props.getVector(name));
			case Properties::EPoint:
				return bp::object(props.getPoint(name));
			case Properties::ETransform:
				return bp::object(props.getTransform(name));
			case Properties::EAnimatedTransform:
				return bp::object(props.getAnimatedTransform(name));
			case Properties::ESpectrum:
				return bp::object(props.getSpectrum(name));
			default:
				SLog(EError, "Properties: type of keyword \"%s\" is not supported!", name.c_str());
				return bp::object();
		}
	}

	static void set(Properties &props, const std::string name, bp::object value) {
		bp::extract<std::string> extractString(value);
		bp::extract<bool> extractBoolean(value);
		bp::extract<int> extractInteger(value);
		bp::extract<Float> extractFloat(value);
		bp::extract<Vector> extractVector(value);
		bp::extract<Point> extractPoint(value);
		bp::extract<Transform> extractTransform(value);
		bp::extract<Spectrum> extractSpectrum(value);
		bp::extract<AnimatedTransform *> extractAnimatedTransform(value);

		if (extractString.check()){
			props.setString(name, extractString(), false);
		} else if (extractBoolean.check() && PyObject_IsInstance(value.ptr(), (PyObject *) &PyBool_Type)) {
			props.setBoolean(name, extractBoolean(), false);
		} else if (extractInteger.check()) {
			props.setInteger(name, extractInteger(), false);
		} else if (extractFloat.check()) {
			props.setFloat(name, extractFloat(), false);
		} else if (extractPoint.check()) {
			props.setPoint(name, extractPoint(), false);
		} else if (extractVector.check()) {
			props.setVector(name, extractVector(), false);
		} else if (extractTransform.check()) {
			props.setTransform(name, extractTransform(), false);
		} else if (extractAnimatedTransform.check()) {
			props.setAnimatedTransform(name, extractAnimatedTransform(), false);
		} else if (extractSpectrum.check()) {
			props.setSpectrum(name, extractSpectrum(), false);
		} else {
			SLog(EError, "Properties: type of keyword \"%s\" is not supported!", name.c_str());
		}
	}
};

struct path_to_python_str {
	static PyObject* convert(fs::path const& path) {
		return bp::incref(
			bp::object(path.string()).ptr());
	}
};

struct TSpectrum_to_Spectrum {
	static PyObject* convert(const TSpectrum<Float, SPECTRUM_SAMPLES> &spectrum) {
		return bp::incref(bp::object(Spectrum(spectrum)).ptr());
	}
};

static void scheduler_cancel(Scheduler *scheduler, ParallelProcess *proc) {
	ReleaseGIL gil;
	scheduler->cancel(proc);
}

static void scheduler_wait(Scheduler *scheduler, const ParallelProcess *proc) {
	ReleaseGIL gil;
	scheduler->wait(proc);
}

static Matrix4x4 *Matrix4x4_fromList(bp::list list) {
	if (bp::len(list) == 4) {
		Float buf[4][4];
		for (int i=0; i<4; ++i) {
			bp::list subList = bp::extract<bp::list>(list[i]);
			if (bp::len(subList) != 4)
				SLog(EError, "Matrix4x4 list constructor: invalid argument");
			for (int j=0; j<4; ++j)
				buf[i][j] = bp::extract<Float>(subList[j]);
		}
		return new Matrix4x4(buf);
	} else if (bp::len(list) == 16) {
		Float buf[16];
		for (int i=0; i<16; ++i)
			buf[i] = bp::extract<Float>(list[i]);
		return new Matrix4x4(buf);
	} else {
		SLog(EError, "Matrix4x4 list constructor: invalid argument");
		return NULL;
	}
}

static void Matrix4x4_setItem(Matrix4x4 *matrix, bp::tuple tuple, Float value) {
	if (bp::len(tuple) != 2)
		SLog(EError, "Invalid matrix indexing operation, required a tuple of length 2");
	int i = bp::extract<int>(tuple[0]);
	int j = bp::extract<int>(tuple[1]);

	if (i < 0 || j < 0 || i >= 4 || j >= 4)
		SLog(EError, "Index (%i, %i) is out of bounds!", i, j);

	matrix->operator()(i, j) = value;
}

static Float Matrix4x4_getItem(Matrix4x4 *matrix, bp::tuple tuple) {
	if (bp::len(tuple) != 2)
		SLog(EError, "Invalid matrix indexing operation, required a tuple of length 2");
	int i = bp::extract<int>(tuple[0]);
	int j = bp::extract<int>(tuple[1]);

	if (i < 0 || j < 0 || i >= 4 || j >= 4)
		SLog(EError, "Index (%i, %i) is out of bounds!", i, j);

	return matrix->operator()(i, j);
}

static Matrix4x4 Matrix4x4_transpose(Matrix4x4 *matrix) {
	Matrix4x4 result;
	matrix->transpose(result);
	return result;
}

static Matrix4x4 Matrix4x4_invert(Matrix4x4 *matrix) {
	Matrix4x4 result;
	matrix->invert(result);
	return result;
}

static bp::tuple Matrix4x4_symEig(Matrix4x4 *matrix) {
	Matrix4x4 Q;
	Float d[4];
	matrix->symEig(Q, d);

	bp::list list;
	for (int i=0; i<4; ++i)
		list.append(d[i]);

	return bp::make_tuple(Q, list);
}

static bp::tuple Matrix4x4_lu(Matrix4x4 *matrix) {
	Matrix4x4 LU;
	int piv[4];
	int pivsign;
	matrix->lu(LU, piv, pivsign);

	bp::list list;
	for (int i=0; i<4; ++i)
		list.append(piv[i]);

	return bp::make_tuple(LU, list, pivsign);
}

static Vector4 Matrix4x4_cholSolve(Matrix4x4 *matrix, Vector B) {
	typedef Matrix<4, 1, Float> Matrix4x1;
	Vector4 X;
	Matrix4x1 & aliasedB MTS_MAY_ALIAS = reinterpret_cast<Matrix4x1 &>(B);
	Matrix4x1 & aliasedX MTS_MAY_ALIAS = reinterpret_cast<Matrix4x1 &>(X);
	matrix->cholSolve<1>(aliasedB, aliasedX);
	return X;
}

static Vector4 Matrix4x4_luSolve(Matrix4x4 *matrix, Vector B, bp::list pivList) {
	typedef Matrix<4, 1, Float> Matrix4x1;
	Vector4 X;
	int piv[4];

	if (bp::len(pivList) != 4)
		SLog(EError, "Matrix4x4 list constructor: invalid argument");
	for (int i=0; i<4; ++i)
		piv[i] = bp::extract<Float>(pivList[i]);

	Matrix4x1 & aliasedB MTS_MAY_ALIAS = reinterpret_cast<Matrix4x1 &>(B);
	Matrix4x1 & aliasedX MTS_MAY_ALIAS = reinterpret_cast<Matrix4x1 &>(X);
	matrix->luSolve<1>(aliasedB, aliasedX, piv);
	return X;
}

static void SHVector_setItem(SHVector *v, bp::tuple tuple, Float value) {
	if (bp::len(tuple) != 2)
		SLog(EError, "Invalid v indexing operation, required a tuple of length 2");
	int i = bp::extract<int>(tuple[0]);
	int j = bp::extract<int>(tuple[1]);

	if (i < 0 || i >= v->getBands() || i < -i || j > i)
		SLog(EError, "Index (%i, %i) is out of bounds!", i, j);

	v->operator()(i, j) = value;
}

static Float SHVector_getItem(SHVector *v, bp::tuple tuple) {
	if (bp::len(tuple) != 2)
		SLog(EError, "Invalid v indexing operation, required a tuple of length 2");
	int i = bp::extract<int>(tuple[0]);
	int j = bp::extract<int>(tuple[1]);

	if (i < 0 || i >= v->getBands() || i < -i || j > i)
		SLog(EError, "Index (%i, %i) is out of bounds!", i, j);

	return v->operator()(i, j);
}

static Class *object_getClass(Object *object) {
	return const_cast<Class *>(object->getClass());
}

static Class *class_forName(const char *name) {
	return const_cast<Class *>(Class::forName(name));
}

static Class *class_getSuperClass(Class *theClass) {
	return const_cast<Class *>(theClass->getSuperClass());
}

static ref<SerializableObject> instance_manager_getinstance(InstanceManager *manager, Stream *stream) {
	return manager->getInstance(stream);
}

void appender_logProgress(Appender *appender, Float progress, const std::string &name,
	const std::string &formatted, const std::string &eta) {
	appender->logProgress(progress, name, formatted, eta, NULL);
}

static void logger_logProgress(Logger *logger, Float progress, const std::string &name,
	const std::string &formatted, const std::string &eta) {
	logger->logProgress(progress, name, formatted, eta, NULL);
}

static void mts_log(ELogLevel level, const std::string &msg) {
	bp::object traceback(bp::import("traceback"));
	bp::object extract_stack(traceback.attr("extract_stack"));
	bp::object stack = extract_stack();
	bp::object top(stack[bp::len(stack)-1]);
	std::string module = bp::extract<std::string>(top[2]);
	Thread::getThread()->getLogger()->log(level,
		NULL, bp::extract<const char *>(top[0]),
		bp::extract<int>(top[1]), "%s%s: %s",
		module.c_str(), module[0] == '<' ? "" : "()",
		msg.c_str());
}

class FormatterWrapper : public Formatter {
public:
	FormatterWrapper(PyObject *self) : m_self(self), m_locked(false) { Py_INCREF(m_self); }

	std::string format(ELogLevel logLevel, const Class *theClass,
			const Thread *thread, const std::string &text,
			const char *file, int line) {
        if (m_locked)
            return "";
        AcquireGIL gil;
        TrivialScopedLock lock(m_locked);
		return bp::call_method<std::string>(m_self, "format", logLevel,
				bp::ptr(const_cast<Class *>(theClass)),
				bp::ptr(const_cast<Thread *>(thread)), text, file, line);
	}

	virtual ~FormatterWrapper() {
		Py_DECREF(m_self);
	}
private:
	PyObject *m_self;
    bool m_locked;
};

class AppenderWrapper : public Appender {
public:
	AppenderWrapper(PyObject *self) : m_self(self), m_locked(false) { Py_INCREF(m_self); }

	void append(ELogLevel level, const std::string &text) {
        CALLBACK_SYNC_GIL();
		bp::call_method<void>(m_self, "append", level, text);
	}

	void logProgress(Float progress, const std::string &name,
		const std::string &formatted, const std::string &eta,
		const void *ptr) {
        CALLBACK_SYNC_GIL();
		bp::call_method<void>(m_self, "logProgress", progress, name, formatted, eta);
	}

	virtual ~AppenderWrapper() {
		Py_DECREF(m_self);
	}
private:
	PyObject *m_self;
    bool m_locked;
};

static Spectrum *spectrum_array_constructor(bp::list list) {
	Float spec[SPECTRUM_SAMPLES];
	if (bp::len(list) != SPECTRUM_SAMPLES)
		SLog(EError, "Spectrum: expected %i arguments", SPECTRUM_SAMPLES);

	for (int i=0; i<bp::len(list); ++i)
		spec[i] = bp::extract<Float>(list[i]);

	return new Spectrum(spec);
}

static Point ray_eval(const Ray &ray, Float t) {
	return ray(t);
}

static bp::tuple spectrum_toLinearRGB(const Spectrum &s) {
	Float r, g, b;
	s.toLinearRGB(r, g, b);
	return bp::make_tuple(r, g, b);
}

static bp::tuple spectrum_toSRGB(const Spectrum &s) {
	Float r, g, b;
	s.toSRGB(r, g, b);
	return bp::make_tuple(r, g, b);
}

static bp::tuple spectrum_toXYZ(const Spectrum &s) {
	Float x, y, z;
	s.toXYZ(x, y, z);
	return bp::make_tuple(x, y, z);
}

static bp::tuple spectrum_toIPT(const Spectrum &s) {
	Float I, P, T;
	s.toIPT(I, P, T);
	return bp::make_tuple(I, P, T);
}

static bp::object bsphere_rayIntersect(BSphere *bsphere, const Ray &ray) {
	Float nearT, farT;
	if (bsphere->rayIntersect(ray, nearT, farT))
		return bp::make_tuple(nearT, farT);
	else
		return bp::object();

}

static bp::object aabb_rayIntersect(AABB *aabb, const Ray &ray) {
	Float nearT, farT;
	if (aabb->rayIntersect(ray, nearT, farT))
		return bp::make_tuple(nearT, farT);
	else
		return bp::object();

}

static bp::object logger_readLog(Logger *logger) {
	std::string string;
	if (logger->readLog(string))
		return bp::object(string);
	else
		return bp::object();
}

static bp::object aabb_rayIntersect2(AABB *aabb, const Ray &ray, Float nearT, Float farT) {
	Point nearP, farP;
	if (aabb->rayIntersect(ray, nearT, farT, nearP, farP))
		return bp::make_tuple(nearT, farT, nearP, farP);
	else
		return bp::object();
}

static void thread_join(Thread *thread) {
	ReleaseGIL gil;
	thread->join();
}

static Vector transform_mul_vector(Transform *transform, const Vector &vector) { return transform->operator()(vector); }
static Vector4 transform_mul_vector4(Transform *transform, const Vector4 &vector) { return transform->operator()(vector); }
static Normal transform_mul_normal(Transform *transform, const Normal &normal) { return transform->operator()(normal); }
static Point transform_mul_point(Transform *transform, const Point &point) { return transform->operator()(point); }
static Ray transform_mul_ray(Transform *transform, const Ray &ray) { return transform->operator()(ray); }
static Transform transform_mul_transform(Transform *transform, const Transform &other) { return *transform * other; }

bp::object cast(ConfigurableObject *obj) {
	if (obj == NULL)
		return bp::object();
	const Class *cls = obj->getClass();
	#define TryCast(ClassName) if (cls->derivesFrom(MTS_CLASS(ClassName))) \
		return bp::object(ref<ClassName>(static_cast<ClassName *>(obj)))
	TryCast(BSDF);
	TryCast(TriMesh);
	TryCast(Shape);
	TryCast(PhaseFunction);
	TryCast(Integrator);
	TryCast(Texture);
	TryCast(Medium);
	TryCast(VolumeDataSource);
	TryCast(Film);
	TryCast(PerspectiveCamera);
	TryCast(ProjectiveCamera);
	TryCast(Sensor);
	TryCast(Emitter);
	TryCast(Sampler);
	TryCast(ReconstructionFilter);
	TryCast(ConfigurableObject);
	#undef TryCast
	SLog(EError, "Internal error in cast()!");
	return bp::object();
}

static bp::object pluginmgr_createobject_1(PluginManager *mgr, const Properties &props) {
	return cast(mgr->createObject(props));
}

static bp::object pluginmgr_createobject_2(PluginManager *mgr, const Class *cls, const Properties &props) {
	return cast(mgr->createObject(cls, props));
}

static ref<ConfigurableObject> pluginmgr_create_helper(PluginManager *manager, bp::dict dict, std::map<std::string, ConfigurableObject *> &objs) {
	Properties properties;
	bp::list items = dict.items();
	std::map<std::string, ConfigurableObject *> children;

	for (bp::stl_input_iterator<bp::tuple> it(items), end; it!=end; ++it) {
		bp::tuple tuple = *it;
		std::string name = bp::extract<std::string>(tuple[0]);
		bp::extract<bp::dict> extractDict(tuple[1]);
		bp::extract<std::string> extractString(tuple[1]);
		bp::extract<ConfigurableObject *> extractConfigurableObject(tuple[1]);

		if (name == "type") {
			if (!extractString.check())
				SLog(EError, "'type' property must map to a string!");
			else
				properties.setPluginName(extractString());
		} else if (name == "id") {
			if (!extractString.check())
				SLog(EError, "'id' property must map to a string!");
			else
				properties.setID(extractString());
		} else if (extractDict.check()) {
			ref<ConfigurableObject> obj = pluginmgr_create_helper(manager, extractDict(), objs);
			children[name] = obj;
			obj->incRef();
		} else if (extractConfigurableObject.check()) {
			ref<ConfigurableObject> obj = extractConfigurableObject();
			children[name] = obj;
			obj->incRef();
		} else {
			properties_wrapper::set(properties, name, tuple[1]);
		}
	}

	ref<ConfigurableObject> object;
	if (properties.getPluginName() == "scene") {
		object = new Scene(properties);
	} else if (properties.getPluginName() == "ref") {
		std::string id = properties.getID();
		if (id == "unnamed")
			SLog(EError, "id parameter of reference is missing!");
		if (objs.find(id) == objs.end())
			SLog(EError, "Could not find referenced object with id \"%s\"", id.c_str());
		for (std::map<std::string, ConfigurableObject *>::iterator it = children.begin();
			it != children.end(); ++it)
			it->second->decRef();
		return objs[id];
	} else {
		object = manager->createObject(properties);
	}

	if (properties.getID() != "unnamed") {
		if (objs.find(properties.getID()) != objs.end())
			SLog(EError, "Duplicate ID '%s'", properties.getID().c_str());
		objs[properties.getID()] = object;
		object->incRef();
	}

	for (std::map<std::string, ConfigurableObject *>::iterator it = children.begin();
		it != children.end(); ++it) {
		object->addChild(it->first, it->second);
		it->second->setParent(object);
		it->second->decRef();
	}

	object->configure();
	return object;
}

static bp::object pluginmgr_create(PluginManager *manager, bp::dict dict) {
	std::map<std::string, ConfigurableObject *> objs;
	ref<ConfigurableObject> result = pluginmgr_create_helper(manager, dict, objs);
	for (std::map<std::string, ConfigurableObject *>::iterator it = objs.begin();
			it != objs.end(); ++it)
		it->second->decRef();
	return cast(result);
}

static bp::tuple mkCoordinateSystem(const Vector &n) {
	Vector s, t;
	coordinateSystem(n, s, t);

	return bp::make_tuple(s, t);
}

static bp::tuple fresnelDielectricExt1(Float cosThetaI, Float eta) {
	Float cosThetaT;
	Float result = fresnelDielectricExt(cosThetaI, cosThetaT, eta);

	return bp::make_tuple(result, cosThetaT);
}

static Float fresnelDielectricExt2(Float cosThetaI, Float eta) {
	return fresnelDielectricExt(cosThetaI, eta);
}

static Vector refract1(const Vector &wi, const Normal &n, Float eta, Float cosThetaT) {
	return refract(wi, n, eta, cosThetaT);
}

static bp::tuple refract2(const Vector &wi, const Normal &n, Float eta) {
	Float cosThetaT, F;
	Vector result = refract(wi, n, eta, cosThetaT, F);
	return bp::make_tuple(result, cosThetaT, F);
}

static void bitmap_applyMatrix(Bitmap *bitmap, bp::list list) {
	int length = bp::len(list);
	if (length != 9)
		SLog(EError, "Require a color matrix specified as a list with 9 entries!");

	Float matrix[3][3];

	int idx = 0;
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			matrix[i][j] = bp::extract<Float>(list[idx++]);

	bitmap->applyMatrix(matrix);
}

static void bitmap_write1(Bitmap *bitmap, Bitmap::EFileFormat fmt, Stream *stream) {
	bitmap->write(fmt, stream);
}

static void bitmap_write2(Bitmap *bitmap, Bitmap::EFileFormat fmt, Stream *stream, int compression) {
	bitmap->write(fmt, stream, compression);
}

static void bitmap_write3(Bitmap *bitmap, Bitmap::EFileFormat fmt, const fs::path &path) {
	bitmap->write(fmt, path);
}

static void bitmap_write4(Bitmap *bitmap, Bitmap::EFileFormat fmt, const fs::path &path, int compression) {
	bitmap->write(fmt, path, compression);
}

static void bitmap_write5(Bitmap *bitmap, const fs::path &path) {
	bitmap->write(path);
}

static void bitmap_write6(Bitmap *bitmap, const fs::path &path, int compression) {
	bitmap->write(path, compression);
}

static void bitmap_convert_0(Bitmap *bitmap, Bitmap *target) {
    bitmap->convert(target);
}

static void bitmap_convert_1(Bitmap *bitmap, Bitmap *target, Float multiplier) {
    bitmap->convert(target, multiplier);
}

static void bitmap_convert_2(Bitmap *bitmap, Bitmap *target, Float multiplier, Spectrum::EConversionIntent intent) {
	bitmap->convert(target, multiplier, intent);
}

static ref<Bitmap> bitmap_convert_3(Bitmap *bitmap, Bitmap::EPixelFormat pixelFormat, Bitmap::EComponentFormat componentFormat,
		Float gamma, Float multiplier, Spectrum::EConversionIntent intent) {
	return bitmap->convert(pixelFormat, componentFormat, gamma, multiplier, intent);
}

static ref<Bitmap> bitmap_convert_4(Bitmap *bitmap, Bitmap::EPixelFormat pixelFormat, Bitmap::EComponentFormat componentFormat,
		Float gamma, Float multiplier) {
	return bitmap->convert(pixelFormat, componentFormat, gamma, multiplier);
}

static ref<Bitmap> bitmap_convert_5(Bitmap *bitmap, Bitmap::EPixelFormat pixelFormat, Bitmap::EComponentFormat componentFormat,
		Float gamma) {
	return bitmap->convert(pixelFormat, componentFormat, gamma);
}

static ref<Bitmap> bitmap_convert_6(Bitmap *bitmap, Bitmap::EPixelFormat pixelFormat, Bitmap::EComponentFormat componentFormat) {
	return bitmap->convert(pixelFormat, componentFormat);
}

static bp::dict bitmap_split(Bitmap *bitmap) {
	std::map<std::string, Bitmap *> map = bitmap->split();
	bp::dict result;
	for (std::map<std::string, Bitmap *>::iterator it = map.begin(); it != map.end(); ++it)
		result[it->first] = bp::object(ref<Bitmap>(it->second));
	return result;
}

static bp::object bitmap_plot(Bitmap *_bitmap) {
	bp::dict locals;
	locals["bitmap"] = ref<Bitmap>(_bitmap);
	bp::exec("import matplotlib.pyplot as plt\n"
			 "handle = plt.imshow(bitmap.buffer())\n", bp::object(), locals);
	return locals["handle"];
}

static ref<Bitmap> bitmap_extractChannels(Bitmap::EPixelFormat fmt, Bitmap *bitmap, bp::list list) {
	std::vector<int> channels(bp::len(list));
	for (int i=0; i<bp::len(list); ++i)
		channels[i] = bp::extract<int>(list[i]);
	return bitmap->extractChannels(fmt, channels);
}

static void bitmap_fromByteArray(Bitmap *bitmap, bp::object obj) {
	if (PyByteArray_Check(obj.ptr())) {
		uint8_t *ptr = (uint8_t *) PyByteArray_AsString(obj.ptr());
		size_t size = PyByteArray_Size(obj.ptr());
		SAssertEx(size == bitmap->getBufferSize(), "Bitmap::fromByteArray(): buffer sizes don't match!");

		memcpy(bitmap->getData(), ptr, size);
	} else {
		SLog(EError, "Bitmap::fromByteArray(): Invalid argument!");
	}
}

static void bitmap_toByteArray_1(const Bitmap *bitmap, bp::object obj) {
	if (PyByteArray_Check(obj.ptr())) {
		uint8_t *ptr = (uint8_t *) PyByteArray_AsString(obj.ptr());
		size_t size = PyByteArray_Size(obj.ptr());
		SAssertEx(size == bitmap->getBufferSize(), "Bitmap::fromByteArray(): buffer sizes don't match!");

		memcpy(ptr, bitmap->getData(), size);
	} else {
		SLog(EError, "Bitmap::toByteArray(): Invalid argument!");
	}
}

static bp::object bitmap_toByteArray_2(const Bitmap *bitmap) {
	return bp::object(bp::handle<>(PyByteArray_FromStringAndSize(
			(char *) bitmap->getUInt8Data(), bitmap->getBufferSize())));
}

static bp::tuple bitmap_tonemapReinhard(Bitmap *bitmap, Float logAvgLuminance, Float maxLuminance, Float key, Float burn) {
	bitmap->tonemapReinhard(logAvgLuminance, maxLuminance, key, burn);
	return bp::make_tuple(logAvgLuminance, maxLuminance);
}

static bp::object bitmap_join(Bitmap::EPixelFormat fmt, bp::list list) {
	std::vector<Bitmap *> bitmaps(bp::len(list));

	for (int i=0; i<bp::len(list); ++i)
		bitmaps[i] = bp::extract<Bitmap *>(list[i]);

	return bp::object(Bitmap::join(fmt, bitmaps));
}

static Transform transform_glOrthographic1(Float clipNear, Float clipFar) {
	return Transform::glOrthographic(clipNear, clipFar);
}

static Transform transform_glOrthographic2(Float clipLeft, Float clipRight,
		Float clipBottom, Float clipTop, Float clipNear, Float clipFar) {
	return Transform::glOrthographic(clipLeft, clipRight,
		clipBottom, clipTop, clipNear, clipFar);
}

static bp::list fileresolver_resolveAll(const FileResolver *fres, const fs::path &path) {
	bp::list result;
	std::vector<fs::path> paths = fres->resolveAll(path);

	for (size_t i=0; i<paths.size(); ++i)
		result.append(paths[i]);

	return result;
}

static bp::tuple DiscreteDistribution_sample(DiscreteDistribution *d, Float sampleValue) {
	Float pdf;
	size_t index = d->sample(sampleValue, pdf);
	return bp::make_tuple(index, pdf);
}

static bp::tuple DiscreteDistribution_sampleReuse(DiscreteDistribution *d, Float sampleValue) {
	Float pdf;
	size_t index = d->sampleReuse(sampleValue, pdf);
	return bp::make_tuple(index, pdf, sampleValue);
}

static Float DiscreteDistribution_getitem(DiscreteDistribution *d, int i) {
	if (i < 0 || i >= (int) d->size()) {
		SLog(EError, "Index %i is out of range!", i);
		return 0.0f;
	}
	return d->operator[](i);
}

static bp::tuple legendrePD_double(int l, double x) {
	std::pair<double, double> result = legendrePD(l, x);
	return bp::make_tuple(result.first, result.second);
}

static bp::tuple gaussLegendre_(int n) {
	Float *nodes  = new Float[n];
	Float *weights= new Float[n];
	gaussLegendre(n, nodes, weights);

	bp::list nodeList, weightList;
	for (int i=0; i<n; ++i) {
		nodeList.append(nodes[i]);
		weightList.append(weights[i]);
	}

	delete[] nodes;
	delete[] weights;
	return bp::make_tuple(nodeList, weightList);
}

static bp::tuple gaussLobatto_(int n) {
	Float *nodes  = new Float[n];
	Float *weights= new Float[n];
	gaussLobatto(n, nodes, weights);

	bp::list nodeList, weightList;
	for (int i=0; i<n; ++i) {
		nodeList.append(nodes[i]);
		weightList.append(weights[i]);
	}

	delete[] nodes;
	delete[] weights;
	return bp::make_tuple(nodeList, weightList);
}

/**
 * \brief Computes the nodes and weights of a Gauss-Lobatto quadrature
 * rule with the given number of evaluations.
 *
 * Integration is over the interval \f$[-1, 1]\f$. Gauss-Lobatto quadrature
 * is preferable to Gauss-Legendre quadrature whenever the endpoints of the
 * integration domain should explicitly be included. It maximizes the order
 * of exactly integrable polynomials subject to this constraint and achieves
 * this up to degree \f$2n-3\f$ (where \f$n\f$ is the number of function
 * evaluations).
 *
 * This method is numerically well-behaved until about \f$n=200\f$
 * and then becomes progressively less accurate. It is generally not a
 * good idea to go much higher---in any case, a composite or
 * adaptive integration scheme will be superior for large \f$n\f$.
 *
 * \param n
 *     Desired number of evalution points
 * \param nodes
 *     Length-\c n array used to store the nodes of the quadrature rule
 * \param nodes
 *     Length-\c n array used to store the weights of the quadrature rule
 */
extern MTS_EXPORT_CORE void gaussLobatto(int n, Float *nodes, Float *weights);

struct NativeBuffer {
	ref<Object> owner;
	void *ptr;
	Bitmap::EComponentFormat format;
	int ndim;
	Py_ssize_t shape[3], strides[4];
	const char* formatString;

	NativeBuffer(Object *owner, void *ptr, Bitmap::EComponentFormat format, int ndim,
			Py_ssize_t shape[3]) : owner(owner), ptr(ptr), format(format), ndim(ndim) {
		size_t itemSize = 0;
		switch (format) {
			case Bitmap::EUInt8:   formatString = "B"; itemSize = 1; break;
			case Bitmap::EUInt16:  formatString = "H"; itemSize = 2; break;
			case Bitmap::EUInt32:  formatString = "I"; itemSize = 4; break;
			case Bitmap::EFloat16: formatString = "e"; itemSize = 2; break;
			case Bitmap::EFloat32: formatString = "f"; itemSize = 4; break;
			case Bitmap::EFloat64: formatString = "d"; itemSize = 8; break;
			default:
				SLog(EError, "Unsupported bufer format!");
		}
		strides[ndim] = itemSize;

		for (int i=ndim-1; i>=0; --i) {
			this->shape[i] = shape[i];
			strides[i] = strides[i+1] * shape[i];
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "NativeBuffer[ndim=" << ndim << ", shape=[";
		for (int i=0; i<ndim; ++i) {
			oss << shape[i];
			if (i+1 < ndim)
				oss << ", ";
		}
		oss << "], strides=[";
		for (int i=0; i<=ndim; ++i) {
			oss << strides[i];
			if (i+1 <= ndim)
				oss << ", ";
		}
		oss << "], format=" << format << ", size=" << memString(strides[0]) << "]";
		return oss.str();
	}

	static int getbuffer(PyObject *obj, Py_buffer *view, int flags) {
		bp::extract<NativeBuffer&> b(obj);
		if (!b.check()) {
			PyErr_SetString(PyExc_BufferError, "Native buffer is invalid!");
			view->obj = NULL;
			return -1;
		}
		NativeBuffer &buf = b();

		if (!buf.ptr) {
			PyErr_SetString(PyExc_BufferError, "Native buffer does not point anywhere!");
			view->obj = NULL;
			return -1;
		}

		if (view == NULL)
			return 0;

		view->obj = obj;
		if (view->obj)
			Py_INCREF(view->obj);
		buf.owner->incRef();

		view->ndim = 1;
		view->buf = buf.ptr;
		view->format = NULL;
		view->shape = NULL;
		view->suboffsets = NULL;
		view->internal = NULL;
		view->strides = NULL;
		view->len = buf.strides[0];
		view->readonly = false;
		view->itemsize = buf.strides[buf.ndim];

		if ((flags & PyBUF_FORMAT) == PyBUF_FORMAT)
			view->format = const_cast<char *>(buf.formatString);

		if ((flags & PyBUF_STRIDES) == PyBUF_STRIDES)
			view->strides = &buf.strides[1];

		if ((flags & PyBUF_ND) == PyBUF_ND) {
			view->ndim = buf.ndim;
			view->shape = &buf.shape[0];
		}

		return 0;
	}

	static void releasebuffer(PyObject *obj, Py_buffer *view) {
		bp::extract<NativeBuffer&> b(obj);
		if (!b.check()) {
			PyErr_SetString(PyExc_BufferError, "Native buffer is invalid!");
			return;
		}
		NativeBuffer &buf = b();
		buf.owner->decRef();
	}

	static Py_ssize_t len(PyObject *obj) {
		bp::extract<NativeBuffer&> b(obj);
		if (!b.check()) {
			PyErr_SetString(PyExc_BufferError, "Native buffer is invalid!");
			return -1;
		}
		NativeBuffer &buf = b();
		return buf.strides[0] / buf.strides[buf.ndim];
	}

	static PyObject* item(PyObject *obj, Py_ssize_t idx) {
		bp::extract<NativeBuffer&> b(obj);
		if (!b.check()) {
			PyErr_SetString(PyExc_BufferError, "Native buffer is invalid!");
			return 0;
		}
		NativeBuffer &buf = b();

		bp::object result;
		switch (buf.format) {
			case Bitmap::EUInt8:   result = bp::object(((uint8_t *) buf.ptr)[idx]); break;
			case Bitmap::EUInt16:  result = bp::object(((uint16_t *) buf.ptr)[idx]); break;
			case Bitmap::EUInt32:  result = bp::object(((uint32_t *) buf.ptr)[idx]); break;
			case Bitmap::EFloat16: result = bp::object((float) ((half *) buf.ptr)[idx]); break;
			case Bitmap::EFloat32: result = bp::object(((float *) buf.ptr)[idx]); break;
			case Bitmap::EFloat64: result = bp::object(((double *) buf.ptr)[idx]); break;
			default:
				PyErr_SetString(PyExc_BufferError, "Unsupported buffer format!");
				return 0;
		}

		return bp::incref(result.ptr());
	}
};

static NativeBuffer bitmap_buffer(Bitmap *bitmap) {
	int ndim = bitmap->getChannelCount() == 1 ? 2 : 3;
	Py_ssize_t shape[3] = {
		(Py_ssize_t) bitmap->getHeight(),
		(Py_ssize_t) bitmap->getWidth(),
		(Py_ssize_t) bitmap->getChannelCount()
	};

	return NativeBuffer(bitmap, bitmap->getUInt8Data(), bitmap->getComponentFormat(), ndim, shape);
}


static ref<Bitmap> bitmap_array_constructor(bp::object _obj) {
	PyObject *obj = _obj.ptr();
	if (!obj)
		SLog(EError, "Expected a non-NULL argument!");

	bp::extract<fs::path> extractPath(_obj);
	if (extractPath.check())
		return new Bitmap(extractPath());

	Py_buffer buffer;
	if (PyObject_GetBuffer(obj, &buffer, PyBUF_CONTIG_RO | PyBUF_FORMAT))
		SLog(EError, "Could not access supplied object using the buffer protocol!");

	Vector2i size(1);
	Bitmap::EPixelFormat pixelFormat = Bitmap::ELuminance;
	Bitmap::EComponentFormat componentFormat = Bitmap::EUInt8;
	int nChannels = 1;

	if (buffer.ndim == 0 || buffer.ndim > 3)
		SLog(EError, "Invalid number of dimensions!");

	if (buffer.ndim == 1) {
		size.x = buffer.shape[0];
	} else if (buffer.ndim > 1) {
		size.y = buffer.shape[0];
		size.x = buffer.shape[1];
	}
	if (buffer.ndim > 2) {
		nChannels = buffer.shape[2];
		if (nChannels == 1)
			pixelFormat = Bitmap::ELuminance;
		else if (nChannels == 2)
			pixelFormat = Bitmap::ELuminanceAlpha;
		else if (nChannels == 3)
			pixelFormat = Bitmap::ERGB;
		else if (nChannels == 4)
			pixelFormat = Bitmap::ERGBA;
		else
			pixelFormat = Bitmap::EMultiChannel;
	}
	if (strlen(buffer.format) != 1)
		SLog(EError, "Invalid buffer format \"%s\"", buffer.format);

	switch (buffer.format[0]) {
		case 'B': componentFormat = Bitmap::EUInt8; break;
		case 'H': componentFormat = Bitmap::EUInt16; break;
		case 'I': componentFormat = Bitmap::EUInt32; break;
		case 'e': componentFormat = Bitmap::EFloat16; break;
		case 'f': componentFormat = Bitmap::EFloat32; break;
		case 'd': componentFormat = Bitmap::EFloat64; break;
		default:
			SLog(EError, "Invalid buffer format \"%s\"", buffer.format);
	}

	ref<Bitmap> result = new Bitmap(pixelFormat, componentFormat, size, nChannels);
	if ((size_t) buffer.len != result->getBufferSize())
		SLog(EError, "Internal error: Python buffer size and Mitsuba bitmap size disagree: "
			SIZE_T_FMT " vs " SIZE_T_FMT, (size_t) buffer.len, (size_t) result->getBufferSize());

	memcpy(result->getData(), buffer.buf, result->getBufferSize());

	PyBuffer_Release(&buffer);
	return result;
}
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(fromLinearRGB_overloads, fromLinearRGB, 3, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(fromXYZ_overloads, fromXYZ, 3, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(fromIPT_overloads, fromIPT, 3, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(reset_overloads, reset, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(resample1_overloads, resample, 4, 7)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(resample2_overloads, resample, 6, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(filter1_overloads, filter, 4, 7)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(filter2_overloads, filter, 3, 5)

#define IMPLEMENT_ANIMATION_TRACK(Name) \
	BP_CLASS(Name, AbstractAnimationTrack, (bp::init<AbstractAnimationTrack::EType, size_t>())) \
		.def(bp::init<Name *>()) \
		.def("reserve", &Name::reserve) \
		.def("prependTransformation", &Name::prependTransformation) \
		.def("appendTransformation", &Name::appendTransformation) \
		.def("eval", &Name::eval, BP_RETURN_VALUE) \
		.def("setValue", &Name::setValue) \
		.def("getValue", &Name::getValue, BP_RETURN_VALUE) \
		.def("append", &Name::append)

struct PythonIntegrand {
	PythonIntegrand(bp::object integrand) : integrand(integrand) {}

	Float operator()(Float value) {
		bp::object obj = integrand(value);
		bp::extract<Float> extract(obj);
		return (Float) extract();
	}

	bp::object integrand;
};

struct PythonIntegrandFromPythonCallable {
	PythonIntegrandFromPythonCallable() {
		bp::converter::registry::push_back(&convertible,
				&construct, bp::type_id<GaussLobattoIntegrator::Integrand>());
	}

	static void* convertible(PyObject* obj) {
		if(!PyCallable_Check(obj))
			return 0;
		return obj;
	}

	static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data* data) {
		bp::object callable(bp::handle<>(bp::borrowed(obj)));
		void* storage = ((bp::converter::rvalue_from_python_storage<GaussLobattoIntegrator::Integrand>*) data)->storage.bytes;
		new (storage) GaussLobattoIntegrator::Integrand(PythonIntegrand(callable));
		data->convertible = storage;
	}
};

static bp::tuple GaussLobattoIntegrator_integrate(GaussLobattoIntegrator *integrator,
		const GaussLobattoIntegrator::Integrand &integrand, Float a, Float b) {
	size_t nEvals = 0;
	Float result = integrator->integrate(integrand, a, b, &nEvals);
	return bp::make_tuple(result, nEvals);
}

static std::string memString1(size_t size) { return mitsuba::memString(size); }
static std::string memString2(size_t size, bool precise) { return mitsuba::memString(size, precise); }
static std::string timeString1(Float size) { return mitsuba::timeString(size); }
static std::string timeString2(Float size, bool precise) { return mitsuba::timeString(size, precise); }

// Based on 'http://stackoverflow.com/questions/15842126/feeding-a-python-list-into-a-function-taking-in-a-vector-with-boost-python'
struct iterable_converter {
	template <typename Container> iterable_converter& from_python() {
		bp::converter::registry::push_back(
			&iterable_converter::convertible,
			&iterable_converter::construct<Container>,
			bp::type_id<Container>());
		return *this;
	}

	static void* convertible(PyObject* object) {
		return PyObject_GetIter(object) ? object : NULL;
	}

	template <typename Container> static void construct(
		PyObject* object,
		bp::converter::rvalue_from_python_stage1_data* data) {

		// Object is a borrowed reference, so create a handle indicting it is
		// borrowed for proper reference counting.
		bp::handle<> handle(bp::borrowed(object));

		// Obtain a handle to the memory block that the converter has allocated
		// for the C++ type.
		typedef bp::converter::rvalue_from_python_storage<Container> storage_type;
		void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;

		typedef bp::stl_input_iterator<typename Container::value_type> iterator;
		// Allocate the C++ type into the converter's memory block, and assign
		// its handle to the converter's convertible variable.  The C++
		// container is populated by passing the begin and end iterators of
		// the python object to the container's constructor.
		data->convertible = new (storage) Container(
			iterator(bp::object(handle)),
			iterator());
	}
};

void export_core() {
	/* Set up various implicit conversions */
	bp::to_python_converter<fs::path, path_to_python_str>();
	bp::to_python_converter<TSpectrum<Float, SPECTRUM_SAMPLES>, TSpectrum_to_Spectrum>();
	bp::implicitly_convertible<std::string, fs::path>();
	PythonIntegrandFromPythonCallable();

	iterable_converter()
		.from_python<StringVector>();

	bp::object coreModule(
		bp::handle<>(bp::borrowed(PyImport_AddModule("mitsuba.core"))));
	bp::scope().attr("core") = coreModule;
	PyObject *oldScope = bp::detail::current_scope;

	BP_SETSCOPE(coreModule);

	coreModule.attr("__path__") = "mitsuba.core";

	/* Basic STL containers */
	bp::class_<StringVector>("StringVector")
		.def(bp::vector_indexing_suite<StringVector>())
		.def("__init__", bp::make_constructor(StringVector_fromList));

	bp::class_<StringMap>("StringMap")
		.def(bp::map_indexing_suite<StringMap>());

	/* Logging */
	bp::enum_<ELogLevel>("ELogLevel")
		.value("ETrace", ETrace)
		.value("EDebug", EDebug)
		.value("EInfo", EInfo)
		.value("EWarn", EWarn)
		.value("EError", EError)
		.export_values();

	bp::def("Log", &mts_log);

	/* Basic constants */
	coreModule.attr("Epsilon") = Epsilon;
	coreModule.attr("ShadowEpsilon") = ShadowEpsilon;
	coreModule.attr("DeltaEpsilon") = DeltaEpsilon;
	coreModule.attr("SPECTRUM_SAMPLES") = SPECTRUM_SAMPLES;
	coreModule.attr("MTS_VERSION") = MTS_VERSION;
	coreModule.attr("MTS_YEAR") = MTS_YEAR;

	#if defined(SINGLE_PRECISION)
		coreModule.attr("DOUBLE_PRECISION") = 0;
		coreModule.attr("SINGLE_PRECISION") = 1;
	#else
		coreModule.attr("DOUBLE_PRECISION") = 1;
		coreModule.attr("SINGLE_PRECISION") = 0;
	#endif

	bp::class_<Class, boost::noncopyable>("Class", bp::no_init)
		.def("getName", &Class::getName, BP_RETURN_CONSTREF)
		.def("isAbstract", &Class::isAbstract)
		.def("isInstantiable", &Class::isInstantiable)
		.def("isSerializable", &Class::isSerializable)
		.def("derivesFrom", &Class::derivesFrom)
		.def("getSuperClass", &class_getSuperClass, BP_RETURN_INTREF)
		.def("forName", &class_forName, BP_RETURN_INTREF)
		.def("unserialize", &Class::unserialize, BP_RETURN_VALUE)
		.def("instantiate", &Class::instantiate, BP_RETURN_VALUE)
		.staticmethod("forName");
	bp::register_ptr_to_python<Class*>();

	bp::class_<Object, ref<Object>, boost::noncopyable>("Object", bp::no_init)
		.def("getRefCount", &Object::getRefCount)
		.def("getClass", &object_getClass, BP_RETURN_INTREF)
		.def("__repr__", &Object::toString);
	bp::register_ptr_to_python<Object*>();

	BP_CLASS(Stream, Object, bp::no_init)
		.def("setByteOrder", &Stream::setByteOrder)
		.def("getByteOrder", &Stream::getByteOrder)
		.def("getHostByteOrder", &Stream::getHostByteOrder)
		.def("truncate", &Stream::truncate)
		.def("seek", &Stream::seek)
		.def("getPos", &Stream::getPos)
		.def("getSize", &Stream::getSize)
		.def("flush", &Stream::flush)
		.def("canWrite", &Stream::canWrite)
		.def("canRead", &Stream::canRead)
		.def("skip", &Stream::skip)
		.def("copyTo", &Stream::copyTo)
		.def("writeString", &Stream::writeString)
		.def("readString", &Stream::readString)
		.def("writeLine", &Stream::writeLine)
		.def("readLine", &Stream::readLine)
		.def("writeChar", &Stream::writeChar)
		.def("readChar", &Stream::readChar)
		.def("writeUChar", &Stream::writeUChar)
		.def("readUChar", &Stream::readUChar)
		.def("writeShort", &Stream::writeShort)
		.def("readShort", &Stream::readShort)
		.def("writeUShort", &Stream::writeUShort)
		.def("readUShort", &Stream::readUShort)
		.def("writeInt", &Stream::writeInt)
		.def("readInt", &Stream::readInt)
		.def("writeUInt", &Stream::writeUInt)
		.def("readUInt", &Stream::readUInt)
		.def("writeLong", &Stream::writeLong)
		.def("readLong", &Stream::readLong)
		.def("writeULong", &Stream::writeULong)
		.def("readULong", &Stream::readULong)
		.def("writeFloat", &Stream::writeFloat)
		.def("readFloat", &Stream::readFloat)
		.def("writeSingle", &Stream::writeSingle)
		.def("readSingle", &Stream::readSingle)
		.def("writeDouble", &Stream::writeDouble)
		.def("readDouble", &Stream::readDouble);

	BP_SETSCOPE(Stream_class);
	bp::enum_<Stream::EByteOrder>("EByteOrder")
		.value("EBigEndian", Stream::EBigEndian)
		.value("ELittleEndian", Stream::ELittleEndian)
		.value("ENetworkByteOrder", Stream::ENetworkByteOrder)
		.export_values();
	BP_SETSCOPE(coreModule);

	BP_CLASS(FileStream, Stream, bp::init<>())
		.def(bp::init<std::string, FileStream::EFileMode>())
		.def("getPath", &FileStream::getPath,
				BP_RETURN_CONSTREF)
		.def("open", &FileStream::open)
		.def("close", &FileStream::close)
		.def("remove", &FileStream::remove)
		.def("createTemporary", &FileStream::createTemporary)
		.staticmethod("createTemporary");

	BP_CLASS(SocketStream, Stream, (bp::init<std::string, int>()))
		.def("getPeer", &SocketStream::getPeer, BP_RETURN_CONSTREF)
		.def("getReceivedBytes", &SocketStream::getReceivedBytes)
		.def("getSentBytes", &SocketStream::getSentBytes);

	BP_CLASS_DECL(ConsoleStream, Stream, bp::init<>());

	BP_CLASS(SSHStream, Stream, (bp::init<std::string, std::string, const StringVector &>()))
		.def(bp::init<std::string, std::string, const StringVector &, int>())
		.def(bp::init<std::string, std::string, const StringVector &, int, int>())
		.def("getUserName", &SSHStream::getUserName, BP_RETURN_CONSTREF)
		.def("getHostName", &SSHStream::getHostName, BP_RETURN_CONSTREF)
		.def("getReceivedBytes", &SSHStream::getReceivedBytes)
		.def("getSentBytes", &SSHStream::getSentBytes);

	BP_SETSCOPE(FileStream_class);
	bp::enum_<FileStream::EFileMode>("EFileMode")
		.value("EReadOnly", FileStream::EReadOnly)
		.value("EReadWrite", FileStream::EReadWrite)
		.value("ETruncWrite", FileStream::ETruncWrite)
		.value("ETruncReadWrite", FileStream::ETruncReadWrite)
		.value("EAppendWrite", FileStream::EAppendWrite)
		.value("EAppendReadWrite", FileStream::EAppendReadWrite)
		.export_values();
	BP_SETSCOPE(coreModule);

	BP_CLASS(MemoryStream, Stream, (bp::init<bp::optional<size_t> >()))
		.def("reset", &MemoryStream::reset);

	BP_CLASS(SerializableObject, Object, bp::no_init)
		.def("serialize", &SerializableObject::serialize);

	void (ConfigurableObject::*cobject_add_child_1)(ConfigurableObject *) = &ConfigurableObject::addChild;
	void (ConfigurableObject::*cobject_add_child_2)(const std::string &, ConfigurableObject *) = &ConfigurableObject::addChild;

	BP_CLASS(ConfigurableObject, SerializableObject, bp::no_init)
		.def("setParent", &ConfigurableObject::setParent)
		.def("getID", &ConfigurableObject::getID, BP_RETURN_VALUE)
		.def("setID", &ConfigurableObject::setID)
		.def("getProperties", &ConfigurableObject::getProperties, BP_RETURN_VALUE)
		.def("addChild", cobject_add_child_1)
		.def("addChild", cobject_add_child_2)
		.def("configure", &ConfigurableObject::configure);

	BP_CLASS(NetworkedObject, ConfigurableObject, bp::no_init)
		.def("bindUsedResources", &NetworkedObject::bindUsedResources);

	Thread *(Thread::*thread_get_parent)() = &Thread::getParent;
	BP_CLASS(Thread, Object, bp::no_init)
		.def("getID", &Thread::getID)
		.def("setPriority", &Thread::setPriority)
		.def("getPriority", &Thread::getPriority)
		.def("setCoreAffinity", &Thread::setCoreAffinity)
		.def("getCoreAffinity", &Thread::getCoreAffinity)
		.def("setCritical", &Thread::setCritical)
		.def("getCritical", &Thread::getCritical)
		.def("setName", &Thread::setName)
		.def("getName", &Thread::getName, BP_RETURN_CONSTREF)
		.def("getParent", thread_get_parent, BP_RETURN_VALUE)
		.def("setLogger", &Thread::setLogger)
		.def("getLogger", &Thread::getLogger, BP_RETURN_VALUE)
		.def("setFileResolver", &Thread::setFileResolver)
		.def("getFileResolver", &Thread::getFileResolver, BP_RETURN_VALUE)
		.def("getThread", &Thread::getThread, BP_RETURN_VALUE)
		.def("isRunning", &Thread::isRunning)
		.def("registerUnmanagedThread", &Thread::registerUnmanagedThread, BP_RETURN_VALUE)
		.def("sleep", &Thread::sleep)
		.def("detach", &Thread::detach)
		.def("join", thread_join)
		.def("start", &Thread::start)
		.staticmethod("sleep")
		.staticmethod("getThread")
		.staticmethod("registerUnmanagedThread");

	BP_SETSCOPE(Thread_class);
	bp::enum_<Thread::EThreadPriority>("EThreadPriority")
		.value("EIdlePriority", Thread::EIdlePriority)
		.value("ELowestPriority", Thread::ELowestPriority)
		.value("ELowPriority", Thread::ELowPriority)
		.value("ENormalPriority", Thread::ENormalPriority)
		.value("EHighPriority", Thread::EHighPriority)
		.value("EHighestPriority", Thread::EHighestPriority)
		.value("ERealtimePriority", Thread::ERealtimePriority)
		.export_values();
	BP_SETSCOPE(coreModule);

	BP_WRAPPED_CLASS(Appender, AppenderWrapper, Object, bp::init<>())
		.def("append", &Appender::append)
		.def("logProgress", &appender_logProgress);

	BP_CLASS(StreamAppender, Appender, bp::init<std::string>())
		.def("logsToFile", &StreamAppender::logsToFile)
		.def("readLog", &StreamAppender::readLog);

	BP_WRAPPED_CLASS(Formatter, FormatterWrapper, Object, bp::init<>())
		.def("format", &Formatter::format);

	BP_CLASS(DefaultFormatter, Formatter, bp::init<>())
		.def("setHaveDate", &DefaultFormatter::setHaveDate)
		.def("setHaveThread", &DefaultFormatter::setHaveThread)
		.def("setHaveLogLevel", &DefaultFormatter::setHaveLogLevel)
		.def("setHaveClass", &DefaultFormatter::setHaveClass);

	Appender *(Logger::*logger_get_appender)(size_t) = &Logger::getAppender;
	BP_CLASS(Logger, Object, bp::init<ELogLevel>())
		.def("logProgress", logger_logProgress)
		.def("setLogLevel", &Logger::setLogLevel)
		.def("getLogLevel", &Logger::getLogLevel)
		.def("setErrorLevel", &Logger::setErrorLevel)
		.def("getErrorLevel", &Logger::getErrorLevel)
		.def("addAppender", &Logger::addAppender)
		.def("removeAppender", &Logger::removeAppender)
		.def("clearAppenders", &Logger::clearAppenders)
		.def("getAppenderCount", &Logger::getAppenderCount)
		.def("getAppender", logger_get_appender, BP_RETURN_VALUE)
		.def("getFormatter", &Logger::getFormatter, BP_RETURN_VALUE)
		.def("setFormatter", &Logger::setFormatter)
		.def("readLog", &logger_readLog)
		.def("getWarningCount", &Logger::getWarningCount);

	BP_CLASS(InstanceManager, Object, bp::init<>())
		.def("serialize", &InstanceManager::serialize)
		.def("getInstance", &instance_manager_getinstance, BP_RETURN_VALUE);

	bp::class_<ContinuousSpectrum, boost::noncopyable>("ContinuousSpectrum", bp::no_init)
		.def("eval", &ContinuousSpectrum::eval)
		.def("average", &ContinuousSpectrum::average)
		.def("__repr__", &ContinuousSpectrum::toString);

	bp::class_<InterpolatedSpectrum, bp::bases<ContinuousSpectrum>, boost::noncopyable>
			("InterpolatedSpectrum", bp::init<>())
		.def(bp::init<size_t>())
		.def(bp::init<fs::path>())
		.def("append", &InterpolatedSpectrum::append)
		.def("clear", &InterpolatedSpectrum::clear)
		.def("zeroExtend", &InterpolatedSpectrum::zeroExtend);

	bp::class_<BlackBodySpectrum, bp::bases<ContinuousSpectrum>, boost::noncopyable>
			("BlackBodySpectrum", bp::init<Float>());

	void (Bitmap::*accumulate_1)(const Bitmap *bitmap, Point2i sourceOffset, Point2i targetOffset, Vector2i size) = &Bitmap::accumulate;
	void (Bitmap::*accumulate_2)(const Bitmap *bitmap, Point2i targetOffset) = &Bitmap::accumulate;
	void (Bitmap::*accumulate_3)(const Bitmap *bitmap) = &Bitmap::accumulate;
	void (Bitmap::*copyFrom_1)(const Bitmap *bitmap, Point2i sourceOffset, Point2i targetOffset, Vector2i size) = &Bitmap::copyFrom;
	void (Bitmap::*copyFrom_2)(const Bitmap *bitmap, Point2i targetOffset) = &Bitmap::copyFrom;
	void (Bitmap::*copyFrom_3)(const Bitmap *bitmap) = &Bitmap::copyFrom;
	const Properties &(Bitmap::*get_metadata)() const = &Bitmap::getMetadata;

	void (Bitmap::*resample_1)(const ReconstructionFilter *,
		ReconstructionFilter::EBoundaryCondition, ReconstructionFilter::EBoundaryCondition,
		Bitmap *, Bitmap *, Float, Float) const  = &Bitmap::resample;

	ref<Bitmap> (Bitmap::*resample_2)(const ReconstructionFilter *,
		ReconstructionFilter::EBoundaryCondition, ReconstructionFilter::EBoundaryCondition,
		const Vector2i &, Float, Float) const  = &Bitmap::resample;

	void (Bitmap::*filter_1)(const ReconstructionFilter *,
		ReconstructionFilter::EBoundaryCondition, ReconstructionFilter::EBoundaryCondition,
		Bitmap *, Bitmap *, Float, Float) const  = &Bitmap::filter;

	ref<Bitmap> (Bitmap::*filter_2)(const ReconstructionFilter *,
		ReconstructionFilter::EBoundaryCondition, ReconstructionFilter::EBoundaryCondition,
		Float, Float) const  = &Bitmap::filter;

	const std::vector<std::string> & (Bitmap::*getChannelNames_1)() const = &Bitmap::getChannelNames;

	BP_CLASS(Bitmap, Object, (bp::init<Bitmap::EPixelFormat, Bitmap::EComponentFormat, const Vector2i &>()))
		.def(bp::init<Bitmap::EPixelFormat, Bitmap::EComponentFormat, const Vector2i &, int>())
		.def(bp::init<Bitmap::EFileFormat, Stream *, bp::optional<std::string> >())
		.def(bp::init<fs::path, bp::optional<std::string> >())
		.def("__init__", bp::make_constructor(bitmap_array_constructor))
		.def("getPixelFormat", &Bitmap::getPixelFormat)
		.def("getComponentFormat", &Bitmap::getComponentFormat)
		.def("getSize", &Bitmap::getSize, BP_RETURN_VALUE)
		.def("getPixelCount", &Bitmap::getPixelCount)
		.def("average", &Bitmap::average)
		.def("getWidth", &Bitmap::getWidth)
		.def("getHeight", &Bitmap::getHeight)
		.def("getChannelCount", &Bitmap::getChannelCount)
		.def("getChannelName", &Bitmap::getChannelName, BP_RETURN_VALUE)
		.def("getChannelNames", getChannelNames_1, BP_RETURN_VALUE)
		.def("setChannelNames", &Bitmap::setChannelNames)
		.def("isSquare", &Bitmap::isSquare)
		.def("hasAlpha", &Bitmap::hasAlpha)
		.def("hasWeight", &Bitmap::hasWeight)
		.def("isMultiChannel", &Bitmap::isMultiChannel)
		.def("getBitsPerComponent", &Bitmap::getBitsPerComponent)
		.def("getBytesPerComponent", &Bitmap::getBytesPerComponent)
		.def("getBytesPerPixel", &Bitmap::getBytesPerPixel)
		.def("getBufferSize", &Bitmap::getBufferSize)
		.def("getPixel", &Bitmap::getPixel, BP_RETURN_VALUE)
		.def("setPixel", &Bitmap::setPixel)
		.def("drawHLine", &Bitmap::drawHLine)
		.def("drawVLine", &Bitmap::drawVLine)
		.def("clone", &Bitmap::clone, BP_RETURN_VALUE)
		.def("clear", &Bitmap::clear)
		.def("write", &bitmap_write1)
		.def("write", &bitmap_write2)
		.def("write", &bitmap_write3)
		.def("write", &bitmap_write4)
		.def("write", &bitmap_write5)
		.def("write", &bitmap_write6)
		.def("tonemapReinhard", &bitmap_tonemapReinhard)
		.def("expand", &Bitmap::expand, BP_RETURN_VALUE)
		.def("extractChannel", &Bitmap::extractChannel, BP_RETURN_VALUE)
		.def("extractChannels", bitmap_extractChannels, BP_RETURN_VALUE)
		.def("join", &bitmap_join, BP_RETURN_VALUE)
		.def("crop", &Bitmap::crop)
		.def("flipVertically", &Bitmap::flipVertically)
		.def("rotateFlip", &Bitmap::rotateFlip, BP_RETURN_VALUE)
		.def("scale", &Bitmap::scale)
		.def("pow", &Bitmap::pow)
		.def("colorBalance", &Bitmap::colorBalance)
		.def("applyMatrix", &bitmap_applyMatrix)
		.def("accumulate", accumulate_1)
		.def("accumulate", accumulate_2)
		.def("accumulate", accumulate_3)
		.def("copyFrom", copyFrom_1)
		.def("copyFrom", copyFrom_2)
		.def("copyFrom", copyFrom_3)
		.def("convolve", &Bitmap::convolve)
		.def("arithmeticOperation", &Bitmap::arithmeticOperation, BP_RETURN_VALUE)
		.def("resample", resample_1, resample1_overloads())
		.def("resample", resample_2, resample2_overloads()[BP_RETURN_VALUE])
		.def("filter", filter_1, filter1_overloads())
		.def("filter", filter_2, filter2_overloads()[BP_RETURN_VALUE])
		.def("setGamma", &Bitmap::setGamma)
		.def("getGamma", &Bitmap::getGamma)
		.def("setMetadataString", &Bitmap::setMetadataString)
		.def("getMetadataString", &Bitmap::getMetadataString, BP_RETURN_VALUE)
		.def("setMetadata", &Bitmap::setMetadata)
		.def("getMetadata", get_metadata, BP_RETURN_VALUE)
		.def("drawRect", &Bitmap::drawRect)
		.def("fillRect", &Bitmap::fillRect)
		.def("drawWorkUnit", &Bitmap::drawWorkUnit)
		.def("convert", &bitmap_convert_0)
		.def("convert", &bitmap_convert_1)
		.def("convert", &bitmap_convert_2)
		.def("convert", &bitmap_convert_3, BP_RETURN_VALUE)
		.def("convert", &bitmap_convert_4, BP_RETURN_VALUE)
		.def("convert", &bitmap_convert_5, BP_RETURN_VALUE)
		.def("convert", &bitmap_convert_6, BP_RETURN_VALUE)
		.def("fromByteArray", &bitmap_fromByteArray)
		.def("toByteArray", &bitmap_toByteArray_1)
		.def("toByteArray", &bitmap_toByteArray_2)
		.def("buffer", bitmap_buffer)
		.def("split", bitmap_split)
		.def("plot", bitmap_plot)
		.staticmethod("join")
		.staticmethod("arithmeticOperation");

	BP_SETSCOPE(Bitmap_class);
	bp::enum_<Bitmap::EPixelFormat>("EPixelFormat")
		.value("ELuminance", Bitmap::ELuminance)
		.value("ELuminanceAlpha", Bitmap::ELuminanceAlpha)
		.value("ERGB", Bitmap::ERGB)
		.value("ERGBA", Bitmap::ERGBA)
		.value("EXYZ", Bitmap::EXYZ)
		.value("EXYZA", Bitmap::EXYZA)
		.value("ESpectrum", Bitmap::ESpectrum)
		.value("ESpectrumAlpha", Bitmap::ESpectrumAlpha)
		.value("ESpectrumAlphaWeight", Bitmap::ESpectrumAlphaWeight)
		.value("EMultiSpectrumAlphaWeight", Bitmap::EMultiSpectrumAlphaWeight)
		.value("EMultiChannel", Bitmap::EMultiChannel)
		.export_values();

	bp::enum_<Bitmap::EComponentFormat>("EComponentFormat")
		.value("EBitmask", Bitmap::EBitmask)
		.value("EUInt8", Bitmap::EUInt8)
		.value("EUInt16", Bitmap::EUInt16)
		.value("EUInt32", Bitmap::EUInt32)
		.value("EFloat16", Bitmap::EFloat16)
		.value("EFloat32", Bitmap::EFloat32)
		.value("EFloat64", Bitmap::EFloat64)
		.value("EFloat", Bitmap::EFloat)
		.value("EInvalid", Bitmap::EInvalid)
		.export_values();

	bp::enum_<Bitmap::EFileFormat>("EFileFormat")
		.value("EPNG", Bitmap::EPNG)
		.value("EOpenEXR", Bitmap::EOpenEXR)
		.value("ETGA", Bitmap::ETGA)
		.value("EPFM", Bitmap::EPFM)
		.value("EPPM", Bitmap::EPPM)
		.value("ERGBE", Bitmap::ERGBE)
		.value("EBMP", Bitmap::EBMP)
		.value("EJPEG", Bitmap::EJPEG)
		.value("EAuto", Bitmap::EAuto)
		.export_values();

	bp::enum_<Bitmap::EArithmeticOperation>("EArithmeticOperation")
		.value("EAddition", Bitmap::EAddition)
		.value("ESubtraction", Bitmap::ESubtraction)
		.value("EMultiplication", Bitmap::EMultiplication)
		.value("EDivision", Bitmap::EDivision)
		.export_values();

	bp::enum_<Bitmap::ERotateFlipType>("ERotateFlipType")
		.value("ERotateNoneFlipNone", Bitmap::ERotateNoneFlipNone)
		.value("ERotate180FlipXY", Bitmap::ERotate180FlipXY)
		.value("ERotate90FlipNone", Bitmap::ERotate90FlipNone)
		.value("ERotate270FlipXY", Bitmap::ERotate270FlipXY)
		.value("ERotate180FlipNone", Bitmap::ERotate180FlipNone)
		.value("ERotateNoneFlipXY", Bitmap::ERotateNoneFlipXY)
		.value("ERotate270FlipNone", Bitmap::ERotate270FlipNone)
		.value("ERotate90FlipXY", Bitmap::ERotate90FlipXY)
		.value("ERotateNoneFlipX", Bitmap::ERotateNoneFlipX)
		.value("ERotate180FlipY", Bitmap::ERotate180FlipY)
		.value("ERotate90FlipX", Bitmap::ERotate90FlipX)
		.value("ERotate270FlipY", Bitmap::ERotate270FlipY)
		.value("ERotate180FlipX", Bitmap::ERotate180FlipX)
		.value("ERotateNoneFlipY", Bitmap::ERotateNoneFlipY)
		.value("ERotate270FlipX", Bitmap::ERotate270FlipX)
		.value("ERotate90FlipY", Bitmap::ERotate90FlipY)
		.export_values();

	BP_SETSCOPE(coreModule);

	/* Native buffers for bitmaps */
	bp::class_<NativeBuffer>("NativeBuffer", bp::no_init)
		.def("__repr__", &NativeBuffer::toString);

	const bp::converter::registration& fb_reg(
		bp::converter::registry::lookup(bp::type_id<NativeBuffer>()));
	PyTypeObject* fb_type = fb_reg.get_class_object();

	static PyBufferProcs NativeBuffer_buffer_procs = {
		#if PY_MAJOR_VERSION < 3
			NULL, NULL, NULL, NULL,
		#endif
		&NativeBuffer::getbuffer,
		&NativeBuffer::releasebuffer
	};

	// partial sequence protocol support
	static PySequenceMethods NativeBuffer_as_sequence = {
		&NativeBuffer::len,
		NULL,
		NULL,
		&NativeBuffer::item
	};

	fb_type->tp_as_sequence = &NativeBuffer_as_sequence;
	fb_type->tp_as_buffer = &NativeBuffer_buffer_procs;

	#if PY_MAJOR_VERSION < 3
		fb_type->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
	#endif

	BP_CLASS(FileResolver, Object, bp::init<>())
		.def("getPathCount", &FileResolver::getPathCount)
		.def("getPath", &FileResolver::getPath, BP_RETURN_VALUE)
		.def("resolve", &FileResolver::resolve, BP_RETURN_VALUE)
		.def("resolveAll", &fileresolver_resolveAll)
		.def("resolveAbsolute", &FileResolver::resolveAbsolute, BP_RETURN_VALUE)
		.def("clone", &FileResolver::clone, BP_RETURN_VALUE)
		.def("appendPath", &FileResolver::appendPath)
		.def("prependPath", &FileResolver::prependPath)
		.def("clear", &FileResolver::clear);

	void (Random::*random_seed_random)(Random *) = &Random::seed;
	void (Random::*random_seed_uint64_t)(uint64_t) = &Random::seed;
	BP_CLASS(Random, Object, bp::init<>())
		.def(bp::init<uint64_t>())
		.def(bp::init<Random *>())
		.def(bp::init<Stream *, InstanceManager *>())
		.def("set", &Random::set)
		.def("seed", random_seed_random)
		.def("seed", random_seed_uint64_t)
		.def("nextULong", &Random::nextULong)
		.def("nextUInt", &Random::nextUInt)
		.def("nextSize", &Random::nextSize)
		.def("nextFloat", &Random::nextFloat)
		.def("nextStandardNormal", &Random::nextStandardNormal)
		.def("serialize", &Random::serialize);

	BP_CLASS(PluginManager, Object, bp::no_init)
		.def("ensurePluginLoaded", &PluginManager::ensurePluginLoaded)
		.def("getLoadedPlugins", &PluginManager::getLoadedPlugins)
		.def("create", pluginmgr_create, BP_RETURN_VALUE)
		.def("createObject", pluginmgr_createobject_1, BP_RETURN_VALUE)
		.def("createObject", pluginmgr_createobject_2, BP_RETURN_VALUE)
		.def("getInstance", &PluginManager::getInstance, BP_RETURN_VALUE)
		.staticmethod("getInstance");

	BP_CLASS(Statistics, Object, bp::no_init)
		.def("getStats", &Statistics::getStats, BP_RETURN_VALUE)
		.def("resetAll", &Statistics::resetAll)
		.def("printStats", &Statistics::printStats)
		.def("getInstance", &Statistics::getInstance, BP_RETURN_VALUE)
		.staticmethod("getInstance");

	BP_CLASS(WorkUnit, Object, bp::no_init)
		.def("set", &WorkUnit::set)
		.def("load", &WorkUnit::load)
		.def("save", &WorkUnit::save);

	BP_CLASS(WorkResult, Object, bp::no_init)
		.def("load", &WorkResult::load)
		.def("save", &WorkResult::save);

	BP_CLASS(WorkProcessor, SerializableObject, bp::no_init)
		.def("createWorkUnit", &WorkProcessor::createWorkUnit, BP_RETURN_VALUE)
		.def("createWorkResult", &WorkProcessor::createWorkResult, BP_RETURN_VALUE)
		.def("clone", &WorkProcessor::clone, BP_RETURN_VALUE)
		.def("prepare", &WorkProcessor::prepare)
		.def("process", &WorkProcessor::process);

	BP_CLASS(ParallelProcess, Object, bp::no_init)
		.def("generateWork", &ParallelProcess::generateWork)
		.def("processResult", &ParallelProcess::processResult)
		.def("handleCancellation", &ParallelProcess::handleCancellation)
		.def("getReturnStatus", &ParallelProcess::getReturnStatus)
		.def("createWorkProcessor", &ParallelProcess::createWorkProcessor)
		.def("bindResource", &ParallelProcess::bindResource)
		.def("isLocal", &ParallelProcess::isLocal)
		.def("getLogLevel", &ParallelProcess::getLogLevel)
		.def("getRequiredPlugins", &ParallelProcess::getRequiredPlugins, BP_RETURN_VALUE);

	BP_SETSCOPE(ParallelProcess_class);
	bp::enum_<ParallelProcess::EStatus>("EStatus")
		.value("EUnknown", ParallelProcess::EUnknown)
		.value("EPause", ParallelProcess::EPause)
		.value("ESuccess", ParallelProcess::ESuccess)
		.value("EFailure", ParallelProcess::EFailure)
		.export_values();
	BP_SETSCOPE(coreModule);

	BP_CLASS(Worker, Thread, bp::no_init)
		.def("getCoreCount", &Worker::getCoreCount)
		.def("isRemoteWorker", &Worker::isRemoteWorker);

	BP_CLASS(LocalWorker, Worker, (bp::init<int, const std::string>()))
		.def(bp::init<int, const std::string, Thread::EThreadPriority>());

	BP_CLASS(RemoteWorker, Worker, (bp::init<const std::string, Stream *>()))
		.def("getNodeName", &RemoteWorker::getNodeName, BP_RETURN_VALUE);

	bp::class_<SerializableObjectVector>("SerializableObjectVector")
		.def(bp::vector_indexing_suite<SerializableObjectVector>());

	BP_CLASS(Scheduler, Object, bp::no_init)
		.def("schedule", &Scheduler::schedule)
		.def("wait", scheduler_wait)
		.def("cancel", scheduler_cancel)
		.def("registerResource", &Scheduler::registerResource)
		.def("registerMultiResource", &Scheduler::registerMultiResource)
		.def("retainResource", &Scheduler::retainResource)
		.def("unregisterResource", &Scheduler::unregisterResource)
		.def("getResourceID", &Scheduler::getResourceID)
		.def("registerWorker", &Scheduler::registerWorker)
		.def("unregisterWorker", &Scheduler::unregisterWorker)
		.def("getWorkerCount", &Scheduler::getWorkerCount)
		.def("getLocalWorkerCount", &Scheduler::getLocalWorkerCount)
		.def("getWorker", &Scheduler::getWorker, BP_RETURN_VALUE)
		.def("start", &Scheduler::start)
		.def("pause", &Scheduler::pause)
		.def("stop", &Scheduler::stop)
		.def("getCoreCount", &Scheduler::getCoreCount)
		.def("hasLocalWorkers", &Scheduler::hasLocalWorkers)
		.def("hasRemoteWorkers", &Scheduler::hasRemoteWorkers)
		.def("getInstance", &Scheduler::getInstance, BP_RETURN_VALUE)
		.def("isRunning", &Scheduler::isRunning)
		.def("isBusy", &Scheduler::isBusy)
		.staticmethod("getInstance");

	BP_CLASS(AbstractAnimationTrack, Object, bp::no_init)
		.def("getType", &AbstractAnimationTrack::getType)
		.def("setTime", &AbstractAnimationTrack::setTime)
		.def("getTime", &AbstractAnimationTrack::getTime)
		.def("getSize", &AbstractAnimationTrack::getSize)
		.def("clone", &AbstractAnimationTrack::clone, BP_RETURN_VALUE);

	BP_CLASS_DECL(StreamBackend, Thread, (bp::init<const std::string, Scheduler *, const std::string &, Stream *, bool>()));

	IMPLEMENT_ANIMATION_TRACK(FloatTrack);
	IMPLEMENT_ANIMATION_TRACK(VectorTrack);
	IMPLEMENT_ANIMATION_TRACK(PointTrack);
	IMPLEMENT_ANIMATION_TRACK(QuatTrack);

	BP_SETSCOPE(AbstractAnimationTrack_class);
	bp::enum_<AbstractAnimationTrack::EType>("EType")
		.value("EInvalid", AbstractAnimationTrack::EInvalid)
		.value("ETranslationX", AbstractAnimationTrack::ETranslationX)
		.value("ETranslationY", AbstractAnimationTrack::ETranslationY)
		.value("ETranslationZ", AbstractAnimationTrack::ETranslationZ)
		.value("ETranslationXYZ", AbstractAnimationTrack::ETranslationXYZ)
		.value("EScaleX", AbstractAnimationTrack::EScaleX)
		.value("EScaleY", AbstractAnimationTrack::EScaleY)
		.value("EScaleZ", AbstractAnimationTrack::EScaleZ)
		.value("EScaleXYZ", AbstractAnimationTrack::EScaleXYZ)
		.value("ERotationX", AbstractAnimationTrack::ERotationX)
		.value("ERotationY", AbstractAnimationTrack::ERotationY)
		.value("ERotationZ", AbstractAnimationTrack::ERotationZ)
		.value("ERotationQuat", AbstractAnimationTrack::ERotationQuat)
		.export_values();

	BP_SETSCOPE(coreModule);

	AbstractAnimationTrack *(AnimatedTransform::*animatedTransform_getTrack)(size_t) = &AnimatedTransform::getTrack;
	AbstractAnimationTrack *(AnimatedTransform::*animatedTransform_findTrack)(AbstractAnimationTrack::EType) = &AnimatedTransform::findTrack;

	BP_CLASS(AnimatedTransform, Object, (bp::init<Transform>()))
		.def(bp::init<>())
		.def(bp::init<Stream *>())
		.def(bp::init<AnimatedTransform *>())
		.def("getTrackCount", &AnimatedTransform::getTrackCount)
		.def("findTrack", animatedTransform_findTrack, BP_RETURN_VALUE)
		.def("getTrack", animatedTransform_getTrack, BP_RETURN_VALUE)
		.def("addTrack", &AnimatedTransform::addTrack)
		.def("appendTransform", &AnimatedTransform::appendTransform)
		.def("isStatic", &AnimatedTransform::isStatic)
		.def("sortAndSimplify", &AnimatedTransform::sortAndSimplify)
		.def("serialize", &AnimatedTransform::serialize)
		.def("getTranslationBounds", &AnimatedTransform::getTranslationBounds, BP_RETURN_VALUE)
		.def("getSpatialBounds", &AnimatedTransform::getSpatialBounds, BP_RETURN_VALUE)
		.def("eval", &AnimatedTransform::eval, BP_RETURN_VALUE);

	BP_STRUCT(Color3, bp::init<>())
		.def(bp::init<Float>())
		.def(bp::init<Float, Float, Float>())
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self -= bp::self)
		.def(bp::self *= Float())
		.def(bp::self * Float())
		.def(bp::self *= bp::self)
		.def(bp::self * bp::self)
		.def(bp::self / Float())
		.def(bp::self /= Float())
		.def(bp::self /= bp::self)
		.def(bp::self / bp::self)
		.def("isValid", &Color3::isValid)
		.def("isNaN", &Color3::isNaN)
		.def("average", &Color3::average)
		.def("abs", &Color3::abs)
		.def("sqrt", &Color3::sqrt)
		.def("safe_sqrt", &Color3::safe_sqrt)
		.def("exp", &Color3::exp)
		.def("log", &Color3::log)
		.def("pow", &Color3::pow)
		.def("clampNegative", &Color3::clampNegative)
		.def("min", &Color3::min)
		.def("max", &Color3::max)
		.def("isZero", &Color3::isZero)
		.def("getLuminance", &Color3::getLuminance)
		.def("__repr__", &Color3::toString)
		.def("__len__", &SpectrumWrapper<Color3>::len)
		.def("__getitem__", &SpectrumWrapper<Color3>::get)
		.def("__setitem__", &SpectrumWrapper<Color3>::set);

	BP_STRUCT(Spectrum, bp::init<>())
		.def("__init__", bp::make_constructor(spectrum_array_constructor))
		.def(bp::init<Float>())
		.def(bp::init<Stream *>())
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self -= bp::self)
		.def(bp::self *= Float())
		.def(bp::self * Float())
		.def(bp::self *= bp::self)
		.def(bp::self * bp::self)
		.def(bp::self / Float())
		.def(bp::self /= Float())
		.def(bp::self /= bp::self)
		.def(bp::self / bp::self)
		.def("isValid", &Spectrum::isValid)
		.def("isNaN", &Spectrum::isNaN)
		.def("average", &Spectrum::average)
		.def("abs", &Spectrum::abs)
		.def("sqrt", &Spectrum::sqrt)
		.def("safe_sqrt", &Spectrum::safe_sqrt)
		.def("exp", &Spectrum::exp)
		.def("log", &Spectrum::log)
		.def("pow", &Spectrum::pow)
		.def("clampNegative", &Spectrum::clampNegative)
		.def("min", &Spectrum::min)
		.def("max", &Spectrum::max)
		.def("isZero", &Spectrum::isZero)
		.def("eval", &Spectrum::eval)
		.def("getLuminance", &Spectrum::getLuminance)
		.def("fromXYZ", &Spectrum::fromXYZ, fromXYZ_overloads())
		.def("toXYZ", &spectrum_toXYZ)
		.def("fromIPT", &Spectrum::fromIPT, fromIPT_overloads())
		.def("toIPT", &spectrum_toIPT)
		.def("fromLinearRGB", &Spectrum::fromLinearRGB, fromLinearRGB_overloads())
		.def("toLinearRGB", &spectrum_toLinearRGB)
		.def("fromSRGB", &Spectrum::fromSRGB)
		.def("toSRGB", &spectrum_toSRGB)
		.def("fromContinuousSpectrum", &Spectrum::fromContinuousSpectrum)
		.def("serialize", &Spectrum::serialize)
		.def("__repr__", &Spectrum::toString)
		.def("__len__", &SpectrumWrapper<Spectrum>::len)
		.def("__getitem__", &SpectrumWrapper<Spectrum>::get)
		.def("__setitem__", &SpectrumWrapper<Spectrum>::set);

	BP_SETSCOPE(Spectrum_struct);
	bp::enum_<Spectrum::EConversionIntent>("EConversionIntent")
		.value("EReflectance", Spectrum::EReflectance)
		.value("EIlluminant", Spectrum::EIlluminant)
		.export_values();
	BP_SETSCOPE(coreModule);

	bp::class_<Properties> properties("Properties");
	properties
		.def(bp::init<std::string>())
		.def(bp::init<const Properties &>())
		.def("getPluginName", &Properties::getPluginName, BP_RETURN_CONSTREF)
		.def("setPluginName", &Properties::setPluginName)
		.def("getID", &Properties::getID, BP_RETURN_CONSTREF)
		.def("setID", &Properties::setID)
		.def("getType", &Properties::getType)
		.def("getPropertyNames", &Properties::getPropertyNames)
		.def("hasProperty", &Properties::hasProperty)
		.def("wasQueried", &Properties::wasQueried)
		.def("markQueried", &Properties::markQueried)
		.def("removeProperty", &Properties::removeProperty)
		.def("merge", &Properties::merge)
		.def("__getitem__", &properties_wrapper::get)
		.def("__setitem__", &properties_wrapper::set)
		.def("__contains__", &Properties::hasProperty)
		.def("__repr__", &Properties::toString);

	BP_SETSCOPE(properties);
	bp::enum_<Properties::EPropertyType>("EPropertyType")
		.value("EBoolean", Properties::EBoolean)
		.value("EInteger", Properties::EInteger)
		.value("EFloat", Properties::EFloat)
		.value("EPoint", Properties::EPoint)
		.value("ETransform", Properties::ETransform)
		.value("ESpectrum", Properties::ESpectrum)
		.value("EString", Properties::EString)
		.value("EData", Properties::EData)
		.export_values();
	BP_SETSCOPE(coreModule);

	BP_STRUCT(Vector1, bp::init<>())
		.def(bp::init< Float>())
		.def(bp::init<Point1>())
		.def_readwrite("x", &Vector1::x);

	BP_STRUCT(Vector2, bp::init<>())
		.def(bp::init<Float>())
		.def(bp::init<Float, Float>())
		.def(bp::init<Point2>())
		.def_readwrite("x", &Vector2::x)
		.def_readwrite("y", &Vector2::y);

	BP_STRUCT(Vector2i, bp::init<>())
		.def(bp::init<int>())
		.def(bp::init<int, int>())
		.def(bp::init<Point2i>())
		.def_readwrite("x", &Vector2i::x)
		.def_readwrite("y", &Vector2i::y);

	BP_STRUCT(Vector3, bp::init<>())
		.def(bp::init<Float>())
		.def(bp::init<Float, Float, Float>())
		.def(bp::init<Point3>())
		.def(bp::init<Normal>())
		.def_readwrite("x", &Vector3::x)
		.def_readwrite("y", &Vector3::y)
		.def_readwrite("z", &Vector3::z);

	BP_SUBSTRUCT(Normal, Vector3, bp::init<>())
		.def(bp::init<Float, Float, Float>())
		.def(bp::init<Vector>())
		.def_readwrite("x", &Normal::x)
		.def_readwrite("y", &Normal::y)
		.def_readwrite("z", &Normal::z);

	BP_STRUCT(Vector3i, bp::init<>())
		.def(bp::init<int>())
		.def(bp::init<int, int, int>())
		.def(bp::init<Point3i>())
		.def_readwrite("x", &Vector3i::x)
		.def_readwrite("y", &Vector3i::y)
		.def_readwrite("z", &Vector3i::z);

	BP_STRUCT(Vector4, bp::init<>())
		.def(bp::init<Float>())
		.def(bp::init<Float, Float, Float, Float>())
		.def(bp::init<Point4>())
		.def_readwrite("x", &Vector4::x)
		.def_readwrite("y", &Vector4::y)
		.def_readwrite("z", &Vector4::z)
		.def_readwrite("w", &Vector4::w);

	BP_STRUCT(Vector4i, bp::init<>())
		.def(bp::init<int>())
		.def(bp::init<int, int, int, int>())
		.def(bp::init<Point4i>())
		.def_readwrite("x", &Vector4i::x)
		.def_readwrite("y", &Vector4i::y)
		.def_readwrite("z", &Vector4i::z)
		.def_readwrite("w", &Vector4i::w);

	BP_STRUCT(Point1, bp::init<>())
		.def(bp::init<Float>())
		.def(bp::init<Vector1>())
		.def_readwrite("x", &Point1::x);

	BP_STRUCT(Point2, bp::init<>())
		.def(bp::init<Float>())
		.def(bp::init<Float, Float>())
		.def(bp::init<Vector2>())
		.def_readwrite("x", &Point2::x)
		.def_readwrite("y", &Point2::y);

	BP_STRUCT(Point2i, bp::init<>())
		.def(bp::init<int>())
		.def(bp::init<int, int>())
		.def(bp::init<Vector2i>())
		.def_readwrite("x", &Point2i::x)
		.def_readwrite("y", &Point2i::y);

	BP_STRUCT(Point3, bp::init<>())
		.def(bp::init<Float>())
		.def(bp::init<Float, Float, Float>())
		.def(bp::init<Vector3>())
		.def(bp::init<Normal>())
		.def_readwrite("x", &Point3::x)
		.def_readwrite("y", &Point3::y)
		.def_readwrite("z", &Point3::z);

	BP_STRUCT(Point3i, bp::init<>())
		.def(bp::init<int>())
		.def(bp::init<int, int, int>())
		.def(bp::init<Vector3i>())
		.def_readwrite("x", &Point3i::x)
		.def_readwrite("y", &Point3i::y)
		.def_readwrite("z", &Point3i::z);

	BP_STRUCT(Point4, bp::init<>())
		.def(bp::init<Float>())
		.def(bp::init<Float, Float, Float, Float>())
		.def(bp::init<Vector4>())
		.def_readwrite("x", &Point4::x)
		.def_readwrite("y", &Point4::y)
		.def_readwrite("z", &Point4::z)
		.def_readwrite("w", &Point4::w);

	BP_STRUCT(Point4i, bp::init<>())
		.def(bp::init<int>())
		.def(bp::init<int, int, int, int>())
		.def(bp::init<Vector4i>())
		.def_readwrite("x", &Point4i::x)
		.def_readwrite("y", &Point4i::y)
		.def_readwrite("z", &Point4i::z)
		.def_readwrite("w", &Point4i::w);

	BP_IMPLEMENT_VECTOR_OPS(Normal, Float, 3);
	BP_IMPLEMENT_VECTOR_OPS(Vector2i, int, 2);
	BP_IMPLEMENT_VECTOR_OPS(Vector3i, int, 3);
	BP_IMPLEMENT_VECTOR_OPS(Vector4i, int, 3);
	BP_IMPLEMENT_VECTOR_OPS(Vector2, Float, 2);
	BP_IMPLEMENT_VECTOR_OPS(Vector3, Float, 3);
	BP_IMPLEMENT_VECTOR_OPS(Vector4, Float, 3);
	BP_IMPLEMENT_POINT_OPS(Point2i, int, 2);
	BP_IMPLEMENT_POINT_OPS(Point3i, int, 3);
	BP_IMPLEMENT_POINT_OPS(Point4i, int, 3);
	BP_IMPLEMENT_POINT_OPS(Point2, Float, 2);
	BP_IMPLEMENT_POINT_OPS(Point3, Float, 3);
	BP_IMPLEMENT_POINT_OPS(Point4, Float, 3);

	Float (*dot2)(const Vector2 &, const Vector2 &) = &dot;
	Float (*dot3)(const Vector3 &, const Vector3 &) = &dot;
	Float (*dot4)(const Vector4 &, const Vector4 &) = &dot;
	Float (*absDot2)(const Vector2 &, const Vector2 &) = &absDot;
	Float (*absDot3)(const Vector3 &, const Vector3 &) = &absDot;
	Float (*absDot4)(const Vector4 &, const Vector4 &) = &absDot;
	Vector2 (*normalize2)(const Vector2 &) = &normalize;
	Vector3 (*normalize3)(const Vector3 &) = &normalize;
	Vector4 (*normalize4)(const Vector4 &) = &normalize;
	Vector3 (*cross3)(const Vector3 &, const Vector3 &) = &cross;
	bp::def("dot", dot2);
	bp::def("dot", dot3);
	bp::def("dot", dot4);
	bp::def("absDot", absDot2);
	bp::def("absDot", absDot3);
	bp::def("absDot", absDot4);
	bp::def("normalize", normalize2);
	bp::def("normalize", normalize3);
	bp::def("normalize", normalize4);
	bp::def("cross", cross3);

	coreModule.attr("Vector") = coreModule.attr("Vector3");
	coreModule.attr("Point") = coreModule.attr("Point3");

	bp::class_<Matrix4x4>("Matrix4x4", bp::init<>())
		.def(bp::init<Float>())
		.def(bp::init<Stream *>())
		.def(bp::init<Matrix4x4>())
		.def(bp::init<Vector4, Vector4, Vector4, Vector4>())
		.def("__init__", bp::make_constructor(Matrix4x4_fromList))
		.def(bp::init<Stream *>())
		.def("__setitem__", &Matrix4x4_setItem)
		.def("__getitem__", &Matrix4x4_getItem)
		.def("setIdentity", &Matrix4x4::setIdentity)
		.def("setZero", &Matrix4x4::setZero)
		.def("isZero", &Matrix4x4::isZero)
		.def("isIdentity", &Matrix4x4::isIdentity)
		.def("trace", &Matrix4x4::trace)
		.def("frob", &Matrix4x4::frob)
		.def("det", &Matrix4x4::det)
		.def("cholDet", &Matrix4x4::cholDet)
		.def("luDet", &Matrix4x4::luDet)
		.def("chol", &Matrix4x4::chol)
		.def("symEig", &Matrix4x4_symEig)
		.def("lu", &Matrix4x4_lu)
		.def("cholSolve", &Matrix4x4_cholSolve)
		.def("luSolve", &Matrix4x4_luSolve)
		.def("transpose", &Matrix4x4_transpose)
		.def("invert", &Matrix4x4_invert)
		.def("serialize", &Matrix4x4::serialize)
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self + Float())
		.def(bp::self += Float())
		.def(bp::self - bp::self)
		.def(bp::self -= bp::self)
		.def(bp::self - Float())
		.def(bp::self -= Float())
		.def(bp::self * Float())
		.def(Float() * bp::self)
		.def(bp::self *= Float())
		.def(bp::self * bp::self)
		.def(bp::self *= bp::self)
		.def(bp::self / Float())
		.def(bp::self /= Float())
		.def("__repr__", &Matrix4x4::toString);

	bp::class_<DiscreteDistribution>("DiscreteDistribution", bp::init<bp::optional<size_t> >())
		.def("clear", &DiscreteDistribution::clear)
		.def("reserve", &DiscreteDistribution::reserve)
		.def("append", &DiscreteDistribution::append)
		.def("isNormalized", &DiscreteDistribution::isNormalized)
		.def("getSum", &DiscreteDistribution::getSum)
		.def("normalize", &DiscreteDistribution::normalize)
		.def("size", &DiscreteDistribution::size)
		.def("sample", &DiscreteDistribution_sample)
		.def("sampleReuse", &DiscreteDistribution_sampleReuse)
		.def("__getitem__", &DiscreteDistribution_getitem)
		.def("__repr__", &DiscreteDistribution::toString);

	bp::class_<Ray>("Ray", bp::init<>())
		.def(bp::init<Ray &>())
		.def(bp::init<Ray &, Float, Float>())
		.def(bp::init<Point, Vector, Float>())
		.def(bp::init<Point, Vector, Float, Float, Float>())
		.def_readwrite("o", &Ray::o)
		.def_readwrite("d", &Ray::d)
		.def_readwrite("dRcp", &Ray::dRcp)
		.def_readwrite("mint", &Ray::mint)
		.def_readwrite("maxt", &Ray::maxt)
		.def_readwrite("time", &Ray::time)
		.def("setOrigin", &Ray::setOrigin)
		.def("setDirection", &Ray::setDirection)
		.def("setTime", &Ray::setTime)
		.def("eval", &ray_eval, BP_RETURN_VALUE)
		.def("__repr__", &Ray::toString);

	bp::class_<RayDifferential, bp::bases<Ray> >("RayDifferential", bp::init<>())
		.def(bp::init<Ray &>())
		.def(bp::init<RayDifferential &>())
		.def(bp::init<Point, Vector, Float>())
		.def_readwrite("rxOrigin", &RayDifferential::rxOrigin)
		.def_readwrite("ryOrigin", &RayDifferential::ryOrigin)
		.def_readwrite("rxDirection", &RayDifferential::rxDirection)
		.def_readwrite("ryDirection", &RayDifferential::ryDirection)
		.def_readwrite("hasDifferentials", &RayDifferential::hasDifferentials)
		.def("scaleDifferential", &RayDifferential::scaleDifferential)
		.def("__repr__", &RayDifferential::toString);

	bp::class_<BSphere>("BSphere", bp::init<>())
		.def(bp::init<BSphere>())
		.def(bp::init<Point, Float>())
		.def(bp::init<Stream *>())
		.def_readwrite("center", &BSphere::center)
		.def_readwrite("radius", &BSphere::radius)
		.def("isEmpty", &BSphere::isEmpty)
		.def("expandBy", &BSphere::expandBy)
		.def("contains", &BSphere::contains)
		.def(bp::self == bp::self)
		.def(bp::self != bp::self)
		.def("rayIntersect", &bsphere_rayIntersect)
		.def("serialize", &BSphere::serialize)
		.def("__repr__", &BSphere::toString);

	BP_STRUCT_DECL(AABB1, bp::init<>());
	BP_IMPLEMENT_AABB_OPS(AABB1, Point1);

	BP_STRUCT_DECL(AABB2, bp::init<>());
	BP_IMPLEMENT_AABB_OPS(AABB2, Point2);

	typedef TAABB<Point3> AABB3;
	BP_STRUCT_DECL(AABB3, bp::init<>());
	BP_IMPLEMENT_AABB_OPS(AABB3, Point3);

	BP_SUBSTRUCT(AABB, AABB3, bp::init<>())
		.def(bp::init<AABB>())
		.def(bp::init<AABB3>())
		.def(bp::init<Point3>())
		.def(bp::init<Point3, Point3>())
		.def(bp::init<Stream *>())
		.def("rayIntersect", &aabb_rayIntersect)
		.def("rayIntersect", &aabb_rayIntersect2)
		.def("getBSphere", &AABB::getBSphere, BP_RETURN_VALUE);

	BP_STRUCT_DECL(AABB4, bp::init<>());
	BP_IMPLEMENT_AABB_OPS(AABB4, Point4);

	bp::class_<Frame>("Frame", bp::init<>())
		.def(bp::init<Vector, Vector, Normal>())
		.def(bp::init<Vector, Vector, Vector>())
		.def(bp::init<Vector>())
		.def(bp::init<Stream *>())
		.def_readwrite("s", &Frame::s)
		.def_readwrite("t", &Frame::t)
		.def_readwrite("n", &Frame::n)
		.def("serialize", &Frame::serialize)
		.def("toLocal", &Frame::toLocal, BP_RETURN_VALUE)
		.def("toWorld", &Frame::toWorld, BP_RETURN_VALUE)
		.def("cosTheta", &Frame::cosTheta)
		.def("cosTheta2", &Frame::cosTheta2)
		.def("sinTheta", &Frame::sinTheta)
		.def("sinTheta2", &Frame::sinTheta2)
		.def("tanTheta", &Frame::tanTheta)
		.def("tanTheta2", &Frame::tanTheta2)
		.def("sinPhi", &Frame::sinPhi)
		.def("cosPhi", &Frame::cosPhi)
		.def("sinPhi2", &Frame::sinPhi2)
		.def("cosPhi2", &Frame::cosPhi2)
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def("__repr__", &Frame::toString)
		.staticmethod("cosTheta")
		.staticmethod("cosTheta2")
		.staticmethod("sinTheta")
		.staticmethod("sinTheta2")
		.staticmethod("tanTheta")
		.staticmethod("tanTheta2")
		.staticmethod("sinPhi")
		.staticmethod("cosPhi")
		.staticmethod("sinPhi2")
		.staticmethod("cosPhi2");

	bp::class_<Transform>("Transform", bp::init<>())
		.def(bp::init<Stream *>())
		.def(bp::init<const Matrix4x4 &>())
		.def(bp::init<const Matrix4x4 &, const Matrix4x4 &>())
		.def("inverse", &Transform::inverse, BP_RETURN_VALUE)
		.def("getMatrix", &Transform::getMatrix, BP_RETURN_VALUE)
		.def("getInverseMatrix", &Transform::getInverseMatrix, BP_RETURN_VALUE)
		.def("det3x3", &Transform::det3x3)
		.def("hasScale", &Transform::hasScale)
		.def("isIdentity", &Transform::isIdentity)
		.def("serialize", &Transform::serialize)
		.def("__mul__", &transform_mul_transform, BP_RETURN_VALUE)
		.def("__mul__", &transform_mul_point, BP_RETURN_VALUE)
		.def("__mul__", &transform_mul_vector, BP_RETURN_VALUE)
		.def("__mul__", &transform_mul_vector4, BP_RETURN_VALUE)
		.def("__mul__", &transform_mul_normal, BP_RETURN_VALUE)
		.def("__mul__", &transform_mul_ray, BP_RETURN_VALUE)
		.def("__repr__", &Transform::toString)
		.def("translate", &Transform::translate, BP_RETURN_VALUE)
		.def("rotate", &Transform::rotate, BP_RETURN_VALUE)
		.def("scale", &Transform::scale, BP_RETURN_VALUE)
		.def("lookAt", &Transform::lookAt, BP_RETURN_VALUE)
		.def("perspective", &Transform::perspective, BP_RETURN_VALUE)
		.def("orthographic", &Transform::orthographic, BP_RETURN_VALUE)
		.def("glPerspective", &Transform::glPerspective, BP_RETURN_VALUE)
		.def("glFrustum", &Transform::glFrustum, BP_RETURN_VALUE)
		.def("glOrthographic", &transform_glOrthographic1, BP_RETURN_VALUE)
		.def("glOrthographic", &transform_glOrthographic2, BP_RETURN_VALUE)
		.def("fromFrame", &Transform::fromFrame, BP_RETURN_VALUE)
		.staticmethod("translate")
		.staticmethod("rotate")
		.staticmethod("scale")
		.staticmethod("lookAt")
		.staticmethod("perspective")
		.staticmethod("orthographic")
		.staticmethod("glPerspective")
		.staticmethod("glFrustum")
		.staticmethod("glOrthographic")
		.staticmethod("fromFrame");

	Float (*fresnelConductorApprox1)(Float, Float, Float) = &fresnelConductorApprox;
	Float (*fresnelConductorExact1)(Float, Float, Float) = &fresnelConductorExact;
	Spectrum (*fresnelConductorApprox2)(Float, const Spectrum &, const Spectrum &) = &fresnelConductorApprox;
	Spectrum (*fresnelConductorExact2)(Float, const Spectrum &, const Spectrum &) = &fresnelConductorExact;

	/* Functions from util.h */
	bp::def("fresnelDielectric", &fresnelDielectric);
	bp::def("fresnelDielectricExt", &fresnelDielectricExt1);
	bp::def("fresnelDielectricExt", &fresnelDielectricExt2);
	bp::def("fresnelConductorApprox", fresnelConductorApprox1, BP_RETURN_VALUE);
	bp::def("fresnelConductorApprox", fresnelConductorApprox2, BP_RETURN_VALUE);
	bp::def("fresnelConductorExact", fresnelConductorExact1, BP_RETURN_VALUE);
	bp::def("fresnelConductorExact", fresnelConductorExact2, BP_RETURN_VALUE);
	bp::def("fresnelDiffuseReflectance", &fresnelDiffuseReflectance);
	bp::def("reflect", &reflect);
	bp::def("refract", &refract1);
	bp::def("refract", &refract2);
	bp::def("coordinateSystem", &mkCoordinateSystem);
	bp::def("memString", &memString1);
	bp::def("memString", &memString2);
	bp::def("timeString", &timeString1);
	bp::def("timeString", &timeString2);
	bp::def("getCoreCount", &getCoreCount);
	bp::def("getHostName", &getHostName);
	bp::def("getPrivateMemoryUsage", &getPrivateMemoryUsage);
	bp::def("getTotalSystemMemory", &getTotalSystemMemory);
	bp::def("getFQDN", &getFQDN);
	bp::def("rdtsc", &rdtsc);

	/* Functions from qmc.h */
	bp::def("radicalInverse2Single", radicalInverse2Single);
	bp::def("radicalInverse2Double", radicalInverse2Double);
	bp::def("radicalInverse2", radicalInverse2Double);
	bp::def("sobol2Single", sobol2Single);
	bp::def("sobol2Double", sobol2Double);
	bp::def("sobol2", sobol2Double);
	bp::def("sample02Single", sample02Single);
	bp::def("sample02Double", sample02Double);
	bp::def("sample02", sample02);
	bp::def("sampleTEA", sampleTEA);
	bp::def("sampleTEAFloat", sampleTEAFloat);
	bp::def("radicalInverse", radicalInverse);
	bp::def("radicalInverseFast", radicalInverseFast);
	bp::def("radicalInverseIncremental", radicalInverseIncremental);

	/* Functions from quad.h */
	double (*legendreP1)(int, double) = &legendreP;
	double (*legendreP2)(int, int, double) = &legendreP;

	bp::def("legendreP", legendreP1);
	bp::def("legendreP", legendreP2);
	bp::def("legendrePD", legendrePD_double);
	bp::def("gaussLegendre", gaussLegendre_);
	bp::def("gaussLobatto", gaussLobatto_);

	bp::class_<GaussLobattoIntegrator>("GaussLobattoIntegrator", (bp::init<size_t, bp::optional<Float, Float, bool, bool> >()))
		.def("integrate", GaussLobattoIntegrator_integrate);

	BP_STRUCT(Quaternion, bp::init<>())
		.def(bp::init<Vector, Float>())
		.def(bp::init<Stream *>())
		.def_readwrite("v", &Quaternion::v)
		.def_readwrite("w", &Quaternion::w)
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self -= bp::self)
		.def(bp::self *= Float())
		.def(bp::self * Float())
		.def(bp::self *= bp::self)
		.def(bp::self * bp::self)
		.def(bp::self *= bp::self)
		.def(bp::self * bp::self)
		.def(bp::self / Float())
		.def(bp::self /= Float())
		.def("isIdentity", &Quaternion::isIdentity)
		.def("axis", &Quaternion::axis)
		.def("angle", &Quaternion::angle)
		.def("exp", &Quaternion::exp)
		.def("log", &Quaternion::log)
		.def("toTransform", &Quaternion::toTransform)
		.def("serialize", &Quaternion::serialize)
		.def("fromAxisAngle", &Quaternion::fromAxisAngle)
		.def("fromTransform", &Quaternion::fromTransform)
		.def("fromDirectionPair", &Quaternion::fromDirectionPair)
		.def("fromMatrix", &Quaternion::fromMatrix)
		.def("fromEulerAngles", &Quaternion::fromEulerAngles)
		.def("__repr__", &Quaternion::toString)
		.staticmethod("fromAxisAngle")
		.staticmethod("fromDirectionPair")
		.staticmethod("fromTransform")
		.staticmethod("fromMatrix")
		.staticmethod("fromEulerAngles");

	BP_SETSCOPE(Quaternion_struct);

	bp::enum_<Quaternion::EEulerAngleConvention>("EEulerAngleConvention")
		.value("EEulerXYZ", Quaternion::EEulerXYZ)
		.value("EEulerXZY", Quaternion::EEulerXZY)
		.value("EEulerYXZ", Quaternion::EEulerYXZ)
		.value("EEulerYZX", Quaternion::EEulerYZX)
		.value("EEulerZXY", Quaternion::EEulerZXY)
		.value("EEulerZYX", Quaternion::EEulerZYX)
		.export_values();
	BP_SETSCOPE(coreModule);

	Float (*dotQ)(const Quaternion &, const Quaternion &) = &dot;
	Quaternion (*normalizeQ)(const Quaternion &) = &normalize;
	Quaternion (*slerpQ)(const Quaternion &, const Quaternion &, Float) = &slerp;

	bp::def("dot", dotQ);
	bp::def("normalize", normalizeQ);
	bp::def("slerp", slerpQ);

	BP_CLASS(ReconstructionFilter, ConfigurableObject, bp::no_init)
		.def("eval", &ReconstructionFilter::eval)
		.def("evalDiscretized", &ReconstructionFilter::evalDiscretized)
		.def("getRadius", &ReconstructionFilter::getRadius)
		.def("getBorderSize", &ReconstructionFilter::getBorderSize);

	BP_SETSCOPE(ReconstructionFilter_class);
	bp::enum_<ReconstructionFilter::EBoundaryCondition>("EBoundaryCondition")
		.value("EClamp", ReconstructionFilter::EClamp)
		.value("ERepeat", ReconstructionFilter::ERepeat)
		.value("EMirror", ReconstructionFilter::EMirror)
		.value("EZero", ReconstructionFilter::EZero)
		.value("EOne", ReconstructionFilter::EOne)
		.export_values();
	BP_SETSCOPE(coreModule);

	Float (SHVector::*shvector_eval1)(Float, Float) const = &SHVector::eval;
	Float (SHVector::*shvector_eval2)(const Vector &) const = &SHVector::eval;
	Float (SHVector::*shvector_evalAzimuthallyInvariant1)(Float, Float) const = &SHVector::evalAzimuthallyInvariant;
	Float (SHVector::*shvector_evalAzimuthallyInvariant2)(const Vector &) const = &SHVector::evalAzimuthallyInvariant;

	BP_STRUCT(SHVector, bp::init<>())
		.def(bp::init<int>())
		.def(bp::init<const SHVector &>())
		.def(bp::init<Stream *>())
		.def(bp::self != bp::self)
		.def(bp::self == bp::self)
		.def(-bp::self)
		.def(bp::self + bp::self)
		.def(bp::self += bp::self)
		.def(bp::self - bp::self)
		.def(bp::self -= bp::self)
		.def(bp::self *= Float())
		.def(bp::self * Float())
		.def(bp::self / Float())
		.def(bp::self /= Float())
		.def("getBands", &SHVector::getBands)
		.def("serialize", &SHVector::serialize)
		.def("energy", &SHVector::energy)
		.def("eval", shvector_eval1)
		.def("eval", shvector_eval2)
		.def("evalAzimuthallyInvariant", shvector_evalAzimuthallyInvariant1)
		.def("evalAzimuthallyInvariant", shvector_evalAzimuthallyInvariant2)
		.def("normalize", &SHVector::normalize)
		.def("mu2", &SHVector::mu2)
		.def("findMinimum", &SHVector::findMinimum)
		.def("addOffset", &SHVector::addOffset)
		.def("convolve", &SHVector::convolve)
		.def("__repr__", &SHVector::toString)
		.def("__getitem__", SHVector_getItem)
		.def("__setitem__", SHVector_setItem)
		.def("rotation", &SHVector::rotation)
		.staticmethod("rotation");

	Float (*dotSH)(const SHVector &, const SHVector &) = &mitsuba::dot;
	bp::def("dot", dotSH);

	BP_CLASS(SHSampler, Object, (bp::init<int, int>()))
		.def("warp", &SHSampler::warp);

	BP_STRUCT(SHRotation, bp::init<int>())
		.def("__call__", &SHRotation::operator());

	BP_CLASS(Timer, Object, bp::init<bp::optional<bool> >())
		.def("start", &Timer::start)
		.def("reset", &Timer::reset, reset_overloads())
		.def("stop", &Timer::stop)
		.def("lap", &Timer::lap)
		.def("getNanoseconds", &Timer::getNanoseconds)
		.def("getMicroseconds", &Timer::getMicroseconds)
		.def("getMilliseconds", &Timer::getMilliseconds)
		.def("getSeconds", &Timer::getSeconds)
		.def("getNanosecondsSinceStart", &Timer::getNanosecondsSinceStart)
		.def("getMicrosecondsSinceStart", &Timer::getMicrosecondsSinceStart)
		.def("getMillisecondsSinceStart", &Timer::getMillisecondsSinceStart)
		.def("getSecondsSinceStart", &Timer::getSecondsSinceStart);

	BP_STRUCT(VonMisesFisherDistr, bp::init<Float>())
		.def("getKappa", &VonMisesFisherDistr::getKappa)
		.def("setKappa", &VonMisesFisherDistr::setKappa)
		.def("eval", &VonMisesFisherDistr::eval)
		.def("getMeanCosine", &VonMisesFisherDistr::getMeanCosine)
		.def("sample", &VonMisesFisherDistr::sample, BP_RETURN_VALUE)
		.def("forMeanCosine", &VonMisesFisherDistr::forMeanCosine)
		.def("forMeanLength", &VonMisesFisherDistr::forMeanLength)
		.def("forPeakValue", &VonMisesFisherDistr::forPeakValue)
		.def("convolve", &VonMisesFisherDistr::convolve)
		.def("__repr__", &VonMisesFisherDistr::toString)
		.staticmethod("forMeanCosine")
		.staticmethod("forPeakValue")
		.staticmethod("convolve");

	bp::detail::current_scope = oldScope;
}

BOOST_PYTHON_MODULE(mitsuba) {
	PyEval_InitThreads();

	bp::object package = bp::scope();
	package.attr("__path__") = "mitsuba";

	/* Automatically take care of the framework
	   initialization / shutdown */
	initializeFramework();
	atexit(shutdownFramework);

	export_core();
	export_render();
}

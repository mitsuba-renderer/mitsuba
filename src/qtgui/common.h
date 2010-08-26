#ifndef QTGUI_COMMON_H
#define QTGUI_COMMON_H

#include <mitsuba/core/platform.h>
#include <QtGui>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/vpl.h>
#include <mitsuba/core/bitmap.h>

using namespace mitsuba;

enum EConnectionType {
	ESSHConnection,
	EDirectConnection
};

namespace mitsuba {
	class RemoteWorker;
};

struct ServerConnection {
	EConnectionType type;
	QString hostName, userName, instDir;
	int port;
	RemoteWorker *worker;
	bool isRegistered;

	inline ServerConnection() : worker(NULL), isRegistered(false) { }

	inline bool operator==(const ServerConnection &c) const {
		return type == c.type && hostName == c.hostName 
			&& userName == c.userName && instDir == c.instDir
			&& port == c.port && worker == c.worker
			&& isRegistered == c.isRegistered;
	}

	inline void fromVariant(QList<QVariant> list) {
		type = (EConnectionType) list[0].toInt();
		hostName = list[1].toString();
		port = list[2].toInt();
		if (type == ESSHConnection) {
			userName = list[3].toString();
			instDir = list[4].toString();
		}
	}

	inline QList<QVariant> toVariant() const {
		QList<QVariant> result;
		result.append(type);
		result.append(hostName);
		result.append(port);
		if (type == ESSHConnection) {
			result.append(userName);
			result.append(instDir);
		}
		return result;
	}

	inline QByteArray toByteArray() const {
		QByteArray a;
		QDataStream stream(&a, QIODevice::WriteOnly);
		stream << toVariant();
		return a;
	}
	
	inline void fromByteArray(QByteArray a) {
		QDataStream stream(a);
		QList<QVariant> variant;
		stream >> variant;
		fromVariant(variant);
	}

	bool createWorker(QWidget *parent);
	
	QString toString() const;
};

enum EMode {
	EPreview = 0,
	ERender
};

enum EPreviewMethod {
	EOpenGL = 0,
	EOpenGLSinglePass,
	ERayTraceCoherent,
	ERayTrace
};

enum EToneMappingMethod {
	EGamma = 0,
	EReinhard
};

namespace mitsuba {
	class GPUSync;
	class GPUTexture;
};

struct PreviewQueueEntry {
	int id, vplSampleOffset;
	GPUTexture *buffer;
	GPUSync *sync;

	inline PreviewQueueEntry(int id = 0) 
		: id(id), vplSampleOffset(0), buffer(NULL), sync(NULL) {
	}
};


struct SceneContext {
	/* Scene-related */
	ref<Scene> scene;
	int sceneResID;
	QString fileName;
	QString shortName;
	Float movementScale;
	Vector up;

	/* Rendering/Preview-related */
	RenderJob *renderJob;
	bool cancelled;
	float progress;
	QString eta, progressName;
	ref<Bitmap> framebuffer;
	EMode mode;
	Float gamma, exposure, clamping;
	bool srgb;
	int pathLength, shadowMapResolution;
	EPreviewMethod previewMethod;
	EToneMappingMethod toneMappingMethod;
	QSize windowSize, sizeIncrease;
	Vector2i scrollOffset;
	Float reinhardKey, reinhardBurn;
	bool allowNonDiffuseVPLs;

	/* Preview state */
	std::deque<VPL> vpls;
	PreviewQueueEntry previewBuffer;

	SceneContext() : scene(NULL), sceneResID(-1), renderJob(NULL) {
	}

	/// Detect the path length
	int detectPathLength() const;

	/* Clone a scene */
	SceneContext(SceneContext *ctx);
	~SceneContext();
};

#endif // QTGUI_COMMON_H

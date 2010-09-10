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

#ifndef QTGUI_COMMON_H
#define QTGUI_COMMON_H

#include <mitsuba/core/platform.h>
#include <QtGui>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/vpl.h>
#include <mitsuba/core/bitmap.h>

using namespace mitsuba;

enum EConnectionType {
	ESSHConnection = 0,
	EDirectConnection
};

enum ENavigationMode {
	EFlythroughFixedYaw = 0,
	EFlythrough
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
	EMode mode, cancelMode;
	Float gamma, exposure, clamping;
	bool srgb;
	int pathLength, shadowMapResolution;
	EPreviewMethod previewMethod;
	EToneMappingMethod toneMappingMethod;
	QSize windowSize, sizeIncrease;
	Vector2i scrollOffset;
	Float reinhardKey, reinhardBurn;
	bool diffuseSources;
	bool diffuseReceivers;

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

class NonClosableDialog : public QDialog {
public:
	NonClosableDialog(QWidget *parent) : QDialog(parent) {
	}

	void closeEvent(QCloseEvent *e) {
		e->ignore();
	}
};

#endif // QTGUI_COMMON_H

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

#ifndef QTGUI_COMMON_H
#define QTGUI_COMMON_H

#include <mitsuba/core/platform.h>
#include <QtGui>
#include <QtXml>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/vpl.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/version.h>

using namespace mitsuba;

enum EConnectionType {
	ESSHConnection = 0,
	EDirectConnection
};

enum ENavigationMode {
	EStandard = 0,
	EFlythrough
};

enum ESelectionMode {
	ENothing = 0,
	EShape,
	EScene
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
	EDisabled = 0,
	EOpenGL,
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
	int id;
	size_t vplSampleOffset;
	GPUTexture *buffer;
	GPUSync *sync;

	inline PreviewQueueEntry(int id = 0) 
		: id(id), vplSampleOffset(0), buffer(NULL), sync(NULL) {
	}
};

struct VisualWorkUnit {
	Point2i offset;
	Vector2i size;
	int worker;
};

struct block_comparator : std::binary_function<VisualWorkUnit, VisualWorkUnit, bool> {
	static int compare(const VisualWorkUnit &v1, const VisualWorkUnit &v2) {
		if (v1.offset.x < v2.offset.x) return -1;
		else if (v1.offset.x > v2.offset.x) return 1;
		if (v1.offset.y < v2.offset.y) return -1;
		else if (v1.offset.y > v2.offset.y) return 1;
		if (v1.size.x < v2.size.x) return -1;
		else if (v1.size.x > v2.size.x) return 1;
		if (v1.size.y < v2.size.y) return -1;
		else if (v1.size.y > v2.size.y) return 1;
		return 0;
	}

	bool operator()(const VisualWorkUnit &v1, const VisualWorkUnit &v2) const {
		return compare(v1, v2) < 0;
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
	std::set<VisualWorkUnit, block_comparator> workUnits;
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
	bool showKDTree;
	int shownKDTreeLevel;
	ESelectionMode selectionMode;
	const Shape *selectedShape;
	QDomDocument doc;

	/* Preview state */
	std::deque<VPL> vpls;
	PreviewQueueEntry previewBuffer;

	SceneContext() : scene(NULL), sceneResID(-1), 
		renderJob(NULL), selectionMode(ENothing),
		selectedShape(NULL) { }

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

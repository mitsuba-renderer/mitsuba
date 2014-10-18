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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "common.h"
#include <QtNetwork/QtNetwork>
#include <mitsuba/render/renderjob.h>

#define MAX_RECENT_FILES 10

// amount of time that a rendering thread will wait for the GUI to accept a
// render view refresh request to show intermediate progress (in ms)
#define REFRESH_TIMEOUT 50

namespace Ui {
    class MainWindow;
}

class QConsoleAppender;
class ServerWidget;
class GLWidget;
class PreviewSettingsDlg;

/**
 * Captures progress notifications from mitsuba and helps to transfer
 * them to the Qt event loop.
 */
class QRenderListener : public QObject, public RenderListener {
	Q_OBJECT
public:

	QRenderListener() {
		m_mutex = new Mutex();
		m_cond = new ConditionVariable(m_mutex);
		m_refreshRequest = NULL;
	}

	/// Called when work has begun in a rectangular image region
	inline void workBeginEvent(const RenderJob *job, const RectangularWorkUnit *wu, int worker) {
		emit workBegin(job, wu, worker);
	}

	/// Called when work has finished in a rectangular image region
	inline void workEndEvent(const RenderJob *job, const ImageBlock *wr, bool cancelled) {
		emit workEnd(job, wr);
	}

	/// Called when work has been canceled in a rectangular image region
	inline void workCanceledEvent(const RenderJob *job, const Point2i &offset, const Vector2i &size) {
		emit workCanceled(job, offset, size);
	}

	/// Called when the whole target image has been altered in some way
	inline void refreshEvent(const RenderJob *job) {
		LockGuard lock(m_mutex);

		/* Potentially overwrite a previous refresh request */
		m_refreshRequest = job;

		/* Asynchronously signal the GUI */
		emit refresh();

		/* Wait for the GUI to draw the image (or another
		   thread to overwrite the refresh request) */
		m_cond->wait(REFRESH_TIMEOUT);

		m_refreshRequest = NULL;
	}

	/// Called when a render job has completed successfully or unsuccessfully
	inline void finishJobEvent(const RenderJob *job, bool cancelled) {
		emit jobFinished(job, cancelled);
	}

	inline const RenderJob *acquireRefreshRequest() { m_mutex->lock(); return m_refreshRequest; }
	inline void releaseRefreshRequest() { m_refreshRequest = NULL; m_cond->signal(); m_mutex->unlock(); }

	MTS_DECLARE_CLASS()
signals:
	void workBegin(const RenderJob *job, const RectangularWorkUnit *wu, int worker);
	void workEnd(const RenderJob *job, const ImageBlock *wr);
	void workCanceled(const RenderJob *job, const Point2i &offset, const Vector2i &size);
	void jobFinished(const RenderJob *job, bool cancelled);
	void refresh();

protected:
	virtual ~QRenderListener() { }

private:
	ref<Mutex> m_mutex;
	ref<ConditionVariable> m_cond;
	const RenderJob *m_refreshRequest;
};

class PreviewSettingsDialog;
class LogWidget;

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();
	void loadFile(QString filename, const QString &destFile = "");
	void adjustSize();
	bool isActive();
	bool initWorkersProcessArgv();

protected:
	SceneContext *loadScene(const QString &filename, const QString &destFile = "");
	void resizeEvent(QResizeEvent *event);
    void changeEvent(QEvent *e);
	void updateRecentFileActions();
	void addRecentFile(QString fileName);
	void closeEvent(QCloseEvent *event);
	SceneContext *getContext(const RenderJob *job, bool failIfNotFound = true);
	void drawHLine(SceneContext *ctx, int x1, int y, int x2, const float *color);
	void drawVLine(SceneContext *ctx, int x, int y1, int y2, const float *color);
	void drawVisualWorkUnit(SceneContext *context, const VisualWorkUnit &block);
	void checkForUpdates(bool notifyIfNone = false);
	void saveAs(SceneContext *ctx, const QString &targetFile);
	void refresh(const RenderJob *job);
	void updateCameraMenu();
	QSize sizeHint() const;

signals:
	void updateView();

private slots:
	void on_actionRenderSettings_triggered();
	void on_actionPreviewSettings_triggered();
	void on_actionOpen_triggered();
	void on_actionExit_triggered();
	void on_actionAbout_triggered();
	void on_actionRefresh_triggered();
	void on_actionRender_triggered();
	void on_actionClose_triggered();
	void on_actionStop_triggered();
	void on_actionMagnify_triggered();
	void on_actionCrop_triggered();
	void on_actionResetView_triggered();
	void on_actionShowLog_triggered();
	void on_actionSettings_triggered();
	void on_actionUpdateCheck_triggered();
    void on_actionStartServer_triggered();
	void on_actionSave_triggered();
	void on_actionSaveAs_triggered();
	void on_actionExportImage_triggered();
	void on_actionReferenceManual_triggered();
	void on_actionImport_triggered();
	void on_actionDuplicateTab_triggered();
	void on_actionAdjustSize_triggered();
	void on_actionReportBug_triggered();
	void on_actionFeedback_triggered();
	void on_actionShowKDTree_triggered();
	void on_actionFocusSelected_triggered();
	void on_actionFocusAll_triggered();
	void on_actionSceneDescription_triggered();
	void on_actionEnableCommandLine_triggered();
	void on_actionCopyImage_triggered();
	void on_tabBar_currentChanged(int index);
	bool on_tabBar_tabCloseRequested(int index);
	void on_tabBar_tabMoved(int from, int to);
	void on_tabBar_customContextMenuRequested(const QPoint &pt);
	void on_glView_loadFileRequest(const QString &string);
	void onOpenRecent();
	void onClearRecent();
	void onWorkBegin(const RenderJob *job, const RectangularWorkUnit *wu, int worker);
	void onWorkEnd(const RenderJob *job, const ImageBlock *wr);
	void onWorkCanceled(const RenderJob *job, const Point2i &offset, const Vector2i &size);
	void onRefresh();
	void onJobFinished(const RenderJob *job, bool cancelled);
	void onProgressMessage(const RenderJob *job, const QString &name,
		float progress, const QString &eta);
	void onStatusMessage(const QString &status);
	void onNetworkFinished(QNetworkReply *reply);
	void onServerClosed();
	void updateUI();
	void updateStatus();
	void onPreviewSettingsClose();
	void onOpenDialogClose(int reason);
	void onExportDialogClose(int reason);
	void onSaveAsDialogClose(int reason);
	void onRenderSettingsClose(int reason);
	void onImportDialogClose(int reason);
	void onSceneInformationClose(int reason);
	void onActivateCamera();
	void on_glView_crop(int type, int x=0, int y=0,
		int width=0, int height=0);
	void onSelectionChanged();
	void onSwitchTab(int rel);

private:
	void exportImage(const QString &fileName);
	void saveSceneAs(const QString &fileName);

    Ui::MainWindow *ui;
	QAction *m_actRecent[MAX_RECENT_FILES];
	QAction *m_clearRecent;
	QList<SceneContext *> m_context;
	ref<RenderQueue> m_renderQueue;
	ref<QRenderListener> m_renderListener;
	Thread::EThreadPriority m_workerPriority;
	QList<ServerConnection> m_connections;
	QMutex m_contextMutex;
	PreviewSettingsDialog *m_previewSettingsDialog;
	LogWidget *m_logWidget;
	ServerWidget *m_serverWidget;
	ref<QConsoleAppender> m_consoleAppender;
	QNetworkAccessManager *m_networkManager;
	QNetworkReply *m_networkReply;
	QProgressBar *m_progress;
	QLabel *m_progressLabel;
	QWidget *m_progressWidget;
	QString m_statusMessage;
	QStringList m_searchPaths;
	QString m_nodeName;
	int m_blockSize, m_listenPort;
	bool m_checkForUpdates, m_manualUpdateCheck;
	bool m_activeWindowHack;
	int m_contextIndex;
	SceneContext *m_lastTab;
	std::map<std::string, std::string, SimpleStringOrdering> m_parameters;
#if defined(__OSX__)
	PreviewSettingsDlg *m_previewSettings;
#endif
	QWidget *m_currentChild;
};

#endif // MAINWINDOW_H

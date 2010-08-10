#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "common.h"
#include <QtNetwork>
#include <mitsuba/render/renderjob.h>

#define MAX_RECENT_FILES 10

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
	/// Called when work has begun in a rectangular image region
	inline void workBeginEvent(const RenderJob *job, const RectangularWorkUnit *wu, int worker) {
		emit workBegin(job, wu, worker);
	}

	/// Called when work has finished in a rectangular image region
	inline void workEndEvent(const RenderJob *job, const ImageBlock *wr) {
		emit workEnd(job, wr);
	}

	/// Called when the whole target image has been altered in some way
	inline void refreshEvent(const RenderJob *job) {
		emit refresh(job);
	}

	/// Called when a render job has completed successfully or unsuccessfully
	inline void finishJobEvent(const RenderJob *job, bool cancelled) {
		emit jobFinished(job, cancelled);
	}

	MTS_DECLARE_CLASS()
signals:
	void workBegin(const RenderJob *job, const RectangularWorkUnit *wu, int worker);
	void workEnd(const RenderJob *job, const ImageBlock *wr);
	void refresh(const RenderJob *job);
	void jobFinished(const RenderJob *job, bool cancelled);
protected:
	virtual ~QRenderListener() { }
};

class PreviewSettingsDialog;
class LogWidget;

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();
	void loadFile(QString filename);
	void adjustSize();
	bool isActive();
	void initWorkers();

protected:
	SceneContext *loadScene(const QString &filename);
	void resizeEvent(QResizeEvent *event);
    void changeEvent(QEvent *e);
	void updateRecentFileActions();
	void addRecentFile(QString fileName);
	void closeEvent(QCloseEvent *event);
	SceneContext *getContext(const RenderJob *job, bool failIfNotFound = true);
	void drawHLine(SceneContext *ctx, int x1, int y, int x2, const float *color);
	void drawVLine(SceneContext *ctx, int x, int y1, int y2, const float *color);
	void checkForUpdates(bool notifyIfNone = false);
	void saveAs(SceneContext *ctx, const QString &targetFile);
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
	void on_actionShowLog_triggered();
	void on_actionSettings_triggered();
	void on_actionUpdateCheck_triggered();
    void on_actionStartServer_triggered();
	void on_actionSave_triggered();
	void on_actionSaveAs_triggered();
	void on_actionExportImage_triggered();
	void on_actionNavigationControls_triggered();
	void on_actionReferenceManual_triggered();
	void on_actionImport_triggered();
	void on_actionDuplicateTab_triggered();
	void on_actionAdjustSize_triggered();
	void on_tabBar_currentChanged(int index);
	bool on_tabBar_tabCloseRequested(int index);
	void on_tabBar_tabMoved(int from, int to);
	void on_tabBar_customContextMenuRequested(const QPoint &pt);
	void on_glView_loadFileRequest(const QString &string);
	void onOpenRecent();
	void onClearRecent();
	void onWorkBegin(const RenderJob *job, const RectangularWorkUnit *wu, int worker);
	void onWorkEnd(const RenderJob *job, const ImageBlock *wr);
	void onRefresh(const RenderJob *job);
	void onJobFinished(const RenderJob *job, bool cancelled);
	void onProgressMessage(const RenderJob *job, const QString &name, 
		float progress, const QString &eta);
	void onStatusMessage(const QString &status);
	void onNetworkFinished(QNetworkReply *reply);
	void onServerClosed();
	void onBugReportError();
	void onBugReportSubmitted();
	void updateUI();
	void updateStatus();
	void onPreviewSettingsClose();
	void onOpenDialogClose(int reason);
	void onSaveAsDialogClose(int reason);
	void onRenderSettingsClose(int reason);

private:
    Ui::MainWindow *ui;
	QAction *m_actRecent[MAX_RECENT_FILES];
	QAction *m_clearRecent;
	QList<SceneContext *> m_context;
	ref<RenderQueue> m_renderQueue;
	ref<QRenderListener> m_renderListener;
	QList<ServerConnection> m_connections;
	QMutex m_contextMutex;
	PreviewSettingsDialog *m_previewSettingsDialog;
	LogWidget *m_logWidget;
	ServerWidget *m_serverWidget;
	ref<QConsoleAppender> m_consoleAppender;
	QNetworkAccessManager *m_networkManager;
	QProgressBar *m_progress;
	QLabel *m_progressLabel;
	QWidget *m_progressWidget;
	QString m_statusMessage;
	QStringList m_searchPaths;
	QString m_nodeName;
	int m_blockSize, m_listenPort;
	bool m_checkForUpdates, m_notifyIfNoUpdates;
	bool m_activeWindowHack;
	int m_bugStatus, m_contextIndex;
	SceneContext *m_lastTab;
#if defined(__OSX__)
	PreviewSettingsDlg *m_previewSettings;
#endif
	QWidget *m_currentChild;
};

#endif // MAINWINDOW_H

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

#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "rendersettingsdlg.h"
#include "previewsettingsdlg.h"
#include "programsettingsdlg.h"
#include "sceneinfodlg.h"
#include "sceneloader.h"
#include "logwidget.h"
#include "aboutdlg.h"
#include "importdlg.h"
#include "loaddlg.h"
#include "updatedlg.h"
#include "server.h"
#include "save.h"
#include "upgrade.h"
#include <QtNetwork/QtNetwork>
#include <mitsuba/core/sched_remote.h>
#include <mitsuba/core/sstream.h>
#include <mitsuba/core/sshstream.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>

#if !defined(WIN32)
#include <pwd.h>
#endif

#if defined(__OSX__)
#include "previewsettingsdlg_cocoa.h"
#else
#include "previewsettingsdlg.h"
#endif

extern bool create_symlinks();

static int remoteWorkerCtr = 0;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent), ui(new Ui::MainWindow),
    m_networkReply(NULL), m_activeWindowHack(false) {
    Logger *logger = Thread::getThread()->getLogger();

#if defined(__OSX__)
    m_previewSettings = NULL;
#endif

    QSettings settings;
    logger->setLogLevel((ELogLevel) settings.value("verbosity", EDebug).toInt());

    m_logWidget = new LogWidget(NULL);

    m_consoleAppender = new QConsoleAppender();
    logger->addAppender(m_consoleAppender);
    connect(m_consoleAppender, SIGNAL(textMessage(ELogLevel, const QString &)),
        m_logWidget, SLOT(onTextMessage(ELogLevel, const QString &)), Qt::QueuedConnection);
    connect(m_consoleAppender, SIGNAL(criticalError(const QString &)),
        m_logWidget, SLOT(onCriticalError(const QString &)), Qt::QueuedConnection);

    SLog(EInfo, "Mitsuba version %s, Copyright (c) " MTS_YEAR " Wenzel Jakob",
        Version(MTS_VERSION).toStringComplete().c_str());

#if SPECTRUM_SAMPLES == 3
    SLog(EInfo, "Configured for RGB-based rendering");
#else
    std::ostringstream oss;
    oss << std::fixed;
    oss.precision(2);
    for (int i=0; i<SPECTRUM_SAMPLES; i++) {
        std::pair<Float, Float> bin = Spectrum::getBinCoverage(i);
        if (i < SPECTRUM_SAMPLES) {
            oss << bin.first << "-" << bin.second << " nm";
            if (i+1<SPECTRUM_SAMPLES)
                oss << ", ";
        }
    }

    SLog(EInfo, "Configured for spectral rendering using "
        "the wavelength discretization %s", oss.str().c_str());
#endif

    m_currentChild = NULL;
    ui->setupUi(this);
    isActive();
    m_lastTab = NULL;
    ui->glView->setScrollBars(ui->hScrollBar, ui->vScrollBar);
    ui->hScrollBar->setVisible(false);
    ui->vScrollBar->setVisible(false);
    ui->actionUpdateCheck->setMenuRole(QAction::ApplicationSpecificRole);
    ui->actionFeedback->setMenuRole(QAction::ApplicationSpecificRole);
    ui->actionReportBug->setMenuRole(QAction::ApplicationSpecificRole);
    ui->actionEnableCommandLine->setMenuRole(QAction::ApplicationSpecificRole);

#if !defined(__OSX__)
    ui->actionEnableCommandLine->setVisible(false);
#endif

    m_progressWidget = new QWidget(centralWidget());
    m_progressLabel = new QLabel(m_progressWidget);
    m_progress = new QProgressBar(m_progressWidget);
    m_progressWidget->setObjectName("progressWidget");
    QHBoxLayout *hlayout = new QHBoxLayout();
    hlayout->addWidget(m_progressLabel);
    hlayout->addWidget(m_progress);
    m_progressWidget->setLayout(hlayout);

#if defined(__OSX__)
    m_progressWidget->setStyleSheet("QWidget#progressWidget {background: qlineargradient(x1:0, y1:0, x2:0, y2:1, "
        "stop:0 palette(dark), stop: 1 palette(mid)); border-top: 1px solid palette(mid); margin: 0px; spacing: 15px; }");
    m_progress->setMaximumSize(SHRT_MAX,18);
    m_progress->setMinimumSize(10,18);
    hlayout->setContentsMargins(10, 8, 10, 3);
    m_progress->setAttribute(Qt::WA_MacSmallSize, true);
    m_progressLabel->setAttribute(Qt::WA_MacSmallSize, true);
#endif

    m_serverWidget = NULL;
    m_contextIndex = -1;

    for (int i = 0; i < MAX_RECENT_FILES; ++i) {
        m_actRecent[i] = new QAction(this);
        m_actRecent[i]->setVisible(false);
        connect(m_actRecent[i], SIGNAL(triggered()),
            this, SLOT(onOpenRecent()));
        ui->menuOpen_Recent->addAction(m_actRecent[i]);
    }

    m_clearRecent = new QAction(tr("Clear Menu"), this);
    connect(m_clearRecent, SIGNAL(triggered()),
        this, SLOT(onClearRecent()));
    ui->menuOpen_Recent->addAction(m_clearRecent);

#if defined(__OSX__)
    /* Make this the default menu bar */
    ui->menuBar->setParent(NULL);
#endif

    updateRecentFileActions();

    ui->tabBar->setDocumentMode(true);
    ui->tabBar->setElideMode(Qt::ElideRight);
    ui->tabBar->setUsesScrollButtons(true);
    ui->tabBar->setTabsClosable(true);
    ui->tabBar->setMovable(true);
    ui->tabBar->setContextMenuPolicy(Qt::CustomContextMenu);
    on_tabBar_currentChanged(-1);

    connect(ui->glView, SIGNAL(quit()), this, SLOT(on_actionExit_triggered()));
    connect(ui->glView, SIGNAL(beginRendering()), this, SLOT(on_actionRender_triggered()));
    connect(ui->glView, SIGNAL(stopRendering()), this, SLOT(on_actionStop_triggered()));
    connect(ui->glView, SIGNAL(statusMessage(const QString &)), this, SLOT(onStatusMessage(const QString &)),
        Qt::QueuedConnection);
    connect(ui->glView, SIGNAL(selectionChanged()), this, SLOT(onSelectionChanged()));
    connect(ui->glView, SIGNAL(switchTab(int)), this, SLOT(onSwitchTab(int)));

    /* Load defaults from app settings file */
    ui->glView->setInvertMouse(settings.value("invertMouse", false).toBool());
    ui->glView->setMouseSensitivity(settings.value("mouseSensitivity", 3).toInt());
    ui->glView->setNavigationMode((ENavigationMode) settings.value("navMode",
        EStandard).toInt());
    m_searchPaths = settings.value("searchPaths", QStringList()).toStringList();
    m_blockSize = settings.value("blockSize", 32).toInt();
    m_listenPort = settings.value("listenPort", MTS_DEFAULT_PORT).toInt();
    m_nodeName = settings.value("nodeName", getFQDN().c_str()).toString();

    m_renderQueue = new RenderQueue();
    m_renderListener = new QRenderListener();
    m_renderQueue->registerListener(m_renderListener);

    connect(m_renderListener, SIGNAL(jobFinished(const RenderJob *, bool)),
        this, SLOT(onJobFinished(const RenderJob *, bool)), Qt::QueuedConnection);
    connect(m_renderListener, SIGNAL(workEnd(const RenderJob *, const ImageBlock *)),
        this, SLOT(onWorkEnd(const RenderJob *, const ImageBlock *)), Qt::DirectConnection);
    connect(m_renderListener, SIGNAL(workCanceled(const RenderJob *, const Point2i &, const Vector2i &)),
        this, SLOT(onWorkCanceled(const RenderJob *, const Point2i &, const Vector2i &)), Qt::DirectConnection);
    connect(m_renderListener, SIGNAL(workBegin(const RenderJob *, const RectangularWorkUnit *, int)),
        this, SLOT(onWorkBegin(const RenderJob *, const RectangularWorkUnit *, int)), Qt::DirectConnection);
    connect(m_renderListener, SIGNAL(refresh()), this, SLOT(onRefresh()), Qt::QueuedConnection);
    connect(m_consoleAppender,
        SIGNAL(progressMessage(const RenderJob *, const QString &, float, const QString &)),
        this, SLOT(onProgressMessage(const RenderJob *, const QString &, float, const QString &)),
        Qt::QueuedConnection);
    connect(this, SIGNAL(updateView()), ui->glView, SLOT(onUpdateView()));

    QPoint windowPos;
    if (settings.contains("pos")) {
        windowPos = settings.value("pos").toPoint();
    } else {
        QDesktopWidget *desktop = QApplication::desktop();
        QRect geo = desktop->screenGeometry();
        windowPos = QPoint(
            geo.left() + (geo.width() - width()) / 2,
            geo.top() + (geo.height() - height())/2
        );
    }

#if defined(__OSX__)
    ui->toolBar->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
    QToolButton *previewButton = static_cast<QToolButton *>(ui->toolBar->widgetForAction(ui->actionPreviewSettings));
    previewButton->setStyleSheet("margin-left: -5px; margin-right:-5px");

    /* Weird Qt/OSX bug -- moving while a window while it is invisible causes
       it to appear move up by 65 pixels (this is related to the unified toolbar) */
    move(windowPos + QPoint(0, 65));
#else
    move(windowPos);
#endif
    show();
    /* Move again just to be sure ... */
    move(windowPos);

    updateUI();
    adjustSize();

    m_networkManager = new QNetworkAccessManager(this);

    if (ui->glView->isUsingSoftwareFallback() &&
        !ui->glView->getErrorString().isEmpty()) {
        QMessageBox::warning(this, tr("Insufficient OpenGL capabilities"),
            ui->glView->getErrorString(), QMessageBox::Ok);
    }

    connect(m_networkManager, SIGNAL(finished(QNetworkReply *)), this, SLOT(onNetworkFinished(QNetworkReply *)));
    m_checkForUpdates = settings.value("checkForUpdates", true).toBool();
    if (m_checkForUpdates)
        checkForUpdates(false);
}

MainWindow::~MainWindow() {
    m_renderQueue->unregisterListener(m_renderListener);
    ref<Scheduler> scheduler = Scheduler::getInstance();
    if (scheduler->isRunning())
        scheduler->pause();
    for (int i=0; i<m_connections.size(); ++i) {
        ServerConnection &c = m_connections[i];
        scheduler->unregisterWorker(c.worker);
    }
    if (m_networkReply)
        m_networkReply->abort();
#if defined(__OSX__)
    delete ui->menuBar;
#endif
    delete ui;
}

bool MainWindow::initWorkersProcessArgv() {
    QSettings settings;
    QStringList args = qApp->arguments();
    std::vector<std::pair<QString, QString> > toBeLoaded;
    int localWorkerCount = settings.value("localWorkers", getCoreCount()).toInt();
    ref<Scheduler> scheduler = Scheduler::getInstance();
    QString outputFile("");
    bool beginRendering = false;

    for (int i=1; i<args.count(); ++i) {
        if (args[i].startsWith("-D")) {
            QString value = args[i].mid(2);
            if (value.length() == 0 && i+1<args.count())
                value = args[++i];
            QStringList list = value.split("=");
            if (list.length() == 2)
                m_parameters[list[0].toStdString()] = list[1].toStdString();
        } else if (args[i].startsWith("-a")) {
            ref<FileResolver> resolver = Thread::getThread()->getFileResolver();
            QString value = args[i].mid(2);
            if (value.length() == 0 && i+1<args.count())
                value = args[++i];
            resolver->appendPath(value.toStdString());
        } else if (args[i].startsWith("-o")) {
            outputFile = args[i].mid(2);
            if (outputFile.length() == 0 && i+1<args.count())
                outputFile = args[++i];
        } else if (args[i].startsWith("-p")) {
            QString value = args[i].mid(2);
            if (value.length() == 0 && i+1<args.count())
                value = args[++i];
            bool ok = false;
            unsigned int numCores = value.toUInt(&ok, 10);
            if (ok)
                localWorkerCount = numCores;
        } else if (args[i].startsWith("-r")) {
            beginRendering = true;
        } else if (args[i].startsWith("-L")) {
            Logger *logger = Thread::getThread()->getLogger();
            QString value = args[i].mid(2);
            if (value.length() == 0 && i+1<args.count())
                value = args[++i];
            value = value.toLower();
            ELogLevel logLevel = EDebug;
            if (value == "trace")
                logLevel = ETrace;
            else if (value == "debug")
                logLevel = EDebug;
            else if (value == "info")
                logLevel = EInfo;
            else if (value == "warn")
                logLevel = EWarn;
            else if (value == "error")
                logLevel = EError;
            else
                SLog(EError, "Invalid log level!");
            logger->setLogLevel(logLevel);
            settings.setValue("verbosity", logLevel);
        } else if (args[i] == "-h") {
            cout <<  "Mitsuba version " << Version(MTS_VERSION).toStringComplete()
                 << ", Copyright (c) " MTS_YEAR " Wenzel Jakob" << endl;
            cout <<  "Usage: mtsgui [options] <One or more scene XML files>" << endl;
            cout <<  "Options/Arguments:" << endl;
            cout <<  "   -h          Display this help text" << endl << endl;
            cout <<  "   -D key=val  Define a constant, which can referenced as \"$key\" in the scene" << endl << endl;
            cout <<  "   -o fname    Write the output image to the file denoted by \"fname\"" << endl << endl;
            cout <<  "   -a p1;p2;.. Add one or more entries to the resource search path" << endl << endl;
            cout <<  "   -p count    Override the detected number of processors." << endl << endl;
            cout <<  "   -L level    Explicitly specify the log level (trace/debug/info/warn/error)" << endl << endl;
            cout <<  "   -r          Immediately begin rendering" << endl << endl;
            cout <<  " For documentation, please refer to http://www.mitsuba-renderer.org/docs.html" << endl;
            return false;
        } else if (args[i].startsWith("-")) {
            cerr << "Unknown option \"" << args[i].toStdString() << "\"" << endl;
        } else {
            toBeLoaded.push_back(std::make_pair(args[i], outputFile));
            outputFile = QString("");
        }
    }

    m_workerPriority = (Thread::EThreadPriority)
        settings.value("workerPriority", (int) Thread::ELowPriority).toInt();
    bool useCoreAffinity = localWorkerCount == getCoreCount();
    for (int i=0; i<localWorkerCount; ++i)
        scheduler->registerWorker(new LocalWorker(useCoreAffinity ? i : -1,
            formatString("wrk%i", i), m_workerPriority));

    int networkConnections = 0;
    QList<QVariant> connectionData = settings.value("connections").toList();
    if (connectionData.size() > 0) {
        QDialog *dialog = new NonClosableDialog(this);
        dialog->setWindowModality(Qt::WindowModal);
        QVBoxLayout *layout = new QVBoxLayout(dialog);
        QLabel *label = new QLabel(tr("Establishing network connections ..."), dialog);
        label->setAlignment(Qt::AlignCenter);
        layout->addWidget(label);
        QProgressBar *progressBar = new QProgressBar(dialog);
        progressBar->setTextVisible(false);
        dialog->resize(200, 50);
        layout->addWidget(progressBar);
        progressBar->setTextVisible(false);
        // weird, Qt/Win needs this to get a busy indicator
        progressBar->setValue(1);
        progressBar->setRange(0, 0);
        dialog->show();

        for (int i=0; i<connectionData.size(); ++i) {
            ServerConnection conn;
            conn.fromByteArray(connectionData[i].toByteArray());
            for (int j=0; j<10; ++j)
                qApp->processEvents();
            if (conn.createWorker(this)) {
                ++networkConnections;
                conn.isRegistered = true;
                scheduler->registerWorker(conn.worker);
                m_connections.append(conn);
            }
        }

        /* Update, this removes connection failures */
        connectionData.clear();
        for (int i=0; i<m_connections.size(); ++i)
            connectionData.append(m_connections[i].toByteArray());
        settings.setValue("connections", connectionData);

        dialog->hide();
        delete dialog;
    }

    if (networkConnections + localWorkerCount == 0) {
        QMessageBox::warning(this, tr("Scheduler warning"),
            tr("There must be at least one worker thread -- forcing creation of one."),
            QMessageBox::Ok);
        scheduler->registerWorker(new LocalWorker(-1, formatString("wrk%i", 0), m_workerPriority));
    }

    for (size_t i=0; i<toBeLoaded.size(); ++i)
        loadFile(toBeLoaded[i].first, toBeLoaded[i].second);

    scheduler->start();

    if (beginRendering) {
        for (int i=0; i<ui->tabBar->count(); ++i) {
            ui->tabBar->setCurrentIndex(i);
            on_actionRender_triggered();
        }
    }

    raise();
    return true;
}


void MainWindow::adjustSize() {
    /* Like QWidget::adjustSize, but clamps the size to 4/5ths of the screen */
    ensurePolished();
    QSize s = sizeHint();
    QRect screen = QApplication::desktop()->screenGeometry(pos());

    s.setWidth(qMin(s.width(), screen.width()*4/5));
    s.setHeight(qMin(s.height(), screen.height()*4/5));
    layout()->activate();

    if (s.isValid())
        resize(s);
}

void MainWindow::checkForUpdates(bool manualUpdateCheck) {
    m_manualUpdateCheck = manualUpdateCheck;
    m_networkReply = m_networkManager->get(QNetworkRequest(QUrl("http://www.mitsuba-renderer.org/version")));
}

void MainWindow::onNetworkFinished(QNetworkReply *reply) {
    if (reply->error() == QNetworkReply::NoError) {
        try {
            QSettings settings;
            Version remote(QString(reply->readAll()).toStdString());
            Version ignoredVersion(settings.value("ignoredVersion", "0.0.0").toString().toStdString());
            Version local(MTS_VERSION);

            if (local < remote) {
                if (!m_manualUpdateCheck && remote == ignoredVersion)
                    return;
                UpdateDialog *dialog = new UpdateDialog(this, local, remote);
                dialog->setAttribute(Qt::WA_DeleteOnClose);
                dialog->setWindowModality(Qt::WindowModal);
                dialog->show();
            } else if (m_manualUpdateCheck) {
                QMessageBox::information(this, tr("Installed version is current"),
                    QString("<p>You're up to date!</p>"
                        "<p>Mitsuba <b>%1</b> is still the newest version available.</p>")
                        .arg(local.toString().c_str()), QMessageBox::Ok);
            }
        } catch (const std::exception &e) {
            /* Got something weird and couldn't parse the version string --
               very likely due to some Internet cafe login page. Ignore unless
               a manual update check was done. */
            if (m_manualUpdateCheck)
                throw e;
        }
    } else {
        if (m_manualUpdateCheck)
            QMessageBox::warning(this, tr("Unable to determine the newest version"),
                QString("<p>Unable to determine the newest Mitsuba version.</p><p>"
                    "Perhaps you are not connected to the Internet?</p>"),
                    QMessageBox::Ok);
    }
}

void MainWindow::on_actionImport_triggered() {
#if defined(MTS_HAS_COLLADA)
    ref<FileResolver> resolver = Thread::getThread()->getFileResolver();
    ref<FileResolver> newResolver = resolver->clone();
    for (int i=(int) m_searchPaths.size()-1; i>=0; --i)
        newResolver->prependPath(toFsPath(m_searchPaths[i]));

    ImportDialog *dialog = new ImportDialog(this, newResolver);
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    connect(dialog, SIGNAL(finished(int)), this, SLOT(onImportDialogClose(int)));
    m_currentChild = dialog;
    // prevent a tab drawing artifact on Qt/OSX
    m_activeWindowHack = true;
    dialog->show();
    qApp->processEvents();
    m_activeWindowHack = false;
#else
    QMessageBox::critical(this, tr("Importer disabled"),
        tr("The importer is disabled in this build. To use it, you will need "
        "to install COLLADA-DOM and recompile Mitsuba -- please see the "
        "documentation for more details."),
        QMessageBox::Ok);
#endif
}

void MainWindow::onImportDialogClose(int reason) {
    m_currentChild = NULL;
}

void MainWindow::on_actionDuplicateTab_triggered() {
    int currentIndex = ui->tabBar->currentIndex();
    if (m_contextIndex != -1)
        currentIndex = m_contextIndex;
    SceneContext *currentContext = m_context[currentIndex];
    SceneContext *newContext = new SceneContext(currentContext);
    if (currentContext->renderJob)
        newContext->windowSize -= currentContext->sizeIncrease;

    m_contextMutex.lock();
    m_context.append(newContext);
    m_contextMutex.unlock();
    ui->tabBar->addTab(newContext->shortName);

    int index = ui->tabBar->count()-1;
    if (ui->tabBar->currentIndex() != index)
        ui->tabBar->setCurrentIndex(index);
}

void MainWindow::on_actionUpdateCheck_triggered() {
    checkForUpdates(true);
}

void MainWindow::on_actionShowKDTree_triggered() {
    int currentIndex = ui->tabBar->currentIndex();
    if (currentIndex == -1)
        return;
    SceneContext *currentContext = m_context[currentIndex];
    bool checked = ui->actionShowKDTree->isChecked();
    currentContext->showKDTree = checked;
    if (currentContext->previewMethod != EOpenGL)
        ui->glView->setPreviewMethod(EOpenGL);
    else
        ui->glView->resetPreview();
}

void MainWindow::on_actionFocusSelected_triggered() {
    QKeyEvent event(QEvent::KeyPress, Qt::Key_F, Qt::NoModifier);
    ui->glView->keyPressEvent(&event);
}

void MainWindow::onSwitchTab(int rel) {
    int index = ui->tabBar->currentIndex() + rel;
    if (index >= 0 && index < ui->tabBar->count())
        ui->tabBar->setCurrentIndex(index);
}

void MainWindow::on_actionFocusAll_triggered() {
    QKeyEvent event(QEvent::KeyPress, Qt::Key_A, Qt::NoModifier);
    ui->glView->keyPressEvent(&event);
}

void MainWindow::on_actionSceneDescription_triggered() {
    int currentIndex = ui->tabBar->currentIndex();
    if (currentIndex == -1)
        return;
    SceneContext *context= m_context[currentIndex];
    SceneInformationDialog *dialog = new SceneInformationDialog(this,
        context->scene);

    /* Center the dialog */
    QDesktopWidget *desktop = QApplication::desktop();
    QRect geo = desktop->screenGeometry(geometry().center());
    QPoint windowPos(
        geo.left() + (geo.width() - dialog->width()) / 2,
        geo.top() + (geo.height() - dialog->height())/2
    );
    dialog->move(windowPos);

    connect(dialog, SIGNAL(finished(int)), this, SLOT(onSceneInformationClose(int)));
    m_currentChild = dialog;
    // prevent a tab drawing artifact on Qt/OSX
    m_activeWindowHack = true;
    dialog->show();
    qApp->processEvents();
    m_activeWindowHack = false;
}

void MainWindow::onSceneInformationClose(int reason) {
    m_currentChild = NULL;
}

void MainWindow::changeEvent(QEvent *e) {
    QMainWindow::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

SceneContext *MainWindow::getContext(const RenderJob *job, bool failIfNotFound) {
    m_contextMutex.lock();
    for (int i=0; i<m_context.size(); ++i) {
        if (m_context[i]->renderJob == job) {
            m_contextMutex.unlock();
            return m_context[i];
        }
    }
    m_contextMutex.unlock();
    if (failIfNotFound)
        SLog(EError, "Internal error: could not find render context!");
    return NULL;
}

void MainWindow::onProgressMessage(const RenderJob *job, const QString &name,
    float progress, const QString &eta) {
    SceneContext *context = getContext(job, false);
    if (context == NULL)
        return;
    context->eta = eta;
    context->progress = progress;
    context->progressName = name + ": ";
    updateUI();
}

#if MTSGUI_STATIC_QFILEDIALOG

void MainWindow::on_actionOpen_triggered() {
    QSettings settings;
    QStringList fileNames = QFileDialog::getOpenFileNames(this, QString(),
        settings.value("fileDir").toString(),
        tr("All supported formats (*.xml *.exr *.rgbe *.hdr *.pfm *.ppm *.png *.jpg *.jpeg);;"
           "Mitsuba scenes (*.xml);;High dynamic-range images (*.exr *.rgbe *.hdr *.pfm *.ppm);;"
           "Low dynamic-range images (*.png *.jpg *.jpeg)"));

    QStringList::ConstIterator it = fileNames.constBegin();
    if (it != fileNames.constEnd()) {
        QFileInfo info(*it);
        settings.setValue("fileDir", info.absolutePath());
    }
    for ( ; it != fileNames.constEnd(); ++it) {
        loadFile(*it);
    }
}

void MainWindow::onOpenDialogClose(int reason) { /* unused */ }

#else // MTSGUI_STATIC_QFILEDIALOG

void MainWindow::on_actionOpen_triggered() {
    QFileDialog *dialog = new QFileDialog(this, Qt::Sheet);
    dialog->setNameFilter(tr("All supported formats (*.xml *.exr *.rgbe *.hdr *.pfm *.ppm *.png *.jpg *.jpeg);;"
            "Mitsuba scenes (*.xml);;High dynamic-range images (*.exr *.rgbe *.hdr *.pfm *.ppm);;Low "
            "dynamic-range images (*.png *.jpg *.jpeg)"));
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    dialog->setAcceptMode(QFileDialog::AcceptOpen);
    dialog->setViewMode(QFileDialog::Detail);
    dialog->setFileMode(QFileDialog::ExistingFiles);
    dialog->setWindowModality(Qt::WindowModal);
    QSettings settings;
    dialog->restoreState(settings.value("fileDialogState").toByteArray());
    connect(dialog, SIGNAL(finished(int)), this, SLOT(onOpenDialogClose(int)));
    m_currentChild = dialog;
    // prevent a tab drawing artifact on Qt/OSX
    m_activeWindowHack = true;
    dialog->show();
    qApp->processEvents();
    m_activeWindowHack = false;
}

void MainWindow::onOpenDialogClose(int reason) {
    QSettings settings;
    QFileDialog *dialog = static_cast<QFileDialog *>(sender());
    m_currentChild = NULL;
    if (reason == QDialog::Accepted) {
        QStringList fileNames = dialog->selectedFiles();
        settings.setValue("fileDialogState", dialog->saveState());
        for (int i=0; i<fileNames.size(); ++i)
            loadFile(fileNames[i]);
    }
}

#endif // MTSGUI_STATIC_QFILEDIALOG

void MainWindow::on_actionExit_triggered() {
    qApp->closeAllWindows();
}

void MainWindow::onOpenRecent() {
    QAction *action = qobject_cast<QAction *>(sender());
    if (action)
        loadFile(action->data().toString());
}

void MainWindow::onClearRecent() {
    QSettings settings;
    for (int j = 0; j < MAX_RECENT_FILES; ++j)
        m_actRecent[j]->setVisible(false);
    settings.setValue("recentFileList", QStringList());
}

SceneContext *MainWindow::loadScene(const QString &qFileName, const QString &destFile) {
    ref<FileResolver> resolver = Thread::getThread()->getFileResolver();
    ref<Logger> logger = Thread::getThread()->getLogger();
    fs::path filename = resolver->resolve(toFsPath(qFileName));
    fs::path filePath = fs::absolute(filename).parent_path();
    ref<FileResolver> newResolver = resolver->clone();
    for (int i=(int) m_searchPaths.size()-1; i>=0; --i)
        newResolver->prependPath(toFsPath(m_searchPaths[i]));
    newResolver->prependPath(filePath);
    LoadDialog *loaddlg = new LoadDialog(this);
    SceneContext *result = NULL;
    ref<SceneLoader> loadingThread;
    loaddlg->setAttribute(Qt::WA_DeleteOnClose);
    loaddlg->setWindowModality(Qt::ApplicationModal);
    loaddlg->setWindowTitle("Loading ...");
    loaddlg->show();

retry:
    loadingThread = new SceneLoader(newResolver, filename, toFsPath(destFile), m_parameters);
    loadingThread->start();

    while (loadingThread->isRunning()) {
        QCoreApplication::processEvents();
        loadingThread->wait(20);
    }
    loadingThread->join();

    result = loadingThread->getResult();
    if (result == NULL) {
        if (loadingThread->isVersionError()) {
            Version version = loadingThread->getVersion();
            int ret;
            if (version.isValid()) {
                ret = QMessageBox::question(this, tr("Version mismatch -- update scene file?"),
                    QString("The requested scene file is from an older version of Mitsuba "
                        "(%1). To work with version %2, it will need to be updated. If you "
                        "continue, Mitsuba will perform a fully automated upgrade (a "
                        "backup copy will be made).\n\nProceed?")
                        .arg(version.toString().c_str())
                        .arg(MTS_VERSION), QMessageBox::Yes | QMessageBox::Cancel);
            } else {
                QMessageBox box(QMessageBox::Question, tr("Version mismatch -- update scene file?"),
                    (loadingThread->getError() + "\n\nAlternatively, if this file is from version 0.2.1 "
                     "(the last release without explicit version numbers), you can perform a fully "
                     "automated upgrade from this version. A backup copy will be made in this case.").c_str());
                QPushButton *version021Button = box.addButton(tr("Assume version 0.2.1?"), QMessageBox::YesRole);
                box.addButton(tr("Cancel"), QMessageBox::RejectRole);
                ret = box.exec();
                if (box.clickedButton() == version021Button) {
                    version = Version(0, 2, 1);
                    ret = QMessageBox::Yes;
                }
            }
            if (ret == QMessageBox::Yes) {
                loaddlg->expand();

                UpgradeManager upgradeMgr(newResolver);
                try {
                    upgradeMgr.performUpgrade(fromFsPath(filename), version);
                    goto retry;
                } catch (const std::exception &ex) {
                    QMessageBox::critical(this, tr("Unable to update %1").arg(qFileName),
                        QString::fromUtf8(ex.what()), QMessageBox::Ok);
                }
            } else {
                QMessageBox::critical(this, tr("Unable to load %1").arg(qFileName),
                    QString("No upgrade was performed -- giving up."), QMessageBox::Ok);
            }
        } else {
            QMessageBox::critical(this, tr("Unable to load %1").arg(qFileName),
                QString::fromUtf8(loadingThread->getError().c_str()),
                QMessageBox::Ok);
        }
    }
    loaddlg->close();

    return result;
}

void MainWindow::loadFile(QString filename, const QString &destFile) {
    QFileInfo fInfo(filename);
    fInfo.makeAbsolute();

    SceneContext *context = loadScene(filename, destFile);
    if (context == NULL)
        return;
    m_contextMutex.lock();
    m_context.append(context);
    m_contextMutex.unlock();

    addRecentFile(context->fileName);

    ui->tabBar->addTab(context->shortName);

    /* Select the newly loaded scene */
    int index = ui->tabBar->count()-1;
    if (ui->tabBar->currentIndex() != index)
        ui->tabBar->setCurrentIndex(index);
    updateUI();
    adjustSize();
}

void MainWindow::addRecentFile(QString fileName) {
    /* Update recent files list */
    QSettings settings;
    QStringList files = settings.value("recentFileList").toStringList();
    files.removeAll(fileName);
    files.prepend(fileName);

    while (files.size() > MAX_RECENT_FILES)
        files.removeLast();

    settings.setValue("recentFileList", files);
    foreach (QWidget *widget, QApplication::topLevelWidgets()) {
        MainWindow *mainWin = qobject_cast<MainWindow *>(widget);
        if (mainWin)
            mainWin->updateRecentFileActions();
    }
}

void MainWindow::updateRecentFileActions() {
    QSettings settings;
    QStringList files = settings.value("recentFileList").toStringList();

    int numRecentFiles = qMin(files.size(), MAX_RECENT_FILES);

    for (int i = 0; i < numRecentFiles; ++i) {
        m_actRecent[i]->setText(QFileInfo(files[i]).fileName());
        m_actRecent[i]->setData(files[i]);
        m_actRecent[i]->setVisible(true);
    }

    for (int j = numRecentFiles; j < MAX_RECENT_FILES; ++j)
        m_actRecent[j]->setVisible(false);
}

void MainWindow::onSelectionChanged() {
    int index = ui->tabBar->currentIndex();
    SceneContext *context = index != -1 ? m_context[index] : NULL;
    bool isInactiveScene = (context && context->scene) ? context->renderJob == NULL : false;
    ui->actionFocusSelected->setEnabled(isInactiveScene && ui->glView->hasSelection());
}

void MainWindow::updateUI() {
    int index = ui->tabBar->currentIndex();
    bool hasTab = (index != -1);

    SceneContext *context = hasTab ? m_context[index] : NULL;
    bool isRendering = hasTab ? context->renderJob != NULL : false;
    bool isShowingRendering = hasTab ? context->mode == ERender : false;
    bool hasScene = hasTab && context->scene != NULL;
    bool isInactive = hasTab ? context->renderJob == NULL : false;
    bool isInactiveScene = (hasTab && hasScene) ? context->renderJob == NULL : false;
    bool slowFallback = ui->glView->isUsingSoftwareFallback() &&
        !ui->glView->hasFastSoftwareFallback();

    ui->actionStop->setEnabled(isShowingRendering);
    if (isShowingRendering && !isRendering)
        ui->actionStop->setToolTip(tr("Return to the realtime preview"));
    else
        ui->actionStop->setToolTip(tr("Stop rendering"));

    ui->actionRender->setEnabled(isInactiveScene);
    ui->actionRefresh->setEnabled(isInactive);
    ui->actionRenderSettings->setEnabled(isInactiveScene);
    ui->actionSave->setEnabled(hasScene);
    ui->actionSaveAs->setEnabled(hasScene);
    ui->actionExportImage->setEnabled(hasTab);
    ui->actionClose->setEnabled(hasTab);
    ui->actionDuplicateTab->setEnabled(hasTab);
    ui->actionAdjustSize->setEnabled(hasTab);
    ui->actionCopyImage->setEnabled(hasTab);
    ui->actionShowKDTree->setEnabled(hasScene);
    ui->actionShowKDTree->setChecked(hasScene && context->showKDTree);
    ui->actionSceneDescription->setEnabled(hasScene);
    ui->actionResetView->setEnabled(hasScene);
    ui->actionMagnify->setEnabled(hasScene);
    ui->actionCrop->setEnabled(hasScene);
    ui->menuCamera->setEnabled(hasScene);
    ui->actionFocusAll->setEnabled(isInactiveScene);
    ui->actionFocusSelected->setEnabled(isInactiveScene && ui->glView->hasSelection());

#if !defined(__OSX__)
    ui->actionPreviewSettings->setEnabled(!slowFallback && hasTab);
#else
    bool isVisible = m_previewSettings != NULL && m_previewSettings->isVisible();
    ui->actionPreviewSettings->setEnabled(hasTab && !isVisible && !slowFallback);
#endif

    if (isRendering) {
        if (!m_progress->isVisible()) {
            static_cast<QGridLayout *>(centralWidget()->layout())->
                addWidget(m_progressWidget, 3, 0, 1, 3);
            m_progressWidget->show();
        }
        m_progress->setValue(context->progress);
        QString etaString;

        if (context->eta != "")
            etaString = QString(" (ETA: ") + context->eta  + QString(")");

        if (context->progressName != "")
            m_progressLabel->setText(QString("<b>%1</b> %2%%3")
                .arg(context->progressName)
                .arg(context->progress, 0, 'f', 1)
                .arg(etaString));
        else
            m_progressLabel->setText("");
    } else {
        if (m_progress->isVisible()) {
            QLayout *centralLayout = centralWidget()->layout();
            m_progressWidget->hide();
            centralLayout->removeWidget(m_progressWidget);
        }
    }
    centralWidget()->updateGeometry();
    layout()->activate();
}

void MainWindow::resizeEvent(QResizeEvent *event) {
    QMainWindow::resizeEvent(event);
    ui->glView->updateScrollBars();
}

void MainWindow::on_tabBar_customContextMenuRequested(const QPoint &pt) {
    if (pt.isNull())
        return;
    int tabIndex = ui->tabBar->tabAt(pt);
    if (tabIndex == -1)
        return;
    m_contextIndex = tabIndex;
    QMenu menu(this);
    menu.addAction(ui->actionDuplicateTab);
    if (tabIndex == ui->tabBar->currentIndex())
        menu.addAction(ui->actionAdjustSize);
    menu.addAction(ui->actionCopyImage);
    menu.addAction(ui->actionClose);
    menu.exec(ui->tabBar->mapToGlobal(pt));
    m_contextIndex = -1;
}

void MainWindow::on_tabBar_tabMoved(int from, int to) {
    m_contextMutex.lock();
    m_context.move(from, to);
    m_contextMutex.unlock();
}

void MainWindow::on_tabBar_currentChanged(int index) {
    if (m_lastTab != NULL) {
        m_lastTab->windowSize = size();
        m_lastTab->wasRendering = m_lastTab->renderJob != NULL;
    }

    ui->glView->ignoreResizeEvents(true);
    if (ui->tabBar->currentIndex() != -1)
        ui->tabBar->show();
    else
        ui->tabBar->hide();
    ui->glView->ignoreResizeEvents(false);

    if (index != -1) {
        if (m_context[index] != m_lastTab) {
            ui->glView->setScene(m_context[index]);
#if defined(__OSX__)
            if (m_previewSettings)
                m_previewSettings->setContext(m_context[index]);
#endif
        }
    } else {
        ui->glView->setScene(NULL);
#if defined(__OSX__)
        if (m_previewSettings && m_previewSettings->isVisible())
            m_previewSettings->hide();
#endif
    }
    m_statusMessage = "";
    updateStatus();
    updateUI();
    updateCameraMenu();

    if (index != -1) {
        SceneContext *context = m_context[index];
        QSize windowSize = context->windowSize;

        if (windowSize.isValid()) {
            if (context->wasRendering && !context->renderJob)
                windowSize -= context->sizeIncrease;
#if defined(__LINUX__)
            int error = (sizeHint()-windowSize).height();
            if (error > 0 && error <= 5)
                resize(windowSize + QSize(0, error));
            else
                resize(windowSize);
#else
            resize(windowSize);
#endif
        } else {
            adjustSize();
        }

        m_lastTab = context;

    } else {
        adjustSize();
        m_lastTab = NULL;
    }
    ui->glView->setFocus();
}

void MainWindow::updateCameraMenu() {
    ui->menuCamera->clear();
    int currentIndex = ui->tabBar->currentIndex();
    if (currentIndex < 0)
        return;
    SceneContext *context = m_context[currentIndex];
    const Scene *scene = context->scene;
    if (!scene)
        return;

    bool isRendering = context->renderJob != NULL;

    const ref_vector<Sensor> &sensors = scene->getSensors();
    for (size_t i = 0; i < sensors.size(); ++i) {
        const Sensor *sensor = sensors[i].get();
        const Properties &props = sensor->getProperties();
        QAction *act = new QAction(ui->menuCamera);
        act->setCheckable(true);
        act->setChecked(sensor == scene->getSensor());
        act->setEnabled(!isRendering);
        act->setText(formatString("\"%s\" (%s)",
            props.getID().c_str(), props.getPluginName().c_str()).c_str());
        act->setData(qVariantFromValue((void *) sensor));
        ui->menuCamera->addAction(act);
        connect(act, SIGNAL(triggered()), this, SLOT(onActivateCamera()));
    }
}

void MainWindow::on_actionClose_triggered() {
    int index = ui->tabBar->currentIndex();
    if (m_contextIndex != -1)
        index = m_contextIndex;
    on_tabBar_tabCloseRequested(index);
}

void MainWindow::on_actionAdjustSize_triggered() {
    adjustSize();
}

bool MainWindow::on_tabBar_tabCloseRequested(int index) {
    SceneContext *context = m_context[index];
    if (context->renderJob != NULL) {
        QMessageBox box(QMessageBox::Question, tr("Really close?"),
            tr("Rendering of scene \"%1\" is unfinished - all progress "
                "will be lost. Are you sure you want to continue?").arg(context->shortName),
            QMessageBox::Yes | QMessageBox::No, this,
            Qt::Dialog | Qt::MSWindowsFixedSizeDialogHint | Qt::Sheet);
        box.setWindowModality(Qt::WindowModal);
        if (box.exec() == QMessageBox::No)
            return false;
    }
    m_contextMutex.lock();
    m_context.removeAt(index);
    m_contextMutex.unlock();
    ui->tabBar->removeTab(index);
    if (context->renderJob != NULL) {
        context->renderJob->cancel();
        context->renderJob->join();
    }
    ui->glView->makeCurrent();
    delete context;
    return true;
}

void MainWindow::closeEvent(QCloseEvent *event) {
    QSettings settings;
    settings.setValue("pos", pos());
    for (int i=ui->tabBar->count()-1; i>=0; --i) {
        if (!on_tabBar_tabCloseRequested(i)) {
            event->ignore();
            return;
        }
    }
    ui->glView->shutdown();
    QMainWindow::closeEvent(event);
    m_logWidget->hide();
    Logger *logger = Thread::getThread()->getLogger();
    logger->removeAppender(m_consoleAppender);
    event->accept();
}

void MainWindow::on_actionCrop_triggered() {
    ui->glView->startCrop(GLWidget::ECrop);
}

void MainWindow::on_actionMagnify_triggered() {
    ui->glView->startCrop(GLWidget::ECropAndMagnify);
}

void MainWindow::on_glView_crop(int type, int x, int y, int width, int height) {
    int currentIndex = ui->tabBar->currentIndex();
    if (currentIndex < 0)
        return;
    SceneContext *context = m_context[currentIndex];
    Scene *scene = context->scene;

    ref<Sensor> oldSensor = scene->getSensor();
    ref<Film> oldFilm = scene->getFilm();

    Properties sensorProps = oldSensor->getProperties();
    Properties filmProps = oldSensor->getFilm()->getProperties();
    Vector2i oldSize = oldFilm->getSize();
    Point2i oldCropOffset = oldFilm->getCropOffset();
    Vector2i oldCropSize = oldFilm->getCropSize();

    if (type == GLWidget::ECrop) {
        filmProps.setInteger("cropOffsetX", x + oldCropOffset.x, false);
        filmProps.setInteger("cropOffsetY", y + oldCropOffset.y, false);
        filmProps.setInteger("cropWidth", width, false);
        filmProps.setInteger("cropHeight", height, false);
    } else if (type == GLWidget::ECropAndMagnify) {
        int maxOldSize = std::max(oldCropSize.x, oldCropSize.y);
        int maxNewSize = std::max(width, height);
        Float magnification = maxOldSize / (Float) maxNewSize;

        width = math::floorToInt(width*magnification);
        height = math::floorToInt(height*magnification);

        filmProps.setInteger("cropOffsetX", math::floorToInt(magnification
                * (x + oldCropOffset.x)), false);
        filmProps.setInteger("cropOffsetY", math::floorToInt(magnification
                * (y + oldCropOffset.y)), false);
        filmProps.setInteger("cropWidth", width);
        filmProps.setInteger("cropHeight", height);
        filmProps.setInteger("width", math::ceilToInt(oldSize.x*magnification), false);
        filmProps.setInteger("height", math::ceilToInt(oldSize.y*magnification), false);
    } else {
        width = context->originalSize.x;
        height = context->originalSize.y;
        filmProps.removeProperty("cropOffsetX");
        filmProps.removeProperty("cropOffsetY");
        filmProps.removeProperty("cropWidth");
        filmProps.removeProperty("cropHeight");
        filmProps.setInteger("width", width, false);
        filmProps.setInteger("height", height, false);
    }

    if (oldCropOffset == Point2i(0) &&
        Vector2i(width, height) == oldCropSize)
        return;

    on_tabBar_currentChanged(-1);
    ref<PluginManager> pluginMgr = PluginManager::getInstance();
    ref<Sensor> newSensor = static_cast<Sensor *>
        (pluginMgr->createObject(MTS_CLASS(Sensor), sensorProps));
    ref<Film> newFilm = static_cast<Film *> (
        pluginMgr->createObject(MTS_CLASS(Film), filmProps));

    newFilm->addChild(oldFilm->getReconstructionFilter());
    newFilm->configure();
    newSensor->addChild(oldSensor->getSampler());
    newSensor->addChild(newFilm);
    newSensor->setMedium(oldSensor->getMedium());
    newSensor->configure();

    scene->removeSensor(oldSensor);
    scene->addSensor(newSensor);
    scene->setSensor(newSensor);

    context->framebuffer = new Bitmap(Bitmap::ERGBA,
        Bitmap::EFloat32, Vector2i(width, height));
    context->framebuffer->clear();
    context->scrollOffset = Vector2i(0, 0);
    context->mode = EPreview;
    context->windowSize = QSize();

    if (context->previewBuffer.buffer) {
        context->previewBuffer.buffer->decRef();
        context->previewBuffer.buffer = NULL;
        context->previewBuffer.vplSampleOffset = 0;
    }

    qApp->processEvents();
    on_tabBar_currentChanged(currentIndex);
    adjustSize();
}

void MainWindow::on_actionResetView_triggered() {
    on_glView_crop(GLWidget::ENone);
}

bool MainWindow::isActive() {
    if (isActiveWindow() || m_activeWindowHack)
        return true;
    else if (m_currentChild != NULL && m_currentChild->isActiveWindow())
        return true;
#if defined(__OSX__)
    if (m_previewSettings != NULL && m_previewSettings->isActiveWindow())
        return true;
#endif
    return false;
}

void MainWindow::drawHLine(SceneContext *ctx, int x1, int y, int x2, const float *color) {
    y = std::max(0, std::min(y, ctx->framebuffer->getHeight()-1));
    x1 = std::max(0, x1); x2 = std::min(x2, ctx->framebuffer->getWidth()-1);
    float *framebuffer = ctx->framebuffer->getFloat32Data();
    int fbOffset = (x1 + y*ctx->framebuffer->getWidth())*4;
    for (int x=x1; x<=x2; x++) {
        framebuffer[fbOffset] = color[0];
        framebuffer[fbOffset+1] = color[1];
        framebuffer[fbOffset+2] = color[2];
        fbOffset+=4;
    }
}

void MainWindow::drawVLine(SceneContext *ctx, int x, int y1, int y2, const float *color) {
    x = std::max(0, std::min(x, ctx->framebuffer->getWidth()-1));
    y1 = std::max(0, y1); y2 = std::min(y2, ctx->framebuffer->getHeight()-1);
    float *framebuffer = ctx->framebuffer->getFloat32Data();
    int width = ctx->framebuffer->getWidth(), fbOffset = (x + y1*width)*4;
    for (int y=y1; y<=y2; y++) {
        framebuffer[fbOffset] = color[0];
        framebuffer[fbOffset+1] = color[1];
        framebuffer[fbOffset+2] = color[2];
        fbOffset += width * 4;
    }
}

void MainWindow::on_actionRenderSettings_triggered() {
    int currentIndex = ui->tabBar->currentIndex();
    SceneContext *context = m_context[currentIndex];
    RenderSettingsDialog *dialog = new RenderSettingsDialog(this);
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    dialog->setWindowModality(Qt::WindowModal);
    dialog->load(context);
    connect(dialog, SIGNAL(finished(int)), this, SLOT(onRenderSettingsClose(int)));
    m_currentChild = dialog;
    // prevent a tab drawing artifact on Qt/OSX
    m_activeWindowHack = true;
    dialog->show();
    qApp->processEvents();
    m_activeWindowHack = false;
}

void MainWindow::onRenderSettingsClose(int reason) {
    RenderSettingsDialog *dialog = static_cast<RenderSettingsDialog *>(sender());
    int currentIndex = ui->tabBar->currentIndex();
    SceneContext *context = m_context[currentIndex];
    m_currentChild = NULL;
    if (reason == QDialog::Accepted) {
        Scene *scene = context->scene;
        scene->incRef();
        dialog->apply(context);
        ui->glView->setPathLength(context->detectPathLength());
        if (dialog->resolutionHasChanged()) {
            ui->glView->refreshScene();
            updateUI();
            adjustSize();
        }
        scene->decRef();
#if defined(__OSX__)
        if (m_previewSettings)
            m_previewSettings->setContext(context);
#endif
    }
}

void MainWindow::on_actionPreviewSettings_triggered() {
#if !defined(__OSX__)
    SceneContext *context = m_context[ui->tabBar->currentIndex()];
    PreviewSettingsDialog d(this, context, ui->glView->getRendererCapabilities());
    connect(&d, SIGNAL(pathLengthChanged(int)), ui->glView, SLOT(setPathLength(int)));
    connect(&d, SIGNAL(clampingChanged(Float)), ui->glView, SLOT(setClamping(Float)));
    connect(&d, SIGNAL(shadowMapResolutionChanged(int)), ui->glView, SLOT(setShadowMapResolution(int)));
    connect(&d, SIGNAL(gammaChanged(bool, Float)), ui->glView, SLOT(setGamma(bool, Float)));
    connect(&d, SIGNAL(exposureChanged(Float)), ui->glView, SLOT(setExposure(Float)));
    connect(&d, SIGNAL(reinhardKeyChanged(Float)), ui->glView, SLOT(setReinhardKey(Float)));
    connect(&d, SIGNAL(reinhardBurnChanged(Float)), ui->glView, SLOT(setReinhardBurn(Float)));
    connect(&d, SIGNAL(previewMethodChanged(EPreviewMethod)), ui->glView, SLOT(setPreviewMethod(EPreviewMethod)));
    connect(&d, SIGNAL(toneMappingMethodChanged(EToneMappingMethod)), ui->glView, SLOT(setToneMappingMethod(EToneMappingMethod)));
    connect(&d, SIGNAL(diffuseReceiversChanged(bool)), ui->glView, SLOT(setDiffuseReceivers(bool)));
    connect(&d, SIGNAL(diffuseSourcesChanged(bool)), ui->glView, SLOT(setDiffuseSources(bool)));
    d.setPreviewEnabled(!ui->glView->isUsingSoftwareFallback());
    d.setMaximumSize(d.minimumSize());
    d.exec();
    QSettings settings;
    settings.setValue("preview_sRGB", context->srgb);
    settings.setValue("preview_gamma", context->gamma);
    settings.setValue("preview_exposure", context->exposure);
    settings.setValue("preview_shadowMapResolution", context->shadowMapResolution);
    settings.setValue("preview_reinhardKey", context->reinhardKey);
    settings.setValue("preview_reinhardBurn", context->reinhardBurn);
    settings.setValue("preview_clamping", context->clamping);
    settings.setValue("preview_method", context->previewMethod);
    settings.setValue("preview_toneMappingMethod", context->toneMappingMethod);
    settings.setValue("preview_diffuseReceivers", context->diffuseReceivers);
    settings.setValue("preview_diffuseSources", context->diffuseSources);
#else
    if (!m_previewSettings) {
        m_previewSettings = new PreviewSettingsDlg(this);
        connect(m_previewSettings, SIGNAL(pathLengthChanged(int)), ui->glView, SLOT(setPathLength(int)));
        connect(m_previewSettings, SIGNAL(clampingChanged(Float)), ui->glView, SLOT(setClamping(Float)));
        connect(m_previewSettings, SIGNAL(shadowMapResolutionChanged(int)), ui->glView, SLOT(setShadowMapResolution(int)));
        connect(m_previewSettings, SIGNAL(gammaChanged(bool, Float)), ui->glView, SLOT(setGamma(bool, Float)));
        connect(m_previewSettings, SIGNAL(exposureChanged(Float)), ui->glView, SLOT(setExposure(Float)));
        connect(m_previewSettings, SIGNAL(reinhardKeyChanged(Float)), ui->glView, SLOT(setReinhardKey(Float)));
        connect(m_previewSettings, SIGNAL(reinhardBurnChanged(Float)), ui->glView, SLOT(setReinhardBurn(Float)));
        connect(m_previewSettings, SIGNAL(previewMethodChanged(EPreviewMethod)), ui->glView, SLOT(setPreviewMethod(EPreviewMethod)));
        connect(m_previewSettings, SIGNAL(toneMappingMethodChanged(EToneMappingMethod)), ui->glView, SLOT(setToneMappingMethod(EToneMappingMethod)));
        connect(m_previewSettings, SIGNAL(close()), this, SLOT(onPreviewSettingsClose()));
        connect(m_previewSettings, SIGNAL(diffuseReceiversChanged(bool)), ui->glView, SLOT(setDiffuseReceivers(bool)));
        connect(m_previewSettings, SIGNAL(diffuseSourcesChanged(bool)), ui->glView, SLOT(setDiffuseSources(bool)));
        m_previewSettings->setPreviewEnabled(!ui->glView->isUsingSoftwareFallback());
    }
    SceneContext *ctx = NULL;
    if (ui->tabBar->currentIndex() != -1)
        ctx = m_context[ui->tabBar->currentIndex()];
    m_previewSettings->setContext(ctx);
    m_previewSettings->show();
    ui->actionPreviewSettings->setEnabled(false);
#endif
}

void MainWindow::onPreviewSettingsClose() {
    int index = ui->tabBar->currentIndex();
    bool hasTab = index != -1;
    ui->actionPreviewSettings->setEnabled(hasTab);
    if (hasTab) {
        SceneContext *context = m_context[index];
        QSettings settings;
        settings.setValue("preview_sRGB", context->srgb);
        settings.setValue("preview_gamma", context->gamma);
        settings.setValue("preview_exposure", context->exposure);
        settings.setValue("preview_shadowMapResolution", context->shadowMapResolution);
        settings.setValue("preview_reinhardKey", context->reinhardKey);
        settings.setValue("preview_reinhardBurn", context->reinhardBurn);
        settings.setValue("preview_clamping", context->clamping);
        settings.setValue("preview_method", context->previewMethod);
        settings.setValue("preview_toneMappingMethod", context->toneMappingMethod);
        settings.setValue("preview_diffuseSources", context->diffuseSources);
        settings.setValue("preview_diffuseReceivers", context->diffuseReceivers);
    }
}

void MainWindow::on_actionSettings_triggered() {
    Logger *logger = Thread::getThread()->getLogger();
    ref<Scheduler> sched = Scheduler::getInstance();
    std::vector<Worker *> localWorkers;

    if (m_renderQueue->getJobCount() != 0) {
        QMessageBox::warning(this, tr("Rendering in progress"),
            tr("The program settings cannot be changed while a rendering is in progress."),
            QMessageBox::Ok);
        return;
    }

    size_t workerCount = sched->getWorkerCount();
    for (size_t i=0; i<workerCount; ++i) {
        Worker *worker = sched->getWorker((int) i);
        if (worker->getClass()->derivesFrom(MTS_CLASS(LocalWorker)))
            localWorkers.push_back(worker);
    }

    ProgramSettingsDialog d(this);
    d.setWindowModality(Qt::ApplicationModal);
    d.setLogLevel(logger->getLogLevel());
    d.setWorkerPriority(m_workerPriority);
    d.setInvertMouse(ui->glView->getInvertMouse());
    d.setNavigationMode(ui->glView->getNavigationMode());
    d.setMouseSensitivity(ui->glView->getMouseSensitivity());
    d.setBlockSize(m_blockSize);
    d.setSearchPaths(m_searchPaths);
    d.setLocalWorkerCount(localWorkers.size());
    d.setConnections(m_connections);
    d.setCheckForUpdates(m_checkForUpdates);
    d.setListenPort(m_listenPort);
    d.setNodeName(m_nodeName);

    if (d.exec()) {
        QList<ServerConnection> &connections = d.getConnections();
        QList<QVariant> connectionData;
        for (int i=0; i<connections.size(); ++i)
            connectionData.append(connections[i].toByteArray());

        QSettings settings;
        settings.setValue("verbosity", d.getLogLevel());
        settings.setValue("invertMouse", d.getInvertMouse());
        settings.setValue("blockSize", d.getBlockSize());
        settings.setValue("searchPaths", d.getSearchPaths());
        settings.setValue("localWorkers", d.getLocalWorkerCount());
        settings.setValue("connections", connectionData);
        settings.setValue("checkForUpdates", d.getCheckForUpdates());
        settings.setValue("mouseSensitivity", d.getMouseSensitivity());
        settings.setValue("listenPort", d.getListenPort());
        settings.setValue("nodeName", d.getNodeName());
        settings.setValue("navMode", (int) d.getNavigationMode());
        settings.setValue("workerPriority", d.getWorkerPriority());

        logger->setLogLevel(d.getLogLevel());
        ui->glView->setInvertMouse(d.getInvertMouse());
        ui->glView->setMouseSensitivity(d.getMouseSensitivity());
        ui->glView->setNavigationMode(d.getNavigationMode());
        m_blockSize = d.getBlockSize();
        m_searchPaths = d.getSearchPaths();
        m_checkForUpdates = d.getCheckForUpdates();
        m_listenPort = d.getListenPort();
        m_nodeName = d.getNodeName();
        m_workerPriority = d.getWorkerPriority();

        bool localWorkersChanged = (int) localWorkers.size() != d.getLocalWorkerCount();

        if (localWorkersChanged || m_connections != d.getConnections()) {
            ref<Scheduler> sched = Scheduler::getInstance();
            sched->pause();

            if ((int) localWorkers.size() != d.getLocalWorkerCount()) {
                /* Completely remove old workers so that CPU affinities can be reassigned */
                while (!localWorkers.empty()) {
                    Worker *worker = localWorkers.back();
                    sched->unregisterWorker(worker);
                    localWorkers.pop_back();
                }
                int workerCount = std::max(1, d.getLocalWorkerCount());
                bool useCoreAffinity = workerCount == getCoreCount();
                for (int i=0; i<workerCount; ++i) {
                    LocalWorker *worker = new LocalWorker(useCoreAffinity ? i : -1,
                        formatString("wrk%i", i), m_workerPriority);
                    sched->registerWorker(worker);
                    localWorkers.push_back(worker);
                }
            }

            QList<ServerConnection> removeList,
                &newConnections = d.getConnections();
            for (int i=0; i<m_connections.size(); ++i) {
                ServerConnection &c = m_connections[i];
                if (!newConnections.contains(c))
                    removeList.append(c);
            }
            for (int i=0; i<newConnections.size(); ++i) {
                ServerConnection &c = newConnections[i];
                if (!m_connections.contains(c)) {
                    sched->registerWorker(c.worker);
                    c.isRegistered = true;
                    c.worker->decRef();
                    m_connections.append(c);
                }
            }
            for (int i=0; i<removeList.size(); ++i) {
                ServerConnection &c = removeList[i];
                sched->unregisterWorker(c.worker);
                m_connections.removeAll(c);
            }
            sched->start();
        }
    }
}

void MainWindow::on_actionStop_triggered() {
    SceneContext *context = m_context[ui->tabBar->currentIndex()];
    m_contextMutex.lock();
    context->workUnits.clear();
    if (context->renderJob) {
        context->cancelled = true;
        ref<RenderJob> renderJob = context->renderJob;
        m_contextMutex.unlock();
        renderJob->cancel();
    } else if (context->mode == ERender) {
        context->mode = EPreview;
        m_contextMutex.unlock();
    } else {
        m_contextMutex.unlock();
    }
    ui->glView->resumePreview();

    updateUI();
}

void MainWindow::on_actionShowLog_triggered() {
    m_logWidget->show();
}

void MainWindow::on_actionRender_triggered() {
    int index = ui->tabBar->currentIndex();
    if (index == -1)
        return;

    SceneContext *context = m_context[index];
    Scene *scene = context->scene;

    if (context->renderJob != NULL || scene == NULL)
        return;

    scene->setBlockSize(m_blockSize);

    if (m_renderQueue->getJobCount() == 0)
        Statistics::getInstance()->resetAll();

    context->renderJob = new RenderJob("rend", scene, m_renderQueue,
        context->sceneResID, -1, -1, false, true);
    context->cancelMode = ERender;
    if (context->mode != ERender)
        ui->glView->downloadFramebuffer();
    context->cancelled = false;
    context->progress = 0;
    context->eta = "";
    context->progressName = "";
    context->mode = ERender;
    m_statusMessage = "";
    updateUI();
    updateCameraMenu();
    context->sizeIncrease = QSize(0, m_progressWidget->sizeHint().height());

#if defined(__LINUX__)
    /* Workaround: on Linux, a few pixels get lost somehow, which
       otherwise causes an unnecessary scrollbar to appear */
    int error = (sizeHint()-size()).height() - context->sizeIncrease.height();
    if (error > 0 && error <= 5)
        resize(size() + context->sizeIncrease + QSize(0, error));
    else
        resize(size() + context->sizeIncrease);
#else
    resize(size() + context->sizeIncrease);
#endif

    updateStatus();

    context->renderJob->start();
}

void MainWindow::on_actionRefresh_triggered() {
    int index = ui->tabBar->currentIndex();
    if (index < 0)
        return;
    SceneContext *context = m_context[index];
    on_tabBar_currentChanged(-1);
    SceneContext *newContext = loadScene(context->fileName);
    if (newContext == NULL)
        return;
    delete context;
    m_contextMutex.lock();
    m_context[index] = newContext;
    m_contextMutex.unlock();
    qApp->processEvents();
    on_tabBar_currentChanged(index);
}

inline float toSRGB(float value) {
    if (value < 0.0031308f)
        return 12.92f * value;
    return 1.055f * std::pow(value, 0.41666f) - 0.055f;
}

void MainWindow::on_actionCopyImage_triggered() {
    exportImage("__clipboard__");
}

#if MTSGUI_STATIC_QFILEDIALOG

void MainWindow::on_actionExportImage_triggered() {
    QSettings settings;
    const QString fileName = QFileDialog::getSaveFileName(this,
        tr("Export image ..."), settings.value("exportFileDir").toString(),
        tr("All supported formats (*.exr *.hdr *.rgbe *.pfm *.png *.jpg *.jpeg);;"
            "High dynamic range OpenEXR image (*.exr);;"
            "High dynamic range Radiance RGBE image (*.rgbe *.hdr);;"
            "High dynamic range Portable Float Map image (*.pfm);;"
            "High dynamic range Portable Pixel Map image (*.ppm);;"
            "Tonemapped low dynamic range image (*.png *.jpg *.jpeg)"));
    if (!fileName.isEmpty()) {
        QSettings settings;
        settings.setValue("exportFileDir", QFileInfo(fileName).absolutePath());
        exportImage(fileName);
    }
}

void MainWindow::onExportDialogClose(int reason) { /* unused */ }

#else // MTSGUI_STATIC_QFILEDIALOG

void MainWindow::on_actionExportImage_triggered() {
    QFileDialog *dialog = new QFileDialog(this, tr("Export image ..."),
        "", tr("All supported formats (*.exr *.hdr *.rgbe *.pfm *.png *.jpg *.jpeg);;"
               "High dynamic range OpenEXR image (*.exr);;"
               "High dynamic range Radiance RGBE image (*.rgbe *.hdr);;"
               "High dynamic range Portable Float Map image (*.pfm);;"
            "High dynamic range Portable Pixel Map image (*.ppm);;"
               "Tonemapped low dynamic range image (*.png *.jpg *.jpeg)"));

    QSettings settings;
    dialog->setViewMode(QFileDialog::Detail);
    dialog->setAcceptMode(QFileDialog::AcceptSave);

#if defined(__OSX__)
    dialog->setOption(QFileDialog::DontUseNativeDialog, true);
#endif

    dialog->restoreState(settings.value("fileDialogState").toByteArray());
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    dialog->setWindowModality(Qt::WindowModal);
    connect(dialog, SIGNAL(finished(int)), this, SLOT(onExportDialogClose(int)));
    m_currentChild = dialog;
    // prevent a tab drawing artifact on Qt/OSX
    m_activeWindowHack = true;
    dialog->show();
    qApp->processEvents();
    m_activeWindowHack = false;
}

void MainWindow::onExportDialogClose(int reason) {
    QSettings settings;
    QFileDialog *dialog = static_cast<QFileDialog *>(sender());
    m_currentChild = NULL;

    if (reason == QDialog::Accepted) {
        QString fileName = dialog->selectedFiles().value(0);
        settings.setValue("fileDialogState", dialog->saveState());
        exportImage(fileName);
    }
}

#endif // MTSGUI_STATIC_QFILEDIALOG

void MainWindow::exportImage(const QString &fileName) {
    if (!fileName.isEmpty()) {
        Bitmap::EComponentFormat compFormat = Bitmap::EInvalid;
        Bitmap::EFileFormat format;

        if (fileName.endsWith(".exr")) {
            format = Bitmap::EOpenEXR;
        } else if (fileName.endsWith(".png") || fileName == "__clipboard__") {
            format = Bitmap::EPNG;
            compFormat = Bitmap::EUInt8;
        } else if (fileName.endsWith(".hdr") || fileName.endsWith(".rgbe")) {
            format = Bitmap::ERGBE;
        } else if (fileName.endsWith(".jpg") || fileName.endsWith(".jpeg")) {
            format = Bitmap::EJPEG;
            compFormat = Bitmap::EUInt8;
        } else if (fileName.endsWith(".pfm")) {
            format = Bitmap::EPFM;
        } else if (fileName.endsWith(".ppm")) {
            format = Bitmap::EPPM;
            compFormat = Bitmap::EUInt16;
        } else {
            SLog(EError, "Unknown file type -- the filename must end in either"
                " .exr, .rgbe, .hdr, .pfm, .ppm, .png, .jpg, or .jpeg");
            return;
        }

        const int currentIndex = ui->tabBar->currentIndex();
        const SceneContext *ctx = m_context[currentIndex];

        if (ctx->mode == EPreview)
            ui->glView->downloadFramebuffer();

        ref<Bitmap> bitmap = ctx->framebuffer;
        if (compFormat == Bitmap::EUInt8 || compFormat == Bitmap::EUInt16) {
            /* Tonemap the image */
            if (ctx->toneMappingMethod == EReinhard) {
                Float logAvgLuminance = 0, maxLuminance = 0; /* Unused */

                Float burn = (ctx->reinhardBurn + 10) / 20.0f;
                bitmap = bitmap->clone();
                bitmap->tonemapReinhard(logAvgLuminance, maxLuminance,
                    ctx->reinhardKey, burn);
            }

            bitmap = bitmap->convert(Bitmap::ERGB, compFormat,
                ctx->srgb ? (Float) -1 : ctx->gamma,
                ctx->toneMappingMethod == EReinhard
                ? (Float) 1.0f : std::pow((Float) 2, ctx->exposure));
        }

        if (fileName == "__clipboard__") {
            QImage image(bitmap->getWidth(), bitmap->getHeight(), QImage::Format_RGB888);
            size_t scanlineSize = (size_t) bitmap->getWidth() * 3 * sizeof(uint8_t);
            for (int i=0; i<bitmap->getHeight(); ++i)
                memcpy(image.scanLine(i), bitmap->getUInt8Data() + scanlineSize * i, scanlineSize);
            QClipboard *clipboard = QApplication::clipboard();
            clipboard->setPixmap(QPixmap::fromImage(image));
        } else {
            ref<FileStream> fs = new FileStream(toFsPath(fileName),
                FileStream::ETruncReadWrite);
            bitmap->write(format, fs);
        }
    }
}

void MainWindow::on_actionSave_triggered() {
    SceneContext *context = m_context[ui->tabBar->currentIndex()];
    saveScene(this, context, context->fileName);
}

#if MTSGUI_STATIC_QFILEDIALOG

void MainWindow::on_actionSaveAs_triggered() {
    SceneContext *context = m_context[ui->tabBar->currentIndex()];
    QString fileDir = QFileInfo(context->fileName).absolutePath();
    QString fileName = QFileDialog::getSaveFileName(this,
        tr("Save scene as ..."), fileDir, tr("Mitsuba scenes (*.xml)"));
    if (!fileName.isEmpty()) {
        QSettings settings;
        settings.setValue("fileDir", QFileInfo(fileName).absolutePath());
        saveSceneAs(fileName);
    }
}

void MainWindow::onSaveAsDialogClose(int reason) { /* unused */ }

#else // MTSGUI_STATIC_QFILEDIALOG

void MainWindow::on_actionSaveAs_triggered() {
    QFileDialog *dialog = new QFileDialog(this, tr("Save scene as ..."),
        "", tr("Mitsuba scenes (*.xml)"));

    m_currentChild = dialog;
    QSettings settings;
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    dialog->setViewMode(QFileDialog::Detail);
    dialog->setAcceptMode(QFileDialog::AcceptSave);
    dialog->restoreState(settings.value("fileDialogState").toByteArray());
    dialog->setWindowModality(Qt::WindowModal);
    connect(dialog, SIGNAL(finished(int)), this, SLOT(onSaveAsDialogClose(int)));
    m_currentChild = dialog;
    // prevent a tab drawing artifact on Qt/OSX
    m_activeWindowHack = true;
    dialog->show();
    qApp->processEvents();
    m_activeWindowHack = false;
}

void MainWindow::onSaveAsDialogClose(int reason) {
    QSettings settings;
    QFileDialog *dialog = static_cast<QFileDialog *>(sender());
    m_currentChild = NULL;
    if (reason == QDialog::Accepted) {
        QString fileName = dialog->selectedFiles().value(0);
        settings.setValue("fileDialogState", dialog->saveState());
        saveSceneAs(fileName);
    }
}

#endif // MTSGUI_STATIC_QFILEDIALOG

void MainWindow::saveSceneAs(const QString &fileName) {
    if (!fileName.isEmpty()) {
        int currentIndex = ui->tabBar->currentIndex();
        SceneContext *context = m_context[currentIndex];

        saveScene(this, context, fileName);
        fs::path pathName = toFsPath(fileName),
                 complete = fs::absolute(pathName),
                 baseName = pathName.stem();
        context->fileName = fileName;
        context->shortName = QFileInfo(fileName).fileName();
        context->scene->setSourceFile(pathName);
        context->scene->setDestinationFile(baseName);
        ui->tabBar->setTabText(currentIndex, context->shortName);
        addRecentFile(fromFsPath(complete));
    }
}

void MainWindow::on_actionReferenceManual_triggered() {
    QDesktopServices::openUrl(QUrl("http://www.mitsuba-renderer.org/docs.html"));
}

void MainWindow::on_actionAbout_triggered() {
    AboutDialog about(this);
    about.exec();
}

void MainWindow::onJobFinished(const RenderJob *job, bool cancelled) {
    SceneContext *context = getContext(job, false);
    if (context == NULL)
        return;
    m_renderQueue->join();
    context->workUnits.clear();
    if (cancelled) {
        if (!context->cancelled) {
            QMessageBox::critical(this, tr("Error while rendering"),
                tr("The rendering job did not complete successfully. Please check the log."),
                QMessageBox::Ok);
        } else {
            context->mode = context->cancelMode;
            if (ui->tabBar->currentIndex() != -1 &&
                m_context[ui->tabBar->currentIndex()] == context)
                ui->glView->resumePreview();
        }
    }
    refresh(job);
    context->renderJob = NULL;
    updateUI();
    updateCameraMenu();
    if (ui->tabBar->currentIndex() != -1 &&
        m_context[ui->tabBar->currentIndex()] == context)
        resize(size() - context->sizeIncrease);
}

void MainWindow::onStatusMessage(const QString &status) {
    m_statusMessage = status;
    updateStatus();
}

void MainWindow::updateStatus() {
    if (m_statusMessage == "")
        setWindowTitle(tr("Mitsuba renderer"));
    else
        setWindowTitle(tr("Mitsuba renderer [%1]").arg(m_statusMessage));
}

void MainWindow::on_actionStartServer_triggered() {
    m_serverWidget = new ServerWidget(NULL, m_nodeName, m_listenPort);
    ui->actionStartServer->setEnabled(false);
    connect(m_serverWidget, SIGNAL(closed()), this, SLOT(onServerClosed()));
    m_serverWidget->show();
}

void MainWindow::on_actionEnableCommandLine_triggered() {
    if (QMessageBox::question(this, tr("Enable command line access"),
        tr("<p>If you proceed, Mitsuba will create symbolic links in <tt>/usr/bin</tt> and <tt>/Library/Python/{2.6,2.7}/site-packages</tt>, as well as an entry in .bashrc, "
            "which enable command line and Python usage. Note that you will have to "
            "repeat this process every time the Mitsuba application is moved.</p>"
            "<p>Create links?</p>"),
            QMessageBox::Yes | QMessageBox::No) != QMessageBox::Yes)
        return;
    if (!create_symlinks())
        QMessageBox::critical(this, tr("Authentication error"),
            tr("Unable to create the symbolic links!"), QMessageBox::Ok);
}

void MainWindow::on_actionReportBug_triggered() {
    QDesktopServices::openUrl(QUrl("https://www.mitsuba-renderer.org/bugtracker/projects/mitsuba"));
}

void MainWindow::on_actionFeedback_triggered() {
    QDesktopServices::openUrl(QUrl("mailto:Wenzel%20Jakob%20<wenzel@cs.cornell.edu>?subject=Feedback%20on%20Mitsuba"));
}

void MainWindow::onServerClosed() {
    delete m_serverWidget;
    m_serverWidget = NULL;
    ui->actionStartServer->setEnabled(true);
}

void MainWindow::onWorkBegin(const RenderJob *job, const RectangularWorkUnit *wu, int worker) {
    SceneContext *context = getContext(job, false);
    if (context == NULL)
        return;
    VisualWorkUnit vwu(wu->getOffset(), wu->getSize(), worker);
    /* This is not executed in the event loop -- take some precautions */
    m_contextMutex.lock();
    context->workUnits.insert(vwu);
    drawVisualWorkUnit(context, vwu);
    bool isCurrentView = ui->tabBar->currentIndex() < m_context.size() &&
        m_context[ui->tabBar->currentIndex()] == context;
    m_contextMutex.unlock();
    if (isCurrentView)
        emit updateView();
}

void MainWindow::drawVisualWorkUnit(SceneContext *context, const VisualWorkUnit &vwu) {
    int ox = vwu.offset.x, oy = vwu.offset.y,
        ex = ox + vwu.size.x, ey = oy + vwu.size.y;
    if (vwu.size.x < 3 || vwu.size.y < 3)
        return;

    const float *color = NULL;

    /* Use desaturated colors to highlight the host
       responsible for rendering the current image vwu */
    const float white[]     = { 1.0f, 1.0f, 1.0f };
    const float red[]       = { 1.0f, 0.3f, 0.3f };
    const float green[]     = { 0.3f, 1.0f, 0.3f };
    const float blue[]      = { 0.3f, 0.3f, 1.0f };
    const float gray[]      = { 0.5f, 0.5f, 0.5f };
    const float yellow[]    = { 1.0f, 1.0f, 0.0f };
    const float magenta[]   = { 1.0f, 0.3f, 1.0f };
    const float turquoise[] = { 0.3f, 1.0f, 1.0f };

    switch (vwu.worker % 8) {
        case 1: color = green; break;
        case 2: color = yellow; break;
        case 3: color = blue; break;
        case 4: color = gray; break;
        case 5: color = red; break;
        case 6: color = magenta; break;
        case 7: color = turquoise; break;
        case 0:
        default:
            color = white;
            break;
    }

    float scale = .7f * (color[0] * 0.212671f + color[1] * 0.715160f + color[2] * 0.072169f);

    float ncol[3] = {
        color[0]*scale,
        color[1]*scale,
        color[2]*scale
    };

    drawHLine(context, ox, oy, ox + 3, ncol);
    drawHLine(context, ex - 4, oy, ex - 1, ncol);
    drawHLine(context, ox, ey - 1, ox + 3, ncol);
    drawHLine(context, ex - 4, ey - 1, ex - 1, ncol);
    drawVLine(context, ox, oy, oy + 3, ncol);
    drawVLine(context, ex - 1, oy, oy + 3, ncol);
    drawVLine(context, ex - 1, ey - 4, ey - 1, ncol);
    drawVLine(context, ox, ey - 4, ey - 1, ncol);
}

void MainWindow::onWorkCanceled(const RenderJob *job, const Point2i &offset, const Vector2i &size) {
    SceneContext *context = getContext(job, false);
    if (context == NULL)
        return;
    VisualWorkUnit vwu(offset, size);
    m_contextMutex.lock();
    if (context->workUnits.find(vwu) != context->workUnits.end())
        context->workUnits.erase(vwu);
    m_contextMutex.unlock();
}

void MainWindow::onWorkEnd(const RenderJob *job, const ImageBlock *block) {
    SceneContext *context = getContext(job, false);
    if (context == NULL)
        return;

    // int border = block->getBorderSize();
    int border = 0; // produces a more appealing preview

    Point2i offset(
        std::max(0, block->getOffset().x - border),
        std::max(0, block->getOffset().y - border));

    Bitmap *target = context->framebuffer;
    Vector2i size(
        std::min(target->getWidth(),  offset.x + block->getSize().x + 2 * border)-offset.x,
        std::min(target->getHeight(), offset.y + block->getSize().y + 2 * border)-offset.y);

    context->scene->getFilm()->develop(offset, size, offset, target);

    /* This is executed by worker threads -- take some precautions */
    m_contextMutex.lock();
    bool isCurrentView = ui->tabBar->currentIndex() < m_context.size() &&
        m_context[ui->tabBar->currentIndex()] == context;
    VisualWorkUnit vwu(block->getOffset(), block->getSize());
    if (context->workUnits.find(vwu) != context->workUnits.end()) {
        context->workUnits.erase(vwu);
    } else if (!context->cancelled) {
        SLog(EWarn, "Internal error: unable to find previously scheduled"
                " rectangular work unit.");
    }
    for (std::set<VisualWorkUnit, block_comparator>::const_iterator it =
        context->workUnits.begin(); it != context->workUnits.end(); ++it)
        drawVisualWorkUnit(context, *it);
    m_contextMutex.unlock();
    if (isCurrentView)
        emit updateView();
}

void MainWindow::onRefresh() {
    const RenderJob *req = m_renderListener->acquireRefreshRequest();
    if (req)
        refresh(req);
    m_renderListener->releaseRefreshRequest();
}

void MainWindow::refresh(const RenderJob *job) {
    SceneContext *context = getContext(job, false);
    if (context == NULL)
        return;

    Bitmap *target = context->framebuffer;
    context->scene->getFilm()->develop(Point2i(0, 0),
        target->getSize(), Point2i(0, 0), target);

    /* This is executed by worker threads -- take some precautions */
    m_contextMutex.lock();
    bool isCurrentView = ui->tabBar->currentIndex() < m_context.size() &&
        m_context[ui->tabBar->currentIndex()] == context;
    for (std::set<VisualWorkUnit, block_comparator>::const_iterator it =
        context->workUnits.begin(); it != context->workUnits.end(); ++it)
        drawVisualWorkUnit(context, *it);
    m_contextMutex.unlock();
    if (isCurrentView)
        emit updateView();
}

void MainWindow::onActivateCamera() {
    int index = ui->tabBar->currentIndex();
    if (index == -1 || m_context[index]->scene == NULL ||
        m_context[index]->renderJob != NULL)
        return;
    QAction *action = qobject_cast<QAction *>(sender());
    Sensor *sensor = static_cast<Sensor *>(action->data().value<void *>());
    if (!sensor)
        return;
    SceneContext *context = m_context[index];
    on_tabBar_currentChanged(-1);
    Vector2i size = sensor->getFilm()->getSize();
    context->scene->setSensor(sensor);

    context->framebuffer = new Bitmap(Bitmap::ERGBA,
        Bitmap::EFloat32, size);
    context->framebuffer->clear();
    context->scrollOffset = Vector2i(0, 0);
    context->mode = EPreview;
    context->windowSize = QSize();

    if (context->previewBuffer.buffer) {
        context->previewBuffer.buffer->decRef();
        context->previewBuffer.buffer = NULL;
        context->previewBuffer.vplSampleOffset = 0;
    }

    on_tabBar_currentChanged(index);
    updateUI();
    adjustSize();
    updateCameraMenu();
}

void MainWindow::on_glView_loadFileRequest(const QString &string) {
    loadFile(string);
}

bool ServerConnection::createWorker(QWidget *parent) {
    ref<Stream> stream;
    try {
        if (type == EDirectConnection) {
            stream = new SocketStream(hostName.toStdString(), port);
        } else {
            std::vector<std::string> cmdLine;
            cmdLine.push_back(formatString("bash -c 'cd %s; . setpath.sh; mtssrv -ls'", instDir.toLatin1().constData()));
            stream = new SSHStream(userName.toStdString(),
                hostName.toStdString(), cmdLine, port);
        }
        worker = new RemoteWorker(formatString("net%i", remoteWorkerCtr++), stream);
        return true;
    } catch (const std::exception &e) {
        QString extra;
        if (type == ESSHConnection && !QString(e.what()).contains("configuration mismatch"))
            extra = parent->tr(" - Please make sure that you can log in manually using the "
                "command line and that <a href=\"http://www.debian-administration.org/articles/152\">"
                "passwordless authentication</a> is active.");
        QMessageBox::critical(parent, parent->tr("Unable to connect to %1").arg(hostName),
            QString("Unable to create a connection to \"%1\": %2%3")
            .arg(hostName).arg(e.what()).arg(extra), QMessageBox::Ok);
        return false;
    }
}

QSize MainWindow::sizeHint() const {
    QSize hint = QMainWindow::sizeHint();
    /* Don't include scroll bars in the size hint */
    if (ui->hScrollBar->isVisible())
        hint -= QSize(0, ui->hScrollBar->sizeHint().height());
    if (ui->vScrollBar->isVisible())
        hint -= QSize(ui->vScrollBar->sizeHint().width(), 0);
    return hint;
}

QString ServerConnection::toString() const {
    return QString("%1 (%2, %3 cores)")
        .arg(worker->getNodeName().c_str())
        .arg(type == EDirectConnection ? "direct" : "ssh")
        .arg(worker->getCoreCount());
}

static ref<Integrator> recursivelyCloneIntegrator(const Integrator *oldIntegrator) {
    ref<Integrator> integrator = static_cast<Integrator *> (PluginManager::getInstance()->createObject(
        MTS_CLASS(Integrator), oldIntegrator->getProperties()));
    int idx = 0;
    while (oldIntegrator->getSubIntegrator(idx) != 0) {
        ref<Integrator> child = recursivelyCloneIntegrator(oldIntegrator->getSubIntegrator(idx));
        integrator->addChild(child);
        child->setParent(integrator);
        idx++;
    }
    integrator->configure();
    return integrator;
}

SceneContext::SceneContext(SceneContext *ctx) {
    if (ctx->scene) {
        /* Temporarily set up a new file resolver */
        ref<Thread> thread = Thread::getThread();
        ref<FileResolver> oldResolver = thread->getFileResolver();
        ref<FileResolver> newResolver = oldResolver->clone();
        newResolver->prependPath(fs::absolute(ctx->scene->getSourceFile()).parent_path());
        thread->setFileResolver(newResolver);

        scene = new Scene(ctx->scene);
        ref<PluginManager> pluginMgr = PluginManager::getInstance();
        ref<Sensor> oldSensor = static_cast<Sensor *>(ctx->scene->getSensor());
        ref<Sensor> newSensor = static_cast<Sensor *>
            (pluginMgr->createObject(MTS_CLASS(Sensor), oldSensor->getProperties()));
        ref<Sampler> sampler = static_cast<Sampler *>
            (pluginMgr->createObject(MTS_CLASS(Sampler), ctx->scene->getSampler()->getProperties()));
        ref<Film> film = static_cast<Film *>
            (pluginMgr->createObject(MTS_CLASS(Film), oldSensor->getFilm()->getProperties()));
        const Integrator *oldIntegrator = ctx->scene->getIntegrator();
        scene->setIntegrator(recursivelyCloneIntegrator(oldIntegrator));

        ref<ReconstructionFilter> rfilter = static_cast<ReconstructionFilter *>
            (pluginMgr->createObject(MTS_CLASS(ReconstructionFilter), oldSensor->getFilm()->
                getReconstructionFilter()->getProperties()));

        rfilter->configure();
        film->addChild(rfilter);
        film->configure();
        sampler->configure();
        newSensor->addChild(sampler);
        newSensor->setMedium(oldSensor->getMedium());
        newSensor->addChild(film);
        newSensor->configure();
        scene->removeSensor(oldSensor);
        scene->addSensor(newSensor);
        scene->setSensor(newSensor);
        scene->setSampler(sampler);
        scene->configure();
        sceneResID = ctx->sceneResID;
        Scheduler::getInstance()->retainResource(sceneResID);
        thread->setFileResolver(oldResolver);
    } else {
        sceneResID = -1;
        renderJob = NULL;
    }
    fileName = ctx->fileName;
    shortName = ctx->shortName;
    movementScale = ctx->movementScale;
    up = ctx->up;
    renderJob = NULL;
    wasRendering = false;
    cancelled = false;
    progress = 0.0f;
    framebuffer = ctx->framebuffer->clone();
    mode = ctx->renderJob ? EPreview : ctx->mode;
    gamma = ctx->gamma;
    exposure = ctx->exposure;
    clamping = ctx->clamping;
    srgb = ctx->srgb;
    pathLength = ctx->pathLength;
    shadowMapResolution = ctx->shadowMapResolution;
    previewMethod = ctx->previewMethod;
    toneMappingMethod = ctx->toneMappingMethod;
    windowSize = ctx->windowSize;
    sizeIncrease = ctx->sizeIncrease;
    scrollOffset = ctx->scrollOffset;
    reinhardKey = ctx->reinhardKey;
    reinhardBurn = ctx->reinhardBurn;
    diffuseReceivers = ctx->diffuseReceivers;
    diffuseSources = ctx->diffuseSources;
    showKDTree = ctx->showKDTree;
    shownKDTreeLevel = ctx->shownKDTreeLevel;
    selectedShape = ctx->selectedShape;
    selectionMode = ctx->selectionMode;
    originalSize = ctx->originalSize;
    doc = ctx->doc.cloneNode(true).toDocument();
    layers = ctx->layers;
    for (size_t i=0; i<layers.size(); ++i)
        layers[i].second->incRef();
    currentLayer = ctx->currentLayer;
}

SceneContext::~SceneContext() {
    if (scene && sceneResID != -1)
        Scheduler::getInstance()->unregisterResource(sceneResID);
    if (previewBuffer.buffer)
        previewBuffer.buffer->decRef();
    if (previewBuffer.sync)
        previewBuffer.sync->decRef();
    for (size_t i=0; i<layers.size(); ++i)
        layers[i].second->decRef();
}

int SceneContext::detectPathLength() const {
    if (!scene)
        return 2;

    const Integrator *integrator = scene->getIntegrator();
    int extraDepth = 0;

    while (integrator->getSubIntegrator(0) != NULL) {
        if (integrator->getClass()->getName() == "IrradianceCacheIntegrator")
            extraDepth = 1;
        integrator = integrator->getSubIntegrator(0);
    }
    const Properties &integratorProps = integrator->getProperties();
    int maxDepth = -1;

    if (integratorProps.hasProperty("maxDepth"))
        maxDepth = integratorProps.getInteger("maxDepth");

    if (maxDepth == -1) {
        if (integratorProps.getPluginName() == "direct")
            maxDepth = 2;
        else
            maxDepth = 5;
    }

    return std::max(2, std::min(maxDepth + extraDepth, 6));
}

MTS_IMPLEMENT_CLASS(QRenderListener, false, RenderListener)

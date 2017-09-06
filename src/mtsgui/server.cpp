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

#include "server.h"
#include <mitsuba/core/sstream.h>
#include <mitsuba/core/sched_remote.h>

ServerThread::ServerThread(Logger *logger, int listenPort, const std::string &nodeName)
    : Thread("serv"), m_listenPort(listenPort), m_nodeName(nodeName) {
    setLogger(logger);
}

ServerThread::~ServerThread() {
}

void ServerThread::shutdown() {
    m_active = false;
    join();
}

void ServerThread::run() {
    SLog(EInfo, "Mitsuba version %s, Copyright (c) " MTS_YEAR " Wenzel Jakob",
            Version(MTS_VERSION).toStringComplete().c_str());
    /* Allocate a socket of the proper type (IPv4/IPv6) */
    ref<Scheduler> scheduler = Scheduler::getInstance();
    struct addrinfo hints, *servinfo, *p = NULL;
    memset(&hints, 0, sizeof(struct addrinfo));
    hints.ai_family = AF_UNSPEC;
    hints.ai_flags = AI_PASSIVE;
    hints.ai_socktype = SOCK_STREAM;
    char portName[8];
    int rv, one = 1;
    m_socket = INVALID_SOCKET;
    m_active = true;
    bool hostNameSet = m_nodeName != getFQDN();

    snprintf(portName, sizeof(portName), "%i", m_listenPort);
    if ((rv = getaddrinfo(hostNameSet ? m_nodeName.c_str() : NULL, portName, &hints, &servinfo)) != 0)
        SLog(EError, "Error in getaddrinfo(%s:%i): %s", m_nodeName.c_str(), m_listenPort, gai_strerror(rv));

    for (p = servinfo; p != NULL; p = p->ai_next) {
        /* Allocate a socket */
        m_socket = socket(p->ai_family, p->ai_socktype, p->ai_protocol);
        if (m_socket == -1)
            SocketStream::handleError("none", "socket");

        /* Avoid "bind: socket already in use" */
        if (setsockopt(m_socket, SOL_SOCKET, SO_REUSEADDR, (char *) &one, sizeof(int)) < 0)
            SocketStream::handleError("none", "setsockopt");

        /* Bind the socket to the port number */
        if (bind(m_socket, p->ai_addr, (socklen_t) p->ai_addrlen) == -1) {
            SocketStream::handleError("none", formatString("bind(%s:%i)", m_nodeName.c_str(), m_listenPort), EError);
#if defined(__WINDOWS__)
            closesocket(m_socket);
#else
            close(m_socket);
#endif
            continue;
        }
        break;
    }
    if (p == NULL)
        SLog(EError, "Failed to bind to port %i!", m_listenPort);

    freeaddrinfo(servinfo);

    if (listen(m_socket, CONN_BACKLOG) == -1)
        SocketStream::handleError("none", "bind");

    SLog(EInfo, "%s: Listening on port %i.. Close the window to stop.",
        m_nodeName.c_str(), m_listenPort);

    fd_set fds;
    struct timeval tv;
    memset(&tv, 0, sizeof(struct timeval));
    tv.tv_usec = 500 * 1000;

    int connectionIndex = 0;
    /* Wait for connections */
    while (m_active) {
        FD_ZERO(&fds);
        FD_SET(m_socket, &fds);
        int ret = select(m_socket+1, &fds, NULL, NULL, &tv);
        if (ret <= 0 || !FD_ISSET(m_socket, &fds))
            continue;

        socklen_t addrlen = sizeof(sockaddr_storage);
        struct sockaddr_storage sockaddr;
        memset(&sockaddr, 0, addrlen);

        SOCKET newSocket = accept(m_socket, (struct sockaddr *) &sockaddr, &addrlen);

        if (newSocket == INVALID_SOCKET) {
#if defined(__WINDOWS__)
            if (!m_active)
                break;
#else
            if (errno == EINTR)
                continue;
#endif
            SocketStream::handleError("none", "accept", EWarn);
            continue;
        }

        ref<StreamBackend> backend = new StreamBackend(formatString("con%i", connectionIndex++),
                scheduler, m_nodeName, new SocketStream(newSocket), true);
        backend->start();
    }

#if defined(__WINDOWS__)
    closesocket(m_socket);
#else
    close(m_socket);
#endif
}

ServerWidget::ServerWidget(QWidget *parent,
        const QString &nodeName, int listenPort)
        : QMainWindow(parent) {
    m_contents = new QTextEdit(this);
    QFont font = QFontDatabase::systemFont(QFontDatabase::FixedFont);
#if defined(__OSX__)
    font.setPointSize(14);
#endif
    m_contents->setFont(font);
    m_contents->setReadOnly(true);
    QPalette palette;
    palette.setColor(QPalette::Base, Qt::black);
    m_contents->setPalette(palette);
    setPalette(palette);
    QToolBar *toolBar = new QToolBar(this);
    toolBar->setMovable(false);
    toolBar->setAllowedAreas(Qt::TopToolBarArea);
    toolBar->setIconSize(QSize(32, 32));
    toolBar->setToolButtonStyle(Qt::ToolButtonIconOnly);
    toolBar->setFloatable(false);

    QAction *actionClear = new QAction(this);
    QIcon clearIcon;
    clearIcon.addFile(QString::fromUtf8(":/resources/clear.png"), QSize(), QIcon::Normal, QIcon::Off);
    actionClear->setIcon(clearIcon);
    actionClear->setToolTip(tr("Clear"));
    actionClear->setText(tr("Clear"));
    connect(actionClear, SIGNAL(triggered()), m_contents, SLOT(clear()));
    toolBar->addAction(actionClear);

#if defined(__OSX__)
    toolBar->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
#endif

    addToolBar(Qt::TopToolBarArea, toolBar);
    setCentralWidget(m_contents);
    setUnifiedTitleAndToolBarOnMac(true);
    setMenuBar(NULL);
    setWindowTitle(tr("Server Log"));
    resize(QSize(1000, 500));

    QConsoleAppender *consoleAppender = new QConsoleAppender();

    m_logger = new Logger(Thread::getThread()->getLogger()->getLogLevel());
    m_logger->addAppender(consoleAppender);
    m_logger->setFormatter(new DefaultFormatter());
    connect(consoleAppender, SIGNAL(textMessage(ELogLevel, const QString &)),
        this, SLOT(onTextMessage(ELogLevel, const QString &)), Qt::QueuedConnection);

    m_thread = new ServerThread(m_logger, listenPort, nodeName.toStdString());
    m_thread->start();
}

ServerWidget::~ServerWidget() {
}

void ServerWidget::onTextMessage(ELogLevel level, const QString &message) {
    QColor color;
    switch (level) {
        case ETrace:
        case EDebug:
        case EInfo:
            color = Qt::gray;
            break;
        case EError:
            color = Qt::red;
            break;
        case EWarn:
        default:
            color = QColor(255,180,180);
            break;
    }
    m_contents->setTextColor(color);
    m_contents->append(message);
    QTextCursor c =  m_contents->textCursor();
    c.movePosition(QTextCursor::End);
    m_contents->setTextCursor(c);
    m_contents->ensureCursorVisible();
}

void ServerWidget::show() {
    if (isVisible()) {
        raise();
    } else {
        /* Center the dialog */
        QDesktopWidget *desktop = QApplication::desktop();
        QRect geo = desktop->screenGeometry();
        QPoint windowPos(
            geo.left() + (geo.width() - width()) / 2,
            geo.top() + (geo.height() - height())/2
        );
        move(windowPos);
        QMainWindow::show();
    }
}

void ServerWidget::closeEvent(QCloseEvent *event) {
    QMainWindow::closeEvent(event);
    m_thread->shutdown();
    event->accept();
    emit closed();
}


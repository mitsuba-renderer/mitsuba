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

#if !defined(__SERVERWIDGET_H)
#define __SERVERWIDGET_H

#include "logwidget.h"

#if defined(__WINDOWS__)
#include <io.h>
#include <ws2tcpip.h>
#else
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>

#define INVALID_SOCKET -1
#define SOCKET int
#endif

/* How many clients are allowed to wait for a connection at a time */
#define CONN_BACKLOG 5

namespace mitsuba {
    class RenderJob;
};

class ServerThread : public Thread {
public:
    ServerThread(Logger *logger, int listenPort,
        const std::string &nodeName);
    virtual ~ServerThread();

    void run();
    void shutdown();
private:
    SOCKET m_socket;
    int m_listenPort;
    bool m_active;
    std::string m_nodeName;
};

class ServerWidget : public QMainWindow {
    Q_OBJECT
public:
    ServerWidget(QWidget *parent,
        const QString &nodeName, int listenPort);
    virtual ~ServerWidget();
    void show();
signals:
    void closed();
protected slots:
    void onTextMessage(ELogLevel level, const QString &message);
protected:
    void closeEvent(QCloseEvent *event);
private:
    QTextEdit *m_contents;
    ref<Logger> m_logger;
    ref<ServerThread> m_thread;
};

#endif /* __SERVERWIDGET_H */

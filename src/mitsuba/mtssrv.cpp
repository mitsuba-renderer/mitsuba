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

#include <mitsuba/core/sched_remote.h>
#include <mitsuba/core/cstream.h>
#include <mitsuba/core/sstream.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/sshstream.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/appender.h>
#include <mitsuba/core/version.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <fstream>
#include <stdexcept>
#include <boost/algorithm/string.hpp>

#ifdef __WINDOWS__
#include <io.h>
#include <ws2tcpip.h>
#include <mitsuba/core/getopt.h>
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

using namespace mitsuba;

static bool running = true;
static SOCKET sock = INVALID_SOCKET;

#if defined(__WINDOWS__)
BOOL CtrlHandler(DWORD type) {
    switch (type) {
        case CTRL_C_EVENT:
            running = false;
            if (sock) {
                closesocket(sock);
                sock = INVALID_SOCKET;
            }
            return TRUE;
        default:
            return FALSE;
    }
}

#else
/* Catch Ctrl+C, SIGKILL */
void sigterm_handler(int) {
    SLog(EInfo, "Caught signal - shutting down..");
    running = false;

    /* The next signal will immediately terminate the
       program. (a precaution for hung processes) */
    signal(SIGTERM, SIG_DFL);
    signal(SIGINT, SIG_DFL);
}

/* Collect zombie processes */
void collect_zombies(int s) {
    while (waitpid(-1, NULL, WNOHANG) > 0);
}
#endif

int mtssrv(int argc, char **argv) {
    int optchar;
    char *end_ptr = NULL;

    try {
        /* Default settings */
        int nprocs = getCoreCount(),
            listenPort = MTS_DEFAULT_PORT;
        std::string nodeName = getHostName(),
                    networkHosts = "";
        bool quietMode = false;
        ELogLevel logLevel = EInfo;
        std::string hostName = getFQDN();
        FileResolver *fileResolver = Thread::getThread()->getFileResolver();
        bool hostNameSet = false;

        optind = 1;
        /* Parse command-line arguments */
        while ((optchar = getopt(argc, argv, "a:c:s:n:p:i:l:L:qhv")) != -1) {
            switch (optchar) {
                case 'a': {
                        std::vector<std::string> paths = tokenize(optarg, ";");
                        for (int i=(int)paths.size()-1; i>=0; --i)
                            fileResolver->prependPath(paths[i]);
                    }
                    break;
                case 'c':
                    networkHosts = networkHosts + std::string(";") + std::string(optarg);
                    break;
                case 'i':
                    hostName = optarg;
                    hostNameSet = true;
                    break;
                case 's': {
                        std::ifstream is(optarg);
                        if (is.fail())
                            SLog(EError, "Could not open host file!");
                        std::string host;
                        while (is >> host) {
                            if (host.length() < 1 || host.c_str()[0] == '#')
                                continue;
                            networkHosts = networkHosts + std::string(";") + host;
                        }
                    }
                    break;
                case 'n':
                    nodeName = optarg;
                    break;
                case 'p':
                    nprocs = strtol(optarg, &end_ptr, 10);
                    if (*end_ptr != '\0')
                        SLog(EError, "Could not parse the processor count!");
                    break;
                case 'v':
                    logLevel = EDebug;
                    break;
                case 'L': {
                        std::string arg = boost::to_lower_copy(std::string(optarg));
                        if (arg == "trace")
                            logLevel = ETrace;
                        else if (arg == "debug")
                            logLevel = EDebug;
                        else if (arg == "info")
                            logLevel = EInfo;
                        else if (arg == "warn")
                            logLevel = EWarn;
                        else if (arg == "error")
                            logLevel = EError;
                        else
                            SLog(EError, "Invalid log level!");
                    }
                    break;
                case 'l':
                    if (!strcmp("s", optarg)) {
                        listenPort = -1;
                        quietMode = true;
                    } else {
                        listenPort = strtol(optarg, &end_ptr, 10);
                        if (*end_ptr != '\0')
                            SLog(EError, "Could not parse the port number");
                    }
                    break;
                case 'q':
                    quietMode = true;
                    break;
                case 'h':
                default:
                    cout <<  "Mitsuba version " << Version(MTS_VERSION).toStringComplete()
                        << ", Copyright (c) " MTS_YEAR " Wenzel Jakob" << endl;
                    cout <<  "Usage: mtssrv [options]" << endl;
                    cout <<  "Options/Arguments:" << endl;
                    cout <<  "   -h          Display this help text" << endl << endl;
                    cout <<  "   -a p1;p2;.. Add one or more entries to the resource search path" << endl << endl;
                    cout <<  "   -p count    Override the detected number of processors. Useful for reducing" << endl;
                    cout <<  "               the load or creating scheduling-only nodes in conjunction with"  << endl;
                    cout <<  "               the -c and -s parameters, e.g. -p 0 -c host1;host2;host3,..." << endl << endl;
                    cout <<  "   -q          Quiet mode - do not print any log messages to stdout" << endl << endl;
                    cout <<  "   -c hosts    Nesting: connect to additional mtssrv instances over a network." << endl;
                    cout <<  "               Requires a semicolon-separated list of host names of the form" << endl;
                    cout <<  "                       host.domain[:port] for a direct connection" << endl;
                    cout <<  "                 or" << endl;
                    cout <<  "                       user@host.domain[:/path] for a SSH connection (where" << endl;
                    cout <<  "                       'path' denotes the place where Mitsuba is checked" << endl;
                    cout <<  "                       out -- by default, \"~/mitsuba\" is used)" << endl << endl;
                    cout <<  "   -s file     Connect to additional Mitsuba servers specified in a file" << endl;
                    cout <<  "               with one name per line (same format as in -c)" << endl<< endl;
                    cout <<  "   -i name     IP address / host name on which to listen for connections" << endl << endl;
                    cout <<  "   -l port     Listen for connections on a certain port (Default: " << MTS_DEFAULT_PORT << ")." << endl;
                    cout <<  "               To listen on stdin, specify \"-ls\" (implies -q)" << endl << endl;
                    cout <<  "   -n name     Assign a node name to this instance (Default: host name)" << endl << endl;
                    cout <<  "   -v          Be more verbose (can be specified twice)" << endl << endl;
                    cout <<  "   -L level    Explicitly specify the log level (trace/debug/info/warn/error)" << endl << endl;
                    cout <<  " For documentation, please refer to http://www.mitsuba-renderer.org/docs.html" << endl;
                    return 0;
            }
        }

        /* Configure the logging subsystem */
        ref<Logger> log = Thread::getThread()->getLogger();
        log->setLogLevel(logLevel);

        /* Initialize OpenMP */
        Thread::initializeOpenMP(nprocs);

        /* Disable the default appenders */
        for (size_t i=0; i<log->getAppenderCount(); ++i) {
            Appender *appender = log->getAppender(i);
            if (appender->getClass()->derivesFrom(MTS_CLASS(StreamAppender)))
                log->removeAppender(appender);
        }

        log->addAppender(new StreamAppender(formatString("mtssrv.%s.log", nodeName.c_str())));
        if (!quietMode)
            log->addAppender(new StreamAppender(&std::cout));

        SLog(EInfo, "Mitsuba version %s, Copyright (c) " MTS_YEAR " Wenzel Jakob",
            Version(MTS_VERSION).toStringComplete().c_str());

#if defined(__WINDOWS__)
        /* Custom handler for Ctrl-C signals */
        SetConsoleCtrlHandler((PHANDLER_ROUTINE) CtrlHandler, TRUE);
#endif

        /* Configure the scheduling subsystem */
        Scheduler *scheduler = Scheduler::getInstance();
        for (int i=0; i<nprocs; ++i)
            scheduler->registerWorker(new LocalWorker(i, formatString("wrk%i", i)));
        std::vector<std::string> hosts = tokenize(networkHosts, ";");

        /* Establish network connections to nested servers */
        for (size_t i=0; i<hosts.size(); ++i) {
            const std::string &hostName = hosts[i];
            ref<Stream> stream;

            if (hostName.find("@") == std::string::npos) {
                int port = MTS_DEFAULT_PORT;
                std::vector<std::string> tokens = tokenize(hostName, ":");
                if (tokens.size() == 0 || tokens.size() > 2) {
                    SLog(EError, "Invalid host specification '%s'!", hostName.c_str());
                } else if (tokens.size() == 2) {
                    port = strtol(tokens[1].c_str(), &end_ptr, 10);
                    if (*end_ptr != '\0')
                        SLog(EError, "Invalid host specification '%s'!", hostName.c_str());
                }
                stream = new SocketStream(tokens[0], port);
            } else {
                std::string path = "~/mitsuba";
                std::vector<std::string> tokens = tokenize(hostName, "@/:");
                if (tokens.size() < 2 || tokens.size() > 3) {
                    SLog(EError, "Invalid host specification '%s'!", hostName.c_str());
                } else if (tokens.size() == 3) {
                    path = tokens[2];
                }
                std::vector<std::string> cmdLine;
                cmdLine.push_back(formatString("bash -c 'cd %s; . setpath.sh; mtssrv -ls'", path.c_str()));
                stream = new SSHStream(tokens[0], tokens[1], cmdLine);
            }
            try {
                scheduler->registerWorker(new RemoteWorker(formatString("net%i", i), stream));
            } catch (std::runtime_error &e) {
                if (hostName.find("@") != std::string::npos) {
#if defined(__WINDOWS__)
                    SLog(EWarn, "Please ensure that passwordless authentication "
                        "using plink.exe and pageant.exe is enabled (see the documentation for more information)");
#else
                    SLog(EWarn, "Please ensure that passwordless authentication "
                        "is enabled (e.g. using ssh-agent - see the documentation for more information)");
#endif
                }
                throw e;
            }
        }
        scheduler->start();

        if (listenPort == -1) {
            ref<StreamBackend> backend = new StreamBackend("con0",
                    scheduler, nodeName, new ConsoleStream(), false);
            backend->start();
            backend->join();
            return 0;
        }

        /* Allocate a socket of the proper type (IPv4/IPv6) */
        struct addrinfo hints, *servinfo, *p = NULL;
        memset(&hints, 0, sizeof(struct addrinfo));
        hints.ai_family = AF_UNSPEC;
        hints.ai_flags = AI_PASSIVE;
        hints.ai_socktype = SOCK_STREAM;
        char portName[8];
        int rv, one = 1;
        sock = INVALID_SOCKET;

        snprintf(portName, sizeof(portName), "%i", listenPort);
        if ((rv = getaddrinfo(hostNameSet ? hostName.c_str() : NULL, portName, &hints, &servinfo)) != 0)
            SLog(EError, "Error in getaddrinfo(%s:%i): %s", hostName.c_str(), listenPort, gai_strerror(rv));

        for (p = servinfo; p != NULL; p = p->ai_next) {
            /* Allocate a socket */
            sock = socket(p->ai_family, p->ai_socktype, p->ai_protocol);
            if (sock == -1)
                SocketStream::handleError("none", "socket");

            /* Avoid "bind: socket already in use" */
            if (setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, (char *) &one, sizeof(int)) < 0)
                SocketStream::handleError("none", "setsockopt");

            /* Bind the socket to the port number */
            if (bind(sock, p->ai_addr, (socklen_t) p->ai_addrlen) == -1) {
                SocketStream::handleError("none", formatString("bind(%s:%i)", hostName.c_str(), listenPort), EError);
#if defined(__WINDOWS__)
                closesocket(sock);
#else
                close(sock);
#endif
                continue;
            }
            break;
        }

        if (p == NULL)
            SLog(EError, "Failed to bind to port %i!", listenPort);
        freeaddrinfo(servinfo);

        if (listen(sock, CONN_BACKLOG) == -1)
            SocketStream::handleError("none", "bind");
        SLog(EInfo, "Enter mtssrv -h for more options");

#if defined(__WINDOWS__)
        SLog(EInfo, "%s: Listening on port %i.. Send Ctrl-C to stop.", hostName.c_str(), listenPort);
#else
        /* Avoid zombies processes */
        struct sigaction sa;
        sa.sa_handler = collect_zombies;
        sigemptyset(&sa.sa_mask);
        sa.sa_flags = SA_RESTART;

        if (sigaction(SIGCHLD, &sa, NULL) == -1)
            SLog(EError, "Error in sigaction(): %s!", strerror(errno));

        sa.sa_handler = sigterm_handler;
        sa.sa_flags = 0; // we want SIGINT/SIGTERM to interrupt accept()

        if (sigaction(SIGTERM, &sa, NULL) == -1)
            SLog(EError, "Error in sigaction(): %s!", strerror(errno));
        if (sigaction(SIGINT, &sa, NULL) == -1)
            SLog(EError, "Error in sigaction(): %s!", strerror(errno));

        /* Ignore SIGPIPE */
        signal(SIGPIPE, SIG_IGN);

        SLog(EInfo, "%s: Listening on port %i.. Send Ctrl-C or SIGTERM to stop.", hostName.c_str(), listenPort);
#endif

        int connectionIndex = 0;
        /* Wait for connections */
        while (running) {
            socklen_t addrlen = sizeof(sockaddr_storage);
            struct sockaddr_storage sockaddr;
            memset(&sockaddr, 0, addrlen);

            SOCKET newSocket = accept(sock, (struct sockaddr *) &sockaddr, &addrlen);
            if (newSocket == INVALID_SOCKET) {
#if defined(__WINDOWS__)
                if (!running)
                    break;
#else
                if (errno == EINTR)
                    continue;
#endif
                SocketStream::handleError("none", "accept", EWarn);
                continue;
            }

            ref<StreamBackend> backend = new StreamBackend(formatString("con%i", connectionIndex++),
                scheduler, nodeName, new SocketStream(newSocket), true);
            backend->start();
        }
#if defined(__WINDOWS__)
        SLog(EInfo, "Caught signal - shutting down..");
#else
        close(sock);
#endif
    } catch (const std::exception &e) {
        std::cerr << "Caught a critical exception: " << e.what() << endl;
    } catch (...) {
        std::cerr << "Caught a critical exception of unknown type!" << endl;
    }

    /* Shutdown */
    Statistics::getInstance()->printStats();

    return 0;
}

int mts_main(int argc, char **argv) {
    /* Initialize the core framework */
    Class::staticInitialization();
    Object::staticInitialization();
    PluginManager::staticInitialization();
    Statistics::staticInitialization();
    Thread::staticInitialization();
    Logger::staticInitialization();
    FileStream::staticInitialization();
    Spectrum::staticInitialization();
    Bitmap::staticInitialization();
    Scheduler::staticInitialization();
    SHVector::staticInitialization();

#if defined(__WINDOWS__)
    /* Initialize WINSOCK2 */
    WSADATA wsaData;
    if (WSAStartup(MAKEWORD(2,2), &wsaData))
        SLog(EError, "Could not initialize WinSock2!");
    if (LOBYTE(wsaData.wVersion) != 2 || HIBYTE(wsaData.wVersion) != 2)
        SLog(EError, "Could not find the required version of winsock.dll!");
#endif

#if defined(__LINUX__) || defined(__OSX__)
    setlocale(LC_NUMERIC, "C");
#endif

    int retval = mtssrv(argc, argv);

    /* Shutdown the core framework */
    SHVector::staticShutdown();
    Scheduler::staticShutdown();
    Bitmap::staticShutdown();
    Spectrum::staticShutdown();
    FileStream::staticShutdown();
    Logger::staticShutdown();
    Thread::staticShutdown();
    Statistics::staticShutdown();
    PluginManager::staticShutdown();
    Object::staticShutdown();
    Class::staticShutdown();


#if defined(__WINDOWS__)
    /* Shut down WINSOCK2 */
    WSACleanup();
#endif

    return retval;
}

#if !defined(__WINDOWS__)
int main(int argc, char **argv) {
    return mts_main(argc, argv);
}
#endif

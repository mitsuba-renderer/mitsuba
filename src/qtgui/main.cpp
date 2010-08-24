#include <xercesc/parsers/SAXParser.hpp>
#include <QtGui/QtGui>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/plugin.h>
#if defined(__OSX__)
#include <ApplicationServices/ApplicationServices.h>
#endif
#include "mainwindow.h"

#if defined(__LINUX__)
#include <X11/Xlib.h>
#endif

#if !defined(WIN32)
#include <signal.h>
#include <sys/wait.h>
#include <errno.h>
#endif

XERCES_CPP_NAMESPACE_USE

using namespace mitsuba;

MainWindow *mainWindow = NULL;

class MitsubaApplication : public QApplication {
public:
	MitsubaApplication(int &argc, char **argv) : QApplication(argc, argv) {
	}

	bool event(QEvent *event) {
		switch (event->type()) {
#if defined(__OSX__)
			case QEvent::Quit:
				quit();
				return true;
#endif
			case QEvent::FileOpen:
				if (mainWindow != NULL)
					mainWindow->loadFile(static_cast<QFileOpenEvent *>(event)->file());
				return true;
			default:
				return QApplication::event(event);
		}
	}

	bool notify(QObject *rec, QEvent *e) {
		try {
			return QApplication::notify(rec, e);
		} catch (const std::exception &e) {
			SLog(EWarn, "Caught exception: %s", e.what());
			QMessageBox::critical(NULL, tr("Critical exception"),
				e.what(), QMessageBox::Ok);
			return false;
		}
	}
};

/* Collect zombie processes */
#if !defined(WIN32)
void collect_zombies(int s) {
	while (waitpid(-1, NULL, WNOHANG) > 0);
}
#endif

int main(int argc, char *argv[]) {
	int retval;

	/* Initialize Xerces-C */
	try {
		XMLPlatformUtils::Initialize();
	} catch(const XMLException &toCatch) {
		fprintf(stderr, "Error during Xerces initialization: %s",
			XMLString::transcode(toCatch.getMessage()));
		return -1;
	}

	/* Initialize the core framework */
	Class::staticInitialization();
	PluginManager::staticInitialization();
	Statistics::staticInitialization();
	Thread::staticInitialization();
	Logger::staticInitialization();
	Spectrum::staticInitialization();
	Scheduler::staticInitialization();
	SHVector::staticInitialization();

#if defined(__LINUX__)
	XInitThreads();
	FileResolver *resolver = FileResolver::getInstance();
	resolver->addPath("/usr/share/mitsuba");
#endif

#if defined(__OSX__)
	MTS_AUTORELEASE_BEGIN()
	FileResolver *resolver = FileResolver::getInstance();
	resolver->addPath(__ubi_bundlepath());
	/* Required for the mouse relocation in GLWidget */
	CGSetLocalEventsSuppressionInterval(0.0f);
	MTS_AUTORELEASE_END() 
#endif

#ifdef WIN32
	char lpFilename[1024];
	if (GetModuleFileNameA(NULL,
		lpFilename, sizeof(lpFilename))) {
		FileResolver *resolver = FileResolver::getInstance();
		resolver->addPathFromFile(lpFilename);
	} else {
		SLog(EWarn, "Could not determine the executable path");
	}

	/* Initialize WINSOCK2 */
	WSADATA wsaData;
	if (WSAStartup(MAKEWORD(2,2), &wsaData)) 
		SLog(EError, "Could not initialize WinSock2!");
	if (LOBYTE(wsaData.wVersion) != 2 || HIBYTE(wsaData.wVersion) != 2)
		SLog(EError, "Could not find the required version of winsock.dll!");
#endif

#ifdef __LINUX__
	char exePath[PATH_MAX];
	memset(exePath, 0, PATH_MAX);
	if (readlink("/proc/self/exe", exePath, PATH_MAX) != -1) {
		FileResolver *resolver = FileResolver::getInstance();
		resolver->addPathFromFile(exePath);
	} else {
		SLog(EWarn, "Could not determine the executable path");
	}
#endif

#if !defined(WIN32)
	/* Avoid zombies processes when running the server */
	struct sigaction sa;
	sa.sa_handler = collect_zombies;
	sigemptyset(&sa.sa_mask);
	sa.sa_flags = SA_RESTART;

	if (sigaction(SIGCHLD, &sa, NULL) == -1)
		SLog(EWarn, "Error in sigaction(): %s!", strerror(errno));
#endif

	qRegisterMetaType<ELogLevel>("ELogLevel");

	MitsubaApplication app(argc, argv);
	try {
		QFile stylesheet(":/resources/stylesheet.css");

		if (!stylesheet.open(QFile::ReadOnly)) {
			QMessageBox::critical(NULL, "Internal error", "Could not open stylesheet!");
			exit(-1);
		}
		app.setStyleSheet(QTextStream(&stylesheet).readAll().toAscii());

#if defined(__OSX__)
		app.setAttribute(Qt::AA_DontShowIconsInMenus); 
#endif
		/* Disable the default appenders */
		ref<Logger> logger = Thread::getThread()->getLogger();
		for (size_t i=0; i<logger->getAppenderCount(); ++i) {
			Appender *appender = logger->getAppender(i);
			if (appender->getClass()->derivesFrom(StreamAppender::m_theClass))
				logger->removeAppender(appender);
		}

#if defined(__OSX__)
		/* Create a log file inside the application bundle */
		MTS_AUTORELEASE_BEGIN() 
		logger->addAppender(new StreamAppender(formatString("%s/mitsuba.%s.log", 
			__ubi_bundlepath().c_str(), getHostName().c_str())));
		MTS_AUTORELEASE_END() 
#else
		/* Create a log file inside the current working directory */
		logger->addAppender(new StreamAppender(formatString("mitsuba.%s.log", getHostName().c_str())));
#endif

#if !defined(WIN32)
		/* Correct number parsing on some locales (e.g. ru_RU) */
		setlocale(LC_NUMERIC, "C");
#endif

		mainWindow = new MainWindow();
		mainWindow->initWorkers();
		retval = app.exec();
		delete mainWindow;
	} catch (const std::exception &e) {
		SLog(EWarn, "Critical exception during startup: %s", e.what());
		QMessageBox::critical(NULL, QString("Critical exception"),
			e.what(), QMessageBox::Ok);
		retval = -1;
	}
	Statistics::getInstance()->printStats();

	XMLPlatformUtils::Terminate();

#ifdef WIN32
	/* Shut down WINSOCK2 */
	WSACleanup();
#endif

	/* Shutdown the core framework */
	SHVector::staticShutdown();
	Scheduler::staticShutdown();
	Spectrum::staticShutdown();
	Logger::staticShutdown();
	Thread::staticShutdown();
	Statistics::staticShutdown();
	PluginManager::staticShutdown();
	Class::staticShutdown();

	return retval;
}

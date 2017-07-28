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

#include <mitsuba/core/platform.h>
#include <QtGui/QtGui>
#include <QtOpenGL/QGLFormat>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/appender.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/scenehandler.h>

#if defined(__OSX__)
#include <ApplicationServices/ApplicationServices.h>
#endif
#include "mainwindow.h"

#if defined(__LINUX__)
#include <X11/Xlib.h>
#endif

#if !defined(__WINDOWS__)
#include <signal.h>
#include <sys/wait.h>
#include <errno.h>
#else
#include <winsock2.h>
#endif

#if defined(__WINDOWS__)
// Always use the High Performance GPU on machines using NVIDIA Optimus
// http://stackoverflow.com/questions/10535950/forcing-nvidia-gpu-programmatically-in-optimus-laptops [January 2014]
extern "C" {
	_declspec(dllexport) DWORD NvOptimusEnablement = 0x00000001;
}
#endif

using namespace mitsuba;
MainWindow *mainWindow = NULL;


/// ================================================================
///  Handle application crashes when compiled with MTS_HAS_BREAKPAD
/// ================================================================
#if defined(MTS_HAS_BREAKPAD)
#if defined(__OSX__)

extern void *__mts_init_breakpad_osx();
extern void  __mts_destroy_breakpad_osx(void *);

#elif defined(__WINDOWS__)

#include <breakpad/client/windows/sender/crash_report_sender.h>
#include <breakpad/client/windows/handler/exception_handler.h>

static bool minidumpCallbackWindows(const wchar_t *dump_path,
		const wchar_t *minidump_id, void *context,
		EXCEPTION_POINTERS *exinfo, MDRawAssertionInfo *assertion,
		bool succeeded) {
	if (!dump_path || !minidump_id)
		return false;

	std::wstring filename = std::wstring(dump_path) + std::wstring(L"\\")
		+ std::wstring(minidump_id) + std::wstring(L".dmp");

	google_breakpad::CrashReportSender sender(L"");
	std::map<std::wstring, std::wstring> parameters;

	#if defined(WIN64)
		parameters[L"prod"] = L"Mitsuba/Win64";
	#else
		parameters[L"prod"] = L"Mitsuba/Win32";
	#endif

	std::string version = MTS_VERSION;
	std::wstring wVersion;
	wVersion.assign(version.begin(), version.end());
	parameters[L"ver"] = wVersion;

	std::wstring resultString;

	if (MessageBox(NULL, TEXT("Mitsuba crashed due to an internal error. If you agree below, a brief "
			"report describing the failure will be submitted to the developers. If you would like to "
			"accelerate the debugging process further, please also create a ticket with information on "
			"the steps that led to the problem on https://www.mitsuba-renderer.org/tracker -- thank you!"),
			TEXT("Error"), MB_OKCANCEL | MB_ICONERROR) != IDOK)
		return false;

	google_breakpad::ReportResult result =
		sender.SendCrashReport(L"http://www.mitsuba-renderer.org/bugreport.php",
		parameters, filename, &resultString);

	if (result != google_breakpad::RESULT_SUCCEEDED) {
		MessageBox(NULL, TEXT("The error report could not be submitted due to a lack of "
			"internet connectivity!"), TEXT("Error"), MB_OK | MB_ICONERROR);
		return false;
	}

	return succeeded;
}

#endif
#endif

/// ================================================================

class MitsubaApplication : public QApplication {
public:
	MitsubaApplication(int &argc, char **argv) : QApplication(argc, argv) {
		setOrganizationName("mitsuba-renderer.org");
		setOrganizationDomain("mitsuba-renderer.org");
		setApplicationName("mtsgui");
		setApplicationVersion(MTS_VERSION);
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

#if defined(__OSX__)
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

int main(int argc, char *argv[]) {
	int retval = -1;

	/* Initialize the core framework */
	Class::staticInitialization();
	Object::staticInitialization();
	PluginManager::staticInitialization();
	Statistics::staticInitialization();
	Thread::staticInitialization();
	Logger::staticInitialization();
	FileStream::staticInitialization();
	Thread::initializeOpenMP(getCoreCount());
	Spectrum::staticInitialization();
	Bitmap::staticInitialization();
	Scheduler::staticInitialization();
	SHVector::staticInitialization();
	SceneHandler::staticInitialization();

#if defined(__LINUX__)
	XInitThreads();
#endif

#if defined(__OSX__)
	MTS_AUTORELEASE_BEGIN()
	/* Required for the mouse relocation in GLWidget */
	CGEventSourceRef evsrc =
		CGEventSourceCreate(kCGEventSourceStateCombinedSessionState);
	CGEventSourceSetLocalEventsSuppressionInterval(evsrc, 0.0);
	CFRelease(evsrc);
	MTS_AUTORELEASE_END()
#endif

#if defined(__WINDOWS__)
	/* Initialize WINSOCK2 */
	WSADATA wsaData;
	if (WSAStartup(MAKEWORD(2,2), &wsaData))
		SLog(EError, "Could not initialize WinSock2!");
	if (LOBYTE(wsaData.wVersion) != 2 || HIBYTE(wsaData.wVersion) != 2)
		SLog(EError, "Could not find the required version of winsock.dll!");
#endif

#if !defined(__WINDOWS__)
	/* Avoid zombies processes when running the server */
	if (signal(SIGCHLD, SIG_IGN) == SIG_ERR)
		SLog(EWarn, "Error in signal(): %s!", strerror(errno));
#endif

	qRegisterMetaType<ELogLevel>("ELogLevel");
	qRegisterMetaType<fs::path>("fs::path");

	MitsubaApplication app(argc, argv);
	try {
		QFile stylesheet(":/resources/stylesheet.css");

		if (!stylesheet.open(QFile::ReadOnly)) {
			QMessageBox::critical(NULL, "Internal error", "Could not open stylesheet!");
			exit(-1);
		}
		app.setStyleSheet(QTextStream(&stylesheet).readAll().toLatin1());

#if defined(__OSX__)
		app.setAttribute(Qt::AA_DontShowIconsInMenus);
#endif
		/* Disable the default appenders */
		ref<Logger> logger = Thread::getThread()->getLogger();
		for (size_t i=0; i<logger->getAppenderCount(); ++i) {
			Appender *appender = logger->getAppender(i);
			if (appender->getClass()->derivesFrom(MTS_CLASS(StreamAppender)))
				logger->removeAppender(appender);
		}

#if defined(__OSX__)
		/* Create a log file inside the application bundle */
		MTS_AUTORELEASE_BEGIN()
		logger->addAppender(new StreamAppender(formatString("%s/mitsuba.%s.log",
			__mts_bundlepath().c_str(), getHostName().c_str())));
		MTS_AUTORELEASE_END()

		/* Set application defaults (disable OSX synchronization feature) */
		__mts_set_appdefaults();
#else
		/* Create a log file inside the current working directory */
		logger->addAppender(new StreamAppender(formatString("mitsuba.%s.log", getHostName().c_str())));
#endif

#if !defined(__WINDOWS__)
		/* Correct number parsing on some locales (e.g. ru_RU) */
		setlocale(LC_NUMERIC, "C");
#endif

#if defined(MTS_HAS_BREAKPAD)
	#if defined(__OSX__)
		void *breakpad = __mts_init_breakpad_osx();
	#elif defined(__WINDOWS__)
		_CrtSetReportMode(_CRT_ASSERT, 0);
		std::wstring dump_path;
		dump_path.resize(1024);
		GetTempPathW(1024, &dump_path[0]);

		google_breakpad::ExceptionHandler *breakpad =
			new google_breakpad::ExceptionHandler(
				dump_path, NULL, minidumpCallbackWindows, NULL,
				google_breakpad::ExceptionHandler::HANDLER_ALL,
				MiniDumpNormal, NULL, NULL);
	#endif
#endif

#if !defined(MTS_GUI_SOFTWARE_FALLBACK)
		/* Be a bit more picky about the rendering target */
		QGLFormat fmt;
		fmt.setDepth(true);
		fmt.setStencil(false);
		fmt.setAlpha(true);
		fmt.setStereo(false);
		fmt.setDoubleBuffer(true);
		QGLFormat::setDefaultFormat(fmt);
#endif

		mainWindow = new MainWindow();
		if (mainWindow->initWorkersProcessArgv())
			retval = app.exec();

#if defined(MTS_HAS_BREAKPAD)
	#if defined(__OSX__)
		__mts_destroy_breakpad_osx(breakpad);
	#elif defined(__WINDOWS__)
		delete breakpad;
	#endif
#endif
		delete mainWindow;
	} catch (const std::exception &e) {
		SLog(EWarn, "Critical exception during startup: %s", e.what());
		QMessageBox::critical(NULL, QString("Critical exception"),
			e.what(), QMessageBox::Ok);
	}
	Statistics::getInstance()->printStats();


#if defined(__WINDOWS__)
	/* Shut down WINSOCK2 */
	WSACleanup();
#endif

	/* Shutdown the core framework */
	SceneHandler::staticShutdown();
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

	return retval;
}

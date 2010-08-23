#include <mitsuba/core/platform.h>
#include <xercesc/parsers/SAXParser.hpp>
#include <mitsuba/core/sched_remote.h>
#include <mitsuba/core/sstream.h>
#include <mitsuba/core/sshstream.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/renderjob.h>
#include <fstream>
#include <stdexcept>
#include <sys/types.h>
#include "shandler.h"

#if !defined(WIN32)
#include <dlfcn.h>
#include <dirent.h>
#endif

#if defined(WIN32)
#include "getopt.h"
#else
#include <signal.h>
#endif

using namespace mitsuba;

class UtilityServicesImpl : public UtilityServices {
public:
	ref<Scene> loadScene(const std::string &filename) {
		/* Prepare for parsing scene descriptions */
		FileResolver *resolver = FileResolver::getInstance();
		SAXParser* parser = new SAXParser();
		std::string schemaPath = resolver->resolveAbsolute("schema/scene.xsd");

		/* Check against the 'scene.xsd' XML Schema */
		parser->setDoSchema(true);
		parser->setValidationSchemaFullChecking(true);
		parser->setValidationScheme(SAXParser::Val_Always);
		parser->setExternalNoNamespaceSchemaLocation(schemaPath.c_str());

		std::map<std::string, std::string> parameters;
		SceneHandler *handler = new SceneHandler(parameters);
		parser->setDoNamespaces(true);
		parser->setDocumentHandler(handler);
		parser->setErrorHandler(handler);
			
		parser->parse(filename.c_str());
		ref<Scene> scene = handler->getScene();

		delete parser;
		delete handler;

		return scene;
	}
};

void help() {
	cout <<  "Mitsuba version " MTS_VERSION ", Copyright (c) " MTS_YEAR " Wenzel Jakob" << endl;
	cout <<  "Usage: mtsutil [mtsutil options] <utility name> [arguments]" << endl;
	cout <<  "Options/Arguments:" << endl;
	cout <<  "   -h          Display this help text" << endl << endl;
	cout <<  "   -a p1;p2;.. Add one or more entries to the resource search path" << endl << endl;
	cout <<  "   -p count    Override the detected number of processors. Useful for reducing" << endl;
	cout <<  "               the load or creating scheduling-only nodes in conjunction with"  << endl;
	cout <<  "               the -c and -s parameters, e.g. -p 0 -c host1;host2;host3,..." << endl << endl;
	cout <<  "   -q          Quiet mode - do not print any log messages to stdout" << endl << endl;
	cout <<  "   -c hosts    Network processing: connect to mtssrv instances over a network." << endl;
	cout <<  "               Requires a semicolon-separated list of host names of the form" << endl;
	cout <<  "                       host.domain[:port] for a direct connection" << endl;
	cout <<  "                 or" << endl;
	cout <<  "                       user@host.domain[:path] for a SSH connection (where" << endl;
	cout <<  "                       \"path\" denotes the place where Mitsuba is checked" << endl;
	cout <<  "                       out -- by default, \"~/mitsuba\" is used)" << endl << endl;
	cout <<  "   -s file     Connect to additional Mitsuba servers specified in a file" << endl;
	cout <<  "               with one name per line (same format as in -c)" << endl<< endl;
	cout <<  "   -n name     Assign a node name to this instance (Default: host name)" << endl << endl;
	cout <<  "   -v          Be more verbose" << endl << endl;

	FileResolver *resolver = FileResolver::getInstance();
	cout << "The following utilities are available:" << endl << endl;

	std::vector<std::string> dirPaths = resolver->resolveAllAbsolute("plugins");
	std::set<std::string> seen;

	for (size_t i=0; i<dirPaths.size(); ++i) {
		std::string dirPath = dirPaths[i];

#if !defined(WIN32)
		DIR *directory;
		struct dirent *dirinfo;

		if ((directory = opendir(dirPath.c_str())) == NULL)
			SLog(EInfo, "Could not open plugin directory");

		while ((dirinfo = readdir(directory)) != NULL) {
			std::string fname(dirinfo->d_name);
			if (!endsWith(fname, ".dylib") && !endsWith(fname, ".so"))
				continue;
			std::string fullName = dirPath + "/" + fname;
#else
		HANDLE hFind;
		WIN32_FIND_DATA findFileData;

		if ((hFind = FindFirstFile((dirPath + "\\*.dll").c_str(), &findFileData)) == INVALID_HANDLE_VALUE)
			SLog(EInfo, "Could not open plugin directory");
		
		do {
			std::string fname = findFileData.cFileName;
			std::string fullName = dirPath + "\\" + fname;
#endif
			std::string shortName = fname.substr(0, strrchr(fname.c_str(), '.') - fname.c_str());
			if (seen.find(shortName) != seen.end())
				continue;
			seen.insert(shortName);
			Plugin utility(shortName, fullName);
			if (!utility.isUtility())
				continue;
			cout << "\t" << shortName;
			for (int i=0; i<22-(int) shortName.length(); ++i)
				cout << ' ';
			cout  << utility.getDescription() << endl;
#if !defined(WIN32)	
		}
#else
		} while (FindNextFile(hFind, &findFileData));
		FindClose(hFind);
#endif
	}
}


int ubi_main(int argc, char **argv) {
	char optchar, *end_ptr = NULL;

	try {
		/* Default settings */
		int nprocs = getProcessorCount();
		std::string nodeName = getHostName(),
					networkHosts = "", destFile="";
		bool quietMode = false;
		ELogLevel logLevel = EInfo;
		FileResolver *resolver = FileResolver::getInstance();

		if (argc < 2) {
			help();
			return 0;
		}

		optind = 1;
		/* Parse command-line arguments */
		while ((optchar = getopt(argc, argv, "a:c:s:n:p:qhv")) != -1) {
			switch (optchar) {
				case 'a': {
						std::vector<std::string> paths = tokenize(optarg, ";");
						for (unsigned int i=0; i<paths.size(); ++i) 
							resolver->addPath(paths[i]);
					}
					break;
				case 'c':
					networkHosts = networkHosts + std::string(";") + std::string(optarg);
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
				case 'v':
					logLevel = EDebug;
					break;
				case 'p':
					nprocs = strtol(optarg, &end_ptr, 10);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the processor count!");
					break;
				case 'q':
					quietMode = true;
					break;
				case 'h':
				default:
					help();
					return 0;
			}
		}

		/* Configure the logging subsystem */
		ref<Logger> log = Thread::getThread()->getLogger();
		log->setLogLevel(logLevel);
		
		/* Disable the default appenders */
		for (size_t i=0; i<log->getAppenderCount(); ++i) {
			Appender *appender = log->getAppender(i);
			if (appender->getClass()->derivesFrom(StreamAppender::m_theClass))
				log->removeAppender(appender);
		}

		log->addAppender(new StreamAppender(formatString("mitsuba.%s.log", nodeName.c_str())));
		if (!quietMode)
			log->addAppender(new StreamAppender(&std::cout));

		SLog(EInfo, "Mitsuba version " MTS_VERSION ", Copyright (c) " MTS_YEAR " Wenzel Jakob");

		/* Configure the scheduling subsystem */
		Scheduler *scheduler = Scheduler::getInstance();
		for (int i=0; i<nprocs; ++i)
			scheduler->registerWorker(new LocalWorker(formatString("wrk%i", i)));
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
				std::string path = "~/mitsuba"; // default path if not specified
				std::vector<std::string> tokens = tokenize(hostName, "@:");
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
#if defined(WIN32)
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

		if (argc <= optind) {
			std::cerr << "A utility name must be supplied!" << endl;
			return -1;
		}

		/* Build the full plugin file name */
#if defined(WIN32)
		std::string shortName = std::string("plugins/") + argv[optind] + std::string(".dll");
#elif defined(__OSX__)
		std::string shortName = std::string("plugins/") + argv[optind] + std::string(".dylib");
#else
		std::string shortName = std::string("plugins/") + argv[optind] + std::string(".so");
#endif
		std::string fullName = resolver->resolve(shortName);

		if (!FileStream::exists(fullName)) {
			/* Plugin not found! */
			SLog(EError, "Utility \"%s\" not found (run \"mtsutil\" without arguments to "
				"see a list of available utilities)", fullName.c_str());
		}
		SLog(EInfo, "Loading utility \"%s\" ..", argv[optind]);
		Plugin *plugin = new Plugin(argv[optind], fullName);
		if (!plugin->isUtility())
			SLog(EError, "This plugin does not implement the 'Utility' interface!");
		Statistics::getInstance()->logPlugin(argv[optind], plugin->getDescription());
	
		UtilityServices *utilityServices = new UtilityServicesImpl();
		Utility *utility = plugin->createUtility(utilityServices);

		return utility->run(argc-optind, argv+optind);
	} catch (const std::exception &e) {
		std::cerr << "Caught a critical exeption: " << e.what() << std::endl;
	} catch (...) {
		std::cerr << "Caught a critical exeption of unknown type!" << endl;
	}

	return 0;
}

int main(int argc, char **argv) {
	/* Initialize the core framework */
	Class::staticInitialization();
	Statistics::staticInitialization();
	Thread::staticInitialization();
	Logger::staticInitialization();
	Spectrum::staticInitialization();
	Scheduler::staticInitialization();
	SHVector::staticInitialization();

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
	FileResolver *resolver = FileResolver::getInstance();
	if (readlink("/proc/self/exe", exePath, PATH_MAX) != -1) {
		resolver->addPathFromFile(exePath);
	} else {
		SLog(EWarn, "Could not determine the executable path");
	}
	resolver->addPath("/usr/share/mitsuba");
#endif

#if defined(__OSX__)
	MTS_AUTORELEASE_BEGIN()
	FileResolver *resolver = FileResolver::getInstance();
	resolver->addPath(__ubi_bundlepath());
	MTS_AUTORELEASE_END() 
#endif
		
#if !defined(WIN32)
	/* Correct number parsing on some locales (e.g. ru_RU) */
	setlocale(LC_NUMERIC, "C");
#endif

	/* Initialize Xerces-C */
	try {
		XMLPlatformUtils::Initialize();
	} catch(const XMLException &toCatch) {
		SLog(EError, "Error during Xerces initialization: %s",
			XMLString::transcode(toCatch.getMessage()));
		return -1;
	}
	
	int retval = ubi_main(argc, argv);

	XMLPlatformUtils::Terminate();

	/* Shutdown the core framework */
	SHVector::staticShutdown();
	Scheduler::staticShutdown();
	Spectrum::staticShutdown();
	Logger::staticShutdown();
	Thread::staticShutdown();
	Statistics::staticShutdown();
	Class::staticShutdown();
	
#ifdef WIN32
	/* Shut down WINSOCK2 */
	WSACleanup();
#endif

	return retval;
}

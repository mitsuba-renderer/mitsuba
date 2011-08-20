#include <mitsuba/core/fresolver.h>

MTS_NAMESPACE_BEGIN

FileResolver::FileResolver() {
	m_paths.push_back(fs::current_path());
#if defined(__LINUX__)
	char exePath[PATH_MAX];
	memset(exePath, 0, PATH_MAX);
	if (readlink("/proc/self/exe", exePath, PATH_MAX) != -1) {
		const fs::path exeParentPath = fs::path(exePath).parent_path();
		addPath(exeParentPath);
		// Handle local installs: ~/local/bin/:~/local/share/mitsuba/*
		fs::path sharedDir = exeParentPath.parent_path();
		sharedDir /= fs::path("share/mitsuba");
		if (fs::exists(sharedDir)) ;
			addPath(sharedDir);
	} else {
		Log(EError, "Could not detect the executable path!");
	}
#elif defined(__OSX__)
	MTS_AUTORELEASE_BEGIN()
	addPath(__mts_bundlepath());
	MTS_AUTORELEASE_END() 
#elif defined(WIN32)
	char lpFilename[1024];
	if (GetModuleFileNameA(NULL, lpFilename, sizeof(lpFilename)))
		addPath(fs::path(lpFilename).parent_path());
	else
		Log(EError, "Could not detect the executable path!");
#endif
}

FileResolver *FileResolver::clone() const {
	FileResolver *cloned = new FileResolver();
	cloned->m_paths = m_paths;
	return cloned;
}

void FileResolver::clear() {
	m_paths.clear();
}

void FileResolver::addPath(const fs::path &path) {
	bool found = false;
	for (size_t i=0; i<m_paths.size(); ++i) {
		if (m_paths[i] == path) {
			found = true;
			break;
		}
	}
	if (!found)
		m_paths.push_back(path);
}

fs::path FileResolver::resolve(const fs::path &path) const {
	if (!fs::exists(path)) {
		for (unsigned int i=0; i<m_paths.size(); i++) {
			fs::path newPath = m_paths[i] / path;
			if (fs::exists(newPath))
				return newPath;
		}
	}
	return path;
}

std::vector<fs::path> FileResolver::resolveAll(const fs::path &path) const {
	std::vector<fs::path> results;

	if (fs::exists(path)) 
		results.push_back(path);
	
	for (unsigned int i=0; i<m_paths.size(); i++) {
		fs::path newPath = m_paths[i] / path;
		if (fs::exists(newPath))
			results.push_back(newPath);
	}
	return results;
}

fs::path FileResolver::resolveAbsolute(const fs::path &path) const {
	return fs::complete(resolve(path));
}

std::string FileResolver::toString() const {
	std::ostringstream oss;
	oss << "FileResolver[" << endl
		<< "  paths = {" << endl;
	for (size_t i=0; i<m_paths.size(); ++i) {
		oss << "    \"" << m_paths[i].file_string() << "\"";
		if (i+1 < m_paths.size())
			oss << ",";
		oss << endl;
	}
	oss << "  }" << endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(FileResolver, false, Object)
MTS_NAMESPACE_END

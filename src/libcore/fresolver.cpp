#include <mitsuba/core/fresolver.h>

#if defined(__WINDOWS__)
# include <windows.h>
#endif

MTS_NAMESPACE_BEGIN

FileResolver::FileResolver() {
	m_paths.push_back(fs::current_path());
#if defined(__LINUX__)
	char exePath[PATH_MAX];
	memset(exePath, 0, PATH_MAX);
	if (readlink("/proc/self/exe", exePath, PATH_MAX) != -1) {
		const fs::path exeParentPath = fs::path(exePath).parent_path();
		prependPath(exeParentPath);
		// Handle local installs: ~/local/bin/:~/local/share/mitsuba/*
		fs::path sharedDir = exeParentPath.parent_path();
		sharedDir /= fs::path("share/mitsuba");
		if (fs::exists(sharedDir)) {
			prependPath(sharedDir);
		}
	} else {
		Log(EError, "Could not detect the executable path!");
	}
#elif defined(__OSX__)
	MTS_AUTORELEASE_BEGIN()
	prependPath(__mts_bundlepath());
	MTS_AUTORELEASE_END() 
#elif defined(__WINDOWS__)
	WCHAR lpFilename[MAX_PATH];
	const DWORD nSize = static_cast<DWORD>(sizeof(lpFilename)/sizeof(WCHAR));
	if (GetModuleFileNameW(NULL, lpFilename, nSize) != 0 &&
			GetLastError() == ERROR_SUCCESS)
		prependPath(fs::path(lpFilename).parent_path());
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

void FileResolver::prependPath(const fs::path &path) {
	for (size_t i=0; i<m_paths.size(); ++i) {
		if (m_paths[i] == path)
			return;
	}
	m_paths.push_front(path);
}

void FileResolver::appendPath(const fs::path &path) {
	for (size_t i=0; i<m_paths.size(); ++i) {
		if (m_paths[i] == path)
			return;
	}
	m_paths.push_back(path);
}

fs::path FileResolver::resolve(const fs::path &path) const {
	for (size_t i=0; i<m_paths.size(); i++) {
		fs::path newPath = m_paths[i] / path;
		if (fs::exists(newPath))
			return newPath;
	}
	return path;
}

std::vector<fs::path> FileResolver::resolveAll(const fs::path &path) const {
	std::vector<fs::path> results;

	for (size_t i=0; i<m_paths.size(); i++) {
		fs::path newPath = m_paths[i] / path;
		if (fs::exists(newPath))
			results.push_back(newPath);
	}

	return results;
}

fs::path FileResolver::resolveAbsolute(const fs::path &path) const {
	return fs::absolute(resolve(path));
}

std::string FileResolver::toString() const {
	std::ostringstream oss;
	oss << "FileResolver[" << endl
		<< "  paths = {" << endl;
	for (size_t i=0; i<m_paths.size(); ++i) {
		oss << "    \"" << m_paths[i].string() << "\"";
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

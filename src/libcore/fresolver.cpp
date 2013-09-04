#include <mitsuba/core/fresolver.h>
#include <boost/algorithm/string.hpp>

#if defined(__LINUX__)
# if !defined(_GNU_SOURCE)
#  define _GNU_SOURCE
# endif
# include <dlfcn.h>
#elif defined(__OSX__)
# include <mach-o/dyld.h>
#elif defined(__WINDOWS__)
# include <windows.h>
# include <vector>
#endif



MTS_NAMESPACE_BEGIN

#if defined(__WINDOWS__) || defined(__LINUX__)
	namespace {
		void dummySymbol() { }
	}
#endif

FileResolver::FileResolver() {
	/* Try to detect the base path of the Mitsuba installation */
	fs::path basePath;
#if defined(__LINUX__)
	Dl_info info;

	dladdr((const void *) &dummySymbol, &info);
	if (info.dli_fname) {
		/* Try to detect a few default setups */
		if (boost::starts_with(info.dli_fname, "/usr/lib") ||
			boost::starts_with(info.dli_fname, "/lib")) {
			basePath = fs::path("/usr/share/mitsuba");
		} else if (boost::starts_with(info.dli_fname, "/usr/local/lib")) {
			basePath = fs::path("/usr/local/share/mitsuba");
		} else {
			/* This is a locally-compiled repository */
			basePath = fs::path(info.dli_fname).parent_path();
		}
	}
#elif defined(__OSX__)
	MTS_AUTORELEASE_BEGIN()
	uint32_t imageCount = _dyld_image_count();
	for (uint32_t i=0; i<imageCount; ++i) {
		const char *imageName = _dyld_get_image_name(i);
		if (boost::ends_with(imageName, "libmitsuba-core.dylib")) {
			basePath = fs::canonical(imageName).parent_path().parent_path().parent_path();
			break;
		}
	}
	MTS_AUTORELEASE_END()
	if (basePath.empty())
		Log(EError, "Could not detect the executable path!");
#elif defined(__WINDOWS__)
	std::vector<WCHAR> lpFilename(MAX_PATH);

	// Module handle to this DLL. If the function fails it sets handle to NULL.
	// In that case GetModuleFileName will get the name of the executable which
	// is acceptable soft-failure behavior.
	HMODULE handle;
	GetModuleHandleExW(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS
	                 | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
          reinterpret_cast<LPCWSTR>(&dummySymbol), &handle);

	// Try to get the path with the default MAX_PATH length (260 chars)
	DWORD nSize = GetModuleFileNameW(handle, &lpFilename[0], MAX_PATH);

	// Adjust the buffer size in case if was too short
	while (nSize != 0 && nSize == lpFilename.size()) {
		lpFilename.resize(nSize * 2);
		nSize = GetModuleFileNameW(handle, &lpFilename[0],
			static_cast<DWORD>(lpFilename.size()));
	}

	// There is an error if and only if the function returns 0
	if (nSize != 0)
		basePath = fs::path(lpFilename).parent_path();
	else
		Log(EError, "Could not detect the executable path! (%s)", lastErrorText().c_str());
#endif
	#if BOOST_VERSION >= 104800
		m_paths.push_back(fs::canonical(basePath));
	#else
		m_paths.push_back(fs::absolute(basePath));
	#endif
	m_paths.push_back(fs::current_path());
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
	/* First, try to resolve in case-sensitive mode */
	for (size_t i=0; i<m_paths.size(); i++) {
		fs::path newPath = m_paths[i] / path;
		if (fs::exists(newPath))
			return newPath;
	}

	#if defined(__LINUX__)
		/* On Linux, also try case-insensitive mode if the above failed */
		fs::path parentPath = path.parent_path();
		std::string filename = boost::to_lower_copy(path.filename().string());

		for (size_t i=0; i<m_paths.size(); i++) {
			fs::path path = m_paths[i] / parentPath;

			if (!fs::is_directory(path))
				continue;

			fs::directory_iterator end, it(path);
			for (; it != end; ++it) {
				if (boost::algorithm::to_lower_copy(it->path().filename().string()) == filename)
					return it->path();
			}
		}
	#endif

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

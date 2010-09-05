/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__FRESOLVER_H)
#define __FRESOLVER_H

#include <mitsuba/core/fstream.h>
#include <mitsuba/core/tls.h>
#if defined(WIN32)
#include <direct.h>
#define PATH_MAX 1024
#else
#include <unistd.h>
#endif

#define PATHSEP '/'

MTS_NAMESPACE_BEGIN

/**
 * Utility class that searchs for files within a set of user-defined
 * directories. Each thread has its own, unique file resolver, which
 * makes it possible to have per-thread current working directories.
 */
class MTS_EXPORT_CORE FileResolver : public Object {
public:
	/// Return the current thread's resolver instance
	inline static FileResolver *getInstance() {
		FileResolver *resolver = m_tls.get();
		if (!resolver) {
			resolver = new FileResolver();
			m_tls.set(resolver);
		}
		return resolver;
	}

	/// Set the current thread's resolver instance
	inline static void setInstance(FileResolver *resolver) {
		m_tls.set(resolver);
	}

	/// Clone a file resolver
	inline FileResolver *clone() {
		FileResolver *cloned = new FileResolver();
		cloned->m_paths = m_paths;
		cloned->m_currentDirectory = m_currentDirectory;
		return cloned;
	}

	/// Reset the search path
	inline void clear() {
		m_paths.clear();
	}

	/// Add a search path
	inline void addPath(const std::string &path) {
		m_paths.push_back(path);
	}

	/// Check if the list of search paths contains a certain entry
	inline bool contains(const std::string &path) const {
		return std::find(m_paths.begin(), m_paths.end(), path) != m_paths.end();
	}

	/// Set the resolver's current working directory
	inline void setCurrentDirectory(const std::string &cwd) {
		m_currentDirectory = cwd;
	}

	/// Return the resolver's current working directory
	inline const std::string &getCurrentDirectory() const {
		return m_currentDirectory;
	}
	
	/// Set the current directory by stripping away the last component
	inline void setCurrentDirectoryFromFile(const std::string &cwd) {
		m_currentDirectory = getParentDirectory(cwd);
	}

	/// Adds a path while stripping away the last component
	inline void addPathFromFile(const std::string &path) {
		addPath(getParentDirectory(path));
	}

	/// Strip the last component from a path
	inline std::string getParentDirectory(const std::string &path) const {
		std::vector<std::string> components = tokenize(path, "/\\");
		if (components.size() == 1)
			return ".";
		std::string newPath = components[0];
		for (unsigned int i=1; i<components.size() - 1; i++) 
			newPath += PATHSEP + components[i];
		if (path[0] == PATHSEP)
			newPath = PATHSEP + newPath;
		else if (path[0] == '\\')
			newPath = "\\" + newPath;
		return newPath;
	}

	/// Return the last component of a path
	inline std::string getChild(const std::string &path) const {
		std::vector<std::string> components = tokenize(path, "/\\");
		return components[components.size()-1];
	}

	/// Try to resolve a path to an existing file
	inline std::string resolve(const std::string &path) const {
		if (!FileStream::exists(path)) {
			for (unsigned int i=0; i<m_paths.size(); i++) {
				std::string testPath = m_paths[i] + "/" + path;
				if (FileStream::exists(testPath))
					return testPath;
			}
		}
		return path;
	}
	
	/// Try to resolve all paths to an existing file
	inline std::vector<std::string> resolveAll(const std::string &path) const {
		std::vector<std::string> result;

		if (FileStream::exists(path)) 
			result.push_back(path);

		for (unsigned int i=0; i<m_paths.size(); i++) {
			std::string testPath = m_paths[i] + "/" + path;
			if (FileStream::exists(testPath))
				result.push_back(testPath);
		}
		return result;
	}

	/**
	 * Try to resolve a path to an existing file (returns an
	 * absolute path
	 */
	inline std::string resolveAbsolute(const std::string &path) const {
		return makeAbsolute(resolve(path));
	}
	
	/// Try to resolve all paths to an existing file (returns absolute paths)
	inline std::vector<std::string> resolveAllAbsolute(const std::string &path) const {
		std::vector<std::string> result;

		if (FileStream::exists(path)) 
			result.push_back(makeAbsolute(path));

		for (unsigned int i=0; i<m_paths.size(); i++) {
			std::string testPath = m_paths[i] + "/" + path;
			if (FileStream::exists(testPath))
				result.push_back(makeAbsolute(testPath));
		}
		return result;
	}

	/**
	 * Create a path for a new file. If the path is relative, it
	 * is appended to the thread's current working directory
	 */
	inline std::string resolveDest(const std::string &path) const {
		if (isAbsolute(path))
			return path;
		else
			return m_currentDirectory + "/" + path;
	}

	/**
	 * Return the filename (e.g. last component) of an absolute
	 * or relative path
	 */
	inline static std::string getFilename(const std::string &path) {
		std::vector<std::string> components = tokenize(path, "/\\");
		if (components.size() == 0)
			return "";
		return components[components.size()-1];
	}
	
	/**
	 * Return the filename of an absolute or relative path. This
	 * version also removes the file extension (everything after
	 * the first period sign)
	 */
	inline static std::string getFilenameWithoutExtension(const std::string &path) {
		std::string filename = getFilename(path);
		size_t pos = filename.find('.');
		if (pos != std::string::npos)
			filename = filename.substr(0, pos);
		return filename;
	}
	 
	inline ~FileResolver() { }

	inline std::string toString() const {
		std::ostringstream oss;
		oss << "FileResolver[" << endl
			<< "  cwd = \"" << m_currentDirectory << "\"," << endl
			<< "  paths = {" << endl;
		for (size_t i=0; i<m_paths.size(); ++i)
			oss << "    \"" << m_paths[i] << "\"," << endl;
		oss << "  }" << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
	inline FileResolver() { }
		
	/// Check whether a given path is absolute
	bool isAbsolute(const std::string &path) const {
#if defined(__OSX__) || defined(__LINUX__)
		return (path.length() > 0 &&
			(path[0] == '/' || path[0] == '~'));
#else
		return (strchr(path.c_str(), ':') != NULL) ||
			(path.length() > 0 && path[0] == '\\');
#endif
	}

	/// Turn a path into an absolute path
	std::string makeAbsolute(const std::string &path) const {
		if (isAbsolute(path))
			return path;
		char cwd[PATH_MAX];
#if !defined(WIN32)
		if (getcwd(cwd, PATH_MAX) == NULL)
#else
		if (_getcwd(cwd, PATH_MAX) == NULL)
#endif
			Log(EError, "Cannot detect current path!");
		return std::string(cwd) + PATHSEP + path;
	}
private:
	static ThreadLocal<FileResolver> m_tls;
	std::vector<std::string> m_paths;
	std::string m_currentDirectory;
};

MTS_NAMESPACE_END

#endif /* __FRESOLVER_H */

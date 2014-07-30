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

#pragma once
#if !defined(__MITSUBA_CORE_VERSION_H_)
#define __MITSUBA_CORE_VERSION_H_

MTS_NAMESPACE_BEGIN

/**
 * \brief Current release of Mitsuba
 * \ingroup libcore
 */
#define MTS_VERSION "0.5.0"

/**
 * \brief Year of the current release
 * \ingroup libcore
 */
#define MTS_YEAR "2014"

/**
 * \brief A simple data structure for representing and
 * comparing Mitsuba version strings
 *
 * \ingroup libcore
 */
struct MTS_EXPORT_CORE Version {
public:
	/// Default constructor: initialize to an invalid version (0.0.0)
	inline Version() : m_major(0), m_minor(0), m_release(0) { }

	/// Initialize with the specified version number
	inline Version(int major, int minor, int release)
		: m_major(major), m_minor(minor), m_release(release) { }

	/**
	 * \brief Parse a version string of the form "major.minor.release"
	 * and turn it into a \ref Version structure
	 */
	Version(const std::string &versionString);

	/// Check if this program version is \a older than \c other
	inline bool operator<(const Version &other) const {
		if (m_major < other.m_major)
			return true;
		else if (m_major > other.m_major)
			return false;
		else if (m_minor < other.m_minor)
			return true;
		else if (m_minor > other.m_minor)
			return false;
		else if (m_release < other.m_release)
			return true;
		else
			return false;
	}

	/// Check if this program version is \a older than or equal to \c other
	inline bool operator<=(const Version &other) const {
		return *this < other || *this == other;
	}

	/// Check if two program versions match
	inline bool operator==(const Version &other) const {
		return m_major == other.m_major
			&& m_minor == other.m_minor
			&& m_release == other.m_release;
	}

	/// Is this a valid version number?
	inline bool isValid() {
		return m_major != 0 || m_minor != 0 || m_release != 0;
	}

	/// Are the following two versions compatible?
	inline bool isCompatible(const Version &other) const {
		return m_major == other.m_major &&
			m_minor == other.m_minor;
	}

	/// Turn into a string of the form "major.minor.release"
	std::string toString() const;

	/// Turn into a string of the form "major.minor.release (Architecture)"
	std::string toStringComplete() const;

	/// Return the major version
	inline int getMajorVersion() const { return m_major; }

	/// Return the minor version
	inline int getMinorVersion() const { return m_minor; }

	/// Return the release
	inline int getRelease() const { return m_release; }
private:
	int m_major;
	int m_minor;
	int m_release;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_VERSION_H_ */

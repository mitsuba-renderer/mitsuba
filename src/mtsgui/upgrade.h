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

#if !defined(__UPGRADEMGR_H)
#define __UPGRADEMGR_H

#include "common.h"
#include <boost/filesystem.hpp>

class UpgradeManager : public QObject {
    Q_OBJECT
public:
	UpgradeManager(const FileResolver *resolver);

	void performUpgrade(const QString &path, const Version &version);
private:
	const FileResolver *m_resolver;
	std::vector<std::pair<Version, fs::path> > m_transformations;
};

#endif // __UPGRADEMGR_H

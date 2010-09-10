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

#if !defined(__SCENELOADER_H)
#define __SCENELOADER_H

#include <QtGui>
#include <mitsuba/core/lock.h>

using namespace mitsuba;

class SceneImporter : public Thread {
public:
	SceneImporter(QWidget *parent, FileResolver *resolver, 
		const std::string &sourceFile, const std::string &directory,
		const std::string &targetScene, const std::string &adjustmentFile,
		bool sRGB);

	void run();

	inline void wait(int ms) { m_wait->wait(ms); }

	inline const std::string &getResult() const { return m_result; }
protected:
	virtual ~SceneImporter();
private:
	QWidget *m_parent;
	ref<FileResolver> m_resolver;
	ref<WaitFlag> m_wait;
	std::string m_sourceFile;
	std::string m_directory;
	std::string m_targetScene;
	std::string m_adjustmentFile;
	std::string m_result;
	bool m_srgb;
};

#endif // __SCENELOADER_H


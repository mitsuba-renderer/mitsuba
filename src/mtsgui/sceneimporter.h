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

#if !defined(__SCENELOADER_H)
#define __SCENELOADER_H

#include <QtGui/QtGui>
#include <mitsuba/core/lock.h>
#include <mitsuba/core/fresolver.h>
#include "../converter/converter.h"

using namespace mitsuba;

class GUIGeometryConverter : public QObject, public GeometryConverter {
    Q_OBJECT
public:
    inline GUIGeometryConverter() { }
protected:
    fs::path locateResource(const fs::path &resource);
signals:
    void locateResource(const fs::path &resource, fs::path *result);
};

class SceneImporter : public Thread {
public:
    SceneImporter(FileResolver *resolver,
        const fs::path &sourceFile, const fs::path &directory,
        const fs::path &targetScene, const fs::path &adjustmentFile,
        bool sRGB);

    void run();

    inline void wait(int ms) { m_wait->wait(ms); }

    inline GUIGeometryConverter *getConverter() { return &m_converter; }

    inline const fs::path &getResult() const { return m_result; }
protected:
    virtual ~SceneImporter();

private:
    GUIGeometryConverter m_converter;
    ref<FileResolver> m_resolver;
    ref<WaitFlag> m_wait;
    fs::path m_sourceFile;
    fs::path m_directory;
    fs::path m_targetScene;
    fs::path m_adjustmentFile;
    fs::path m_result;
    bool m_srgb;
};

#endif // __SCENELOADER_H


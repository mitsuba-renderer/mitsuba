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

#include "sceneimporter.h"
#include "../converter/converter.h"
#include "locateresourcedlg.h"

class GUIGeometryConverter : public GeometryConverter {
public:
	inline GUIGeometryConverter(QWidget *parent) : m_parent(parent) {
	}

	std::string locateResource(const std::string &resource) {
		LocateResourceDialog locateResource(m_parent, resource.c_str());
		locateResource.setWindowModality(Qt::ApplicationModal);
		if (locateResource.exec()) 
			return locateResource.getFilename().toStdString();

		return "";
	}
private:
	QWidget *m_parent;
};

SceneImporter::SceneImporter(QWidget *parent, FileResolver *resolver, 
		const std::string &sourceFile, const std::string &directory,
		const std::string &targetScene, const std::string &adjustmentFile,
		bool sRGB)
	: Thread("impt"), m_parent(parent), m_resolver(resolver), 
	  m_sourceFile(sourceFile), m_directory(directory), 
	  m_targetScene(targetScene), m_adjustmentFile(adjustmentFile), m_srgb(sRGB) {
	m_wait = new WaitFlag();
}

SceneImporter::~SceneImporter() {
}

void SceneImporter::run() {
	Thread::getThread()->setFileResolver(m_resolver);
#if defined(MTS_HAS_COLLADA)
	try {
		GUIGeometryConverter cvt(m_parent);
		cvt.setSRGB(m_srgb);
		cvt.convert(m_sourceFile, m_directory, m_targetScene, m_adjustmentFile);
		m_result = cvt.getFilename();
	} catch (const std::exception &ex) {
		SLog(EWarn, "Conversion failed: %s", ex.what());
	} catch (...) {
		SLog(EWarn, "An unknown type of error occurred!");
	}
#else
	SLog(EWarn, "The importer was disabled in this build!");
#endif
	m_wait->set(true);
}


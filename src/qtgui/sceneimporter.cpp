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
	FileResolver::setInstance(m_resolver);
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
	m_wait->set(true);
}


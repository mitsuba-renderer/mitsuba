#if !defined(__SCENELOADER_H)
#define __SCENELOADER_H

#include <QtGui>
#include <mitsuba/core/fresolver.h>
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


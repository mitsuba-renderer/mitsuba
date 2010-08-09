#if !defined(__SCENELOADER_H)
#define __SCENELOADER_H

#include <QtGui>
#include <mitsuba/mitsuba.h>
#include <mitsuba/core/fresolver.h>

class SceneContext;

using namespace mitsuba;

class SceneLoader : public Thread {
public:
	SceneLoader(FileResolver *resolver, const std::string &filename); 
	void run();

	inline void wait(int ms) { m_wait->wait(ms); }

	inline SceneContext *getResult() { return m_result; }
	inline const std::string &getError() const { return m_error; }
protected:
	virtual ~SceneLoader();
private:
	ref<FileResolver> m_resolver;
	ref<WaitFlag> m_wait;
	SceneContext *m_result;
	std::string m_error, m_filename;
};

#endif // __SCENELOADER_H


#if !defined(__LOGWIDGET_H)
#define __LOGWIDGET_H

#include "common.h"
#include <mitsuba/mitsuba.h>

namespace mitsuba {
	class RenderJob;
};


class QConsoleAppender : public QObject, public Appender {
	Q_OBJECT
public:
    void append(ELogLevel level, const std::string &message) {
		emit textMessage(level, QString(message.c_str()));
	}
    void logProgress(Float progress, const std::string &name, 
		const std::string &formatted, const std::string &eta,
		const void *ptr) {
		emit progressMessage((RenderJob *) ptr, 
			QString(name.c_str()), (float) progress, QString(eta.c_str()));
	}
signals:
	void textMessage(ELogLevel level, const QString &message);
	void progressMessage(const RenderJob *job, const QString &name, 
		float progress, const QString &eta);
};

class LogWidget : public QMainWindow {
	Q_OBJECT
public:
	LogWidget(QWidget *parent);
	virtual ~LogWidget();
	void show();
public slots:
	void onTextMessage(ELogLevel level, const QString &message);
	void onShowStats();
private:
	QTextEdit *m_contents;
};

#endif /* __LOGWIDGET_H */

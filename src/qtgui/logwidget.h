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

#if !defined(__LOGWIDGET_H)
#define __LOGWIDGET_H

#include "common.h"
#include <mitsuba/core/appender.h>

namespace mitsuba {
	class RenderJob;
};


class QConsoleAppender : public QObject, public Appender {
	Q_OBJECT
public:
    void append(ELogLevel level, const std::string &message) {
		emit textMessage(level, QString::fromLatin1(message.c_str()));
	}
    void logProgress(Float progress, const std::string &name, 
		const std::string &formatted, const std::string &eta,
		const void *ptr) {
		emit progressMessage((RenderJob *) ptr, 
			QString::fromLatin1(name.c_str()), (float) progress,
			QString::fromLatin1(eta.c_str()));
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

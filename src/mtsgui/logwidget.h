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
	QConsoleAppender() {
		m_timer = new Timer();
		m_ignoreMessages = false;
		m_messageCount = 0;
	}

    void append(ELogLevel level, const std::string &message) {
		if (!m_ignoreMessages) {
			emit textMessage(level, QString::fromUtf8(message.c_str()));
			if (level >= EWarn)
				floodCheck();
		}
	}

    void logProgress(Float progress, const std::string &name,
		const std::string &formatted, const std::string &eta,
		const void *ptr) {
		emit progressMessage((RenderJob *) ptr,
			QString::fromUtf8(name.c_str()), (float) progress,
			QString::fromUtf8(eta.c_str()));
	}

	void floodCheck() {
		++m_messageCount;
		int ms = m_timer->getMilliseconds();
		if (ms > 1000) {
			Float messagesPerSecond = m_messageCount / (ms / 1000.0f);
			if (messagesPerSecond > 1000) {
				emit textMessage(EError,
						QString("Flood alert: received %1 messages in %2 ms! Ignoring "
							"future messages to prevent the user interface from freezing. "
							"Note: this only concerns the UI, all messages will still be written to the log file.")
						.arg(m_messageCount).arg(m_timer->getMilliseconds()));
				m_ignoreMessages = true;
				emit criticalError("The console is being flooded with messages and will now shut down to "
						"prevent the user interface from freezing. Please refer to the console messages to "
						"find out what happened.");
			}
			m_messageCount = 0;
			m_timer->reset();
		}
	}
signals:
	void criticalError(const QString &message);
	void textMessage(ELogLevel level, const QString &message);
	void progressMessage(const RenderJob *job, const QString &name,
		float progress, const QString &eta);
private:
	ref<Timer> m_timer;
	bool m_ignoreMessages;
	size_t m_messageCount;
};

class LogWidget : public QMainWindow {
	Q_OBJECT
public:
	LogWidget(QWidget *parent);
	virtual ~LogWidget();
	void show();
public slots:
	void onTextMessage(ELogLevel level, const QString &message);
	void onCriticalError(const QString &message);
	void onShowStats();
private:
	QTextEdit *m_contents;
};

#endif /* __LOGWIDGET_H */

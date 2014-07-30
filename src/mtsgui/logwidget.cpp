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

#include "logwidget.h"
#include <mitsuba/core/statistics.h>

LogWidget::LogWidget(QWidget *parent)
 : QMainWindow(parent) {
	m_contents = new QTextEdit(this);
	QFont font("Monospace");
	font.setStyleHint(QFont::TypeWriter);
	m_contents->setFont(font);
	m_contents->setReadOnly(true);
	QPalette palette;
	palette.setColor(QPalette::Base, Qt::black);
	m_contents->setPalette(palette);
	setPalette(palette);
	QToolBar *toolBar = new QToolBar(this);
	toolBar->setMovable(false);
	toolBar->setAllowedAreas(Qt::TopToolBarArea);
	toolBar->setIconSize(QSize(32, 32));
	toolBar->setToolButtonStyle(Qt::ToolButtonIconOnly);
	toolBar->setFloatable(false);

	QAction *actionShowStats = new QAction(this);
	QIcon showStatsIcon;
	showStatsIcon.addFile(QString::fromUtf8(":/resources/showStats.png"), QSize(), QIcon::Normal, QIcon::Off);
	actionShowStats->setIcon(showStatsIcon);
    actionShowStats->setToolTip(tr("Show statistics"));
    actionShowStats->setText(tr("Show Statistics"));
	connect(actionShowStats, SIGNAL(triggered()), this, SLOT(onShowStats()));
	toolBar->addAction(actionShowStats);
	QAction *actionClear = new QAction(this);
	QIcon clearIcon;
	clearIcon.addFile(QString::fromUtf8(":/resources/clear.png"), QSize(), QIcon::Normal, QIcon::Off);
	actionClear->setIcon(clearIcon);
    actionClear->setToolTip(tr("Clear"));
    actionClear->setText(tr("Clear"));
	connect(actionClear, SIGNAL(triggered()), m_contents, SLOT(clear()));
	toolBar->addAction(actionClear);
#if defined(__OSX__)
	toolBar->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
#endif

    addToolBar(Qt::TopToolBarArea, toolBar);
	setCentralWidget(m_contents);
    setUnifiedTitleAndToolBarOnMac(true);
	setWindowTitle(tr("Log"));
	resize(QSize(1000, 500));
}

LogWidget::~LogWidget() {
}

void LogWidget::onCriticalError(const QString &message) {
	QMessageBox::critical(this, tr("Critical error"),
			message, QMessageBox::Ok);
}

void LogWidget::onTextMessage(ELogLevel level, const QString &message) {
	QColor color;
	switch (level) {
		case ETrace:
		case EDebug:
		case EInfo:
			color = Qt::gray;
			break;
		case EError:
			color = Qt::red;
			break;
		case EWarn:
		default:
			color = QColor(255,180,180);
			break;
	}
	m_contents->setTextColor(color);
	m_contents->append(message);
	QTextCursor c =  m_contents->textCursor();
	c.movePosition(QTextCursor::End);
	m_contents->setTextCursor(c);
	m_contents->ensureCursorVisible();
}

void LogWidget::show() {
	if (isVisible()) {
		raise();
	} else {
		/* Center the dialog */
		QDesktopWidget *desktop = QApplication::desktop();
		QRect geo = desktop->screenGeometry();
		QPoint windowPos(
			geo.left() + (geo.width() - width()) / 2,
			geo.top() + (geo.height() - height())/2
		);
		move(windowPos);
		QMainWindow::show();
	}
}

void LogWidget::onShowStats() {
	Statistics::getInstance()->printStats();
}

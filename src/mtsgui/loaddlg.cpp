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

#include "ui_loaddlg.h"
#include "loaddlg.h"

LoadDialog::LoadDialog(QWidget *parent) :
        QDialog(parent),
    ui(new Ui::LoadDialog) {
    ui->setupUi(this);
    ui->progressBar->setTextVisible(false);
    // weird, Qt/Win needs this to get a busy indicator
    ui->progressBar->setValue(1);
    ui->progressBar->setRange(0, 0);
    ui->console->hide();
    resize(400,180);
    m_consoleAppender = new QConsoleAppender();
    Logger *logger = Thread::getThread()->getLogger();
    m_oldLogLevel = logger->getLogLevel();
    //logger->setLogLevel(EDebug);
    logger->addAppender(m_consoleAppender);
    connect(m_consoleAppender, SIGNAL(textMessage(ELogLevel, const QString &)),
        this, SLOT(onTextMessage(ELogLevel, const QString &)), Qt::QueuedConnection);
    QFont font("Monospace");
    font.setStyleHint(QFont::TypeWriter);
#if defined(__OSX__)
    font.setPointSize(10);
    QFont headingFont = ui->heading->font();
    headingFont.setPointSize(14);
    ui->heading->setFont(headingFont);
#else
    font.setPointSize(8);
#endif
    ui->console->setFont(font);
    ui->console->setReadOnly(true);
    QPalette palette;
    palette.setColor(QPalette::Base, Qt::black);
    ui->console->setPalette(palette);
}

void LoadDialog::on_toggleButton_clicked() {
    if (ui->console->isVisible()) {
        ui->console->hide();
        updateGeometry();
        adjustSize();
        resize(400, 180);
        ui->toggleButton->setText("+");
    } else {
        ui->console->show();
        updateGeometry();
        resize(std::max(width(), 800), 500);
        ui->toggleButton->setText("-");
    }
}

void LoadDialog::expand() {
    if (!ui->console->isVisible()) {
        ui->console->show();
        updateGeometry();
        resize(std::max(width(), 800), 500);
        ui->toggleButton->setText("-");
    }
}

void LoadDialog::close() {
    Logger *logger = Thread::getThread()->getLogger();
    logger->setLogLevel(m_oldLogLevel);
    logger->removeAppender(m_consoleAppender);
    disconnect(m_consoleAppender, SIGNAL(textMessage(ELogLevel, const QString &)),
        this, SLOT(onTextMessage(ELogLevel, const QString &)));
    m_consoleAppender = NULL;
    if (ui->console->isVisible()) {
        ui->statusLabel->setText(tr("Done."));
        ui->closeButton->setEnabled(true);
        ui->progressBar->setRange(1, 100);
        ui->progressBar->setValue(100);
    } else {
        hide();
    }
}

void LoadDialog::onTextMessage(ELogLevel level, const QString &message) {
    QColor color;
    int idx = message.indexOf("] ");
    if (idx != -1 && ui->progressBar->value() != 100 && message.indexOf("\n") == -1) {
        ui->statusLabel->setText(message.mid(idx+2));
    }
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
    ui->console->setTextColor(color);
    ui->console->append(message);
    QTextCursor c =  ui->console->textCursor();
    c.movePosition(QTextCursor::End);
    ui->console->setTextCursor(c);
    ui->console->ensureCursorVisible();
}

LoadDialog::~LoadDialog() {
    delete ui;
    delete m_consoleAppender;
}

void LoadDialog::changeEvent(QEvent *e) {
    QDialog::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

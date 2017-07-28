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

#if !defined(__LOADDLG_H)
#define __LOADDLG_H

#include "common.h"
#include "logwidget.h"

namespace Ui {
    class LoadDialog;
}

class LoadDialog : public QDialog {
    Q_OBJECT
public:
    LoadDialog(QWidget *parent);
    ~LoadDialog();

    void close();
    void closeEvent(QCloseEvent *e) {
        e->ignore();
    }
    void expand();
protected slots:
    void onTextMessage(ELogLevel level, const QString &message);
    void on_toggleButton_clicked();
protected:
    void changeEvent(QEvent *e);
private:
    Ui::LoadDialog *ui;
    ref<QConsoleAppender> m_consoleAppender;
    ELogLevel m_oldLogLevel;
};

#endif // __LOADDLG_H

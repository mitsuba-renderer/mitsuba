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

#if !defined(__ADDSERVERDLG_H)
#define __ADDSERVERDLG_H

#include "common.h"

namespace Ui {
    class AddServerDialog;
}

class AddServerDialog : public QDialog {
    Q_OBJECT
public:
    AddServerDialog(QWidget *parent);
    ~AddServerDialog();

    ServerConnection getConnection() const;
protected slots:
    void on_sshConnection_toggled();
    void on_directConnection_toggled();
    void validate();
protected:
    void changeEvent(QEvent *e);
private:
    Ui::AddServerDialog *ui;
};

#endif // __ADDSERVERDLG_H

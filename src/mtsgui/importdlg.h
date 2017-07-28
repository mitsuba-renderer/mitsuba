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

#if !defined(__IMPORTDLG_H)
#define __IMPORTDLG_H

#include "common.h"

namespace Ui {
    class ImportDialog;
}

class ImportDialog : public QDialog {
    Q_OBJECT
public:
    ImportDialog(QWidget *parent, FileResolver *resolver);
    ~ImportDialog();
public slots:
    void accept();
protected slots:
    void on_inputBrowse_clicked(bool checked);
    void on_directoryBrowse_clicked(bool checked);
    void on_adjustmentBrowse_clicked(bool checked);
    void onLocateResource(const fs::path &path, fs::path *target);
    void refresh();
protected:
    void changeEvent(QEvent *e);
private:
    Ui::ImportDialog *ui;
    ref<FileResolver> m_resolver;
};

#endif // __IMPORTDLG_H

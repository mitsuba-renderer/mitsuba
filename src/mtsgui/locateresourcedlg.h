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

#if !defined(__LOCATERESOURCEDLG_H)
#define __LOCATERESOURCEDLG_H

#include "common.h"

namespace Ui {
	class LocateResourceDialog;
}

class LocateResourceDialog : public QDialog {
    Q_OBJECT
public:
	LocateResourceDialog(QWidget *parent, const QString &resourceName);
	~LocateResourceDialog();

	inline const QString &getFilename() const {
		return m_filename;
	}
protected slots:
	void on_pathBrowse_clicked();
protected:
    void changeEvent(QEvent *e);
private:
	Ui::LocateResourceDialog *ui;
	QString m_filename;
};

#endif // __LOCATERESOURCEDLG_H

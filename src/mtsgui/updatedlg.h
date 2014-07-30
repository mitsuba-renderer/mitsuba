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

#if !defined(__UPDATEDLG_H)
#define __UPDATEDLG_H

#include "common.h"
#include <QtNetwork/QtNetwork>

namespace Ui {
	class UpdateDialog;
}

class UpdateDialog : public QDialog {
    Q_OBJECT
public:
	UpdateDialog(QWidget *parent, const Version &local,
			const Version &remote);
	~UpdateDialog();
protected slots:
	void on_skipButton_clicked();
	void on_remindButton_clicked();
	void on_downloadButton_clicked();
	void onNetworkFinished(QNetworkReply *reply);
protected:
    void changeEvent(QEvent *e);
private:
	Ui::UpdateDialog *ui;
	QNetworkAccessManager *m_networkManager;
	QNetworkReply *m_networkReply;
	QString m_remoteVersion;
};

#endif // __UPDATEDLG_H

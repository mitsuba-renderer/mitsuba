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

#include "ui_locateresourcedlg.h"
#include "locateresourcedlg.h"

LocateResourceDialog::LocateResourceDialog(QWidget *parent, const QString &resName) :
		QDialog(parent),
	ui(new Ui::LocateResourceDialog) {
	ui->setupUi(this);
	ui->resourceName->setText(resName);
	ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(false);
}

void LocateResourceDialog::on_pathBrowse_clicked() {
	QString fname;
#if MTSGUI_STATIC_QFILEDIALOG
	QSettings settings;
	fname = QFileDialog::getOpenFileName(this, QString(),
		settings.value("importDir").toString());
#else
	QFileDialog dialog(this);
	dialog.setAcceptMode(QFileDialog::AcceptOpen);
	dialog.setViewMode(QFileDialog::Detail);
	dialog.setWindowModality(Qt::ApplicationModal);
	if (dialog.exec()) {
		fname = dialog.selectedFiles()[0];
	}
#endif
	if (!fname.isEmpty()) {
		ui->pathEdit->setText(fname);
		m_filename = fname;
		ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(true);
	}
}

LocateResourceDialog::~LocateResourceDialog() {
	delete ui;
}

void LocateResourceDialog::changeEvent(QEvent *e) {
	QDialog::changeEvent(e);
	switch (e->type()) {
	case QEvent::LanguageChange:
		ui->retranslateUi(this);
		break;
	default:
		break;
	}
}

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

#include "ui_addserverdlg.h"
#include "addserverdlg.h"

#if !defined(__WINDOWS__)
#include <pwd.h>
#endif

AddServerDialog::AddServerDialog(QWidget *parent) :
		QDialog(parent),
	ui(new Ui::AddServerDialog) {
	ui->setupUi(this);
#if !defined(__WINDOWS__)
	uid_t uid = getuid();
	struct passwd *info = getpwuid(uid);
	if (info)
		ui->userName->setText(QString(info->pw_name));
#endif
	ui->port->setValidator(new QIntValidator(this));
}

AddServerDialog::~AddServerDialog() {
	delete ui;
}

void AddServerDialog::on_sshConnection_toggled() {
	if (!ui->sshConnection->isChecked())
		return;
	ui->port->setText("22");
	ui->userName->setEnabled(true);
	ui->installDir->setEnabled(true);
}

void AddServerDialog::on_directConnection_toggled() {
	if (!ui->directConnection->isChecked())
		return;
	ui->port->setText("7554");
	ui->userName->setEnabled(false);
	ui->installDir->setEnabled(false);
}

void AddServerDialog::validate() {
	bool valid = true;

	valid &= !ui->hostName->text().isEmpty();
	valid &= !ui->port->text().isEmpty();
	valid &= ui->port->text().toInt() > 1;
	valid &= ui->port->text().toInt() < 65536;

	if (ui->sshConnection->isChecked()) {
		valid &= !ui->userName->text().isEmpty();
		valid &= !ui->installDir->text().isEmpty();
	}

	ui->buttons->button(QDialogButtonBox::Ok)->setEnabled(valid);
}

ServerConnection AddServerDialog::getConnection() const {
	ServerConnection conn;
	conn.hostName = ui->hostName->text();
	conn.userName = ui->userName->text();
	conn.instDir = ui->installDir->text();
	conn.port = ui->port->text().toInt();
	conn.type = ui->directConnection->isChecked()
		? EDirectConnection : ESSHConnection;
	return conn;
}

void AddServerDialog::changeEvent(QEvent *e) {
    QDialog::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
		ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

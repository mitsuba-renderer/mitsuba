#include "ui_addserverdlg.h"
#include "addserverdlg.h"

#if !defined(WIN32)
#include <pwd.h>
#endif

AddServerDialog::AddServerDialog(QWidget *parent) :
		QDialog(parent),
	ui(new Ui::AddServerDialog) {
	ui->setupUi(this);
#if !defined(WIN32)
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

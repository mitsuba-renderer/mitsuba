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
	QFileDialog dialog(this);
	dialog.setAcceptMode(QFileDialog::AcceptOpen);
	dialog.setViewMode(QFileDialog::Detail);
	dialog.setWindowModality(Qt::ApplicationModal);

	if (dialog.exec()) {
		QString fname = dialog.selectedFiles()[0];
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

#include "ui_importdlg.h"
#include "importdlg.h"
#include "acknowledgmentdlg.h"

ImportDialog::ImportDialog(QWidget *parent) :
		QDialog(parent, Qt::Sheet),
	ui(new Ui::ImportDialog) {
	ui->setupUi(this);
}

ImportDialog::~ImportDialog() {
	delete ui;
}
	
void ImportDialog::changeEvent(QEvent *e) {
    QDialog::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
		ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

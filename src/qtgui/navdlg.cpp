#include "ui_navdlg.h"
#include "navdlg.h"

NavigationDialog::NavigationDialog(QWidget *parent) :
		QDialog(parent),
	ui(new Ui::NavigationDialog) {
	ui->setupUi(this);
}

NavigationDialog::~NavigationDialog() {
	delete ui;
}
	
void NavigationDialog::changeEvent(QEvent *e) {
    QDialog::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
		ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

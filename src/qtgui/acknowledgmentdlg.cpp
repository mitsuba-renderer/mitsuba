#include "ui_acknowledgmentdlg.h"
#include "acknowledgmentdlg.h"

AcknowledgmentDialog::AcknowledgmentDialog(QWidget *parent) :
		QDialog(parent),
	ui(new Ui::AcknowledgmentDialog) {
	ui->setupUi(this);
#if defined(__OSX__)
	ui->textBrowser->setHtml(ui->textBrowser->toHtml().replace("font-size:10pt", "font-size:14pt"));
#endif
}

AcknowledgmentDialog::~AcknowledgmentDialog() {
	delete ui;
}

void AcknowledgmentDialog::changeEvent(QEvent *e) {
    QDialog::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
		ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

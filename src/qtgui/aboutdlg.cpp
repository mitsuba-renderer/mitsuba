#include "ui_aboutdlg.h"
#include "aboutdlg.h"
#include "acknowledgmentdlg.h"

AboutDialog::AboutDialog(QWidget *parent) :
		QDialog(parent),
	ui(new Ui::AboutDialog) {
	ui->setupUi(this);
	QString configFlags;

#if defined(SINGLE_PRECISION)
	configFlags += "SINGLE_PRECISION ";
#elif defined(DOUBLE_PRECISION)
	configFlags += "DOUBLE_PRECISION ";
#else 
#error Unknown precision
#endif

#if defined(MTS_DEBUG)
	configFlags += "MTS_DEBUG ";
#endif

#if defined(MTS_DEBUG_FP)
	configFlags += "MTS_DEBUG_FP ";
#endif

#if defined(MTS_SSE)
	configFlags += "MTS_SSE ";
#endif

#if defined(MTS_HAS_COHERENT_RT)
	configFlags += "MTS_HAS_COHERENT_RT ";
#endif

	ui->label1->setText(ui->label1->text().replace("MTS_VERSION", MTS_VERSION));
	ui->label1->setText(ui->label1->text().replace("CONFIG_FLAGS", configFlags));
#if defined(__OSX__)
	ui->label1->setText(ui->label1->text().replace("font-size:10pt", "font-size:14pt"));
#endif
}

AboutDialog::~AboutDialog() {
	delete ui;
}
	
void AboutDialog::onCredits() {
	AcknowledgmentDialog ackdlg(this);
	ackdlg.exec();
}

void AboutDialog::changeEvent(QEvent *e) {
    QDialog::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
		ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

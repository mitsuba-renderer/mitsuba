#include "ui_importdlg.h"
#include "importdlg.h"
#include "acknowledgmentdlg.h"
#include "../converter/converter.h"
#include "mainwindow.h"
#include "locateresourcedlg.h"

class GUIGeometryConverter : public GeometryConverter {
public:
	inline GUIGeometryConverter(QWidget *parent) : m_parent(parent) {
	}

	std::string locateResource(const std::string &resource) {
		LocateResourceDialog locateResource(m_parent, resource.c_str());
		locateResource.setWindowModality(Qt::ApplicationModal);
		if (locateResource.exec()) 
			return locateResource.getFilename().toStdString();

		return "";
	}
private:
	QWidget *m_parent;
};

ImportDialog::ImportDialog(QWidget *parent) :
		QDialog(parent, Qt::Sheet),
	ui(new Ui::ImportDialog) {
	ui->setupUi(this);
	connect(ui->sceneEdit, SIGNAL(textChanged(const QString &)),
		this, SLOT(refresh()));
	refresh();
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

void ImportDialog::on_inputBrowse_clicked(bool checked) {
	QFileDialog dialog(this);
	dialog.setNameFilter(tr("COLLADA 1.4 scenes (*.dae);; Wavefront OBJ scenes (*.obj)"));
	dialog.setAcceptMode(QFileDialog::AcceptOpen);
	dialog.setViewMode(QFileDialog::Detail);
	dialog.setWindowModality(Qt::ApplicationModal);

	if (dialog.exec()) {
		QString fname = dialog.selectedFiles()[0];
		ui->inputEdit->setText(fname);
		QFileInfo info(fname);
		ui->directoryEdit->setText(info.absoluteDir().absolutePath());
		ui->sceneEdit->setText(info.completeBaseName() + ".xml");
		refresh();
	}
}

void ImportDialog::on_directoryBrowse_clicked(bool checked) {
	QFileDialog dialog(this);
	dialog.setAcceptMode(QFileDialog::AcceptOpen);
	dialog.setFileMode(QFileDialog::DirectoryOnly);
	dialog.setWindowModality(Qt::ApplicationModal);
	if (dialog.exec()) {
		QString fname = dialog.selectedFiles()[0];
		ui->directoryEdit->setText(fname);
		refresh();
	}
}

void ImportDialog::on_adjustmentBrowse_clicked(bool checked) {
	QFileDialog dialog(this);
	dialog.setNameFilter(tr("Import adjustment files (*.xml)"));
	dialog.setAcceptMode(QFileDialog::AcceptOpen);
	dialog.setViewMode(QFileDialog::Detail);
	dialog.setWindowModality(Qt::ApplicationModal);
	if (dialog.exec()) {
		QString fname = dialog.selectedFiles()[0];
		ui->adjustmentEdit->setText(fname);
		refresh();
	}
}

void ImportDialog::refresh() {
	bool hasInput = ui->inputEdit->text() != "";
	bool hasOutput = ui->sceneEdit->text().endsWith(".xml");

	ui->directoryBrowse->setEnabled(hasInput);
	ui->directoryEdit->setEnabled(hasInput);
	ui->adjustmentBrowse->setEnabled(hasInput);
	ui->adjustmentEdit->setEnabled(hasInput);
	ui->sceneEdit->setEnabled(hasInput);
	ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(hasInput && hasOutput);
}

void ImportDialog::accept() {
	QDialog::accept();

	QString sourceFile = ui->inputEdit->text();
	QString directory = ui->directoryEdit->text();
	QString targetScene = ui->sceneEdit->text();
	QString adjustmentFile = ui->adjustmentEdit->text();

	QDialog *dialog = new QDialog(this);
	dialog->setWindowModality(Qt::WindowModal);
	dialog->setWindowTitle("Converting ..");
	QVBoxLayout *layout = new QVBoxLayout(dialog);
	QProgressBar *progressBar = new QProgressBar(dialog);
	dialog->resize(200, 50);
	layout->addWidget(progressBar);
	progressBar->setTextVisible(false);
	progressBar->setValue(1);
	progressBar->setRange(0, 0);
	dialog->show();
	progressBar->show();

	for (int i=0; i<10; ++i)
		QCoreApplication::processEvents();

	GUIGeometryConverter cvt(this);
	cvt.setSRGB(ui->sRGBButton->isChecked());

	try {
		cvt.convert(sourceFile.toStdString(), directory.toStdString(),
			targetScene.toStdString(), adjustmentFile.toStdString());
		((MainWindow *) parent())->loadFile(QString(cvt.getFilename().c_str()));
		dialog->hide();
		delete dialog;
	} catch (const std::exception &ex) {
		dialog->hide();
		delete dialog;
		SLog(EWarn, "Conversion failed: %s", ex.what());
		QMessageBox::critical(this, tr("Scene Import"),
			tr("Conversion failed -- please see the log for details."),
			QMessageBox::Ok);
	}
}

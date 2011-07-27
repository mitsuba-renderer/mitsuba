/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#include "ui_importdlg.h"
#include "importdlg.h"
#include "acknowledgmentdlg.h"
#include "mainwindow.h"
#include "sceneimporter.h"
#include "locateresourcedlg.h"

ImportDialog::ImportDialog(QWidget *parent, FileResolver *resolver) :
	QDialog(parent, Qt::Sheet), ui(new Ui::ImportDialog), m_resolver(resolver) {
	ui->setupUi(this);
	connect(ui->sceneEdit, SIGNAL(textChanged(const QString &)),
		this, SLOT(refresh()));
	refresh();

#if defined(__OSX__)
	ui->inputBrowse->setMaximumSize(24, 24);
	ui->directoryBrowse->setMaximumSize(24, 24);
	ui->adjustmentBrowse->setMaximumSize(24, 24);
#endif
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
	dialog.setNameFilter(tr("All supported formats (*.dae *.zae *.obj);;"
		"COLLADA 1.4 scenes (*.dae *.zae);; Wavefront OBJ scenes (*.obj)"));
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

	NonClosableDialog *dialog = new NonClosableDialog(static_cast<QWidget *>(parent()));
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

	fs::path filePath = fs::complete(fs::path(sourceFile.toStdString())).parent_path();
	ref<FileResolver> resolver = m_resolver->clone();
	resolver->addPath(filePath);

	const Logger *logger = Thread::getThread()->getLogger();
	size_t initialWarningCount = logger->getWarningCount();

	ref<SceneImporter> importingThread = new SceneImporter(
		resolver, sourceFile.toStdString(), directory.toStdString(),
		targetScene.toStdString(), adjustmentFile.toStdString(),
		ui->sRGBButton->isChecked());
	importingThread->start();

	connect(importingThread->getConverter(),
			SIGNAL(locateResource(const fs::path &, fs::path *)),
			this, SLOT(onLocateResource(const fs::path &, fs::path *)),
			Qt::BlockingQueuedConnection);

	while (importingThread->isRunning()) {
		QCoreApplication::processEvents();
		importingThread->wait(20);
	}

	importingThread->join();

	dialog->hide();
	delete dialog;
	
	if (!importingThread->getResult().empty()) {
		size_t warningCount = logger->getWarningCount() - initialWarningCount;
		if (warningCount > 0)
			QMessageBox::warning(this, tr("Scene Import"),
				tr("Encountered %1 warnings while importing -- please see "
				"the log for details.").arg(warningCount), QMessageBox::Ok);
		((MainWindow *) parent())->loadFile(QString(importingThread->getResult().file_string().c_str()));
	} else {
		QMessageBox::critical(this, tr("Scene Import"),
			tr("Conversion failed -- please see the log for details."),
			QMessageBox::Ok);
	}
}

void ImportDialog::onLocateResource(const fs::path &path, fs::path *target) {
	LocateResourceDialog locateResource(this, path.file_string().c_str());
	locateResource.setWindowModality(Qt::ApplicationModal);
	if (locateResource.exec())
		*target = fs::path(locateResource.getFilename().toStdString());
}

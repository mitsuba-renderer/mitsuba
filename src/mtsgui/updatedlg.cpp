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

#include "ui_updatedlg.h"
#include "updatedlg.h"

UpdateDialog::UpdateDialog(QWidget *parent, const Version &local,
            const Version &remote) : QDialog(parent),
    ui(new Ui::UpdateDialog) {
    ui->setupUi(this);
    m_remoteVersion = remote.toString().c_str();
    ui->versionLabel->setText(QApplication::translate("UpdateDialog",
            "Version %1 has been released (you are using %2). Would you like to visit the download page?",
            0).arg(m_remoteVersion).arg(local.toString().c_str()));
    ui->changeView->setHtml("Loading change log ..");
    m_networkManager = new QNetworkAccessManager(this);
    connect(m_networkManager, SIGNAL(finished(QNetworkReply *)),
            this, SLOT(onNetworkFinished(QNetworkReply *)));
    m_networkReply = m_networkManager->get(QNetworkRequest(QUrl("http://www.mitsuba-renderer.org/changelog.html")));
}

void UpdateDialog::onNetworkFinished(QNetworkReply *reply) {
    if (reply->error() == QNetworkReply::NoError)
        ui->changeView->setHtml(QString(reply->readAll()));
    else
        ui->changeView->setHtml("Unable to load the change log!");
    m_networkReply = NULL;
}

UpdateDialog::~UpdateDialog() {
    if (m_networkReply)
        m_networkReply->abort();
    delete ui;
}

void UpdateDialog::on_skipButton_clicked() {
    QSettings settings;
    settings.setValue("ignoredVersion", m_remoteVersion);
    accept();
}

void UpdateDialog::on_remindButton_clicked() {
    accept();
}

void UpdateDialog::on_downloadButton_clicked() {
    QDesktopServices::openUrl(QUrl("http://www.mitsuba-renderer.org/download.html"));
    accept();
}

void UpdateDialog::changeEvent(QEvent *e) {
    QDialog::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

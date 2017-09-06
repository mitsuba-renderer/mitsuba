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

#include "ui_sceneinfodlg.h"
#include "sceneinfodlg.h"

SceneInformationDialog::SceneInformationDialog(QWidget *parent, Scene *scene) :
        QDialog(parent),
    ui(new Ui::SceneInformationDialog) {
    ui->setupUi(this);
#if defined(__OSX__)
    QFont font = ui->textEdit->currentFont();
    font.setPointSize(12);
    ui->textEdit->setCurrentFont(font);
#endif
    ui->textEdit->setText(scene->toString().c_str());
}

SceneInformationDialog::~SceneInformationDialog() {
    delete ui;
}

void SceneInformationDialog::changeEvent(QEvent *e) {
    QDialog::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

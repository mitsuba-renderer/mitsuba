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

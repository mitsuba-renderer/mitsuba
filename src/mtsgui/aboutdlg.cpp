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

#if defined(MTS_KD_CONSERVE_MEMORY)
    configFlags += "MTS_KD_CONSERVE_MEMORY ";
#endif

#if defined(MTS_KD_DEBUG)
    configFlags += "MTS_KD_DEBUG ";
#endif

#if defined(MTS_HAS_COHERENT_RT)
    configFlags += "MTS_HAS_COHERENT_RT ";
#endif

#if defined(MTS_HAS_COLLADA)
    configFlags += "MTS_HAS_COLLADA ";
#endif

#if defined(MTS_HAS_LIBJPEG)
    configFlags += "MTS_HAS_LIBJPEG ";
#endif

#if defined(MTS_HAS_LIBPNG)
    configFlags += "MTS_HAS_LIBPNG ";
#endif

#if defined(MTS_HAS_OPENEXR)
    configFlags += "MTS_HAS_OPENEXR ";
#endif

#if defined(MTS_HAS_BREAKPAD)
    configFlags += "MTS_HAS_BREAKPAD ";
#endif

#if defined(MTS_HAS_FFTW)
    configFlags += "MTS_HAS_FFTW ";
#endif

    configFlags += formatString("SPECTRUM_SAMPLES=%i ",
        SPECTRUM_SAMPLES).c_str();

    ui->label->setText(ui->label->text().replace("MTS_VERSION", MTS_VERSION));
    ui->label->setText(ui->label->text().replace("MTS_YEAR", MTS_YEAR));
    ui->label->setText(ui->label->text().replace("CONFIG_FLAGS", configFlags));

#if defined(__OSX__)
    ui->label->setText(ui->label->text().replace("font-size:10pt", "font-size:14pt"));
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

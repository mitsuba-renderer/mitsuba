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

#include "glwidget.h"
#include "ui_previewsettingsdlg.h"
#include "previewsettingsdlg.h"

#define REINHARD_MIN 0.01
#define REINHARD_MAX 0.8
#define REINHARD_RANGE (REINHARD_MAX-REINHARD_MIN)

class MethodModel : public QStringListModel {
public:
	MethodModel(QObject *parent, bool supportsSinglePass)
		: QStringListModel(parent), m_supportsSinglePass(supportsSinglePass) {
		QStringList tmp;
		tmp << "Disable"
			<< "OpenGL";
		setStringList(tmp);
	}

	Qt::ItemFlags flags(const QModelIndex &index) const {
		if (index.row() == 2 && !m_supportsSinglePass)
			return Qt::NoItemFlags;
//#if !defined(MTS_HAS_COHERENT_RT)
//		if (index.row() == 3)
//			return Qt::NoItemFlags;
//#endif
		return Qt::ItemIsSelectable | Qt::ItemIsEnabled;
	}
private:
	bool m_supportsSinglePass;
};

PreviewSettingsDialog::PreviewSettingsDialog(QWidget *parent, SceneContext *ctx, const RendererCapabilities *cap) :
		QDialog(parent),
	ui(new Ui::PreviewSettingsDialog), m_context(ctx) {
	ui->setupUi(this);
	ui->pathLengthSlider->setValue(ctx->pathLength);

	int index;
	switch (ctx->shadowMapResolution) {
		case 64: index = 0; break;
		case 128: index = 1; break;
		case 256: index = 2; break;
		case 512: index = 3; break;
		case 1024: index = 4; break;
		case 2048: index = 5; break;
		default:
			SLog(EError, "Invalid shadow map resolution!");
			return;
	}

	ui->shadowResolutionCombo->setCurrentIndex(index);

	int clamping = (int) (ctx->clamping * 100);
	ui->clampingSlider->setValue(clamping);
	ui->gammaSpinBox->setValue(ctx->gamma);
	ui->sRGBCheckBox->setCheckState(ctx->srgb ? Qt::Checked : Qt::Unchecked);
	ui->diffuseSourcesBox->setCheckState(ctx->diffuseSources ? Qt::Checked : Qt::Unchecked);
	ui->diffuseReceiversBox->setCheckState(ctx->diffuseReceivers ? Qt::Checked : Qt::Unchecked);
	ui->previewMethodCombo->setModel(new MethodModel(this,
		cap->isSupported(RendererCapabilities::EGeometryShaders)));
	ui->previewMethodCombo->setCurrentIndex(ctx->previewMethod);
	ui->toneMappingMethodCombo->setCurrentIndex(ctx->toneMappingMethod);
	m_ignoreEvent = false;
	ui->exposureSlider->setValue((int) ((ctx->toneMappingMethod == EGamma
		? ctx->exposure : ctx->reinhardBurn)*100));
	ui->keySlider->setValue((int) ((ctx->reinhardKey-REINHARD_MIN)/REINHARD_RANGE * 100));
	ui->diffuseReceiversBox->setEnabled(ui->diffuseSourcesBox->isChecked());

	on_previewMethodCombo_activated(ctx->previewMethod);
	on_toneMappingMethodCombo_activated(ctx->toneMappingMethod);

#if defined(__OSX__)
	layout()->setContentsMargins(10,10,10,5);
#endif
}

void PreviewSettingsDialog::on_resetButton_clicked() {
	ui->pathLengthSlider->setValue(3);
	ui->shadowResolutionCombo->setCurrentIndex(2);
	on_shadowResolutionCombo_activated(2);
	ui->previewMethodCombo->setCurrentIndex(EOpenGL);
	on_previewMethodCombo_activated(EOpenGL);
	ui->toneMappingMethodCombo->setCurrentIndex(EGamma);
	on_toneMappingMethodCombo_activated(EGamma);
	ui->clampingSlider->setValue(10);
	ui->gammaSpinBox->setValue(2.2);
	ui->exposureSpinBox->setValue(0.0);
	ui->keySlider->setValue((int) ((0.18-REINHARD_MIN)/REINHARD_RANGE * 100));

	ui->sRGBCheckBox->setCheckState(Qt::Checked);
	ui->diffuseSourcesBox->setCheckState(Qt::Checked);
	ui->diffuseReceiversBox->setCheckState(Qt::Unchecked);
}

void PreviewSettingsDialog::on_keySlider_valueChanged(int value) {
	emit reinhardKeyChanged((Float) ((value / 100.0f) * REINHARD_RANGE + REINHARD_MIN));

}

void PreviewSettingsDialog::on_pathLengthSlider_valueChanged(int value) {
	ui->pathLengthEdit->setText(QString("%1").arg(value));
	emit pathLengthChanged(value);
}

void PreviewSettingsDialog::on_clampingSlider_valueChanged(int value) {
	emit clampingChanged(value / 100.0f);
}

void PreviewSettingsDialog::on_exposureSlider_valueChanged(int value) {
	if (m_ignoreEvent)
		return;
	m_ignoreEvent = true;
	ui->exposureSpinBox->setValue(value / 100.0f);
	m_ignoreEvent = false;

	if (ui->toneMappingMethodCombo->currentIndex() == EGamma)
		emit exposureChanged(value / 100.0f);
	else
		emit reinhardBurnChanged(value / 100.0f);
}

void PreviewSettingsDialog::on_exposureSpinBox_valueChanged(double value) {
	if (m_ignoreEvent)
		return;
	m_ignoreEvent = true;
	ui->exposureSlider->setValue(value * 100.0f);
	m_ignoreEvent = false;
	if (ui->toneMappingMethodCombo->currentIndex() == EGamma)
		emit exposureChanged((Float) value);
	else
		emit reinhardBurnChanged((Float) value);
}

void PreviewSettingsDialog::on_gammaSpinBox_valueChanged(double value) {
	emit gammaChanged(ui->sRGBCheckBox->checkState() == Qt::Checked, (Float) value);
}

void PreviewSettingsDialog::on_sRGBCheckBox_stateChanged(int state) {
	ui->gammaSpinBox->setEnabled(state == Qt::Unchecked);
	emit gammaChanged(state == Qt::Checked, (Float) ui->gammaSpinBox->value());
}

void PreviewSettingsDialog::on_diffuseReceiversBox_stateChanged(int state) {
	emit diffuseReceiversChanged(state == Qt::Checked);
}

void PreviewSettingsDialog::on_diffuseSourcesBox_stateChanged(int state) {
	emit diffuseSourcesChanged(state == Qt::Checked);
	if (state == Qt::Unchecked && ui->diffuseReceiversBox->isChecked()) {
		emit diffuseReceiversChanged(Qt::Unchecked);
		ui->diffuseReceiversBox->setCheckState(Qt::Unchecked);
	}
	ui->diffuseReceiversBox->setEnabled(state == Qt::Checked);
}

void PreviewSettingsDialog::on_previewMethodCombo_activated(int index) {
	bool visible = true; //index != ERayTrace && index != ERayTraceCoherent;
	ui->shadowResolutionCombo->setVisible(visible);
	ui->shadowResolutionLabel->setVisible(visible);
	emit previewMethodChanged((EPreviewMethod) index);
}

void PreviewSettingsDialog::on_toneMappingMethodCombo_activated(int index) {
	bool keyVisible = index == EReinhard;
	if (index == EReinhard) {
		ui->exposureLabel->setText("Burn :");
		ui->exposureSlider->setValue((int) (m_context->reinhardBurn*100));
	} else {
		ui->exposureLabel->setText("E&xposure : 2 ^");
		ui->exposureSlider->setValue((int) (m_context->exposure*100));
	}

	ui->keySlider->setVisible(keyVisible);
	ui->keyLabel->setVisible(keyVisible);
	emit toneMappingMethodChanged((EToneMappingMethod) index);
}

void PreviewSettingsDialog::on_shadowResolutionCombo_activated(int index) {
	int res;
	switch (index) {
		case 0: res = 64; break;
		case 1: res = 128; break;
		case 2: res = 256; break;
		case 3: res = 512; break;
		case 4: res = 1024; break;
		case 5: res = 2048; break;
		default:
			SLog(EError, "Invalid shadow map resolution!");
			return;
	}
	emit shadowMapResolutionChanged(res);
}

PreviewSettingsDialog::~PreviewSettingsDialog() {
	delete ui;
}

void PreviewSettingsDialog::changeEvent(QEvent *e) {
	QDialog::changeEvent(e);
	switch (e->type()) {
	case QEvent::LanguageChange:
		ui->retranslateUi(this);
		break;
	default:
		break;
	}
}

void PreviewSettingsDialog::setPreviewEnabled(bool enabled) {
	ui->previewLabel->setEnabled(enabled);
	ui->previewMethodLabel->setEnabled(enabled);
	ui->previewMethodCombo->setEnabled(enabled);
	ui->pathLengthLabel->setEnabled(enabled);
	ui->pathLengthSlider->setEnabled(enabled);
	ui->pathLengthEdit->setEnabled(enabled);
	ui->clampingLabel->setEnabled(enabled);
	ui->clampingSlider->setEnabled(enabled);
	ui->shadowResolutionLabel->setEnabled(enabled);
	ui->shadowResolutionCombo->setEnabled(enabled);
	ui->diffuseLabel->setEnabled(enabled);
	ui->diffuseSourcesBox->setEnabled(enabled);
	ui->diffuseReceiversBox->setEnabled(enabled);
}

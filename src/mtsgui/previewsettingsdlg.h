#if !defined(__PREVIEWSETTINGSDLG_H)
#define __PREVIEWSETTINGSDLG_H

#include "common.h"
#include <mitsuba/hw/renderer.h>

namespace Ui {
    class PreviewSettingsDialog;
}

struct SceneContext;
class PreviewSettingsDialog : public QDialog {
    Q_OBJECT
public:
    PreviewSettingsDialog(QWidget *parent,
        SceneContext *ctx, const RendererCapabilities *cap);
    virtual ~PreviewSettingsDialog();

signals:
    void pathLengthChanged(int length);
    void shadowMapResolutionChanged(int resolution);
    void clampingChanged(Float clamping);
    void exposureChanged(Float exposure);
    void gammaChanged(bool srgb, Float gamma);
    void previewMethodChanged(EPreviewMethod method);
    void toneMappingMethodChanged(EToneMappingMethod method);
    void reinhardKeyChanged(Float key);
    void reinhardBurnChanged(Float burn);
    void diffuseReceiversChanged(bool);
    void diffuseSourcesChanged(bool);

public slots:
    void setPreviewEnabled(bool);

protected slots:
    void on_pathLengthSlider_valueChanged(int value);
    void on_clampingSlider_valueChanged(int value);
    void on_shadowResolutionCombo_activated(int index);
    void on_exposureSlider_valueChanged(int value);
    void on_exposureSpinBox_valueChanged(double value);
    void on_gammaSpinBox_valueChanged(double value);
    void on_sRGBCheckBox_stateChanged(int state);
    void on_diffuseReceiversBox_stateChanged(int state);
    void on_diffuseSourcesBox_stateChanged(int state);
    void on_resetButton_clicked();
    void on_previewMethodCombo_activated(int index);
    void on_toneMappingMethodCombo_activated(int index);
    void on_keySlider_valueChanged(int value);

protected:
    void changeEvent(QEvent *e);

private:
    Ui::PreviewSettingsDialog *ui;
    bool m_ignoreEvent;
    SceneContext *m_context;
};

#endif // __PREVIEWSETTINGSDLG_H

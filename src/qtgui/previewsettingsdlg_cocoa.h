#if defined(__OBJC__)
#include <Cocoa/Cocoa.h>
#include <BWToolkitFramework/BWToolkitFramework.h>
#endif
#include "common.h"

class PreviewSettingsDlg;

#if defined(__OBJC__)
@interface CocoaRenderSettingsDlg : NSObject
{
	IBOutlet NSPanel *panel;
	IBOutlet BWTransparentPopUpButton *previewMethod;
	IBOutlet BWTransparentPopUpButton *shadowMapResolution;
	IBOutlet BWTransparentPopUpButton *toneMappingMethod;
	IBOutlet BWTransparentSlider *pathLength;
	IBOutlet BWTransparentSlider *clamping;
	IBOutlet BWTransparentSlider *reinhardKey;
	IBOutlet BWTransparentSlider *exposure;
	IBOutlet BWTransparentSlider *gamma;
	IBOutlet BWTransparentCheckbox *sRGB;
	IBOutlet BWTransparentCheckbox *diffuseSourcesBox;
	IBOutlet BWTransparentCheckbox *diffuseReceiversBox;
	IBOutlet NSTextField *shadowMapLabel;
	IBOutlet NSTextField *reinhardKeyLabel;
	IBOutlet NSTextField *gammaLabel;
	IBOutlet NSTextField *gammaValueLabel;
	IBOutlet NSTextField *exposureLabel;
	IBOutlet NSTextField *exposureValueLabel;
	NSView *parent;
	SceneContext *context;
	PreviewSettingsDlg *delegate;
	bool active;
}
- (id) initWithParentView: (NSView *) view andDelegate: (PreviewSettingsDlg *) del;
- (IBAction) previewMethodChanged: (id) sender;
- (IBAction) clampingChanged: (id) sender;
- (IBAction) pathLengthChanged: (id) sender;
- (IBAction) shadowMapResolutionChanged: (id) sender;
- (IBAction) toneMappingMethodChanged: (id) sender;
- (IBAction) reinhardKeyChanged: (id) sender;
- (IBAction) exposureChanged: (id) sender;
- (IBAction) gammaChanged: (id) sender;
- (IBAction) sRGBChanged: (id) sender;
- (IBAction) diffuseSourcesChanged: (id) sender;
- (IBAction) diffuseReceiversChanged: (id) sender;
- (IBAction) reset: (id) sender;
- (void) updateUI;
- (void) showAt: (NSPoint) point;
@end
#else
class CocoaRenderSettingsDlg;
#endif

class PreviewSettingsDlg : public QObject {
	Q_OBJECT
public:
	PreviewSettingsDlg(QWidget *parent);
	virtual ~PreviewSettingsDlg();

	void show();
	void hide();
	bool isVisible();
	bool isActiveWindow() const;
	void setContext(SceneContext *ctx);
public:
	void triggerPreviewMethodChanged(EPreviewMethod method) {
		emit previewMethodChanged(method);
	}
	void triggerToneMappingMethodChanged(EToneMappingMethod method) {
		emit toneMappingMethodChanged(method);
	}
	void triggerClampingChanged(Float clamping) {
		emit clampingChanged(clamping);
	}
	void triggerPathLengthChanged(int pathLength) {
		emit pathLengthChanged(pathLength);
	}
	void triggerShadowMapResolutionChanged(int res) {
		emit shadowMapResolutionChanged(res);
	}
	void triggerGammaChanged(bool srgb, Float gamma) {
		emit gammaChanged(srgb, gamma);
	}
	void triggerExposureChanged(float exposure) {
		emit exposureChanged(exposure);
	}
	void triggerReinhardKeyChanged(float key) {
		emit reinhardKeyChanged(key);
	}
	void triggerReinhardBurnChanged(float burn) {
		emit reinhardBurnChanged(burn);
	}
	void triggerDiffuseReceiversChanged(bool value) {
		emit diffuseReceiversChanged(value);
	}
	void triggerDiffuseSourcesChanged(bool value) {
		emit diffuseSourcesChanged(value);
	}
	void triggerClose() {
		emit close();
	}
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
	void close();
	void diffuseReceiversChanged(bool);
	void diffuseSourcesChanged(bool);
private:
	CocoaRenderSettingsDlg *m_dlg;
};

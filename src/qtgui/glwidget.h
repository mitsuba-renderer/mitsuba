#if !defined(__GLWIDGET_H)
#define __GLWIDGET_H

#include "common.h"
#include <QGLWidget>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/render/vpl.h>
#include <mitsuba/hw/renderer.h>
#include <mitsuba/hw/gputexture.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/hw/gpusync.h>
#include <mitsuba/hw/vpl.h>
#if defined(WIN32)
#include <mitsuba/hw/wgldevice.h>
#endif

class PreviewThread;

class GLWidget : public QGLWidget {
	Q_OBJECT
public:
	GLWidget(QWidget *parent = 0);
	virtual ~GLWidget();
	QSize sizeHint() const;
	void setScene(SceneContext *context);
	void downloadFramebuffer();
	void resumePreview();
	void refreshScene();
	inline const RendererCapabilities *getRendererCapabilities() const {
		return m_renderer->getCapabilities();
	}

	inline bool getInvertMouse() const { return m_invertMouse; }
	void setInvertMouse(bool invert) { m_invertMouse = invert; }
	inline int getMouseSensitivity() const { return m_mouseSensitivity; }
	void setMouseSensitivity(int sensitivity) { m_mouseSensitivity = sensitivity; }
	void setScrollBars(QScrollBar *hScroll, QScrollBar *vScroll);
	inline void ignoreResizeEvents(bool value) { m_ignoreResizeEvents = value; }
	void updateScrollBars();
	inline const QString &getErrorString() const { return m_errorString; }
	inline bool isUsingSoftwareFallback() const { return m_softwareFallback; }

signals:
	void beginRendering();
	void stopRendering();
	void quit();
	void statusMessage(const QString &status);
	void loadFileRequest(const QString &fileName);

public slots:
	void timerImpulse();
	void onUpdateView();
	void setPathLength(int length);
	void setShadowMapResolution(int shadowMapResolution);
	void setPreviewMethod(EPreviewMethod method);
	void setToneMappingMethod(EToneMappingMethod method);
	void setClamping(Float clamping);
	void setGamma(bool srgb, Float gamma);
	void setReinhardKey(Float value);
	void setReinhardBurn(Float value);
	void setExposure(Float exposure);
	void onException(const QString &what);
	void onScroll();

protected:
	void initializeGL();
	void paintGL();
	void resizeGL(int width, int height);
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	void keyPressEvent(QKeyEvent *event);
	void keyReleaseEvent(QKeyEvent *event);
	void focusOutEvent(QFocusEvent *event);
	void resetPreview();
	void resizeEvent(QResizeEvent *event);
	void wheelEvent(QWheelEvent *event);
	void dragEnterEvent(QDragEnterEvent *event);
	void dropEvent(QDropEvent *event);

	/* Masquerade QGLWidget as a GL device for libhw */
#if defined(WIN32)
	class QtDevice : public WGLDevice {
	public:
		QtDevice(QGLWidget *widget) : WGLDevice(NULL), m_widget(widget) { }
		void init(Device *other = NULL) { 
			m_hwnd = m_widget->winId();
			m_hdc = wglGetCurrentDC();
		}
#else
	class QtDevice : public Device {
	public:
		QtDevice(QGLWidget *widget) : Device(NULL), m_widget(widget) { }
		void init(Device *other = NULL) { }
#endif
		void shutdown() { /* Unsupported */ }
		void setVisible(bool) { /* Unsupported */ }
		void makeCurrent(Renderer *) { /* Unsupported */ }
		virtual ~QtDevice() { }
	private:
		QGLWidget *m_widget;
	};
private:
	ref<Renderer> m_renderer;
	ref<PreviewThread> m_preview;
	ref<GPUTexture> m_logoTexture, m_framebuffer, m_luminanceBuffer[2];
	ref<GPUProgram> m_gammaTonemap, m_reinhardTonemap;
	ref<GPUProgram> m_downsamplingProgram, m_luminanceProgram;
	ref<QtDevice> m_device;
	SceneContext *m_context;
	bool m_framebufferChanged, m_mouseButtonDown;
	bool m_leftKeyDown, m_rightKeyDown;
	bool m_upKeyDown, m_downKeyDown;
	Vector2 m_logoSize;
	Vector m_up;
	QPoint m_mousePos, m_ignoreMouseEvent, m_initialMousePos;
	QTimer *m_movementTimer, *m_redrawTimer;
	QScrollBar *m_hScroll, *m_vScroll;
	ref<Timer> m_clock;
	bool m_invertMouse, m_didSetCursor;
	bool m_ignoreScrollEvents, m_ignoreResizeEvents;
	int m_mouseSensitivity, m_softwareFallback;
	ref<Bitmap> m_fallbackBitmap;
	QString m_errorString;
};

#endif /* __GLWIDGET_H */

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

#if !defined(__GLWIDGET_H)
#define __GLWIDGET_H

#include "common.h"
#if MTS_SSE
#include "simdtonemap.h"
#endif
#include <QtOpenGL/QGLWidget>
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
    enum ECropType {
        ENone = 0,
        ECrop,
        ECropAndMagnify
    };

    GLWidget(QWidget *parent = 0);
    virtual ~GLWidget();
    QSize sizeHint() const;
    void setScene(SceneContext *context);
    void downloadFramebuffer();
    void resumePreview();
    void refreshScene();
    void resetPreview();
    void shutdown();
    inline const RendererCapabilities *getRendererCapabilities() const {
        return m_renderer->getCapabilities();
    }

    inline bool getInvertMouse() const { return m_invertMouse; }
    void setInvertMouse(bool invert) { m_invertMouse = invert; }
    inline ENavigationMode getNavigationMode() const { return m_navigationMode; }
    void setNavigationMode(ENavigationMode mode) { m_navigationMode = mode; }
    inline int getMouseSensitivity() const { return m_mouseSensitivity; }
    void setMouseSensitivity(int sensitivity) { m_mouseSensitivity = sensitivity; }
    void setScrollBars(QScrollBar *hScroll, QScrollBar *vScroll);
    inline void ignoreResizeEvents(bool value) { m_ignoreResizeEvents = value; }
    void updateScrollBars();
    void keyPressEvent(QKeyEvent *event);
    inline const QString &getErrorString() const { return m_errorString; }
    inline bool isUsingSoftwareFallback() const { return m_softwareFallback; }
    inline bool hasSelection() const { return m_aabb.isValid(); }
    inline bool hasFastSoftwareFallback() const {
#if MTS_SSE
        return true;
#else
        return false;
#endif
    }

    void startCrop(ECropType type);

signals:
    void beginRendering();
    void stopRendering();
    void quit();
    void statusMessage(const QString &status);
    void loadFileRequest(const QString &fileName);
    void crop(int type, int x, int y, int width, int height);
    void selectionChanged();
    void switchTab(int rel);

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
    void setDiffuseSources(bool value);
    void setDiffuseReceivers(bool value);
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
    void mouseDoubleClickEvent(QMouseEvent *event);
    void keyReleaseEvent(QKeyEvent *event);
    void focusOutEvent(QFocusEvent *event);
    void resizeEvent(QResizeEvent *event);
    void wheelEvent(QWheelEvent *event);
    void dragEnterEvent(QDragEnterEvent *event);
    void dropEvent(QDropEvent *event);
    void oglRenderKDTree(const KDTreeBase<AABB> *kdtree);
    Point2i upperLeft(bool flipY = false) const;
    void reveal(const AABB &aabb);
    Float autoFocus() const;
    bool askReallyCancelRendering();

    inline ProjectiveCamera *getProjectiveCamera() {
        Sensor *sensor = m_context->scene->getSensor();
        return (sensor->getType() & Sensor::EProjectiveCamera) ?
            static_cast<ProjectiveCamera *>(sensor) : NULL;
    }

    inline const ProjectiveCamera *getProjectiveCamera() const {
        const Sensor *sensor = m_context->scene->getSensor();
        return (sensor->getType() & Sensor::EProjectiveCamera) ?
            static_cast<const ProjectiveCamera *>(sensor) : NULL;
    }

    inline PerspectiveCamera *getPerspectiveCamera() {
        Sensor *sensor = m_context->scene->getSensor();
        return (sensor->getType() & Sensor::EPerspectiveCamera) ?
            static_cast<PerspectiveCamera *>(sensor) : NULL;
    }

    inline const PerspectiveCamera *getPerspectiveCamera() const {
        const Sensor *sensor = m_context->scene->getSensor();
        return (sensor->getType() & Sensor::EPerspectiveCamera) ?
            static_cast<const PerspectiveCamera *>(sensor) : NULL;
    }

    inline Transform getWorldTransform() const {
        const ProjectiveCamera *camera = getProjectiveCamera();
        return camera->getWorldTransform(
            camera->getShutterOpen() + 0.5f * camera->getShutterOpenTime()
        );
    }

    inline bool isRightHanded() {
        return getWorldTransform().det3x3() > 0;
    }

    inline void setWorldTransform(const Transform &trafo) {
        /* Preserve the handedness of the current camera transformation */
        if (getWorldTransform().det3x3() * trafo.det3x3() > 0)
            getProjectiveCamera()->setWorldTransform(trafo);
        else
            getProjectiveCamera()->setWorldTransform(trafo *
                Transform::scale(Vector(-1,1,1)));
    }

#if defined(__WINDOWS__)
    /* Masquerade QGLWidget as a GL device for libhw */
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
    static void setSourceFromResource(GPUProgram *program,
        GPUProgram::EType type, const QString &resourceName);

    GPUProgram* createGPUProgram(const std::string &name,
        const QString &vertexResource, const QString &fragmentResource);

    ref<Renderer> m_renderer;
    ref<PreviewThread> m_preview;
    ref<GPUTexture> m_logoTexture, m_framebuffer, m_luminanceBuffer[2];
    ref<GPUProgram> m_gammaTonemap, m_reinhardTonemap;
    ref<GPUProgram> m_downsamplingProgram, m_luminanceProgram;
#if MTS_SSE
    ref<TonemapCPU> m_cpuTonemap;
#endif
    ref<QtDevice> m_device;
    ref<Font> m_font;
    ref<Bitmap> m_fallbackBitmap;
    SceneContext *m_context;
    int m_mouseSensitivity;
    Vector2 m_logoSize;
    Point m_animationOrigin0, m_animationOrigin1;
    Point m_animationTarget0, m_animationTarget1;
    QPoint m_mousePos, m_ignoreMouseEvent, m_initialMousePos;
    QTimer *m_movementTimer, *m_redrawTimer;
    QScrollBar *m_hScroll, *m_vScroll;
    ref<Timer> m_clock, m_wheelTimer, m_animationTimer;
    ref<Timer> m_statusTimer;
    std::string m_statusMessage;
    bool m_framebufferChanged, m_mouseDrag;
    bool m_leftKeyDown, m_rightKeyDown;
    bool m_upKeyDown, m_downKeyDown, m_animation;
    bool m_invertMouse, m_didSetCursor;
    bool m_ignoreScrollEvents, m_ignoreResizeEvents;
    bool m_softwareFallback, m_cropping;
    Point2i m_cropStart, m_cropEnd;
    ENavigationMode m_navigationMode;
    ECropType m_cropType;
    QString m_errorString;
    AABB m_aabb;
};

#endif /* __GLWIDGET_H */

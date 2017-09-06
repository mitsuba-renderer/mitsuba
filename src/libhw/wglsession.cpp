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

#include <mitsuba/hw/wgldevice.h>

MTS_NAMESPACE_BEGIN

WGLSession::WGLSession()
 : Session() {
}

WGLSession::~WGLSession() {
    if (m_initialized)
        shutdown();
}

void WGLSession::init() {
    Session::init();

    Log(EDebug, "Initializing WGL session");

    m_hinstance = GetModuleHandle(NULL);
    WNDCLASS wndclass;
    ZeroMemory(&wndclass, sizeof(WNDCLASS));
    wndclass.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
    wndclass.lpfnWndProc = (WNDPROC) WGLDevice::WndProc;
    wndclass.cbClsExtra = 0;
    wndclass.cbWndExtra = sizeof(void *);
    wndclass.hInstance = m_hinstance;
    wndclass.hIcon = LoadIcon(NULL, IDI_APPLICATION);
    wndclass.hCursor = LoadCursor(NULL, IDC_ARROW);
    wndclass.hbrBackground = (HBRUSH) GetStockObject(NULL_BRUSH);
    wndclass.lpszMenuName = NULL;

    for (int i=0; i<100; ++i) {
        m_wndClassName = formatString("Mitsuba_WGL_%i", i);
        wndclass.lpszClassName = m_wndClassName.c_str();
        if (!RegisterClass(&wndclass)) {
            Log(EWarn, "Unable to register window class '%s'", m_wndClassName.c_str());
        } else {
            m_initialized = true;
            break;
        }
    }

    if (!m_initialized)
        Log(EError, "Unable to register window class!");
}

void WGLSession::shutdown() {
    Session::shutdown();

    Log(EDebug, "Shutting down WGL session");

    if (!UnregisterClass(m_wndClassName.c_str(), m_hinstance))
        Log(EWarn, "Unable to unregister window class: %s", lastErrorText().c_str());

    m_initialized = false;
}

void WGLSession::processEvents() {
    MSG msg;

    while (PeekMessage(&msg, NULL, 0, 0, PM_NOREMOVE)) {
        if (GetMessage(&msg, NULL, 0, 0) > 0)
            DispatchMessage(&msg);
    }
}

void WGLSession::processEventsBlocking(bool &stop) {
    MSG msg;

    while (true) {
        if (!PeekMessage(&msg, NULL, 0, 0, PM_NOREMOVE) && stop)
            break;
        if (GetMessage(&msg, NULL, 0, 0) > 0)
            DispatchMessage(&msg);
    }
}

MTS_IMPLEMENT_CLASS(WGLSession, false, Session)
MTS_NAMESPACE_END

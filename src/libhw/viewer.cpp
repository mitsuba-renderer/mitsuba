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

#include <mitsuba/hw/viewer.h>

MTS_NAMESPACE_BEGIN

Viewer::Viewer() {
    m_session = Session::create();
    m_device = Device::create(m_session);
    m_renderer = Renderer::create(m_session);
    m_device->setSize(Vector2i(768, 576));
}

int Viewer::run(int argc, char **argv) {
    m_session->init();
    m_device->init();
    m_renderer->init(m_device);
    m_device->addCallback(this);
    m_device->setVisible(true);
    m_font = new Font(Font::EBitstreamVeraMono14);
    m_font->init(m_renderer);
    m_quit = false;
    m_leaveEventLoop = true;
    DeviceEvent event(Device::EResizeEvent);
    windowResized(event);

    if (init(argc, argv)) {
        while (true) {
            m_session->processEventsBlocking(m_leaveEventLoop);
            m_leaveEventLoop = false;
            if (m_quit)
                break;
            m_renderer->clear();
            draw();
            m_device->flip();
        }
        shutdown();
    }

    m_font->cleanup();
    m_renderer->shutdown();
    m_device->shutdown();
    m_session->shutdown();
    return 0;
}

void Viewer::drawHUD(const std::string &text) {
    m_renderer->setColor(Spectrum(0.9f));
    m_renderer->drawText(Point2i(10, 10), m_font, text.c_str());
}

bool Viewer::init(int argc, char **argv) {
    return true;
}

void Viewer::shutdown() { }
void Viewer::keyPressed(const DeviceEvent &event) { }
void Viewer::keyReleased(const DeviceEvent &event) { }
void Viewer::mouseButtonPressed(const DeviceEvent &event) { }
void Viewer::mouseButtonReleased(const DeviceEvent &event) { }
void Viewer::mouseMoved(const DeviceEvent &event) { }
void Viewer::mouseBeginDrag(const DeviceEvent &event) { }
void Viewer::mouseDragged(const DeviceEvent &event) { }
void Viewer::mouseEndDrag(const DeviceEvent &event) { }
void Viewer::windowResized(const DeviceEvent &event) { }

bool Viewer::deviceEventOccurred(const DeviceEvent &event) {
    switch (event.getType()) {
        case Device::EKeyDownEvent:
            if (event.getKeyboardKey() == 'q'
                || event.getKeyboardSpecial() == Device::EKeyEscape) {
                m_quit = true;
                m_leaveEventLoop = true;
            } else {
                keyPressed(event);
            }
            break;
        case Device::EKeyUpEvent: keyReleased(event); break;
        case Device::EMouseMotionEvent: mouseMoved(event); break;
        case Device::EMouseDragEvent: mouseDragged(event); break;
        case Device::EMouseButtonDownEvent: mouseButtonPressed(event); break;
        case Device::EMouseButtonUpEvent: mouseButtonReleased(event); break;
        case Device::EMouseBeginDragEvent: mouseBeginDrag(event); break;
        case Device::EMouseEndDragEvent: mouseEndDrag(event); break;
        case Device::EQuitEvent:
            m_quit = true;
            m_leaveEventLoop = true;
            break;
        case Device::EResizeEvent:
            m_renderer->reconfigure(m_device);
            windowResized(event);
            // no break
        case Device::EGainFocusEvent:
            m_leaveEventLoop = true;
            break;
    }

    return true;
}

MTS_IMPLEMENT_CLASS(Viewer, true, Utility);
MTS_NAMESPACE_END

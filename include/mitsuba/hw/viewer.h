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

#pragma once
#if !defined(__MITSUBA_HW_VIEWER_H_)
#define __MITSUBA_HW_VIEWER_H_

#include <mitsuba/render/util.h>
#include <mitsuba/hw/renderer.h>
#include <mitsuba/hw/font.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Template for simple GLUT-style viewer applications
 *
 * This class makes it possible to rapidly prototype simple OpenGL
 * applications, which can be started using the 'mtsutil' launcher.
 *
 * \ingroup libhw
 */
class MTS_EXPORT_HW Viewer : public Utility, public DeviceEventListener {
public:
	/// Construct a new viewer
	Viewer();

	/**
	 * \brief Program entry point
	 *
	 * The viewer is initialized by this method -- to add custom
	 * initialization code, please override \ref init
	 */
	int run(int argc, char **argv);
protected:
	/// Draw a heads-up display
	void drawHUD(const std::string &text);

	/// Request that the draw() routine be called
	inline void redraw() { m_leaveEventLoop = true; }

	/// To be overwritten by the subclass: main drawing routine
	virtual void draw() = 0;

	/**
	 * \brief To be overwritten (optionally): perform any necessary
	 * initializations
	 *
	 * The default implementation does nothing and returns \c true.
	 *
	 * \param argc Number of command line arguments
	 * \param argv List of command line argument strings
	 * \return \c true upon success. Otherwise, the viewer will quit.
	 */
	virtual bool init(int argc, char **argv);

	/// To be overwritten (optionally): perform any necessary cleanups
	virtual void shutdown();

	/// To be overwritten (optionally): handle a key press event
	virtual void keyPressed(const DeviceEvent &event);

	/// To be overwritten (optionally): handle a key release event
	virtual void keyReleased(const DeviceEvent &event);

	/// To be overwritten (optionally): handle a mouse button press event
	virtual void mouseButtonPressed(const DeviceEvent &event);

	/// To be overwritten (optionally): handle a mouse button release event
	virtual void mouseButtonReleased(const DeviceEvent &event);

	/// To be overwritten (optionally): handle a mouse motion event
	virtual void mouseMoved(const DeviceEvent &event);

	/// To be overwritten (optionally): handle a mouse begin drag event
	virtual void mouseBeginDrag(const DeviceEvent &event);

	/// To be overwritten (optionally): handle a mouse drag event
	virtual void mouseDragged(const DeviceEvent &event);

	/// To be overwritten (optionally): handle a mouse end drag event
	virtual void mouseEndDrag(const DeviceEvent &event);

	/// To be overwritten (optionally): handle a window resize event
	virtual void windowResized(const DeviceEvent &event);

private:
	bool deviceEventOccurred(const DeviceEvent &event);

	MTS_DECLARE_UTILITY()
protected:
	ref<Session> m_session;
	ref<Device> m_device;
	ref<Renderer> m_renderer;
	ref<Font> m_font;
	bool m_quit, m_leaveEventLoop;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_VIEWER_H_ */

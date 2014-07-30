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

#include <mitsuba/hw/session.h>
#if  defined(WIN32)
#include <mitsuba/hw/wglsession.h>
#elif defined(__OSX__)
#include <mitsuba/hw/nsglsession.h>
#else
#include <mitsuba/hw/x11session.h>
#endif

MTS_NAMESPACE_BEGIN

Session::Session() {
	m_initialized = false;
}

Session *Session::create() {
#if defined(WIN32)
	return new WGLSession();
#elif defined(__OSX__)
	return new NSGLSession();
#else
	return new X11Session();
#endif
}

void Session::init() {
	Assert(!m_initialized);
}

void Session::shutdown() {
	Assert(m_initialized);
}

MTS_IMPLEMENT_CLASS(Session, true, Object)
MTS_NAMESPACE_END

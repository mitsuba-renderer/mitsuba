/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#include <mitsuba/core/appender.h>
#include <fstream>
#ifdef WIN32
#include <io.h>
#endif

MTS_NAMESPACE_BEGIN

StreamAppender::StreamAppender(std::ostream *stream)
 : m_stream(stream), m_isFile(false) {
	m_lastMessageWasProgress = false;
}

StreamAppender::StreamAppender(const std::string &filename)
 : m_fileName(filename), m_isFile(true) {
	m_stream = new std::ofstream();
	((std::ofstream *) m_stream)->open(filename.c_str());
	m_lastMessageWasProgress = false;
}

void StreamAppender::append(ELogLevel level, const std::string &text) {
	/* Insert a newline if the last message was a progress message */
	if (m_lastMessageWasProgress && !m_isFile)
		std::cout << std::endl;
	(*m_stream) << text << std::endl;
	m_lastMessageWasProgress = false;
}
	
void StreamAppender::logProgress(Float progress, const std::string &name,
	const std::string &formatted, const std::string &eta, const void *ptr) {
	if (!m_isFile) {
		std::cout << formatted;
		std::cout.flush();
	}
	m_lastMessageWasProgress = true;
}

std::string StreamAppender::toString() const {
	std::ostringstream oss;

	oss << "StreamAppender[stream=";

	if (m_isFile) {
		oss << "\"" << m_fileName << "\"";
	} else {
		oss << "<std::ostream>";
	}

	oss << "]";
	
	return oss.str();
}

StreamAppender::~StreamAppender() {
	if (m_isFile) {
		((std::ofstream *) m_stream)->close();
		delete m_stream;
	}
}

UnbufferedAppender::UnbufferedAppender(int fd) 
 : m_fd(fd) {
	m_lastMessageWasProgress = false;
}

void UnbufferedAppender::append(ELogLevel level, const std::string &text) {
	std::string value = text + std::string("\n");
#if defined(WIN32)
	write(m_fd, value.c_str(), (unsigned int) value.length());
#else
	if (write(m_fd, value.c_str(), value.length()) != (ssize_t) value.length())
		Log(EError, "Unsuccessful write!");
#endif
}

void UnbufferedAppender::logProgress(Float progress, const std::string &name,
	const std::string &formatted, const std::string &eta, const void *ptr) {
	/* Ignore */
}

std::string UnbufferedAppender::toString() const {
	return "UnbufferedAppender[]";
}

UnbufferedAppender::~UnbufferedAppender() {
	close(m_fd);
}

MTS_IMPLEMENT_CLASS(Appender, true, Object)
MTS_IMPLEMENT_CLASS(StreamAppender, false, Appender)
MTS_IMPLEMENT_CLASS(UnbufferedAppender, false, Appender)
MTS_NAMESPACE_END

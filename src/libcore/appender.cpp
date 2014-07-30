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

#include <mitsuba/core/appender.h>
#include <fstream>

#if defined(__WINDOWS__)
# include <io.h>
#endif

MTS_NAMESPACE_BEGIN

StreamAppender::StreamAppender(std::ostream *stream)
 : m_stream(stream), m_isFile(false) {
	m_lastMessageWasProgress = false;
}

StreamAppender::StreamAppender(const std::string &filename)
 : m_fileName(filename), m_isFile(true) {
	std::fstream *stream = new std::fstream();
	stream->open(filename.c_str(),
		std::fstream::in | std::fstream::out | std::fstream::trunc);
	m_stream = stream;
	m_lastMessageWasProgress = false;
}

void StreamAppender::readLog(std::string &target) {
	Assert(m_isFile);
	std::fstream &stream = * ((std::fstream *) m_stream);
	if (!stream.good()) {
		target = "";
		return;
	}
	stream.flush();
	stream.seekg(0, std::ios::end);
	std::streamoff size = stream.tellg();
	if (stream.fail() || size == 0) {
		target = "";
		return;
	}
	target.resize((size_t) size);
	stream.seekg(0, std::ios::beg);

	std::istreambuf_iterator<std::string::value_type> it(stream);
	std::istreambuf_iterator<std::string::value_type> it_eof;
	target.insert(target.begin(), it, it_eof);

	stream.seekg(0, std::ios::end);
	Assert(!stream.fail());
}

void StreamAppender::append(ELogLevel level, const std::string &text) {
	/* Insert a newline if the last message was a progress message */
	if (m_lastMessageWasProgress && !m_isFile)
		(*m_stream) << endl;
	(*m_stream) << text << endl;
	m_lastMessageWasProgress = false;
}

void StreamAppender::logProgress(Float progress, const std::string &name,
	const std::string &formatted, const std::string &eta, const void *ptr) {
	if (!m_isFile) {
		(*m_stream) << formatted;
		m_stream->flush();
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
		((std::fstream *) m_stream)->close();
		delete m_stream;
	}
}

UnbufferedAppender::UnbufferedAppender(int fd)
 : m_fd(fd) {
	m_lastMessageWasProgress = false;
}

void UnbufferedAppender::append(ELogLevel level, const std::string &text) {
	std::string value = text + std::string("\n");
#if defined(__WINDOWS__)
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

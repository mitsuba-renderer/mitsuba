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

#include <mitsuba/core/sshstream.h>
#include <mitsuba/core/statistics.h>

#if !defined(WIN32)
#include <unistd.h>
#endif
#include <errno.h>

MTS_NAMESPACE_BEGIN

SSHStream::SSHStream(const std::string &userName,
	const std::string &hostName, const std::vector<std::string> &cmdLine, int port, int timeout)
 : m_userName(userName), m_hostName(hostName), m_port(port), m_timeout(timeout), m_received(0), m_sent(0)  {
	setByteOrder(ENetworkByteOrder);

	Log(EInfo, "Establishing a SSH connection to \"%s@%s\"", 
		userName.c_str(), hostName.c_str());

#if defined(WIN32)
	/* Inherit pipe handles */
	SECURITY_ATTRIBUTES sAttr;
	sAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
	sAttr.bInheritHandle = TRUE;
	sAttr.lpSecurityDescriptor = NULL;

	/* Create stdout pipe */
	if (!CreatePipe(&m_childOutRd, &m_childOutWr, &sAttr, 0))
		Log(EError, "Error in CreatePipe(): %s", lastErrorText().c_str());

	/* Create stdin pipe */
	if (!CreatePipe(&m_childInRd, &m_childInWr, &sAttr, 0))
		Log(EError, "Error in CreatePipe(): %s", lastErrorText().c_str());

	/* Only inherit one side of the pipes */
	if (!SetHandleInformation(m_childOutRd, HANDLE_FLAG_INHERIT, 0))
		Log(EError, "Error in SetHandleInformation(): %s", lastErrorText().c_str());
	if (!SetHandleInformation(m_childInWr, HANDLE_FLAG_INHERIT, 0))
		Log(EError, "Error in SetHandleInformation(): %s", lastErrorText().c_str());

	/* Start the plink process */
	PROCESS_INFORMATION pi;
	STARTUPINFO si;
	ZeroMemory(&pi, sizeof(PROCESS_INFORMATION));
	ZeroMemory(&si, sizeof(STARTUPINFO));
	si.cb = sizeof(STARTUPINFO);
	si.hStdError  = m_childOutWr;
	si.hStdOutput = m_childOutWr;
	si.hStdInput  = m_childInRd;
	si.dwFlags |= STARTF_USESTDHANDLES;

	std::string params = formatString("-batch -P %i -T %s@%s", m_port, m_userName.c_str(), m_hostName.c_str());

	for (size_t i=0; i<cmdLine.size(); ++i)
		params = params + " " + cmdLine[i];

	if (!CreateProcess("plink.exe", (LPSTR) params.c_str(), 
		NULL, NULL, // Process & thread security attributes
		TRUE, // Inherit handles
		CREATE_NEW_PROCESS_GROUP, // Ignore Ctrl-C
		NULL, // Use environment of the calling process
		NULL, // Use CWD of the calling process
		&si, &pi))
		Log(EError, "Are you sure plink.exe exists? Take a look at sshstream.h -- OS error: CreateProcess() failed: %s", lastErrorText().c_str());

	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);
	CloseHandle(m_childOutWr);
	CloseHandle(m_childInRd);
#else
	int infd[2], outfd[2];

	if (pipe(outfd) == -1)
		Log(EError, "Error in pipe(): %s", strerror(errno));

	if (pipe(infd) == -1)
		Log(EError, "Error in pipe(): %s", strerror(errno));

	char **argv = new char*[15+cmdLine.size()];
	int argc = 0;

	argv[argc++] = strdup("ssh");
	if (m_port != 22) {
		argv[argc++] = strdup("-p"); argv[argc++] = strdup(formatString("%i", m_port).c_str());
	}
	argv[argc++] = strdup("-e"); argv[argc++] = strdup("none"); /* No escape character */
	argv[argc++] = strdup("-T"); /* Disable pseudo TTY allocation (need to be able to transfer binary data) */
	argv[argc++] = strdup("-o"); argv[argc++] = strdup("PasswordAuthentication no"); /* Don't proceed if password auth. is required */
	argv[argc++] = strdup("-o"); argv[argc++] = strdup("StrictHostKeyChecking no"); /* Don't ask whether to accept host keys */
	argv[argc++] = strdup("-o"); argv[argc++] = strdup(formatString("ConnectTimeout %i", m_timeout).c_str()); /* Don't get stuck too long */
	argv[argc++] = strdup("-q"); /* Quiet */
	argv[argc++] = strdup(formatString("%s@%s", m_userName.c_str(), m_hostName.c_str()).c_str());
	for (size_t i=0; i<cmdLine.size(); ++i)
		argv[argc++] = strdup(cmdLine[i].c_str());
	argv[argc++] = NULL;

	if (!fork()) {
		if (close(0) == -1)
			Log(EError, "Error in close(): %s!", strerror(errno));
		if (close(1) == -1)
			Log(EError, "Error in close(): %s!", strerror(errno));
		if (dup2(outfd[0], 0) == -1)
			Log(EError, "Error in dup2(): %s!", strerror(errno));
		if (dup2(infd[1], 1) == -1)
			Log(EError, "Error in dup2(): %s!", strerror(errno));
		if (close(outfd[0]) == -1)
			Log(EError, "Error in close(): %s!", strerror(errno));
		if (close(outfd[1]) == -1)
			Log(EError, "Error in close(): %s!", strerror(errno));
		if (close(infd[0]) == -1)
			Log(EError, "Error in close(): %s!", strerror(errno));
		if (close(infd[1]) == -1)
			Log(EError, "Error in close(): %s!", strerror(errno));
		if (execvp("/usr/bin/ssh", (char * const *) argv) == -1)
			Log(EError, "Error in execvp(): %s!", strerror(errno));
	} else {
		close(outfd[0]);
		close(infd[1]);
		m_infd = infd[0];
		m_outfd = outfd[1];
		m_input = fdopen(infd[0], "rb");
		m_output = fdopen(outfd[1], "wb");
	}
	for (int i=0; i<argc-1; ++i)
		free(argv[i]);
#endif
}

SSHStream::~SSHStream() {
	Log(EDebug, "Closing SSH connection");
#if defined(WIN32)
	CloseHandle(m_childInWr);
	CloseHandle(m_childOutRd);
#else
	fclose(m_input);
	fclose(m_output);
#endif
}

std::string SSHStream::toString() const {
	std::ostringstream oss;
	oss << "SSHStream[userName='"<< m_userName << "', hostName='" 
		<< m_hostName << "', sent=" << (m_sent / 1024) << " KB, "
		"received=" << (m_received/1024) << " KB]" << endl;
	return oss.str();
}

void SSHStream::setPos(size_t pos) {
	Log(EError, "Cannot seek within a socket stream!");
}

size_t SSHStream::getPos() const {
	Log(EError, "Cannot determine the position within a socket stream!");
	return 0;
}

size_t SSHStream::getSize() const {
	Log(EError, "Cannot determine the size of a socket stream!");
	return 0;
}

void SSHStream::truncate(size_t size) {
	Log(EError, "Cannot truncate a socket stream!");
}

void SSHStream::flush() {
#if defined(WIN32)
	// No-op
#else
	if (fflush(m_output) == EOF)
		Log(EError, "Error in fflush(): %s!", strerror(errno));
#endif
}

void SSHStream::read(void *ptr, size_t size) {
	static StatsCounter bytesRcvd("Network", "Bytes received (SSH)");
#if defined(WIN32)
	size_t left = size;
	char *data = (char *) ptr;
	while (left > 0) {
		DWORD nRead = 0;
		if (!ReadFile(m_childOutRd, ptr, (DWORD) left, &nRead, NULL))
			Log(EError, "Connection closed while reading: %s", lastErrorText().c_str());
		left -= nRead;
		data += nRead;
	}
#else
	if (fread(ptr, size, 1, m_input) != 1) {
		if (feof(m_input))
			Log(EError, "Error in fread(): end of file!");
		else if (ferror(m_input))
			Log(EError, "Error in fread(): stream error!");
		/* Otherwise, ignore (strange, but seems to be required) */
	}
#endif
	m_received += size;
	bytesRcvd += size;
}

void SSHStream::write(const void *ptr, size_t size) {
	static StatsCounter bytesSent("Network", "Bytes sent (SSH)");
#if defined(WIN32)
	size_t left = size;
	char *data = (char *) ptr;
	while (left > 0) {
		DWORD nWritten = 0;
		if (!WriteFile(m_childInWr, ptr, (DWORD) left, &nWritten, NULL))
			Log(EError, "Connection closed while writing: %s", lastErrorText().c_str());
		left -= nWritten;
		data += nWritten;
	}
#else
	if (fwrite(ptr, size, 1, m_output) != 1) {
		if (feof(m_output))
			Log(EError, "Error in fwrite(): end of file!");
		else if (ferror(m_output))
			Log(EError, "Error in fwrite(): stream error!");
		/* Otherwise, ignore (strange, but seems to be required) */
	}
#endif
	m_sent += size;
	bytesSent += size;
}

bool SSHStream::canRead() const {
	return true;
}

bool SSHStream::canWrite() const {
	return true;
}

MTS_IMPLEMENT_CLASS(SSHStream, false, Stream)
MTS_NAMESPACE_END

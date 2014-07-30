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

#include <mitsuba/core/sshstream.h>
#include <mitsuba/core/statistics.h>

#if !defined(__WINDOWS__)
# include <unistd.h>
#else
# include <windows.h>
#endif
#include <errno.h>

MTS_NAMESPACE_BEGIN

struct SSHStream::SSHStreamPrivate
{
	const std::string userName, hostName;
	const int port, timeout;
	size_t received, sent;
#if defined(__WINDOWS__)
	HANDLE childInRd,  childInWr;
	HANDLE childOutRd, childOutWr;
#else
	int infd, outfd;
	FILE *input, *output;
#endif

	SSHStreamPrivate(const std::string& uname, const std::string& hname,
		int p, int tm) :
	userName(uname), hostName(hname), port(p), timeout(tm), received(0), sent(0),
#if defined(__WINDOWS__)
	childInRd(0), childInWr(0), childOutRd(0), childOutWr(0)
#else
	infd(-1), outfd(-1), input(0), output(0)
#endif
	{}
};

SSHStream::SSHStream(const std::string &userName,
	const std::string &hostName, const std::vector<std::string> &cmdLine, int port, int timeout)
 : d(new SSHStreamPrivate(userName, hostName, port, timeout)) {
	setByteOrder(ENetworkByteOrder);

	Log(EInfo, "Establishing a SSH connection to \"%s@%s\"",
		userName.c_str(), hostName.c_str());

#if defined(__WINDOWS__)
	/* Inherit pipe handles */
	SECURITY_ATTRIBUTES sAttr;
	sAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
	sAttr.bInheritHandle = TRUE;
	sAttr.lpSecurityDescriptor = NULL;

	/* Create stdout pipe */
	if (!CreatePipe(&(d->childOutRd), &(d->childOutWr), &sAttr, 0))
		Log(EError, "Error in CreatePipe(): %s", lastErrorText().c_str());

	/* Create stdin pipe */
	if (!CreatePipe(&(d->childInRd), &(d->childInWr), &sAttr, 0))
		Log(EError, "Error in CreatePipe(): %s", lastErrorText().c_str());

	/* Only inherit one side of the pipes */
	if (!SetHandleInformation(d->childOutRd, HANDLE_FLAG_INHERIT, 0))
		Log(EError, "Error in SetHandleInformation(): %s", lastErrorText().c_str());
	if (!SetHandleInformation(d->childInWr, HANDLE_FLAG_INHERIT, 0))
		Log(EError, "Error in SetHandleInformation(): %s", lastErrorText().c_str());

	/* Start the plink process */
	PROCESS_INFORMATION pi;
	STARTUPINFO si;
	ZeroMemory(&pi, sizeof(PROCESS_INFORMATION));
	ZeroMemory(&si, sizeof(STARTUPINFO));
	si.cb = sizeof(STARTUPINFO);
	si.hStdError  = d->childOutWr;
	si.hStdOutput = d->childOutWr;
	si.hStdInput  = d->childInRd;
	si.dwFlags |= STARTF_USESTDHANDLES;

	std::string params = formatString("-batch -P %i -T %s@%s",
		d->port, d->userName.c_str(), d->hostName.c_str());

	for (size_t i=0; i<cmdLine.size(); ++i)
		params = params + " " + cmdLine[i];

	if (!CreateProcess("plink.exe", (LPSTR) params.c_str(),
		NULL, NULL, // Process & thread security attributes
		TRUE, // Inherit handles
		CREATE_NEW_PROCESS_GROUP, // Ignore Ctrl-C
		NULL, // Use environment of the calling process
		NULL, // Use CWD of the calling process
		&si, &pi))
		Log(EError, "Are you sure plink.exe exists? Take a look at sshstream.h "
			"-- OS error: CreateProcess() failed: %s", lastErrorText().c_str());

	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);
	CloseHandle(d->childOutWr);
	CloseHandle(d->childInRd);
#else
	int infd[2], outfd[2];

	if (pipe(outfd) == -1)
		Log(EError, "Error in pipe(): %s", strerror(errno));

	if (pipe(infd) == -1)
		Log(EError, "Error in pipe(): %s", strerror(errno));

	char **argv = new char*[15+cmdLine.size()];
	int argc = 0;

	argv[argc++] = strdup("ssh");
	if (d->port != 22) {
		argv[argc++] = strdup("-p"); argv[argc++] = strdup(formatString("%i", d->port).c_str());
	}
	argv[argc++] = strdup("-e"); argv[argc++] = strdup("none"); /* No escape character */
	argv[argc++] = strdup("-T"); /* Disable pseudo TTY allocation (need to be able to transfer binary data) */
	argv[argc++] = strdup("-o"); argv[argc++] = strdup("PasswordAuthentication no"); /* Don't proceed if password auth. is required */
	argv[argc++] = strdup("-o"); argv[argc++] = strdup("StrictHostKeyChecking no"); /* Don't ask whether to accept host keys */
	argv[argc++] = strdup("-o"); argv[argc++] = strdup(formatString("ConnectTimeout %i", d->timeout).c_str()); /* Don't get stuck too long */
	argv[argc++] = strdup("-q"); /* Quiet */
	argv[argc++] = strdup(formatString("%s@%s", d->userName.c_str(), d->hostName.c_str()).c_str());
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
		d->infd = infd[0];
		d->outfd = outfd[1];
		d->input = fdopen(infd[0], "rb");
		d->output = fdopen(outfd[1], "wb");
	}
	for (int i=0; i<argc-1; ++i)
		free(argv[i]);
#endif
}

SSHStream::~SSHStream() {
	Log(EDebug, "Closing SSH connection");
#if defined(__WINDOWS__)
	CloseHandle(d->childInWr);
	CloseHandle(d->childOutRd);
#else
	fclose(d->input);
	fclose(d->output);
#endif
}

const std::string& SSHStream::getHostName() const {
	return d->hostName;
}

const std::string& SSHStream::getUserName() const {
	return d->userName;
}

size_t SSHStream::getReceivedBytes() const {
	return d->received;
}

size_t SSHStream::getSentBytes() const {
	return d->sent;
}

std::string SSHStream::toString() const {
	std::ostringstream oss;
	oss << "SSHStream[userName='"<< d->userName << "', hostName='"
		<< d->hostName << "', sent=" << (d->sent / 1024) << " KB, "
		"received=" << (d->received/1024) << " KB]" << endl;
	return oss.str();
}

void SSHStream::seek(size_t pos) {
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
#if defined(__WINDOWS__)
	// No-op
#else
	if (fflush(d->output) == EOF)
		Log(EError, "Error in fflush(): %s!", strerror(errno));
#endif
}

void SSHStream::read(void *ptr, size_t size) {
	static StatsCounter bytesRcvd("Network", "Bytes received (SSH)");
#if defined(__WINDOWS__)
	size_t left = size;
	char *data = (char *) ptr;
	while (left > 0) {
		DWORD nRead = 0;
		if (!ReadFile(d->childOutRd, ptr, (DWORD) left, &nRead, NULL))
			Log(EError, "Connection closed while reading: %s", lastErrorText().c_str());
		left -= nRead;
		data += nRead;
	}
#else
	if (fread(ptr, size, 1, d->input) != 1) {
		if (feof(d->input))
			Log(EError, "Error in fread(): end of file!");
		else if (ferror(d->input))
			Log(EError, "Error in fread(): stream error!");
		/* Otherwise, ignore (strange, but seems to be required) */
	}
#endif
	d->received += size;
	bytesRcvd += size;
}

void SSHStream::write(const void *ptr, size_t size) {
	static StatsCounter bytesSent("Network", "Bytes sent (SSH)");
#if defined(__WINDOWS__)
	size_t left = size;
	char *data = (char *) ptr;
	while (left > 0) {
		DWORD nWritten = 0;
		if (!WriteFile(d->childInWr, ptr, (DWORD) left, &nWritten, NULL))
			Log(EError, "Connection closed while writing: %s", lastErrorText().c_str());
		left -= nWritten;
		data += nWritten;
	}
#else
	if (fwrite(ptr, size, 1, d->output) != 1) {
		if (feof(d->output))
			Log(EError, "Error in fwrite(): end of file!");
		else if (ferror(d->output))
			Log(EError, "Error in fwrite(): stream error!");
		/* Otherwise, ignore (strange, but seems to be required) */
	}
#endif
	d->sent += size;
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

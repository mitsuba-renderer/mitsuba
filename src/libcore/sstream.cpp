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

#include <mitsuba/core/sstream.h>
#include <mitsuba/core/statistics.h>

#if !defined(__WINDOWS__)
# include <unistd.h>
# include <errno.h>
# include <sys/types.h>
# include <sys/socket.h>
# include <netinet/in.h>
# include <netdb.h>
# include <arpa/inet.h>
# include <sys/wait.h>
# include <signal.h>

# define INVALID_SOCKET -1
# define SOCKET_ERROR   -1
#else
# include <winsock2.h>
# include <ws2tcpip.h>
#endif

MTS_NAMESPACE_BEGIN

namespace
{
#if defined(__WINDOWS__)
// This function is natively avaiable since Windows Vista
const char *inet_ntop(int af, const void *src, char *dst, socklen_t len) {
	if (af == AF_INET) {
		struct sockaddr_in in;
		memset(&in, 0, sizeof(in));
		in.sin_family = AF_INET;
		memcpy(&in.sin_addr, src, sizeof(struct in_addr));
		if (getnameinfo((struct sockaddr *)&in,
			sizeof(struct sockaddr_in), dst, len, NULL, 0, NI_NUMERICHOST) != 0)
			return NULL;
		return dst;
	} else if (af == AF_INET6) {
		struct sockaddr_in6 in;
		memset(&in, 0, sizeof(in));
		in.sin6_family = AF_INET6;
		memcpy(&in.sin6_addr, src, sizeof(struct in_addr6));
		if (getnameinfo((struct sockaddr *)&in,
			sizeof(struct sockaddr_in6), dst, len, NULL, 0, NI_NUMERICHOST) != 0)
			return NULL;
		return dst;
	}
	return NULL;
}
#endif

void *get_in_addr(struct sockaddr_storage *sa)
{
	if (sa->ss_family == AF_INET)
		return &(((struct sockaddr_in*)sa)->sin_addr);
	return &(((struct sockaddr_in6*)sa)->sin6_addr);
}

} // namespace

SocketStream::SocketStream(socket_t socket)
 : m_socket(socket), m_received(0), m_sent(0) {
	setByteOrder(ENetworkByteOrder);
	struct sockaddr_storage sockaddr;
	socklen_t addrlen = sizeof(sockaddr);
	char s[INET6_ADDRSTRLEN];

	if (getpeername(m_socket, (struct sockaddr *) &sockaddr, &addrlen) == SOCKET_ERROR)
		handleError("getpeername");

	if (inet_ntop(sockaddr.ss_family, get_in_addr(&sockaddr), s, sizeof(s)) == NULL)
		handleError("inet_ntop");

	m_peer = s;
}

SocketStream::SocketStream(const std::string &host, int port)
 : m_socket(0), m_received(0), m_sent(0) {
	struct addrinfo hints, *servinfo = NULL;
	char portName[8];
	int rv;

	memset(&hints, 0, sizeof(addrinfo));
	hints.ai_family = AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;
	snprintf(portName, sizeof(portName), "%i", port);

	Log(EInfo, "Connecting to \"%s:%i\"", host.c_str(), port);
	if ((rv = getaddrinfo(host.c_str(), portName, &hints, &servinfo)) != 0)
		Log(EError, "Error in getaddrinfo(): %s", gai_strerror(rv));

	if (servinfo == NULL)
		Log(EError, "Error in getaddrinfo(): did not return results!");

	m_socket = socket(servinfo->ai_family, servinfo->ai_socktype, servinfo->ai_protocol);
	if (m_socket == INVALID_SOCKET) {
		freeaddrinfo(servinfo);
		handleError("socket");
	}

	if (connect(m_socket, servinfo->ai_addr, (socklen_t) servinfo->ai_addrlen) == SOCKET_ERROR) {
		freeaddrinfo(servinfo);
		handleError("connect");
	}

	freeaddrinfo(servinfo);

	setByteOrder(ENetworkByteOrder);
	struct sockaddr_storage sockaddr;
	socklen_t addrlen = sizeof(sockaddr);
	char s[INET6_ADDRSTRLEN];

	if (getpeername(m_socket, (struct sockaddr *) &sockaddr, &addrlen) == SOCKET_ERROR)
		handleError("getpeername");

	if (inet_ntop(sockaddr.ss_family, get_in_addr(&sockaddr), s, sizeof(s)) == NULL)
		handleError("inet_ntop");

#if defined(__OSX__)
	int on = 1;
	/* Turn of SIGPIPE when writing to a disconnected socket */
	if (setsockopt(m_socket, SOL_SOCKET, SO_NOSIGPIPE, &on, sizeof(int)) != 0)
		handleError("setsockopt");
#endif

	m_peer = s;
}

SocketStream::~SocketStream() {
#ifdef __WINDOWS__
	if (closesocket(m_socket) == SOCKET_ERROR)
		handleError("closesocket");
#else
	if (close(m_socket))
		handleError("close");
#endif
}

void SocketStream::read(void *ptr, size_t size) {
	static StatsCounter bytesRcvd("Network", "Bytes received");
	const size_t total = size;
	char *data = (char *) ptr;
	while (size > 0) {
#if defined(__WINDOWS__)
		ssize_t n = recv(m_socket, data, (int) size, 0);
#else
		ssize_t n = recv(m_socket, data, size, 0);
#endif
		if (n == 0) {
			throw EOFException("Connection closed while reading!",
					(size_t) (data - (char *) ptr));
		} else if (n == SOCKET_ERROR) {
			if (!handleError("recv", EWarn))
				continue; /* Wasn't an error after all -- continue */

			throw EOFException("Connection closed while reading!",
						(size_t) (data - (char *) ptr));
		}
		size -= n;
		data += n;
	}
	m_received += total;
	bytesRcvd += total;
}

void SocketStream::write(const void *ptr, size_t size) {
	static StatsCounter bytesSent("Network", "Bytes sent");
	const size_t total = size;
	char *data = (char *) ptr;
	while (size > 0) {
#if defined(__LINUX__)
		/* Linux: Don't send the EPIPE signal when the connection breaks */
		ssize_t n = send(m_socket, data, size, MSG_NOSIGNAL);
#elif defined(__WINDOWS__)
		ssize_t n = send(m_socket, data, (int) size, 0);
#else
		ssize_t n = send(m_socket, data, size, 0);
#endif
		if (n == SOCKET_ERROR) {
			if (!handleError("send", EWarn))
				continue; /* Wasn't an error after all -- continue */

			throw EOFException("Connection closed while writing!",
					(size_t) (data - (char *) ptr));
		}
		size -= n;
		data += n;
	}
	m_sent += total;
	bytesSent += total;
}

bool SocketStream::canRead() const {
	return true;
}

bool SocketStream::canWrite() const {
	return true;
}

std::string SocketStream::toString() const {
	std::ostringstream oss;
	oss << "SocketStream[peer='" << m_peer << "', sent=" << (m_sent / 1024) << " KB, "
		"received=" << (m_received/1024) << " KB]" << endl;
	return oss.str();
}

void SocketStream::seek(size_t pos) {
	Log(EError, "Cannot seek within a socket stream!");
}

size_t SocketStream::getPos() const {
	Log(EError, "Cannot determine the position within a socket stream!");
	return 0;
}

size_t SocketStream::getSize() const {
	Log(EError, "Cannot determine the size of a socket stream!");
	return 0;
}

void SocketStream::truncate(size_t size) {
	Log(EError, "Cannot truncate a socket stream!");
}

void SocketStream::flush() {
	/* Ignore */
}

bool SocketStream::handleError(const std::string &peer, const std::string &cmd, ELogLevel level) {
#if !defined(__WINDOWS__)
	if (level == EWarn && errno == EINTR) /* This is not really a warning -- just retry the operation. */
		return false;

	if (cmd.find("(") == std::string::npos)
		Log(level, "[%s] Error in %s(): %s!", peer.c_str(), cmd.c_str(), strerror(errno));
	else
		Log(level, "[%s] Error in %s: %s!", peer.c_str(), cmd.c_str(), strerror(errno));
#else
	std::string err;
	int error = WSAGetLastError();
	switch (error) {
		case WSABASEERR: err = "Internal error (no reason given)"; break;
		case WSAEINTR:
			if (level == EWarn) /* This is not really a warning -- just retry the operation. */
				return false;
			else
			   err = "Interrupted system call";
			break;
		case WSAEBADF: err = "Bad file number"; break;
		case WSAEACCES: err = "Permission denied"; break;
		case WSAEFAULT: err = "Bad address"; break;
		case WSAEINVAL: err = "Invalid argument"; break;
		case WSAEMFILE: err = "Too many open files"; break;
		case WSAEWOULDBLOCK: err = "Operation would block"; break;
		case WSAEINPROGRESS: err = "Operation now in progress"; break;
		case WSAEALREADY: err = "Operation already in progress"; break;
		case WSAENOTSOCK: err = "Socket operation on non-socket"; break;
		case WSAEDESTADDRREQ: err = "Destination address required"; break;
		case WSAEMSGSIZE: err = "Message too long"; break;
		case WSAEPROTOTYPE: err = "Protocol wrong type for socket"; break;
		case WSAENOPROTOOPT: err = "Bad protocol option"; break;
		case WSAEPROTONOSUPPORT: err = "Protocol not supported"; break;
		case WSAESOCKTNOSUPPORT: err = "Socket type not supported"; break;
		case WSAEOPNOTSUPP: err = "Operation not supported on socket"; break;
		case WSAEPFNOSUPPORT: err = "Protocol family not supported"; break;
		case WSAEAFNOSUPPORT: err = "Address family not supported by protocol family"; break;
		case WSAEADDRINUSE: err = "Address already in use"; break;
		case WSAEADDRNOTAVAIL: err = "Can't assign requested address"; break;
		case WSAENETDOWN: err = "Network is down"; break;
		case WSAENETUNREACH: err = "Network is unreachable"; break;
		case WSAENETRESET: err = "Net dropped connection or reset"; break;
		case WSAECONNABORTED: err = "Software caused connection abort"; break;
		case WSAECONNRESET: err = "Connection reset by peer"; break;
		case WSAENOBUFS: err = "No buffer space available"; break;
		case WSAEISCONN: err = "Socket is already connected"; break;
		case WSAENOTCONN: err = "Socket is not connected"; break;
		case WSAESHUTDOWN: err = "Can't send after socket shutdown"; break;
		case WSAETOOMANYREFS: err = "Too many references can't splice"; break;
		case WSAETIMEDOUT: err = "Connection timed out"; break;
		case WSAECONNREFUSED: err = "Connection refused"; break;
		case WSAELOOP: err = "Too many levels of symbolic links"; break;
		case WSAENAMETOOLONG: err = "File name too long"; break;
		case WSAEHOSTDOWN: err = "Host is down"; break;
		case WSAEHOSTUNREACH: err = "No Route to Host"; break;
		case WSAENOTEMPTY: err = "Directory not empty"; break;
		case WSAEPROCLIM: err = "Too many processes"; break;
		case WSAEUSERS: err = "Too many users"; break;
		case WSAEDQUOT: err = "Disc Quota Exceeded"; break;
		case WSAESTALE: err = "Stale NFS file handle"; break;
		case WSAEREMOTE: err = "Too many levels of remote in path"; break;
		case WSASYSNOTREADY: err = "Network SubSystem is unavailable"; break;
		case WSAVERNOTSUPPORTED: err = "WINSOCK DLL Version out of range"; break;
		case WSANOTINITIALISED: err = "Successful WSASTARTUP not yet performed"; break;
		case WSAHOST_NOT_FOUND: err = "Host not found"; break;
		case WSATRY_AGAIN: err = "Non-Authoritative Host not found"; break;
		case WSANO_DATA: err = "Valid name no data record of requested"; break;
		default: err = "Unknown error";
	}
	if (cmd.find("(") == std::string::npos)
		Log(level, "Error %i in %s(): %s!", error, cmd.c_str(), err.c_str());
	else
		Log(level, "Error %i in %s: %s!", error, cmd.c_str(), err.c_str());
#endif
	return true;
}

MTS_IMPLEMENT_CLASS(SocketStream, false, Stream)
MTS_NAMESPACE_END

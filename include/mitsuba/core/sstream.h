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

#if !defined(__SSTREAM_H)
#define __SSTREAM_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/stream.h>

#ifdef WIN32
#include <io.h>
#include <ws2tcpip.h>
#endif

MTS_NAMESPACE_BEGIN

/** \brief Portable %Stream implementation, which encapsulates a socket 
 * for IPv4/IPv6 network communications. 
 *
 * By default, this type of stream is configured to use network byte
 * order (= big endian).
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE SocketStream : public Stream {
public:
	// =============================================================
	//! @{ \name Constructors
	// =============================================================

	/**
	 * \brief Create a stream from an existing socket
	 * \remark This function is not exposed in the Python bindings
	 */
#if defined(WIN32)
	SocketStream(SOCKET socket);
#else
	SocketStream(int socket);
#endif

	/// Connect to the given host/port
	SocketStream(const std::string &host, int port);

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Socket stream-specific features
	// =============================================================

	/// Return the peer's name
	inline const std::string &getPeer() const { return m_peer; }
	
	/// Return the number of received bytes
	inline size_t getReceivedBytes() const { return m_received; }
	
	/// Return the number of sent bytes
	inline size_t getSentBytes() const { return m_sent; }

	/// Return a string representation
	std::string toString() const;

	/// Handle the last socket-specific error (looks up the appropriate OS description)
	static void handleError(const std::string &cmd, ELogLevel level = EError);

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Implementation of the Stream interface
	// =============================================================

	void read(void *ptr, size_t size);
	void write(const void *ptr, size_t size);
	void setPos(size_t pos);
	size_t getPos() const;
	size_t getSize() const;
	void truncate(size_t size);
	void flush();
	bool canWrite() const;
	bool canRead() const;

	//! @}
	// =============================================================

	MTS_DECLARE_CLASS()
protected:
	/** \brief Virtual destructor
	 *
	 * The destructor frees all resources and closes
	 * the socket if it is still open
	 */
	virtual ~SocketStream();
protected:
#if defined(WIN32)
	SOCKET m_socket;
#else
	int m_socket;
#endif
	size_t m_received, m_sent;
	std::string m_peer;
};

#ifdef WIN32
extern MTS_EXPORT_CORE const char *inet_ntop(int af, const void *src, char *dst, socklen_t len);
#endif

MTS_NAMESPACE_END

#endif /* __SSTREAM_H */

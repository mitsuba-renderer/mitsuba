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
#if !defined(__MITSUBA_CORE_SSTREAM_H_)
#define __MITSUBA_CORE_SSTREAM_H_

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/stream.h>

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
	/// Socket typedef. For Windows it is based on the code in WinSock2.h
#if defined(_WIN64)
	typedef uint64_t socket_t;
#elif defined(_WIN32)
	typedef uint32_t socket_t;
#else
	typedef int socket_t;
#endif

	// =============================================================
	//! @{ \name Constructors
	// =============================================================

	/**
	 * \brief Create a stream from an existing socket
	 * \remark This function is not exposed in the Python bindings
	 */
	SocketStream(socket_t socket);

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
	inline bool handleError(const std::string &cmd, ELogLevel level = EError) { return SocketStream::handleError(m_peer, cmd, level); }

	/// Handle the last socket-specific error (looks up the appropriate OS description)
	static bool handleError(const std::string &peer, const std::string &cmd, ELogLevel level = EError);

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Implementation of the Stream interface
	// =============================================================

	void read(void *ptr, size_t size);
	void write(const void *ptr, size_t size);
	void seek(size_t pos);
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
	socket_t m_socket;
	size_t m_received, m_sent;
	std::string m_peer;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_SSTREAM_H_ */

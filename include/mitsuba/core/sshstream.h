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

#if !defined(__SSHSTREAM_H)
#define __SSHSTREAM_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN


/** \brief Stream implementation based on an encrypted SSH tunnel
 * 
 * This class remotely starts a program and exposes its stdin/stdout 
 * streams through an instance of \ref Stream. To make all 
 * of this work, passwordless authentication must be enabled (for 
 * example by using public key authentication in addition to a 
 * running ssh-agent, which stores the decrypted private key).
 *
 * On Windows, things are implemented a bit differently: Instead
 * of OpenSSH, plink.exe (from PUTTY) is used and must be available
 * in $PATH. For passwordless authentication, convert your private
 * key to PuTTY's format (with the help of puttygen.exe). Afterwards, 
 * pageant.exe is required to load and authenticate the key.
 *
 * Note: SSH streams are set to use network byte order by default.
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE SSHStream : public Stream {
public:
	// =============================================================
	//! @{ \name Constructors
	// =============================================================

	/**
	 * \brief Create a new SSH stream. 
	 *
	 * The timeout parameter specifies specifies the maximum amount of 
	 * time that can be spent before failing to create the initial 
	 * connection. This feature is unsupported (and ignored) on Windows.
	 * 
	 * \param userName Username to use for the authentication
	 * \param hostName Destination host name
	 * \param cmdLine  Command (with arguments) to be executed on the remote side
	 * \param port     Destination port
	 * \param timeout  Maximum time to use for the connection attempt (in seconds)
	 */
	SSHStream(const std::string &userName, 
		const std::string &hostName,
		const std::vector<std::string> &cmdLine,
		int port = 22, int timeout = 10
	);

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name SSH stream-specific features
	// =============================================================

	/// Return the destination machine's host name
	inline const std::string &getHostName() const { return m_hostName; }
	
	/// Return the user name used for authentication
	inline const std::string &getUserName() const { return m_userName; }
	
	/// Return the number of received bytes
	inline size_t getReceivedBytes() const { return m_received; }
	
	/// Return the number of sent bytes
	inline size_t getSentBytes() const { return m_sent; }

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


	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/** \brief Virtual destructor
	 *
	 * The destructor frees all resources and closes
	 * the socket if it is still open
	 */
	virtual ~SSHStream();
private:
	std::string m_userName, m_hostName;
	int m_port, m_timeout;
	size_t m_received, m_sent;
#if defined(WIN32)
	HANDLE m_childInRd, m_childInWr;
	HANDLE m_childOutRd, m_childOutWr;
#else
	int m_infd, m_outfd;
	FILE *m_input, *m_output;
#endif
};

MTS_NAMESPACE_END

#endif /* __SSHSTREAM_H */

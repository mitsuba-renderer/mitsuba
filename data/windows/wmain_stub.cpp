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

#include <mitsuba/core/platform.h>
#if defined(__WINDOWS__)

// Stub for generating UTF-8 command line arguments from wmain (UTF-16)
#include <Windows.h>

extern int mts_main(int argc, char **argv);


namespace {

class ArgsUTF8 {
public:
	ArgsUTF8(int argc, wchar_t *wargv[]) :
		m_argc(-1), m_argv(NULL), m_data(NULL)
	{
		if (argc > 0)
			m_argc = argc;
		else
			return;

		m_argv = new char*[argc];
		int total = 0;

		// Pass 1: get the lengths of each converted string an allocate data
		for (int i = 0; i < argc; ++i) {
			const int lenUtf8 = WideCharToMultiByte(CP_UTF8, 0,
				wargv[i], -1, NULL, 0, NULL, NULL);
			if (lenUtf8 != 0) {
				total += lenUtf8;
				m_argv[i] = reinterpret_cast<char*>(lenUtf8);
			} else {
				m_argc = i;
				break;
			}
		}

		if (m_argc < 1)
			return;

		m_data = new char[total];
		int currOffset = 0;

		// Pass 2: perform the conversion
		for (int i = 0; i < m_argc; ++i) {
			int lenUtf8 = reinterpret_cast<int>(m_argv[i]);
			m_argv[i] = m_data + currOffset;
			lenUtf8 = WideCharToMultiByte(CP_UTF8, 0,
				wargv[i], -1, m_argv[i], lenUtf8, NULL, NULL);
			if (lenUtf8 != 0) {
				currOffset += lenUtf8;
			} else {
				m_argc = i;
				return;
			}
		}
	}

	~ArgsUTF8() {
		if (m_argv != NULL) {
			delete [] m_argv;
		}
		if (m_data != NULL) {
			delete [] m_data;
		}
	}

	inline int argc() const {
		return m_argc;
	}

	inline char** argv() const {
		return m_argv;
	}

private:
	int m_argc;
	char** m_argv;
	char*  m_data;  
};

} // namespace


// MSDN Documentation:
//   http://msdn.microsoft.com/en-US/library/fzc2cy7w%28v=vs.110%29.aspx
int wmain(int argc, wchar_t *wargv[], wchar_t *envp[]) {
	ArgsUTF8 argsUTF8(argc, wargv);
	return mts_main(argsUTF8.argc(), argsUTF8.argv());
}

#endif // __WINDOWS__

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

#include <mitsuba/core/random.h>
#include <mitsuba/core/quad.h>
#include <boost/bind.hpp>
#include <stdarg.h>
#include <iomanip>
#include <errno.h>

/* Some of the implementations in this file are based on PBRT */

#if defined(__OSX__)
#include <sys/sysctl.h>
#elif defined(WIN32)
#include <direct.h>
#else
#include <malloc.h>
#endif

#if defined(WIN32)
#include <winsock2.h>
#include <ws2tcpip.h>
#else
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <fenv.h>
#endif

#if !defined(L1_CACHE_LINE_SIZE)
#define L1_CACHE_LINE_SIZE 64
#endif

MTS_NAMESPACE_BEGIN

#ifdef MTS_SSE
static const float pinf = std::numeric_limits<float>::infinity();
static const float flt_max = std::numeric_limits<float>::max();
const MM_ALIGN16 SSEVector SSEConstants::zero	  = SSEVector(-1.0f, 0.0f, 0.0f, 0.0f);
const MM_ALIGN16 SSEVector SSEConstants::one	  = SSEVector(1.0f, 1.0f, 1.0f, 1.0f);
const MM_ALIGN16 SSEVector SSEConstants::max	  = SSEVector(flt_max, flt_max, flt_max, flt_max);
const MM_ALIGN16 SSEVector SSEConstants::eps	  = SSEVector(Epsilon, Epsilon, Epsilon, Epsilon);
const MM_ALIGN16 SSEVector SSEConstants::op_eps   = SSEVector(1+Epsilon, 1+Epsilon, 1+Epsilon, 1+Epsilon);
const MM_ALIGN16 SSEVector SSEConstants::om_eps   = SSEVector(1-Epsilon, 1-Epsilon, 1-Epsilon, 1-Epsilon);
const MM_ALIGN16 SSEVector SSEConstants::p_inf	  = SSEVector(pinf, pinf, pinf, pinf);
const MM_ALIGN16 SSEVector SSEConstants::n_inf	  = SSEVector(-pinf, -pinf, -pinf, -pinf);
const MM_ALIGN16 SSEVector SSEConstants::ffffffff = SSEVector((int32_t) 0xFFFFFFFF, (int32_t) 0xFFFFFFFF, (int32_t) 0xFFFFFFFF, (int32_t) 0xFFFFFFFF);
#endif

const int primeTable[primeTableSize] = {
	2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 
	109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 
	229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 
	353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 
	479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 
	617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 
	757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 
	907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 
	1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 
	1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 
	1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 
	1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 
	1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 
	1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 
	1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 
	1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 
	2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 
	2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 
	2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 
	2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 
	2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 
	2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 
	2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 
	2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 
	3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299, 
	3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 
	3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 
	3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677, 
	3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 
	3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 
	3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 
	4111, 4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243, 
	4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 
	4409, 4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 
	4549, 4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 
	4691, 4703, 4721, 4723, 4729, 4733, 4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 
	4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 
	4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 
	5119, 5147, 5153, 5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, 5281, 
	5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 
	5441, 5443, 5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569, 
	5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701, 5711, 
	5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 
	5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981, 5987, 6007, 6011, 6029, 
	6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133, 6143, 6151, 6163, 
	6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, 
	6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373, 6379, 6389, 6397, 6421, 6427, 6449, 
	6451, 6469, 6473, 6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581, 6599, 6607, 
	6619, 6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 
	6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 
	6911, 6917, 6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027, 7039, 
	7043, 7057, 7069, 7079, 7103, 7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 7211, 7213, 
	7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 
	7411, 7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 
	7547, 7549, 7559, 7561, 7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 
	7681, 7687, 7691, 7699, 7703, 7717, 7723, 7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 
	7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919
};


// -----------------------------------------------------------------------
//  General utility functions
// -----------------------------------------------------------------------

std::vector<std::string> tokenize(const std::string &string, const std::string &delim) {
	std::string::size_type lastPos = string.find_first_not_of(delim, 0);
	std::string::size_type pos = string.find_first_of(delim, lastPos);
	std::vector<std::string> tokens;

	while (std::string::npos != pos || std::string::npos != lastPos) {
		tokens.push_back(string.substr(lastPos, pos - lastPos));
		lastPos = string.find_first_not_of(delim, pos);
		pos = string.find_first_of(delim, lastPos);
	}

	return tokens;
}

std::string trim(const std::string& str) {
	std::string::size_type 
		start = str.find_first_not_of(" \t\r\n"),
		end = str.find_last_not_of(" \t\r\n");

	return str.substr(start == std::string::npos ? 0 : start, 
			end == std::string::npos ? str.length() - 1 : end - start + 1);
}

std::string indent(const std::string &string, int amount) {
	/* This could probably be done faster (is not
	   really speed-critical though) */
	std::istringstream iss(string);
	std::ostringstream oss;
	std::string str;
	bool firstLine = true;
	while (!iss.eof()) {
		std::getline(iss, str);
		if (!firstLine) {
			for (int i=0; i<amount; ++i)
				oss << "  ";
		}
		oss << str;
		if (!iss.eof())
			oss << endl;
		firstLine = false;
	}
	return oss.str();
}

std::string memString(size_t size) {
	Float value = (Float) size;
	const char *prefixes[] = {
		"B", "KiB", "MiB", "GiB", "TiB", "PiB"
	};
	int prefix = 0;
	while (prefix < 5 && value > 1024.0f) {
		value /= 1024.0f; ++prefix;
	}
	return formatString(prefix == 0 ?
			"%.0f %s" : "%.2f %s", value, prefixes[prefix]);
}

void * __restrict allocAligned(size_t size) {
#if defined(WIN32)
	return _aligned_malloc(size, L1_CACHE_LINE_SIZE);
#elif defined(__OSX__)
	/* OSX malloc already returns 16-byte aligned data suitable
	   for AltiVec and SSE computations */
	return malloc(size);
#else
	return memalign(L1_CACHE_LINE_SIZE, size);
#endif
}

void freeAligned(void *ptr) {
#if defined(WIN32)
	_aligned_free(ptr);
#else
	free(ptr);
#endif
}

int getProcessorCount() {
#if defined(WIN32)
	SYSTEM_INFO sys_info;
	GetSystemInfo(&sys_info);
	return sys_info.dwNumberOfProcessors;
#elif defined(__OSX__)
	int nprocs;
	size_t nprocsSize = sizeof(int);
	if (sysctlbyname("hw.activecpu", &nprocs, &nprocsSize, NULL, 0))
		SLog(EError, "Could not detect the number of processors!");
	return (int) nprocs;
#else
	return sysconf(_SC_NPROCESSORS_CONF);
#endif
}

#if defined(WIN32)
std::string lastErrorText() {
	DWORD errCode = GetLastError();
	char *errorText = NULL;
	if (!FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER 
		| FORMAT_MESSAGE_FROM_SYSTEM
		| FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL,
		errCode,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		(LPTSTR) &errorText,
		0,
		NULL)) {
		return "Internal error while looking up an error code";
	}
	std::string result(errorText);
	LocalFree(errorText);
	return result;
}
#endif

bool enableFPExceptions() {
	bool exceptionsWereEnabled = false;
#if defined(WIN32)
	_clearfp();
	uint32_t cw = _controlfp(0, 0);
	exceptionsWereEnabled = ~cw & (_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
	cw &= ~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
	_controlfp(cw, _MCW_EM);
#elif defined(__OSX__)
#if !defined(MTS_SSE)
#warning SSE must be enabled to handle FP exceptions on OSX
#else
	exceptionsWereEnabled = query_fpexcept_sse() != 0;
#endif
#else
	exceptionsWereEnabled = 
		fegetexcept() & (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
	feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
#if defined(MTS_SSE)
	enable_fpexcept_sse();
#endif
	return exceptionsWereEnabled;
}

bool disableFPExceptions() {
	bool exceptionsWereEnabled = false;
#if defined(WIN32)
	_clearfp();
	uint32_t cw = _controlfp(0, 0);
	exceptionsWereEnabled = ~cw & (_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
	cw |= _EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW;
	_controlfp(cw, _MCW_EM);
#elif defined(__OSX__)
#if !defined(MTS_SSE)
#warning SSE must be enabled to handle FP exceptions on OSX
#else
	exceptionsWereEnabled = query_fpexcept_sse() != 0;
#endif
#else
	exceptionsWereEnabled = 
		fegetexcept() & (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
	fedisableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
#if defined(MTS_SSE)
	disable_fpexcept_sse();
#endif
	return exceptionsWereEnabled;
}

void restoreFPExceptions(bool oldState) {
	bool currentState;
#if defined(WIN32)
	uint32_t cw = _controlfp(0, 0);
	currentState = ~cw & (_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
#elif defined(__OSX__)
#if !defined(MTS_SSE)
#warning SSE must be enabled to handle FP exceptions on OSX
#else
	currentState = query_fpexcept_sse() != 0;
#endif
#else
	currentState = fegetexcept() & (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
	if (oldState != currentState) {
		if (oldState)
			enableFPExceptions();
		else
			disableFPExceptions();
	}
}


std::string getHostName() {
	char hostName[128];
	if (gethostname(hostName, sizeof(hostName)) != 0)
#if defined(WIN32)
		SLog(EError, "Could not retrieve the computer's host name: %s!",
			lastErrorText().c_str());
#else
		SLog(EError, "Could not retrieve the computer's host name : %s!",
			strerror(errno));
#endif
	return hostName;
}

std::string getFQDN() {
	struct addrinfo *addrInfo = NULL, hints;
	memset(&hints, 0, sizeof(addrinfo));
	// Only look for IPv4 addresses -- perhaps revisit this later
	hints.ai_family = AF_INET; // AF_UNSPEC
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_protocol = IPPROTO_TCP;

	int retVal = getaddrinfo(getHostName().c_str(), NULL, &hints, &addrInfo);
	if (addrInfo == NULL || retVal != 0) {
		SLog(EWarn, "Could not retrieve the computer's fully "
			"qualified domain name: could not resolve host address \"%s\"!",
			getHostName().c_str());
		return getHostName();
	}

	char fqdn[NI_MAXHOST];
	retVal = getnameinfo(addrInfo->ai_addr, sizeof(struct sockaddr_in), 
		fqdn, NI_MAXHOST, NULL, 0, 0);
	if (retVal != 0) {
		freeaddrinfo(addrInfo);
#if defined(WIN32)
		SLog(EWarn, "Could not retrieve the computer's fully "
			"qualified domain name: error %i!", WSAGetLastError());
#else
		SLog(EWarn, "Could not retrieve the computer's fully "
			"qualified domain name: error %i!", gai_strerror(retVal));
#endif
		return getHostName();
	}
		
	freeaddrinfo(addrInfo);

	return fqdn;
}

Float log2(Float value) {
	const Float invLn2 = (Float) 1.0f / std::fastlog((Float) 2.0f);
	return std::fastlog(value) * invLn2;
}

std::string formatString(const char *fmt, ...) {
	char tmp[512];
	va_list iterator;

#if defined(WIN32)
	va_start(iterator, fmt);
	size_t size = _vscprintf(fmt, iterator) + 1;

	if (size >= sizeof(tmp)) {
		char *dest = new char[size];
		vsnprintf_s(dest, size, size-1, fmt, iterator);
		va_end(iterator);
		std::string result(dest);
		delete[] dest;
		return result;
	}

	vsnprintf_s(tmp, size, size-1, fmt, iterator);
	va_end(iterator);
#else
	va_start(iterator, fmt);
	size_t size = vsnprintf(tmp, sizeof(tmp), fmt, iterator);
	va_end(iterator);

	if (size >= sizeof(tmp)) {
		/* Overflow! -- dynamically allocate memory */
		char *dest = new char[size+1];
		va_start(iterator, fmt);
		vsnprintf(dest, size+1, fmt, iterator);
		va_end(iterator);

		std::string result(dest);
		delete[] dest;
		return result;
	}
#endif

	return std::string(tmp);
}

int log2i(uint32_t value) {
	int r = 0;
	while ((value >> r) != 0)
		r++;
	return r-1;
}

int log2i(uint64_t value) {
	int r = 0;
	while ((value >> r) != 0)
		r++;
	return r-1;
}

int modulo(int a, int b) {
	int result = a - (a/b) * b;
	return (result < 0) ? result+b : result;
}

Float modulo(Float a, Float b) {
	Float result = a - int(a/b) * b;
	return (result < 0) ? result+b : result;
}

/* Fast rounding & power-of-two test algorithms from PBRT */
uint32_t roundToPowerOfTwo(uint32_t i) {
	i--;
	i |= i >> 1; i |= i >> 2;
	i |= i >> 4; i |= i >> 8;
	i |= i >> 16;
	return i+1;
}

uint64_t roundToPowerOfTwo(uint64_t i) {
	i--;
	i |= i >> 1;  i |= i >> 2;
	i |= i >> 4;  i |= i >> 8;
	i |= i >> 16; i |= i >> 32;
	return i+1;
}

// -----------------------------------------------------------------------
//  Numerical utility functions
// -----------------------------------------------------------------------

bool solveQuadratic(Float a, Float b, Float c, Float &x0, Float &x1) {
	/* Linear case */
	if (a == 0) {
		if (b != 0) {
			x0 = x1 = -c / b;
			return true;
		}
		return false;
	}

	Float discrim = b*b - 4.0f*a*c;

	/* Leave if there is no solution */
	if (discrim < 0)
		return false;

	Float temp, sqrtDiscrim = std::sqrt(discrim);

	/* Numerically stable version of (-b (+/-) sqrtDiscrim) / (2 * a)
	 *
	 * Based on the observation that one solution is always
	 * accurate while the other is not. Finds the solution of
	 * greater magnitude which does not suffer from loss of
	 * precision and then uses the identity x1 * x2 = c / a
	 */
	if (b < 0)
		temp = -0.5f * (b - sqrtDiscrim);
	else
		temp = -0.5f * (b + sqrtDiscrim);

	x0 = temp / a;
	x1 = c / temp;

	/* Return the results so that x0 < x1 */
	if (x0 > x1)
		std::swap(x0, x1);

	return true;
}

bool solveQuadraticDouble(double a, double b, double c, double &x0, double &x1) {
	/* Linear case */
	if (a == 0) {
		if (b != 0) {
			x0 = x1 = -c / b;
			return true;
		}
		return false;
	}

	double discrim = b*b - 4.0f*a*c;

	/* Leave if there is no solution */
	if (discrim < 0)
		return false;

	double temp, sqrtDiscrim = std::sqrt(discrim);

	/* Numerically stable version of (-b (+/-) sqrtDiscrim) / (2 * a)
	 *
	 * Based on the observation that one solution is always
	 * accurate while the other is not. Finds the solution of
	 * greater magnitude which does not suffer from loss of
	 * precision and then uses the identity x1 * x2 = c / a
	 */
	if (b < 0)
		temp = -0.5f * (b - sqrtDiscrim);
	else
		temp = -0.5f * (b + sqrtDiscrim);

	x0 = temp / a;
	x1 = c / temp;

	/* Return the results so that x0 < x1 */
	if (x0 > x1)
		std::swap(x0, x1);

	return true;
}

bool solveLinearSystem2x2(const Float a[2][2], const Float b[2], Float x[2]) {
	Float det = a[0][0] * a[1][1] - a[0][1] * a[1][0];

	if (std::abs(det) < Epsilon)
		return false;

	Float inverse = (Float) 1.0f / det;

	x[0] = (a[1][1] * b[0] - a[0][1] * b[1]) * inverse;
	x[1] = (a[0][0] * b[1] - a[1][0] * b[0]) * inverse;

	return true;
}

Float interpCubic1D(Float p, const Float *data, Float min, Float max, size_t size) {
	if (p < min || p > max)
		return 0.0f;

	/* Transform 'p' so that knots lie at integer positions */
	Float t = ((p - min) * (size - 1)) / (max - min);

	/* Index of the left knot in the queried subinterval,
	   be robust to cases where t=b. */
	size_t k = std::min((size_t) t, size - 2);

	Float f0  = data[k], 
		  f1  = data[k+1],
		  d0, d1;

	/* Approximate the derivatives */
	if (k > 0)
		d0 = 0.5f * (data[k+1] - data[k-1]);
	else
		d0 = data[k+1] - data[k];

	if (k + 2 < size)
		d1 = 0.5f * (data[k+2] - data[k]);
	else
		d1 = data[k+1] - data[k];

	/* Compute the relative position within the interval */
	t = t - (Float) k;

	Float t2 = t*t, t3 = t2*t;

	return 
		( 2*t3 - 3*t2 + 1) * f0 +
		(-2*t3 + 3*t2)     * f1 +
		(   t3 - 2*t2 + t) * d0 +
		(   t3 - t2)       * d1;
}


Float interpCubic2D(const Point2 &p, const Float *data, 
		const Point2 &min, const Point2 &max, const Size2 &size) {
	Float knotWeights[2][4];
	Point2 knot;

	/* Compute interpolation weights separately for each dimension */
	for (int dim=0; dim<2; ++dim) {
		Float *weights = knotWeights[dim];
		if (p[dim] < min[dim] || p[dim] > max[dim])
			return 0.0f;

		/* Transform 'p' so that knots lie at integer positions */
		Float t = ((p[dim] - min[dim]) * (size[dim] - 1))
			/ (max[dim]-min[dim]);

		/* Find the index of the left knot in the queried 
		   subinterval, be robust to cases where t=b. */
		knot[dim] = std::min((size_t) t, size[dim] - 2);

		/* Compute the relative position within the interval */
		t = t - (Float) knot[dim];

		/* Compute node weights */
		Float t2 = t*t, t3 = t2*t;
		weights[0] = 0.0f;
		weights[1] = 2*t3 - 3*t2 + 1;
		weights[2] = -2*t3 + 3*t2;
		weights[3] = 0.0f;

		/* Derivative weights */
		Float d0 = t3 - 2*t2 + t,
			  d1 = t3 - t2;

		/* Turn derivative weights into node weights using
		   an appropriate chosen finite differences stencil */
		if (knot[dim] > 0) {
			weights[2] +=  0.5f * d0;
			weights[0] -=  0.5f * d0;
		} else {
			weights[2] += d0;
			weights[1] -= d0;
		}

		if (knot[dim] + 2 < size[dim]) {
			weights[3] += 0.5f * d1;
			weights[1] -= 0.5f * d1;
		} else {
			weights[2] += d1;
			weights[1] -= d1;
		}
	}

	Float result = 0.0f;
	for (int y=-1; y<=2; ++y) {
		Float wy = knotWeights[1][y+1];
		for (int x=-1; x<=2; ++x) {
			Float wxy = knotWeights[0][x+1] * wy;

			if (wxy == 0)
				continue;

			size_t pos = (knot[1] + y) * size[0] + knot[0] + x;

			result += data[pos] * wxy;
		}
	}
	return result;
}

Float interpCubic3D(const Point3 &p, const Float *data, 
		const Point3 &min, const Point3 &max, const Size3 &size) {
	Float knotWeights[3][4];
	Point3 knot;

	/* Compute interpolation weights separately for each dimension */
	for (int dim=0; dim<3; ++dim) {
		Float *weights = knotWeights[dim];
		if (p[dim] < min[dim] || p[dim] > max[dim])
			return 0.0f;

		/* Transform 'p' so that knots lie at integer positions */
		Float t = ((p[dim] - min[dim]) * (size[dim] - 1))
			/ (max[dim]-min[dim]);

		/* Find the index of the left knot in the queried 
		   subinterval, be robust to cases where t=b. */
		knot[dim] = std::min((size_t) t, size[dim] - 2);

		/* Compute the relative position within the interval */
		t = t - (Float) knot[dim];

		/* Compute node weights */
		Float t2 = t*t, t3 = t2*t;
		weights[0] = 0.0f;
		weights[1] = 2*t3 - 3*t2 + 1;
		weights[2] = -2*t3 + 3*t2;
		weights[3] = 0.0f;

		/* Derivative weights */
		Float d0 = t3 - 2*t2 + t,
			  d1 = t3 - t2;

		/* Turn derivative weights into node weights using
		   an appropriate chosen finite differences stencil */
		if (knot[dim] > 0) {
			weights[2] +=  0.5f * d0;
			weights[0] -=  0.5f * d0;
		} else {
			weights[2] += d0;
			weights[1] -= d0;
		}

		if (knot[dim] + 2 < size[dim]) {
			weights[3] += 0.5f * d1;
			weights[1] -= 0.5f * d1;
		} else {
			weights[2] += d1;
			weights[1] -= d1;
		}
	}

	Float result = 0.0f;
	for (int z=-1; z<=2; ++z) {
		Float wz = knotWeights[2][z+1];
		for (int y=-1; y<=2; ++y) {
			Float wyz = knotWeights[1][y+1] * wz;
			for (int x=-1; x<=2; ++x) {
				Float wxyz = knotWeights[0][x+1] * wyz;

				if (wxyz == 0)
					continue;

				size_t pos = ((knot[2] + z) * size[1] + (knot[1] + y))
					* size[0] + knot[0] + x;

				result += data[pos] * wxyz;
			}
		}
	}
	return result;
}

void stratifiedSample1D(Random *random, Float *dest, int count, bool jitter) {
	Float invCount = 1.0f / count;

	for (int i=0; i<count; i++) {
		Float offset = jitter ? random->nextFloat() : 0.5f;
		*dest++ = (i + offset) * invCount;
	}
}

void stratifiedSample2D(Random *random, Point2 *dest, int countX, int countY, bool jitter) {
	Float invCountX = 1.0f / countX;
	Float invCountY = 1.0f / countY;

	for (int x=0; x<countX; x++) {
		for (int y=0; y<countY; y++) {
			Float offsetX = jitter ? random->nextFloat() : 0.5f;
			Float offsetY = jitter ? random->nextFloat() : 0.5f;
			*dest++ = Point2(
				(x + offsetX) * invCountX,
				(y + offsetY) * invCountY
			);
		}
	}
}

void latinHypercube(Random *random, Float *dest, size_t nSamples, size_t nDim) {
	Float delta = 1 / (Float) nSamples;
	for (size_t i = 0; i < nSamples; ++i)
		for (size_t j = 0; j < nDim; ++j)
			dest[nDim * i + j] = (i + random->nextFloat()) * delta;
	for (size_t i = 0; i < nDim; ++i) {
		for (size_t j = 0; j < nSamples; ++j) {
			size_t other = random->nextSize(nSamples);
			std::swap(dest[nDim * j + i], dest[nDim * other + i]);
		}
	}
}

Vector sphericalDirection(Float theta, Float phi) {
	Float sinTheta, cosTheta, sinPhi, cosPhi;

	std::sincos(theta, &sinTheta, &cosTheta);
	std::sincos(phi, &sinPhi, &cosPhi);

	return Vector(
		sinTheta * cosPhi,
		sinTheta * sinPhi,
		cosTheta
	);
}

Vector squareToSphere(const Point2 &sample) {
	Float z = 1.0f - 2.0f * sample.y;
	Float r = std::sqrt(std::max((Float) 0.0f, 1.0f - z*z));
	Float sinPhi, cosPhi;
	std::sincos(2.0f * M_PI * sample.x, &sinPhi, &cosPhi);
	return Vector(r * cosPhi, r * sinPhi, z);
}

Vector squareToHemisphere(const Point2 &sample) {
	Float z = sample.y;
	Float tmp = std::sqrt(std::min((Float) 0, 1-z*z));

	Float sinPhi, cosPhi;
	std::sincos(2.0f * M_PI * sample.x, &sinPhi, &cosPhi);

	return Vector(cosPhi * tmp, sinPhi * tmp, z);
}

Vector squareToHemispherePSA(const Point2 &sample) {
	Point2 p = squareToDiskConcentric(sample);
	Float z = std::sqrt(std::max((Float) 0, 
		1.0f - p.x*p.x - p.y*p.y));

	/* Guard against numerical imprecisions */
	if (EXPECT_NOT_TAKEN(z == 0))
		z = 1e-10f;

	return Vector(p.x, p.y, z);
}

Point2 squareToDisk(const Point2 &sample) {
	Float r = std::sqrt(sample.x);
	Float sinPhi, cosPhi;
	std::sincos(2.0f * M_PI * sample.y, &sinPhi, &cosPhi);

	return Point2(
		cosPhi * r,
		sinPhi * r
	);
}

void coordinateSystem(const Vector &a, Vector &b, Vector &c) {
	if (std::abs(a.x) > std::abs(a.y)) {
		Float invLen = 1.0f / std::sqrt(a.x * a.x + a.z * a.z);
		c = Vector(a.z * invLen, 0.0f, -a.x * invLen);
	} else {
		Float invLen = 1.0f / std::sqrt(a.y * a.y + a.z * a.z);
		c = Vector(0.0f, a.z * invLen, -a.y * invLen);
	}
	b = cross(c, a);
}

Point2 squareToTriangle(const Point2 &sample) {
	Float a = std::sqrt(1.0f - sample.x);
	return Point2(1 - a, a * sample.y);
}

Point2 toSphericalCoordinates(const Vector &v) {
	Point2 result(
		std::acos(v.z),
		std::atan2(v.y, v.x)
	);
	if (result.y < 0)
		result.y += 2*M_PI;
	return result;
}

Point2 squareToDiskConcentric(const Point2 &sample) {
	Float r1 = 2.0f*sample.x - 1.0f;
	Float r2 = 2.0f*sample.y - 1.0f;

	Point2 coords;
	if (r1 == 0 && r2 == 0) {
		coords = Point2(0, 0);
	} else if (r1 > -r2) { /* Regions 1/2 */
		if (r1 > r2)
			coords = Point2(r1, (M_PI/4.0f) * r2/r1);
		else
			coords = Point2(r2, (M_PI/4.0f) * (2.0f - r1/r2));
	} else { /* Regions 3/4 */
		if (r1<r2)
			coords = Point2(-r1, (M_PI/4.0f) * (4.0f + r2/r1));
		else 
			coords = Point2(-r2, (M_PI/4.0f) * (6.0f - r1/r2));
	}

	Point2 result;
	std::sincos(coords.y, &result.y, &result.x);
	return result*coords.x;
}

Point2 diskToSquareConcentric(const Point2 &p) {
	Float r   = std::sqrt(p.x * p.x + p.y * p.y),
		  phi = std::atan2(p.y, p.x),
		  a, b;

	if (phi < -M_PI/4) {
  		/* in range [-pi/4,7pi/4] */
		phi += 2*M_PI;
	}

	if (phi < M_PI/4) { /* region 1 */
		a = r;
		b = phi * a / (M_PI/4);
	} else if (phi < 3*M_PI/4) { /* region 2 */
		b = r;
		a = -(phi - M_PI/2) * b / (M_PI/4);
	} else if (phi < 5*M_PI/4) { /* region 3 */
		a = -r;
		b = (phi - M_PI) * a / (M_PI/4);
	} else { /* region 4 */
		b = -r;
		a = -(phi - 3*M_PI/2) * b / (M_PI/4);
	}

	return Point2(0.5f * (a+1), 0.5f * (b+1));
}

Float squareToConePdf(Float cosCutoff) {
	return 1 / (2 * M_PI * (1 - cosCutoff));
}

Vector squareToCone(Float cosCutoff, const Point2 &sample) {
	Float cosTheta = (1-sample.x) + sample.x * cosCutoff;
	Float sinTheta = std::sqrt(1 - cosTheta * cosTheta);
	Float phi = sample.y * (2 * M_PI);
	return Vector(std::cos(phi) * sinTheta,
		std::sin(phi) * sinTheta, cosTheta);
}

Point2 squareToStdNormal(const Point2 &sample) {
	Float r   = std::sqrt(-2 * std::fastlog(1-sample.x)),
		  phi = 2 * M_PI * sample.y;
	Point2 result;
	std::sincos(phi, &result.y, &result.x);
	return result * r;
}

Float lanczosSinc(Float t, Float tau) {
	t = std::abs(t);
	if (t < Epsilon)
		return 1.0f;
	else if (t > 1.0f)
		return 0.0f;
	t *= M_PI;
	Float sincTerm = std::sin(t*tau)/(t*tau);
	Float windowTerm = std::sin(t)/t;
	return sincTerm * windowTerm;
}

/* The following functions calculate the reflected and refracted 
	directions in addition to the fresnel coefficients. Based on 
	PBRT and the paper "Derivation of Refraction Formulas" 
	by Paul S. Heckbert. */
Float fresnelDielectric(Float cosThetaI, Float cosThetaT, 
						Float etaI, Float etaT) {
	if (etaI == etaT)
		return 0.0f;

	Float Rs = (etaI * cosThetaI - etaT * cosThetaT)
			/ (etaI * cosThetaI + etaT * cosThetaT);
	Float Rp = (etaT * cosThetaI - etaI * cosThetaT)
			/ (etaT * cosThetaI + etaI * cosThetaT);

	return (Rs * Rs + Rp * Rp) / 2.0f;
}

Spectrum fresnelConductor(Float cosThetaI, const Spectrum &eta, const Spectrum &k) {
	Spectrum tmp = (eta*eta + k*k) * (cosThetaI * cosThetaI);

	Spectrum rParl2 = (tmp - (eta * (2.0f * cosThetaI)) + Spectrum(1.0f))
					/ (tmp + (eta * (2.0f * cosThetaI)) + Spectrum(1.0f));

	Spectrum tmpF = eta*eta + k*k;

	Spectrum rPerp2 = (tmpF - (eta * (2.0f * cosThetaI)) + Spectrum(cosThetaI*cosThetaI)) /
					  (tmpF + (eta * (2.0f * cosThetaI)) + Spectrum(cosThetaI*cosThetaI));

	return (rParl2 + rPerp2) / 2.0f;
}

Float fresnel(Float cosThetaI, Float extIOR, Float intIOR) {
	Float etaI = extIOR, etaT = intIOR;

	/* Swap the indices of refraction if the interaction starts
	   at the inside of the object */
	if (cosThetaI < 0.0f)
		std::swap(etaI, etaT);

	/* Using Snell's law, calculate the squared sine of the
	   angle between the normal and the transmitted ray */
	Float eta = etaI / etaT,
		  sinThetaTSqr = eta*eta * (1-cosThetaI*cosThetaI);

	if (sinThetaTSqr > 1.0f)
		return 1.0f;  /* Total internal reflection! */

	Float cosThetaT = std::sqrt(1.0f - sinThetaTSqr);

	/* Finally compute the reflection coefficient */
	return fresnelDielectric(std::abs(cosThetaI),
		cosThetaT, etaI, etaT);
}

static Float fresnelDiffuseIntegrand(Float eta, Float xi) {
	if (eta > 1)
		return fresnel(std::sqrt(xi), 1, eta);
	else
		return fresnel(std::sqrt(xi), 1/eta, 1);
}

Float fresnelDiffuseReflectance(Float eta, bool fast) {
	if (fast) {
		/* Fast mode: the following code approximates the
		 * diffuse Frensel reflectance for the eta<1 and 
		 * eta>1 cases. An evalution of the accuracy led
		 * to the following scheme, which cherry-picks
		 * fits from two papers where they are best.
		 */
		if (eta < 1) {
			/* Fit by Egan and Hilgeman (1973). Works
			   reasonably well for "normal" IOR values (<2).
	
			   Max rel. error in 1.0 - 1.5 : 0.1%
			   Max rel. error in 1.5 - 2   : 0.6%
			   Max rel. error in 2.0 - 5   : 9.5%
			*/
			return -1.4399f * (eta * eta) 
				  + 0.7099f * eta 
				  + 0.6681f 
				  + 0.0636f / eta;
		} else {
			/* Fit by d'Eon and Irving (2011)
			 *
			 * Maintains a good accuracy even for 
			 * unrealistic IOR values.
			 *
			 * Max rel. error in 1.0 - 2.0   : 0.1%
			 * Max rel. error in 2.0 - 10.0  : 0.2%
			 */
			Float invEta = 1.0f / eta,
				  invEta2 = invEta*invEta,
				  invEta3 = invEta2*invEta,
				  invEta4 = invEta3*invEta,
				  invEta5 = invEta4*invEta;

			return 0.919317f - 3.4793f * invEta 
				 + 6.75335f * invEta2
				 - 7.80989f * invEta3 
				 + 4.98554f * invEta4 
				 - 1.36881f * invEta5;
		}
	} else {
		GaussLobattoIntegrator quad(1024, 0, 1e-5f);
		return quad.integrate(
			boost::bind(&fresnelDiffuseIntegrand, eta, _1), 0, 1);
	}

	return 0.0f;
}

Float radicalInverse(int b, size_t i) {
	Float invB = (Float) 1 / (Float) b;
	Float x = 0.0f, f = invB;
	
	while (i) {
		x += f * (Float) (i % b);
		i /= b;
		f *= invB;
	}
	return x;
}

Float radicalInverseIncremental(int b, Float x) {
	Float invB = (Float) 1 / (Float) b;
	Float h, hh, r = 1.0f - x - (Float) 1e-10;

	if (invB < r) {
		x += invB;
	} else {
		h = invB;
		do {
			hh = h;
			h *= invB;
		} while (h >= r);

		x += hh + h - 1.0f;
	}
	return x;
}

std::string timeString(Float time, bool precise) {
	std::ostringstream os;
	char suffix = 's';
#ifdef WIN32
	if (mts_isnan(time) || std::isinf(time)) {
#else
	if (mts_isnan(time) || std::fpclassify(time) == FP_INFINITE) {
#endif
		return "inf";
	}
	os << std::setprecision(precise ? 4 : 1) << std::fixed;
	if (time > 60) {
		time /= 60; suffix = 'm';
		if (time > 60) {
			time /= 60; suffix = 'h';
			if (time > 12) {
				time /= 12; suffix = 'd';
			}
		}
	}
	os << time << suffix;
	return os.str();
}

double normalQuantile(double p) {
	// By Peter J. Acklam
	// http://home.online.no/~pjacklam/notes/invnorm/impl/sprouse/ltqnorm.c
	static const double LOW = 0.02425;
	static const double HIGH = 0.97575;
	double q, r;

	/* Coefficients in rational approximations. */
	static const double a[] = {
		-3.969683028665376e+01,
		 2.209460984245205e+02,
		-2.759285104469687e+02,
		 1.383577518672690e+02,
		-3.066479806614716e+01,
		 2.506628277459239e+00
	};
	static const double b[] = {
		-5.447609879822406e+01,
		 1.615858368580409e+02,
		-1.556989798598866e+02,
		 6.680131188771972e+01,
		-1.328068155288572e+01
	};
	static const double c[] = {
		-7.784894002430293e-03,
		-3.223964580411365e-01,
		-2.400758277161838e+00,
		-2.549732539343734e+00,
		 4.374664141464968e+00,
		 2.938163982698783e+00
	};
	static const double d[] = {
		7.784695709041462e-03,
		3.224671290700398e-01,
		2.445134137142996e+00,
		3.754408661907416e+00
	};

	errno = 0;

	if (p < 0 || p > 1) {
		errno = EDOM;
		return 0.0;
	} else if (p == 0) {
		errno = ERANGE;
		return -HUGE_VAL /* minus "infinity" */;
	} else if (p == 1) {
		errno = ERANGE;
		return HUGE_VAL /* "infinity" */;
	} else if (p < LOW) {
		/* Rational approximation for lower region */
		q = sqrt(-2*log(p));
		return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	} else if (p > HIGH) {
		/* Rational approximation for upper region */
		q  = sqrt(-2*log(1-p));
		return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	} else {
		/* Rational approximation for central region */
		q = p - 0.5;
		r = q*q;
		return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
			(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
	}
}

Float hypot2(Float a, Float b) {
	Float r;
	if (std::abs(a) > std::abs(b)) {
		r = b / a;
		r = std::abs(a) * std::sqrt(1 + r*r);
	} else if (b != 0) {
		r = a / b;
		r = std::abs(b) * std::sqrt(1 + r*r);
	} else {
		r = 0;
	}
	return r;
}
MTS_NAMESPACE_END

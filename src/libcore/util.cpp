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

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/util.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/quad.h>
#include <mitsuba/core/sse.h>
#include <mitsuba/core/frame.h>
#include <boost/bind.hpp>
#include <stdarg.h>
#include <iomanip>
#include <errno.h>

#if defined(__OSX__)
#include <sys/sysctl.h>
#include <mach/mach.h>
#elif defined(__WINDOWS__)
#include <windows.h>
#include <direct.h>
#include <psapi.h>
#else
#include <malloc.h>
#endif

#if defined(__WINDOWS__)
# include <windows.h>
# include <winsock2.h>
# include <ws2tcpip.h>
#else
# include <sys/types.h>
# include <sys/socket.h>
# include <netdb.h>
# include <fenv.h>
#endif

// SSE is not enabled in general when using double precision, however it is
// required in OS X for FP exception handling
#if defined(__OSX__) && !defined(MTS_SSE)
#include <xmmintrin.h>
#undef enable_fpexcept_sse
#undef query_fpexcept_sse
#undef disable_fpexcept_sse

namespace {
inline void enable_fpexcept_sse() {
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() &
        ~(_MM_MASK_INVALID | _MM_MASK_DIV_ZERO));
}
inline unsigned int query_fpexcept_sse() {
    return (~_MM_GET_EXCEPTION_MASK() &
        (_MM_MASK_INVALID | _MM_MASK_DIV_ZERO));
}
inline void disable_fpexcept_sse() {
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() |
        _MM_MASK_INVALID | _MM_MASK_DIV_ZERO);
}
} // namespace

#endif

MTS_NAMESPACE_BEGIN

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

void * __restrict allocAligned(size_t size) {
#if defined(__WINDOWS__)
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
#if defined(__WINDOWS__)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

static int __cached_core_count = 0;

int getCoreCount() {
    // assumes atomic word size memory access
    if (__cached_core_count)
        return __cached_core_count;

#if defined(__WINDOWS__)
    SYSTEM_INFO sys_info;
    GetSystemInfo(&sys_info);
    __cached_core_count = sys_info.dwNumberOfProcessors;
    return sys_info.dwNumberOfProcessors;
#elif defined(__OSX__)
    int nprocs;
    size_t nprocsSize = sizeof(int);
    if (sysctlbyname("hw.activecpu", &nprocs, &nprocsSize, NULL, 0))
        SLog(EError, "Could not detect the number of processors!");
    __cached_core_count = nprocs;
    return nprocs;
#else
    /* Determine the number of present cores */
    int nCores = sysconf(_SC_NPROCESSORS_CONF);

    /* Don't try to query CPU affinity if running inside Valgrind */
    if (getenv("VALGRIND_OPTS") == NULL) {
        /* Some of the cores may not be available to the user
           (e.g. on certain cluster nodes) -- determine the number
           of actual available cores here. */
        int nLogicalCores = nCores;
        size_t size = 0;
        cpu_set_t *cpuset = NULL;
        int retval = 0;

        /* The kernel may expect a larger cpu_set_t than would
           be warranted by the physical core count. Keep querying
           with increasingly larger buffers if the
           pthread_getaffinity_np operation fails */
        for (int i = 0; i<6; ++i) {
            size = CPU_ALLOC_SIZE(nLogicalCores);
            cpuset = CPU_ALLOC(nLogicalCores);
            if (!cpuset) {
                SLog(EWarn, "getCoreCount(): could not allocate cpu_set_t");
                goto done;
            }
            CPU_ZERO_S(size, cpuset);

            int retval = pthread_getaffinity_np(pthread_self(), size, cpuset);
            if (retval == 0)
                break;
            CPU_FREE(cpuset);
            nLogicalCores *= 2;
        }

        if (retval) {
            SLog(EWarn, "getCoreCount(): pthread_getaffinity_np(): could "
                "not read thread affinity map: %s", strerror(retval));
            goto done;
        }

        int availableCores = 0;
        for (int i=0; i<nLogicalCores; ++i)
            availableCores += CPU_ISSET_S(i, size, cpuset) ? 1 : 0;
        nCores = availableCores;
        CPU_FREE(cpuset);
    }

done:
    __cached_core_count = nCores;
    return nCores;
#endif
}

size_t getTotalSystemMemory() {
#if defined(__WINDOWS__)
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return (size_t) status.ullTotalPhys;
#elif defined(__OSX__)
    int mib[2] = { CTL_HW, HW_MEMSIZE };
    uint64_t size;
    size_t len = sizeof(size);
    if (sysctl(mib, 2, &size, &len, NULL, 0) < 0)
        return 0;
    return (size_t) size;
#else
    size_t pages = sysconf(_SC_PHYS_PAGES);
    size_t page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
#endif
}

size_t getPrivateMemoryUsage() {
#if defined(__WINDOWS__)
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS *) &pmc, sizeof(pmc));
    return (size_t) pmc.PrivateUsage; /* Process-private memory usage (RAM + swap) */
#elif defined(__OSX__)
    struct task_basic_info_64 t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_64_COUNT;

    if (task_info(mach_task_self(), TASK_BASIC_INFO_64,
            (task_info_t)&t_info, &t_info_count) != KERN_SUCCESS)
        return 0;

    return (size_t) t_info.resident_size; /* Not exactly what we want -- oh well.. */
#else
    FILE* file = fopen("/proc/self/status", "r");
    if (!file)
        return 0;

    char buffer[128];
    size_t result = 0;
    while (fgets(buffer, sizeof(buffer), file) != NULL) {
        if (strncmp(buffer, "VmRSS:", 6) != 0 && /* Non-swapped physical memory specific to this process */
            strncmp(buffer, "VmSwap:", 7) != 0)  /* Swapped memory specific to this process */
            continue;

        char *line = buffer;
        while (*line < '0' || *line > '9')
            ++line;
        line[strlen(line)-3] = '\0';
        result += (size_t) atoi(line) * 1024;
    }

    fclose(file);
    return result;
#endif
}

#if defined(__WINDOWS__)
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
#if defined(__WINDOWS__)
    _clearfp();
    uint32_t cw = _controlfp(0, 0);
    exceptionsWereEnabled = ~cw & (_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
    cw &= ~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
    _controlfp(cw, _MCW_EM);
#elif defined(__OSX__)
    exceptionsWereEnabled = query_fpexcept_sse() != 0;
#else
    exceptionsWereEnabled =
        fegetexcept() & (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
    feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
    enable_fpexcept_sse();
    return exceptionsWereEnabled;
}

bool disableFPExceptions() {
    bool exceptionsWereEnabled = false;
#if defined(__WINDOWS__)
    _clearfp();
    uint32_t cw = _controlfp(0, 0);
    exceptionsWereEnabled = ~cw & (_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
    cw |= _EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW;
    _controlfp(cw, _MCW_EM);
#elif defined(__OSX__)
    exceptionsWereEnabled = query_fpexcept_sse() != 0;
#else
    exceptionsWereEnabled =
        fegetexcept() & (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
    fedisableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
    disable_fpexcept_sse();
    return exceptionsWereEnabled;
}

void restoreFPExceptions(bool oldState) {
    bool currentState;
#if defined(__WINDOWS__)
    uint32_t cw = _controlfp(0, 0);
    currentState = ~cw & (_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
#elif defined(__OSX__)
    currentState = query_fpexcept_sse() != 0;
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
#if defined(__WINDOWS__)
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
#if defined(__WINDOWS__)
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

std::string formatString(const char *fmt, ...) {
    char tmp[512];
    va_list iterator;

#if defined(__WINDOWS__)
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

    if (std::abs(det) <= RCPOVERFLOW)
        return false;

    Float inverse = (Float) 1.0f / det;

    x[0] = (a[1][1] * b[0] - a[0][1] * b[1]) * inverse;
    x[1] = (a[0][0] * b[1] - a[1][0] * b[0]) * inverse;

    return true;
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

    math::sincos(theta, &sinTheta, &cosTheta);
    math::sincos(phi, &sinPhi, &cosPhi);

    return Vector(
        sinTheta * cosPhi,
        sinTheta * sinPhi,
        cosTheta
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

void computeShadingFrame(const Vector &n, const Vector &dpdu, Frame &frame) {
    frame.n = n;
    frame.s = normalize(dpdu - frame.n
        * dot(frame.n, dpdu));
    frame.t = cross(frame.n, frame.s);
}

void computeShadingFrameDerivative(const Vector &n, const Vector &dpdu, const Vector &dndu, const Vector &dndv, Frame &du, Frame &dv) {
    Vector s = dpdu - n * dot(n, dpdu);
    Float invLen_s = 1.0f / s.length();
    s *= invLen_s;

    du.s = invLen_s * (-dndu * dot(n, dpdu) - n * dot(dndu, dpdu));
    dv.s = invLen_s * (-dndv * dot(n, dpdu) - n * dot(dndv, dpdu));

    du.s -= s * dot(du.s, s);
    dv.s -= s * dot(dv.s, s);

    du.t = cross(dndu, s) + cross(n, du.s);
    dv.t = cross(dndv, s) + cross(n, dv.s);

    du.n = dndu;
    dv.n = dndv;
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

Float fresnelDielectric(Float cosThetaI, Float cosThetaT, Float eta) {
    if (EXPECT_NOT_TAKEN(eta == 1))
        return 0.0f;

    Float Rs = (cosThetaI - eta * cosThetaT)
             / (cosThetaI + eta * cosThetaT);
    Float Rp = (eta * cosThetaI - cosThetaT)
             / (eta * cosThetaI + cosThetaT);

    /* No polarization -- return the unpolarized reflectance */
    return 0.5f * (Rs * Rs + Rp * Rp);
}

Float fresnelDielectricExt(Float cosThetaI_, Float &cosThetaT_, Float eta) {
    if (EXPECT_NOT_TAKEN(eta == 1)) {
        cosThetaT_ = -cosThetaI_;
        return 0.0f;
    }

    /* Using Snell's law, calculate the squared sine of the
       angle between the normal and the transmitted ray */
    Float scale = (cosThetaI_ > 0) ? 1/eta : eta,
          cosThetaTSqr = 1 - (1-cosThetaI_*cosThetaI_) * (scale*scale);

    /* Check for total internal reflection */
    if (cosThetaTSqr <= 0.0f) {
        cosThetaT_ = 0.0f;
        return 1.0f;
    }

    /* Find the absolute cosines of the incident/transmitted rays */
    Float cosThetaI = std::abs(cosThetaI_);
    Float cosThetaT = std::sqrt(cosThetaTSqr);

    Float Rs = (cosThetaI - eta * cosThetaT)
             / (cosThetaI + eta * cosThetaT);
    Float Rp = (eta * cosThetaI - cosThetaT)
             / (eta * cosThetaI + cosThetaT);

    cosThetaT_ = (cosThetaI_ > 0) ? -cosThetaT : cosThetaT;

    /* No polarization -- return the unpolarized reflectance */
    return 0.5f * (Rs * Rs + Rp * Rp);
}

Float fresnelConductorApprox(Float cosThetaI, Float eta, Float k) {
    Float cosThetaI2 = cosThetaI*cosThetaI;

    Float tmp = (eta*eta + k*k) * cosThetaI2;

    Float Rp2 = (tmp - (eta * (2 * cosThetaI)) + 1)
              / (tmp + (eta * (2 * cosThetaI)) + 1);

    Float tmpF = eta*eta + k*k;

    Float Rs2 = (tmpF - (eta * (2 * cosThetaI)) + cosThetaI2) /
                (tmpF + (eta * (2 * cosThetaI)) + cosThetaI2);

    return 0.5f * (Rp2 + Rs2);
}

Spectrum fresnelConductorApprox(Float cosThetaI, const Spectrum &eta, const Spectrum &k) {
    Float cosThetaI2 = cosThetaI*cosThetaI;

    Spectrum tmp = (eta*eta + k*k) * cosThetaI2;

    Spectrum Rp2 = (tmp - (eta * (2 * cosThetaI)) + Spectrum(1.0f))
                 / (tmp + (eta * (2 * cosThetaI)) + Spectrum(1.0f));

    Spectrum tmpF = eta*eta + k*k;

    Spectrum Rs2 = (tmpF - (eta * (2 * cosThetaI)) + Spectrum(cosThetaI2)) /
                   (tmpF + (eta * (2 * cosThetaI)) + Spectrum(cosThetaI2));

    return 0.5f * (Rp2 + Rs2);
}

Float fresnelConductorExact(Float cosThetaI, Float eta, Float k) {
    /* Modified from "Optics" by K.D. Moeller, University Science Books, 1988 */

    Float cosThetaI2 = cosThetaI*cosThetaI,
          sinThetaI2 = 1-cosThetaI2,
          sinThetaI4 = sinThetaI2*sinThetaI2;

    Float temp1 = eta*eta - k*k - sinThetaI2,
          a2pb2 = math::safe_sqrt(temp1*temp1 + 4*k*k*eta*eta),
          a     = math::safe_sqrt(0.5f * (a2pb2 + temp1));

    Float term1 = a2pb2 + cosThetaI2,
          term2 = 2*a*cosThetaI;

    Float Rs2 = (term1 - term2) / (term1 + term2);

    Float term3 = a2pb2*cosThetaI2 + sinThetaI4,
          term4 = term2*sinThetaI2;

    Float Rp2 = Rs2 * (term3 - term4) / (term3 + term4);

    return 0.5f * (Rp2 + Rs2);
}

Spectrum fresnelConductorExact(Float cosThetaI, const Spectrum &eta, const Spectrum &k) {
    /* Modified from "Optics" by K.D. Moeller, University Science Books, 1988 */

    Float cosThetaI2 = cosThetaI*cosThetaI,
          sinThetaI2 = 1-cosThetaI2,
          sinThetaI4 = sinThetaI2*sinThetaI2;

    Spectrum temp1 = eta*eta - k*k - Spectrum(sinThetaI2),
             a2pb2 = (temp1*temp1 + k*k*eta*eta*4).safe_sqrt(),
             a     = ((a2pb2 + temp1) * 0.5f).safe_sqrt();

    Spectrum term1 = a2pb2 + Spectrum(cosThetaI2),
             term2 = a*(2*cosThetaI);

    Spectrum Rs2 = (term1 - term2) / (term1 + term2);

    Spectrum term3 = a2pb2*cosThetaI2 + Spectrum(sinThetaI4),
             term4 = term2*sinThetaI2;

    Spectrum Rp2 = Rs2 * (term3 - term4) / (term3 + term4);

    return 0.5f * (Rp2 + Rs2);
}

Vector reflect(const Vector &wi, const Normal &n) {
    return 2 * dot(wi, n) * Vector(n) - wi;
}

Vector refract(const Vector &wi, const Normal &n, Float eta, Float cosThetaT) {
    if (cosThetaT < 0)
        eta = 1 / eta;

    return n * (dot(wi, n) * eta + cosThetaT) - wi * eta;
}

Vector refract(const Vector &wi, const Normal &n, Float eta) {
    if (EXPECT_NOT_TAKEN(eta == 1))
        return -wi;

    Float cosThetaI = dot(wi, n);
    if (cosThetaI > 0)
        eta = 1 / eta;

    /* Using Snell's law, calculate the squared sine of the
       angle between the normal and the transmitted ray */
    Float cosThetaTSqr = 1 - (1-cosThetaI*cosThetaI) * (eta*eta);

    /* Check for total internal reflection */
    if (cosThetaTSqr <= 0.0f)
        return Vector(0.0f);

    return n * (cosThetaI * eta - math::signum(cosThetaI)
        * std::sqrt(cosThetaTSqr)) - wi * eta;
}

Vector refract(const Vector &wi, const Normal &n, Float eta, Float &cosThetaT, Float &F) {
    Float cosThetaI = dot(wi, n);
    F = fresnelDielectricExt(cosThetaI, cosThetaT, eta);

    if (F == 1.0f) /* Total internal reflection */
        return Vector(0.0f);

    if (cosThetaT < 0)
        eta = 1 / eta;

    return n * (eta * cosThetaI + cosThetaT) - wi * eta;
}

namespace {
    /// Integrand used by fresnelDiffuseReflectance
    inline Float fresnelDiffuseIntegrand(Float eta, Float xi) {
        return fresnelDielectricExt(std::sqrt(xi), eta);
    }
};

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

std::string timeString(Float time, bool precise) {
    if (std::isnan(time) || std::isinf(time))
        return "inf";

    char suffix = 's';
    if (time > 60) {
        time /= 60; suffix = 'm';
        if (time > 60) {
            time /= 60; suffix = 'h';
            if (time > 12) {
                time /= 12; suffix = 'd';
            }
        }
    }

    std::ostringstream os;
    os << std::setprecision(precise ? 4 : 1)
       << std::fixed << time << suffix;

    return os.str();
}

std::string memString(size_t size, bool precise) {
    Float value = (Float) size;
    const char *suffixes[] = {
        "B", "KiB", "MiB", "GiB", "TiB", "PiB"
    };
    int suffix = 0;
    while (suffix < 5 && value > 1024.0f) {
        value /= 1024.0f; ++suffix;
    }

    std::ostringstream os;
    os << std::setprecision(suffix == 0 ? 0 : (precise ? 4 : 1))
       << std::fixed << value << " " << suffixes[suffix];

    return os.str();
}

MTS_NAMESPACE_END

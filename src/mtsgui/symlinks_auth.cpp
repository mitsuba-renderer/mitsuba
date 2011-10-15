#if defined(__APPLE__)
#ifndef __MAC_10_7
	/// Work around a few missing definitions in the MacOS X Lion header files
    #define __MAC_10_7 1070

    #if !defined(__MAC_OS_X_VERSION_MAX_ALLOWED) || (__MAC_OS_X_VERSION_MAX_ALLOWED == __MAC_10_6)
        #ifdef __MAC_OS_X_VERSION_MAX_ALLOWED
            #undef __MAC_OS_X_VERSION_MAX_ALLOWED
        #endif
        #define __MAC_OS_X_VERSION_MAX_ALLOWED __MAC_10_7
    #endif

    #if __MAC_OS_X_VERSION_MAX_ALLOWED < __MAC_10_7
        #define __AVAILABILITY_INTERNAL__MAC_10_7        __AVAILABILITY_INTERNAL_UNAVAILABLE
    #elif __MAC_OS_X_VERSION_MIN_REQUIRED < __MAC_10_7
        #define __AVAILABILITY_INTERNAL__MAC_10_7        __AVAILABILITY_INTERNAL_WEAK_IMPORT
    #else
        #define __AVAILABILITY_INTERNAL__MAC_10_7        __AVAILABILITY_INTERNAL_REGULAR
    #endif
    #if __MAC_OS_X_VERSION_MIN_REQUIRED >= __MAC_10_7
        #define __AVAILABILITY_INTERNAL__MAC_10_0_DEP__MAC_10_7        __AVAILABILITY_INTERNAL_DEPRECATED
        #define __AVAILABILITY_INTERNAL__MAC_10_1_DEP__MAC_10_7        __AVAILABILITY_INTERNAL_DEPRECATED
        #define __AVAILABILITY_INTERNAL__MAC_10_2_DEP__MAC_10_7        __AVAILABILITY_INTERNAL_DEPRECATED
        #define __AVAILABILITY_INTERNAL__MAC_10_3_DEP__MAC_10_7        __AVAILABILITY_INTERNAL_DEPRECATED
        #define __AVAILABILITY_INTERNAL__MAC_10_4_DEP__MAC_10_7        __AVAILABILITY_INTERNAL_DEPRECATED
        #define __AVAILABILITY_INTERNAL__MAC_10_5_DEP__MAC_10_7        __AVAILABILITY_INTERNAL_DEPRECATED
    #else
        #define __AVAILABILITY_INTERNAL__MAC_10_0_DEP__MAC_10_7        __AVAILABILITY_INTERNAL__MAC_10_0
        #define __AVAILABILITY_INTERNAL__MAC_10_1_DEP__MAC_10_7        __AVAILABILITY_INTERNAL__MAC_10_1
        #define __AVAILABILITY_INTERNAL__MAC_10_2_DEP__MAC_10_7        __AVAILABILITY_INTERNAL__MAC_10_2
        #define __AVAILABILITY_INTERNAL__MAC_10_3_DEP__MAC_10_7        __AVAILABILITY_INTERNAL__MAC_10_3
        #define __AVAILABILITY_INTERNAL__MAC_10_4_DEP__MAC_10_7        __AVAILABILITY_INTERNAL__MAC_10_4
        #define __AVAILABILITY_INTERNAL__MAC_10_5_DEP__MAC_10_7        __AVAILABILITY_INTERNAL__MAC_10_5
    #endif
    #define __AVAILABILITY_INTERNAL__MAC_10_7_DEP__MAC_NA             __AVAILABILITY_INTERNAL__MAC_10_7

#endif  /* !defined(__MAC_10_7) */

#include <Authorization.h>
#include <AuthorizationTags.h>
#include <unistd.h>
#include <iostream>

namespace mitsuba {
	extern std::string __mts_bundlepath();
};

bool create_symlinks() {
	AuthorizationFlags flags = kAuthorizationFlagDefaults;
	AuthorizationRef ref;
	OSStatus status = AuthorizationCreate(NULL, kAuthorizationEmptyEnvironment, flags, &ref);
	if (status != errAuthorizationSuccess)
		return false;

	AuthorizationItem items = {kAuthorizationRightExecute, 0, NULL, 0};
	AuthorizationRights rights = { 1, &items };
	flags = kAuthorizationFlagDefaults |
		kAuthorizationFlagInteractionAllowed | 
		kAuthorizationFlagPreAuthorize | 
		kAuthorizationFlagExtendRights;
	status = AuthorizationCopyRights(ref, &rights, NULL, flags, NULL);
	if (status != errAuthorizationSuccess) {
		AuthorizationFree(ref, kAuthorizationFlagDefaults);
		return false;
	}
	std::string bundlePath = mitsuba::__mts_bundlepath();
	std::string path = bundlePath + "/Contents/MacOS/symlinks_install";
	char *args[] = { const_cast<char *>(bundlePath.c_str()), NULL };
	FILE *pipe = NULL;
	flags = kAuthorizationFlagDefaults;
	status = AuthorizationExecuteWithPrivileges(ref, const_cast<char *>(path.c_str()), flags, args, &pipe);
	if (status != errAuthorizationSuccess) {
		AuthorizationFree(ref, kAuthorizationFlagDefaults);
		return false;
	}
	char buffer[128];
	for (;;) {
		int bytesRead = read(fileno(pipe), buffer, sizeof(buffer));
		if (bytesRead<1)
			break;
		write(fileno(stdout), buffer, bytesRead);
	}
	AuthorizationFree(ref, kAuthorizationFlagDefaults);
	return true;
}

#else
bool create_symlinks() { return false; }
#endif


#if defined(__APPLE__)
#include <Carbon/Carbon.h>
#include <Security/Authorization.h>
#include <Security/AuthorizationTags.h>
#include <unistd.h>
#include <iostream>
#include <sstream>

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
    std::ostringstream oss;
    oss << getuid();
    std::string uid = oss.str();
    char *args[] = { const_cast<char *>(bundlePath.c_str()), const_cast<char *>(uid.c_str()), NULL };
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


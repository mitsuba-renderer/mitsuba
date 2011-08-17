#if defined(__APPLE__)
#include <Authorization.h>
#include <AuthorizationTags.h>
#include <unistd.h>

bool create_symlinks() {
	AuthorizationFlags flags = kAuthorizationFlagDefaults;
	AuthorizationRef ref;

	OSStatus status = AuthorizationCreate(NULL, kAuthorizationEmptyEnvironment, flags, &ref);
	if (status != errAuthorizationSuccess) {
		return false;
	}

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
	char *path = "/bin/bash";
	char *args[] = { "data/install-symlinks.sh", NULL };
	FILE *pipe = pipe = NULL;
	flags = kAuthorizationFlagDefaults;
	status = AuthorizationExecuteWithPrivileges(ref, path, flags, args, &pipe);
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


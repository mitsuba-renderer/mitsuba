#include <mitsuba/mitsuba.h>
#import <Cocoa/Cocoa.h>

MTS_NAMESPACE_BEGIN

pthread_key_t __ubi_autorelease_key;
	
void __ubi_autorelease_init() {
	pthread_key_create(&__ubi_autorelease_key, NULL);
}

void __ubi_autorelease_shutdown() {
	pthread_key_delete(__ubi_autorelease_key);
}

void __ubi_autorelease_begin() {
	NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
	pthread_setspecific(__ubi_autorelease_key, pool);
}

void __ubi_autorelease_end() {
	NSAutoreleasePool *pool = 
		static_cast<NSAutoreleasePool *>(pthread_getspecific(__ubi_autorelease_key));
	[pool release]; 
}

void __ubi_chdir_to_bundlepath() {
	chdir([[[NSBundle mainBundle] bundlePath] fileSystemRepresentation]);
}

std::string __ubi_bundlepath() {
	return [[[NSBundle mainBundle] bundlePath] fileSystemRepresentation];
}

MTS_NAMESPACE_END

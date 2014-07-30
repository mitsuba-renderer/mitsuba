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
#import <Cocoa/Cocoa.h>

MTS_NAMESPACE_BEGIN

static pthread_key_t __mts_autorelease_key;
static bool __mts_cocoa_initialized = false;

void __mts_autorelease_init() {
	pthread_key_create(&__mts_autorelease_key, NULL);
}

void __mts_autorelease_shutdown() {
	pthread_key_delete(__mts_autorelease_key);
}

void __mts_autorelease_begin() {
	NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
	pthread_setspecific(__mts_autorelease_key, pool);
}

void __mts_autorelease_end() {
	NSAutoreleasePool *pool =
		static_cast<NSAutoreleasePool *>(pthread_getspecific(__mts_autorelease_key));
	[pool release];
}

void __mts_chdir_to_bundlepath() {
	chdir([[[NSBundle mainBundle] bundlePath] fileSystemRepresentation]);
}

std::string __mts_bundlepath() {
	return [[[NSBundle mainBundle] bundlePath] fileSystemRepresentation];
}

void __mts_init_cocoa() {
	if (!__mts_cocoa_initialized) {
		[NSApplication sharedApplication]; /* Creates a connection to the windowing environment */
		[NSApp activateIgnoringOtherApps:YES]; /* Pop to front */
		__mts_cocoa_initialized = true;
	}
}

void __mts_set_appdefaults() {
	if ([NSApp respondsToSelector:@selector(invalidateRestorableState)])
		[NSApp invalidateRestorableState];
#if 0

	/* Disable annoying OSX synchronization feature. It's not supported by Qt in any case.. */
	NSUserDefaults *prefs = [NSUserDefaults standardUserDefaults];
	if (![prefs boolForKey: @"ApplePersistenceIgnoreState"]) {
		[prefs setBool: YES forKey: @"ApplePersistenceIgnoreState"];
		[prefs synchronize];
	}
#endif
}

MTS_NAMESPACE_END

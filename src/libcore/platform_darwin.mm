/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

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

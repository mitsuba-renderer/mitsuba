#include <Breakpad/Breakpad.h>

void *__mts_init_breakpad_osx() {
  NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
  BreakpadRef breakpad = 0;
  NSDictionary *plist = [[NSBundle mainBundle] infoDictionary];
  if (plist) {
    // Note: version 1.0.0.4 of the framework changed the type of the argument
    // from CFDictionaryRef to NSDictionary * on the next line:
    breakpad = BreakpadCreate(plist);
  }
  [pool release];

  return breakpad;
}

void __mts_destroy_breakpad_osx(void *ptr) {
	BreakpadRelease((BreakpadRef) ptr);
}

#import <Cocoa/Cocoa.h>

extern int ubi_main(int argc, char **argv);

int main(int argc, char **argv) {
//	    return NSApplicationMain(argc,  (const char **) argv);
	NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
	[NSApplication sharedApplication]; /* Creates a connection to the windowing environment */
	int retval = ubi_main(argc, argv);
	[pool release];
	return retval;
}


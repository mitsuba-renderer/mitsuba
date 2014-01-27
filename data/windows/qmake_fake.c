/* Minimal QMake stub to trick CMake into using the Mitsuba Qt distribution */
#include <stdio.h>
#include <string.h>

/*
MSVC Compile flags:
  cl qmake_fake.c /nologo /Os /MT /arch:SSE2 /link /OUT:"qmake.exe"
*/

/* Qt version currently provided in the Mitsuba dependencies */
#define QT_VERSION "4.8.5"

void printUsage(FILE* out)
{
    fprintf(out, "Usage: qmake -query QT_VERSION\n\n");
    fprintf(out, "This is a minimal QMake stub so that CMake may find "
        "the Qt distribution bundled with the Mitsuba dependencies.\n");
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, "***Unknown options\n");
        printUsage(stderr);
        return 1;
    } else if (strcmp(argv[1], "-query") != 0) {
        fprintf(stderr, "***Unknown option %s\n", argv[1]);
        printUsage(stderr);
        return 1;
    }
    
    if (strcmp(argv[2], "QT_VERSION") == 0) {
        puts(QT_VERSION);
        return 0;
    }
    else {
        puts("**Unknown**");
        return 101;
    }
}

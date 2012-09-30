/* Minimal QMake stub to trick CMake into using the Mitsuba Qt distribution */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <libgen.h>
#include <mach-o/dyld.h>

/* Qt version currently provided in the Mitsuba dependencies */
#define QT_VERSION "4.8.2"

void printUsage(FILE* out)
{
    fprintf(out, "Usage: qmake -query <QT_VERSION|QT_INSTALL_HEADERS>\n\n");
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

    /* Derive the dependencies name from the executable */
    char tmpbuf1[PATH_MAX + 1];
    char tmpbuf2[PATH_MAX + 1];
    int32_t bufsize = sizeof(tmpbuf1);
    _NSGetExecutablePath(tmpbuf1, &bufsize);
    const char *exe_dir = dirname(tmpbuf1);
    strcpy(tmpbuf2, exe_dir);
    strcat(tmpbuf2, "/../");
    const char* dependencies_dir = realpath(tmpbuf2, tmpbuf1);
    
    if (strcmp(argv[2], "QT_VERSION") == 0) {
        puts(QT_VERSION);
    }
    else if (strcmp(argv[2], "QT_INSTALL_HEADERS") == 0) {
        printf("%s/include\n", dependencies_dir);
    }
    else {
        puts("**Unknown**");
        return 101;
    }
    return 0;
}

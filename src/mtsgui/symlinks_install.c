#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/param.h>

void installPython(const char *basedir) {
	const char *fname = "/Library/Python/2.6/site-packages/mitsuba.pth";

	FILE *f = fopen(fname, "w");
	if (!f) {
		fprintf(stderr, "Unable to write to file \"%s\"!\n", fname);
		exit(-1);
	}

	if (fprintf(f, "import sys; sys.path.append(\"%s/python\")\n", basedir) < 1) {
		fprintf(stderr, "Unexpected I/O error while "
			"writing to \"%s\"!\n", fname);
		exit(-1);
	}

	fclose(f);
}

void install(const char *basedir, const char *name) {
	char fname[MAXPATHLEN];
	FILE *f;

	snprintf(fname, sizeof(fname), "/usr/bin/%s", name);

	f = fopen(fname, "w");
	if (!f) {
		fprintf(stderr, "Unable to write to file \"%s\"!\n", fname);
		exit(-1);
	}

	if (fprintf(f, "%s/Contents/MacOS/%s $@\n", basedir, name) < 1) {
		fprintf(stderr, "Unexpected I/O error while "
			"writing to \"%s\"!\n", fname);
		exit(-1);
	}

	fclose(f);

	if (chmod(fname, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0) {
		fprintf(stderr, "Unexpected I/O error while setting "
			"the permssions of \"%s\"!\n", fname);
		exit(-1);
	}
}

int main(int argc, char **argv) {
	if (argc != 2) {
		fprintf(stderr, "Incorrect number of arguments!\n");
		return -1;
	}
		
	if (setuid(0) != 0) {
		fprintf(stderr, "setuid(): failed!\n");
		return -1;
	}

	install(argv[1], "mitsuba");
	install(argv[1], "mtsgui");
	install(argv[1], "mtssrv");
	install(argv[1], "mtsutil");
	install(argv[1], "mtsimport");
	installPython(argv[1]);

	return 0;
}


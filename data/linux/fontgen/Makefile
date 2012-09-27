CPPFLAGS = -I ../../../include -I /usr/include/freetype2 -DSINGLE_PRECISION -DSPECTRUM_SAMPLES=3
LDFLAGS = -lfreetype -L../../../dist -lmitsuba-core -lmitsuba-hw -lboost_filesystem

all: bin2c fontgen

bin2c: bin2c.c

fontgen: fontgen.cpp 

clean:
	rm -f fontgen

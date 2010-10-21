CPPFLAGS = -I ../../../include -I /usr/include/freetype2
LDFLAGS = -lfreetype -L ../../../src/libcore -L ../../../src/libhw \
		  -lmitsuba-core -lmitsuba-hw

all: bin2c fontgen

bin2c: bin2c.c

fontgen: fontgen.cpp 

clean:
	rm -f fontgen

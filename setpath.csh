#!/bin/tcsh

set called=($_)
if ("$called" != "") then
    set reldir=`dirname $called[2]`
else if ("$0" != "tcsh") then
    set reldir=`dirname 0`
else
	echo "Unable to detect path!"
	exit 1
endif
set MITSUBA_DIR=`cd $reldir && pwd`

if ("`uname`" == "Darwin") then
	setenv PATH "$MITSUBA_DIR/Mitsuba.app/Contents/MacOS:$PATH"
else
	if (! ($?LD_LIBRARY_PATH) ) then
		setenv LD_LIBRARY_PATH "$MITSUBA_DIR/dist"
	else
		setenv LD_LIBRARY_PATH "$MITSUBA_DIR/dist:$LD_LIBRARY_PATH"
	endif
	setenv PATH "$MITSUBA_DIR/dist:$PATH"
	# Generate core dumps if something goes wrong
	limit coredumpsize 1000000000
endif

unset reldir

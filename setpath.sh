#!/bin/bash
unamestr=`uname`
if [[ "$unamestr" == 'Darwin' ]]; then
	export PATH=$PATH:`pwd`/Mitsuba.app/Contents/MacOS
else
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`/src/libcore:`pwd`/src/librender:`pwd`/src/libhw
	export PATH=$PATH:`pwd`
fi

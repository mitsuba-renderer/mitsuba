#!/bin/bash
unamestr=`uname`
if [[ "$unamestr" == 'Darwin' ]]; then
	export PATH=`pwd`/Mitsuba.app/Contents/MacOS:$PATH
else
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`/src/libcore:`pwd`/src/librender:`pwd`/src/libhw
	export PATH=`pwd`:$PATH
	ulimit -c 1000000000
fi

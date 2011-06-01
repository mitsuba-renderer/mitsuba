#!/bin/bash
unamestr=`uname`
if [[ "$unamestr" == 'Darwin' ]]; then
	export PATH=`pwd`/Mitsuba.app/Contents/MacOS:$PATH
else
	export LD_LIBRARY_PATH=`pwd`/dist
	export PATH=`pwd`/dist:$PATH
	ulimit -c 1000000000
fi

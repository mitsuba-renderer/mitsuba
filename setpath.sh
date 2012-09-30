#!/bin/bash

UNAME=`uname`
DIR=$(dirname "$BASH_SOURCE")
ABSDIR=$(cd "$DIR"; pwd)

if [[ "$UNAME" == 'Darwin' ]]; then
	export PATH="$ABSDIR/Mitsuba.app/Contents/MacOS":$PATH
else
	export LD_LIBRARY_PATH="$ABSDIR/dist":$LD_LIBRARY_PATH
	export PATH="$ABSDIR/dist":$PATH

	# Generate core dumps if something goes wrong
	ulimit -c 1000000000
fi

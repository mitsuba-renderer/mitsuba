#!/bin/sh
export SCONS_LIB_DIR=`pwd`/tools/scons-1.2.0/engine
./tools/scons-1.2.0/script/scons $@

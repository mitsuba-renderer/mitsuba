#!/bin/bash
APP=`pwd`
PY_DSO=/Library/Python/2.6/site-packages/mitsuba.so

install() {
	echo $APP/Contents/MacOS/$1 \$@ > /usr/bin/$1
	chmod +x /usr/bin/$1
}
fix_import() {
	install_name_tool -change @loader_path/../Contents/Frameworks/$1.dylib $APP/Contents/Frameworks/$1.dylib $PY_DSO
}

install mitsuba
install mtsgui
install mtssrv
install mtsimport
install mtsutil

if [ -e $APP/python/mitsuba.so ]
then
	cp $APP/python/mitsuba.so $PY_DSO 
	fix_import libboost_python
	fix_import libboost_system
	fix_import libboost_filesystem
	fix_import libmitsuba-core
	fix_import libmitsuba-render
	fix_import libmitsuba-hw
	fix_import libmitsuba-bidir
	fix_import libiomp5
fi


#!/bin/bash

tar --exclude-vcs -czvf mitsuba-$1.tar.gz ChangeLog README SConstruct \
	config include/ schema/ src/bsdfs/ src/cameras/ src/collada/ src/films/ src/libcore/ src/libhw/ \
	src/librender/ src/luminaires/ src/medium/ src/mitsuba/ src/phase/ src/qtgui/ src/rfilters/ src/samplers/ \
	src/shapes/ src/subsurface/ src/textures/ src/utils/ src/volume/ src/integrators/direct src/integrators/misc \
	src/integrators/path src/integrators/photonmapper src/integrators/vpl tools/boost tools/darwin \
	tools/build.sh tools/scons-1.2.0 tools/windows tools/linux/*.sh tools/linux/mitsuba.desktop tools/qt4.py \
	tools/build.bat doc/license.txt setpath.sh

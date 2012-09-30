rm -Rf symbols32 symbols64
dependencies/bin/symbolstore.py -a "i386" dependencies/bin/dump_syms symbols32 \
	`find Mitsuba.app -type f -and \( -name *.dylib -or -name mtsgui -or -name mitsuba -or -name mtsutil \)` \
	dependencies/frameworks/QtCore.framework/QtCore dependencies/frameworks/QtGui.framework/QtGui dependencies/frameworks/QtNetwork.framework/QtNetwork \
	dependencies/frameworks/QtXml.framework/QtXml dependencies/frameworks/QtXmlPatterns.framework/QtXmlPatterns dependencies/frameworks/QtOpenGL.framework/QtOpenGL

dependencies/bin/symbolstore.py -a "x86_64" dependencies/bin/dump_syms symbols64 \
	`find Mitsuba.app -type f -and \( -name *.dylib -or -name mtsgui -or -name mitsuba -or -name mtsutil \)` \
	dependencies/frameworks/QtCore.framework/QtCore dependencies/frameworks/QtGui.framework/QtGui dependencies/frameworks/QtNetwork.framework/QtNetwork \
	dependencies/frameworks/QtXml.framework/QtXml dependencies/frameworks/QtXmlPatterns.framework/QtXmlPatterns dependencies/frameworks/QtOpenGL.framework/QtOpenGL

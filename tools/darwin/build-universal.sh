#/bin/bash
mkdir i386 x86_64 
cp config/config-darwin-x86.py config.py
rm .sconsign*
tools/build.sh -j 3
cp -R Mitsuba.app i386/Mitsuba.app
cp config/config-darwin-x86_64.py config.py
rm .sconsign*
tools/build.sh -j 3
cp -R Mitsuba.app x86_64/Mitsuba.app

rm Mitsuba.app/Contents/MacOS/*
rm Mitsuba.app/Contents/Frameworks/librender.dylib Mitsuba.app/Contents/Frameworks/libhw.dylib Mitsuba.app/Contents/Frameworks/libcore.dylib
lipo -create -output Mitsuba.app/Contents/MacOS/mitsuba i386/Mitsuba.app/Contents/MacOS/mitsuba x86_64/Mitsuba.app/Contents/MacOS/mitsuba
lipo -create -output Mitsuba.app/Contents/MacOS/collada i386/Mitsuba.app/Contents/MacOS/collada x86_64/Mitsuba.app/Contents/MacOS/collada
lipo -create -output Mitsuba.app/Contents/MacOS/qtgui i386/Mitsuba.app/Contents/MacOS/qtgui x86_64/Mitsuba.app/Contents/MacOS/qtgui
lipo -create -output Mitsuba.app/Contents/MacOS/mtssrv i386/Mitsuba.app/Contents/MacOS/mtssrv x86_64/Mitsuba.app/Contents/MacOS/mtssrv
lipo -create -output Mitsuba.app/Contents/Frameworks/libcore.dylib i386/Mitsuba.app/Contents/Frameworks/libcore.dylib x86_64/Mitsuba.app/Contents/Frameworks/libcore.dylib
lipo -create -output Mitsuba.app/Contents/Frameworks/libhw.dylib i386/Mitsuba.app/Contents/Frameworks/libhw.dylib x86_64/Mitsuba.app/Contents/Frameworks/libhw.dylib
lipo -create -output Mitsuba.app/Contents/Frameworks/librender.dylib i386/Mitsuba.app/Contents/Frameworks/librender.dylib x86_64/Mitsuba.app/Contents/Frameworks/librender.dylib
for i in Mitsuba.app/plugins/*.dylib; do rm $i; lipo -create -output $i i386/$i x86_64/$i; done

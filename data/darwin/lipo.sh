rm -Rf Mitsuba.app
cp -R Mitsuba64.app/ Mitsuba.app
for i in Mitsuba32.app/plugins/*; do
	base=$(basename $i)
	echo "Running lipo on $base"
	lipo -create Mitsuba32.app/plugins/$base Mitsuba64.app/plugins/$base -output Mitsuba.app/plugins/$base
done
for i in Mitsuba32.app/Contents/MacOS/*; do
	base=$(basename $i)
	echo "Running lipo on $base"
	lipo -create Mitsuba32.app/Contents/MacOS/$base Mitsuba64.app/Contents/MacOS/$base -output Mitsuba.app/Contents/MacOS/$base
done
for i in Mitsuba32.app/Contents/Frameworks/libmitsuba*; do
	base=$(basename $i)
	echo "Running lipo on $base"
	lipo -create Mitsuba32.app/Contents/Frameworks/$base Mitsuba64.app/Contents/Frameworks/$base -output Mitsuba.app/Contents/Frameworks/$base
done
lipo -create Mitsuba32.app/python/2.6/mitsuba.so Mitsuba64.app/python/2.6/mitsuba.so -output Mitsuba.app/python/2.6/mitsuba.so
lipo -create Mitsuba32.app/python/2.7/mitsuba.so Mitsuba64.app/python/2.7/mitsuba.so -output Mitsuba.app/python/2.7/mitsuba.so

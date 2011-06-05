@echo off
xcopy dist "Mitsuba %1" /e /i /h
del "Mitsuba %1.zip"
7z u -tzip "Mitsuba %1 %2bit.zip" "Mitsuba %1"
rmdir /s /q "Mitsuba %1"

#!/bin/bash
rm -f Mitsuba.app/*.log
dmgcanvas data/darwin/Mitsuba.dmgCanvas Mitsuba\ $1.dmg -v "Mitsuba $1" -leopard-compatible yes
#echo $1 > /tmp/version
#scp /tmp/version ChangeLog wazlaf@mitsuba-renderer.org:/home/httpd/mitsuba-renderer.org/htdocs
#scp Mitsuba\ $1.dmg wazlaf@mitsuba-renderer.org:/home/httpd/mitsuba-renderer.org/htdocs/releases

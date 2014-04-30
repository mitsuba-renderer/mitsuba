#!/usr/bin/env python

import subprocess, sys, re

if len(sys.argv) == 1:
    print('rendertime.py: Simple utility to extract the rendering time of an EXR file created using Mitsuba')
    print('Syntax: rendertime.py <one or more EXR files>')

for arg in sys.argv[1:]:
    sys.stdout.write("%s: " % arg)

    p = subprocess.Popen(["exrheader", arg], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    matches = re.findall('Render time: ([^\n]*)', p.communicate()[0])
    if len(matches) != 1:
        sys.stdout.write('error reading metadata!\n')
        continue

    value = matches[0]

    units = [
        ['d', 60*60*24],
        ['h', 60*60],
        ['m', 60],
        ['s', 1],
        ['ms', 0.001]
    ];

    for unit in units:
        if value[-len(unit[0]):] == unit[0]:
            value = float(value[0:-len(unit[0])]) * unit[1]
            break

    force = False
    for unit in units:
        if value > unit[1] or force:
            amount = int(value / unit[1])
            value -= amount * unit[1]
            force = True
            sys.stdout.write('%i%s ' % (amount, unit[0]))
    sys.stdout.write('\n')

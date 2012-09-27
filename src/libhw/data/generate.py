#!/usr/bin/python
# This script turns the available GLSL shaders into a header file
# that can be included in the libhw binary.

import sys

def append(fname):
	f = open(fname)
	sys.stdout.write('static const char *' + fname.replace('.', '_') + ' MAYBE_UNUSED = ')
	lineIndex = 0
	for line in f:
		lineIndex = lineIndex + 1
		if lineIndex < 19:
			continue
		sys.stdout.write('\n\t\"' + line[:-1] + '\\n\"')
	sys.stdout.write(';\n\n')
	f.close()

print("/*")
print("    This file is part of Mitsuba, a physically based rendering system.")
print("")
print("    Copyright (c) 20072011 by Wenzel Jakob and others.")
print("")
print("    Mitsuba is free software; you can redistribute it and/or modify")
print("    it under the terms of the GNU General Public License Version 3")
print("    as published by the Free Software Foundation.")
print("")
print("    Mitsuba is distributed in the hope that it will be useful,")
print("    but WITHOUT ANY WARRANTY; without even the implied warranty of")
print("    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the")
print("    GNU General Public License for more details.")
print("")
print("    You should have received a copy of the GNU General Public License")
print("    along with this program. If not, see <http://www.gnu.org/licenses/>.")
print("*/")
print("")

print("#ifdef __GNUC__")
print("\t#define MAYBE_UNUSED __attribute__((used))")
print("#else")
print("\t#define MAYBE_UNUSED")
print("#endif")
print("")

append('sh_paraboloid.vert')
append('sh_paraboloid.geom')
append('sh_paraboloid.frag')
append('sh_directional.vert')
append('sh_directional.frag')
append('sh_cube_6pass.vert')
append('sh_cube_6pass.frag')
append('sh_cube_1pass.vert')
append('sh_cube_1pass.geom')
append('sh_cube_1pass.frag')
append('sh_hemicube_1pass.vert')
append('sh_hemicube_1pass.geom')
append('sh_hemicube_1pass.frag')
append('sh_background.vert')
append('sh_background.frag')
append('sh_unsupported.vert')
append('sh_unsupported.frag')
append('sh_render.vert')
append('sh_render.geom')
append('sh_render.frag')

/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#extension GL_EXT_geometry_shader4 : enable

/* Transformation matrix for each cube map face */
uniform mat4 transform[6];

/* Depth projection axis */
uniform vec4 projDir[6];

/* -> Fragment shader */
varying float depth;

void main() {
	depth = 0.0; // avoid an (incorrect?) compiler warning

	/* Replicate the geometry six times and rasterize to each cube map layer */
	for (int side = 0; side < 6; side++) {
		gl_Layer = side;
		for (int i = 0; i < gl_VerticesIn; i++) {
			gl_Position = transform[side] * gl_PositionIn[i];
			depth = dot(projDir[side], gl_PositionIn[i]);
			EmitVertex();
		}
		EndPrimitive();
	}
}

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

/* From vertex program */
varying in vec3 posInVPLSpace_vertex[3];
varying in vec3 posInWorldSpace_vertex[3];
varying in vec2 uv_vertex[3];

/* To fragment program */
varying out vec3 posInVPLSpace;
varying out vec3 posInWorldSpace;
varying out vec3 normal;
varying out vec2 uv;

#ifdef VERTEX_COLORS
    varying in vec3 vertexColor_vertex[3];
    varying out vec3 vertexColor;
#endif

void main() {
    vec3 edge1 = posInWorldSpace_vertex[0]-posInWorldSpace_vertex[1];
    vec3 edge2 = posInWorldSpace_vertex[0]-posInWorldSpace_vertex[2];

    normal = cross(edge1, edge2);
    for (int i=0; i<gl_VerticesIn; ++i) {
        gl_Position = gl_PositionIn[i];
        posInWorldSpace = posInWorldSpace_vertex[i];
        posInVPLSpace = posInVPLSpace_vertex[i];
        uv = uv_vertex[i];

        #ifdef VERTEX_COLORS
            vertexColor = vertexColor_vertex[i];
        #endif

        EmitVertex();
    }
    EndPrimitive();
}

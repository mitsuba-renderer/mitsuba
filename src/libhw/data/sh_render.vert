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

/* Uniform parameters */
uniform mat4 vplTransform;
uniform mat4 instanceTransform;

#ifndef FACE_NORMALS
    /* -> Fragment program */
    varying vec3 posInVPLSpace;
    varying vec3 posInWorldSpace;
    varying vec2 uv;
    varying vec3 normal;
#else
    /* -> Geometry program */
    varying vec3 posInVPLSpace_vertex;
    varying vec3 posInWorldSpace_vertex;
    varying vec2 uv_vertex;
#endif

#ifdef ANISOTROPIC
    varying vec3 tangent;
#endif

#ifdef VERTEX_COLORS
    #ifndef FACE_NORMALS
        varying vec3 vertexColor;
    #else
        varying vec3 vertexColor_vertex;
    #endif
#endif

void main() {
    vec4 pos = instanceTransform * gl_Vertex;
    gl_Position = gl_ModelViewProjectionMatrix * pos;

    #ifndef FACE_NORMALS
        posInWorldSpace = pos.xyz;
        posInVPLSpace   = (vplTransform * pos).xyz;
        uv = gl_MultiTexCoord0.xy;

        /* Multiply by instanceTransform (only rigid transformations allowed) */
        normal = (instanceTransform * vec4(gl_Normal, 0.0)).xyz;

        #ifdef VERTEX_COLORS
            vertexColor = gl_Color.rgb;
        #endif
    #else
        posInWorldSpace_vertex = pos.xyz;
        posInVPLSpace_vertex   = (vplTransform * pos).xyz;
        uv_vertex = gl_MultiTexCoord0.xy;

        #ifdef VERTEX_COLORS
            vertexColor_vertex = gl_Color.rgb;
        #endif
    #endif

    #ifdef ANISOTROPIC
        tangent = gl_MultiTexCoord1.xyz;
    #endif
}

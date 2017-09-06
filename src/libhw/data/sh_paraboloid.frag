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

/* Triangle data (from geometry shader) */
varying vec3 p0, edge1, edge2;

/* Projection space position */
varying vec2 pos;

/* Depth range */
uniform float minDepth;
uniform float invDepthRange;

void main() {
    float r2 = dot(pos, pos);
    if (r2 > 1.0)
        discard;

    /* Turn into a direction */
    float cosTheta = (1.0-r2) / (1.0+r2);
    float factor = sqrt((1.0-cosTheta*cosTheta) / r2);
    vec3 d = vec3(pos * factor, cosTheta);

    vec3 pvec = cross(d, edge2);
    float det = dot(edge1, pvec);
    if (det == 0.0)
        discard;

    float inv_det = 1.0 / det;
    float u = -dot(p0, pvec) * inv_det;

    if (u < 0.0 || u > 1.0)
        discard;

    vec3 qvec = cross(p0, edge1);
    float v = -dot(d, qvec) * inv_det;
    if (v < 0.0 || u+v > 1.0)
        discard;

    float t = -dot(edge2, qvec) * inv_det - minDepth;
    if (!(t > 0.0)) // catch NaNs as well
        discard;

    float depth = t * invDepthRange;

    float dx = dFdx(depth), dy = dFdy(depth);
    gl_FragDepth = depth + sqrt(dx*dx+dy*dy);
}

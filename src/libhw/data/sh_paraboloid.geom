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

/* -> Fragment shader (triangle data, projected position) */
varying vec3 p0, edge1, edge2;
varying vec2 pos;

vec2 project(vec3 p) {
	return p.xy / (length(p) + p.z);
}

float cosAlpha(vec3 v1, vec3 v2) {
	vec3 v3 = cross(v1, v2);
	return v3.z / length(v3);
}

void main() {
	vec3 cur = gl_PositionIn[0].xyz;
	cur.z = -cur.z;

	bool curIsInside = cur.z > 0.0;
	vec3 vertices[4];
	int vCount = 0;

	/* Clip the geometry against the forward hemisphere */
	for (int i=0; i<3; ++i) {
		vec3 next = gl_PositionIn[i < 2 ? i+1 : 0].xyz;
		next.z = -next.z;
		bool nextIsInside = next.z > 0.0;

		if (curIsInside && nextIsInside) {
			vertices[vCount++] = next;
		} else if (curIsInside != nextIsInside) {
			float t = -cur.z / (next.z-cur.z);

			vertices[vCount++] = 
				vec3((1-t)*cur.xy + t*next.xy, 0.0);

			if (nextIsInside)
				vertices[vCount++] = next;
		}

		cur = next;
		curIsInside = nextIsInside;
	}

	if (vCount < 3)
		return;

	/* Transform into projection space */
	vec2 projected[4];
	for (int i=0; i<vCount; ++i)
		projected[i] = project(vertices[i]);

	/* Find the longest edge and construct an OBB rotation matrix */
	float largest = 0.0;
	mat2 rot;
	rot[0] = vec2(0.0);

	for (int i=0; i<vCount; i++) {
		int next = (i < vCount-1) ? i+1 : 0;
		vec2 p0 = projected[i],
		     p1 = projected[next],
		     d  = p1-p0;

		float l2 = dot(d, d);
		if (l2 > largest) {
			largest = l2;
			rot[0] = d;
		}
	}
	rot[0] /= sqrt(largest);
	rot[1]  = vec2(rot[0].y, -rot[0].x);

	/* Find a bounding rectangle for each edge and expand the OBB */
	vec2 bmax = vec2(-1.0, -1.0), bmin = vec2(1.0, 1.0);
	for (int i=0; i<vCount; i++) {
		int next = (i < vCount-1) ? i+1 : 0;
		vec2 a  = projected[i],
		     b  = projected[next],
		     p0 = a * rot,
		     p1 = b * rot,
			 p2 = p0, p3 = p1;

		float ca = cosAlpha(vertices[i], vertices[next]);
		if (ca != 0.0) {
			vec2 d = vec2(b.y-a.y, a.x-b.x) * rot;
			float l = length(d), ica = 1.0 / ca;
			d *= (ica - sign(ica) * sqrt(ica*ica - 0.25 * l*l)) / l;
			p2 += d; p3 += d;
		}

		bmax = max(bmax, max(max(p0, p1), max(p2, p3)));
		bmin = min(bmin, min(min(p0, p1), min(p2, p3)));
	}

	/* Plumb triangle data to the fragment shader */
	p0 = gl_PositionIn[0].xyz;
	edge1 = gl_PositionIn[1].xyz - gl_PositionIn[0].xyz;
	edge2 = gl_PositionIn[2].xyz - gl_PositionIn[0].xyz;
	p0.z = -p0.z; edge1.z = -edge1.z; edge2.z = -edge2.z;

	/* Emit a quad */
	pos = rot[0] * bmin.x + rot[1] * bmin.y;
	gl_Position = vec4(pos, 0.5, 1.0);
	EmitVertex();
	pos = rot[0] * bmin.x + rot[1] * bmax.y;
	gl_Position = vec4(pos, 0.5, 1.0);
	EmitVertex();
	pos = rot[0] * bmax.x + rot[1] * bmin.y;
	gl_Position = vec4(pos, 0.5, 1.0);
	EmitVertex();
	pos = rot[0] * bmax.x + rot[1] * bmax.y;
	gl_Position = vec4(pos, 0.5, 1.0);
	EmitVertex();
	EndPrimitive();
}

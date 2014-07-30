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

#extension GL_EXT_gpu_shader4 : enable

#define EPSILON 0.001

/* Some helper functions for BSDF implementations */
float cosTheta(vec3 v) { return v.z; }
float sinTheta2(vec3 v) { return 1.0-v.z*v.z; }
float sinTheta(vec3 v) { return sqrt(max(0.0, sinTheta2(v))); }
float tanTheta(vec3 v) { return sinTheta(v)/cosTheta(v); }
float sinPhi(vec3 v) { return v.y/sinTheta(v); }
float cosPhi(vec3 v) { return v.x/sinTheta(v); }
const float pi = 3.141592653589;
const float inv_pi = 0.318309886183791;
const float inv_twopi = 0.159154943091895;
const float inv_fourpi = 0.0795774715459477;

/* From the vertex program */
varying vec3 posInVPLSpace;
varying vec3 posInWorldSpace;
varying vec3 normal;
varying vec2 uv;

#ifdef ANISOTROPIC
	varying vec3 tangent;
#endif

#ifdef VERTEX_COLORS
	varying vec3 vertexColor;
#endif

/* Uniform parameters */
uniform mat3 vplFrame;
uniform vec3 vplPower;
uniform vec2 vplUV;
uniform vec3 vplWi;
uniform float minDistSqr;
uniform float emitterScale;

uniform mat4 vplTransform;

#ifdef CUBEMAP_VPL
	uniform samplerCube shadowMap;
#else
	uniform sampler2D shadowMap;
#endif

#ifdef DIRECTIONAL_VPL
	uniform vec3 vplDirection;
#else
	uniform vec3 vplPosition;
	uniform vec2 depthRange;
#endif

#ifdef DIRECTIONAL_CAMERA
	uniform vec3 camDirection;
#else
	uniform vec3 camPosition;
#endif

#if defined(DIRECTIONAL_VPL)
	bool isShadowed(vec3 p) {
		p = p * 0.5 + 0.5;
		return texture2D(shadowMap, p.xy).r * (1 + EPSILON) < p.z;
	}
#elif defined(PARABOLOIDAL_VPL)
	bool isShadowed(vec3 p) {
		float depth = texture2D(shadowMap, (p.xy / (-p.z + length(p))) * 0.5 + 0.5).r;
		depth = (depth * (depthRange[1]-depthRange[0]) + depthRange[0]) * (1.0+EPSILON);
		return dot(p, p) > depth * depth && p.z < 0.0;
	}
#elif defined(CUBEMAP_VPL)
	bool isShadowed(vec3 d) {
		float depth = textureCube(shadowMap, vec3(d.x, -d.y, d.z)).r;
		float ref_depth = max(max(abs(d.x), abs(d.y)), abs(d.z));
		depth = (depth * (depthRange[1]-depthRange[0]) + depthRange[0]) * (1.0+EPSILON);
		return depth < ref_depth;
	}
#endif

{{ SUPPLEMENTAL CODE }}

void main() {
	/* Set up an ONB */
	vec3 N = normalize(normal);
	mat3 frame;

	#ifdef ANISOTROPIC
		/* Use the per-vertex tangent information to construct a frame */

		frame[0] = normalize(tangent - dot(tangent, N)*N);
	#else
		/* The material is isotropic -- any frame will do */

		if (abs(N.x) > abs(N.y)) {
			float invLen = 1.0 / sqrt(N.x*N.x + N.z*N.z);
			frame[0] = vec3(-N.z * invLen, 0.0, N.x * invLen);
		} else {
			float invLen = 1.0 / sqrt(N.y*N.y + N.z*N.z);
			frame[0] = vec3(0.0, -N.z * invLen, N.y * invLen);
		}
	#endif

	frame[1] = cross(N, frame[0]);
	frame[2] = N;

	/* Compute the incident direction in local coordinates (at the point being rendered) */
	vec3 wiWorld;
	#ifdef DIRECTIONAL_CAMERA
		wiWorld = -camDirection;
	#else
		wiWorld = normalize(camPosition - posInWorldSpace);
	#endif
	vec3 wi = wiWorld * frame;

	vec3 emission;

	#if defined(EMITTER_AREA_EVAL_NAME) && defined(EMITTER_DIR_EVAL_NAME)
		emission = EMITTER_AREA_EVAL_NAME(uv) *
					EMITTER_DIR_EVAL_NAME(wi) * emitterScale;
	#else
		emission = vec3(0.0);
	#endif

	if (isShadowed(posInVPLSpace)) {
		gl_FragColor = vec4(emission, 1.0);
		return;
	}

	/* Compute the outgoing direction in local coordinates (at the point being rendered) */
	vec3 woWorld;

	#ifdef DIRECTIONAL_VPL
		woWorld = -vplDirection;
	#else
		woWorld = vplPosition - posInWorldSpace;
		float distSqr = dot(woWorld, woWorld);
		woWorld /= sqrt(distSqr);
	#endif
	vec3 wo = woWorld * frame;

	/* Compute the outgoing direction in local coordinates (at the VPL) */
	vec3 vplWo = -(woWorld * vplFrame); /* The parentheses are required or incorrect code is generated on OSX .. */

	vec3 result = vplPower;

	#ifdef EMITTER_VPL
		#ifdef VPL_ON_SURFACE
			result *= VPL_EVAL_NAME(vplWo) * cosTheta(vplWo);
		#else
			result *= VPL_EVAL_NAME(vplWo);
		#endif
	#else
		result *= VPL_EVAL_NAME(vplUV, vplWi, vplWo);
	#endif

	result *= BSDF_EVAL_NAME(uv, wi, wo);

	#ifndef DIRECTIONAL_VPL
		result *= (1.0 / max(distSqr, minDistSqr));
	#endif

	gl_FragColor = vec4(result + emission, 1.0);
}

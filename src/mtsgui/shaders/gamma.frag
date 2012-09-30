uniform sampler2D colorSource, depthSource;
uniform float invWhitePoint, invGamma;
uniform bool sRGB, hasDepth;

float toSRGB(float value) {
	if (value < 0.0031308)
		return 12.92 * value;
	return 1.055 * pow(value, 1.0/2.4) - 0.055;
}

void main() {
	vec4 color = texture2D(colorSource, gl_TexCoord[0].xy) * invWhitePoint;
	if (sRGB)
		gl_FragColor = vec4(toSRGB(color.r), toSRGB(color.g), toSRGB(color.b), 1);
	else
		gl_FragColor = vec4(pow(color.rgb, vec3(invGamma)), 1);
	gl_FragDepth = hasDepth ? texture2D(depthSource, gl_TexCoord[0].xy).r : 0.5;
}

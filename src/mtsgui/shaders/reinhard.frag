uniform sampler2D colorSource, depthSource;
uniform float scale, invWp2, invGamma, multiplier;
uniform bool sRGB, hasDepth;

float toSRGB(float value) {
	if (value < 0.0031308)
		return 12.92 * value;
	return 1.055 * pow(value, 1.0 / 2.4) - 0.055;
}

void main() {
	const mat3 RGB2XYZ = mat3(0.412453, 0.212671, 0.019334,
	                          0.357580, 0.715160, 0.119193,
	                          0.180423, 0.072169, 0.950227);

	const mat3 XYZ2RGB = mat3(3.240479, -0.969256,  0.055648,
	                         -1.537150,  1.875991, -0.204043,
	                         -0.498535,  0.041556,  1.057311);

	vec3 XYZ = RGB2XYZ * texture2D(colorSource, gl_TexCoord[0].xy).rgb * multiplier;
	float normalization = 1.0 / (XYZ.x + XYZ.y + XYZ.z);
	float x = XYZ.x*normalization;
	float y = XYZ.y*normalization;
	float Lp = XYZ.y*scale;
	XYZ.y = Lp * (1.0 + Lp*invWp2) / (1.0 + Lp);
	float ratio = XYZ.y / y;
	XYZ.x = ratio * x;
	XYZ.z = ratio * (1.0 - x - y);
	vec3 color = XYZ2RGB * XYZ;
	if (sRGB)
		gl_FragColor = vec4(toSRGB(color.r), toSRGB(color.g), toSRGB(color.b), 1);
	else
		gl_FragColor = vec4(pow(color.rgb, vec3(invGamma)), 1);
	gl_FragDepth = hasDepth ? texture2D(depthSource, gl_TexCoord[0].xy).r : 0.5;
}

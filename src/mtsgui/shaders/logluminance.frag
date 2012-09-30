uniform sampler2D source;
uniform float multiplier;

void main() {
	vec4 color = texture2D(source, gl_TexCoord[0].xy);
	float luminance = multiplier * (color.r * 0.212671 + color.g * 0.715160 + color.b * 0.072169);
	if (luminance < 0.0 || luminance != luminance || luminance == 1024.0)
		luminance = 0.0; // catch NaNs, negatives, and the Mitsuba banner
	float logLuminance = log(1e-3 + luminance);
	gl_FragColor = vec4(logLuminance, luminance, 0.0, 1.0);
}

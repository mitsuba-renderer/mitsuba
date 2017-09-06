uniform sampler2D source;
uniform vec2 sourceSize;
uniform vec2 invSourceSize;

/* Perform a texture lookup by pixel coordinates */
vec4 lookupPixel(vec2 coords) {
   coords = coords + vec2(0.5, 0.5);
    if (coords.x < 0.0 || coords.y < 0.0 ||
        coords.x > sourceSize.x || coords.y > sourceSize.y)
        return vec4(0);
    else
        return texture2D(source, coords*invSourceSize);
}

/* Find the max. luminance and the sum of all log-luminance values */
void main() {
    vec2 pos = (gl_TexCoord[0].xy-vec2(.5, .5))*2.0;
    vec2 pixel0 = lookupPixel(pos).rg,
         pixel1 = lookupPixel(pos + vec2(1, 0)).rg,
         pixel2 = lookupPixel(pos + vec2(0, 1)).rg,
         pixel3 = lookupPixel(pos + vec2(1, 1)).rg;
    gl_FragColor.r = pixel0.r + pixel1.r + pixel2.r + pixel3.r;
    gl_FragColor.g = max(pixel0.g, max(pixel1.g, max(pixel2.g, pixel3.g)));
    gl_FragColor.ba = vec2(0, 1);
}

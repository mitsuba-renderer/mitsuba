uniform vec2 targetSize;
void main() {
    gl_Position = ftransform();
    gl_TexCoord[0].xy = vec2(gl_MultiTexCoord0.x * targetSize.x,
                             gl_MultiTexCoord0.y * targetSize.y);
}

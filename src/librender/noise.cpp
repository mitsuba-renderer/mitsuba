#include <mitsuba/render/noise.h>

MTS_NAMESPACE_BEGIN

#define GRAD_PERLIN 1

// Based on Ken Perlin's improved noise reference implementation
#define NOISE_PERM_SIZE 256
static int NoisePerm[2 * NOISE_PERM_SIZE] = {
    151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140,
    36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190,  6, 148, 247, 120,
    234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
    88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168,  68, 175, 74, 165, 71,
    134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133,
    230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54,  65, 25, 63,
    161, 1, 216, 80, 73, 209, 76, 132, 187, 208,  89, 18, 169, 200, 196, 135, 130,
    116, 188, 159, 86, 164, 100, 109, 198, 173, 186,  3, 64, 52, 217, 226, 250,
    124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47,
    16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119, 248, 152,  2, 44, 154,
    163,  70, 221, 153, 101, 155, 167,  43, 172, 9, 129, 22, 39, 253,  19, 98,
    108, 110, 79, 113, 224, 232, 178, 185,  112, 104, 218, 246, 97, 228, 251, 34,
    242, 193, 238, 210, 144, 12, 191, 179, 162, 241,  81, 51, 145, 235, 249, 14,
    239, 107, 49, 192, 214,  31, 181, 199, 106, 157, 184,  84, 204, 176, 115, 121,
    50, 45, 127,  4, 150, 254, 138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243,
    141, 128, 195, 78, 66, 215, 61, 156, 180, 151, 160, 137, 91, 90, 15, 131, 13,
    201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240,
    21, 10, 23, 190,  6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219,
    203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136,
    171, 168,  68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83,
    111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102,
    143, 54,  65, 25, 63, 161,  1, 216, 80, 73, 209, 76, 132, 187, 208,  89, 18,
    169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186,  3,
    64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212,
    207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119,
    248, 152,  2, 44, 154, 163,  70, 221, 153, 101, 155, 167,  43, 172, 9, 129,
    22, 39, 253,  19, 98, 108, 110, 79, 113, 224, 232, 178, 185,  112, 104, 218,
    246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241,  81,
    51, 145, 235, 249, 14, 239, 107, 49, 192, 214,  31, 181, 199, 106, 157, 184,
    84, 204, 176, 115, 121, 50, 45, 127,  4, 150, 254, 138, 236, 205, 93, 222,
    114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
};

inline static Float grad(int x, int y, int z, Float dx, Float dy, Float dz) {
    int h = NoisePerm[NoisePerm[NoisePerm[x]+y]+z];
    h &= 15;
#if defined(GRAD_PERLIN)
    /* Based on Ken Perlin's improved Noise reference implementation */
    Float u = h<8 ? dx : dy;
    Float v = h<4 ? dy : h==12 || h==14 ? dx : dz;
#elif defined(GRAD_PBRT)
    /* PBRT's implementation uses the hashes somewhat
       differently. Possibly, this is just a typo */
    Float u = h<8 || h==12 || h==13 ? dx : dy;
    Float v = h<4 || h==12 || h==13 ? dy : dz;
#endif
    return ((h&1) ? -u : u) + ((h&2) ? -v : v);
}

inline static Float noiseWeight(Float t) {
    Float t3 = t*t*t, t4 = t3*t, t5 = t4*t;
    return 6.0f*t5 - 15.0f*t4 + 10.0f*t3;
}

Float Noise::perlinNoise(const Point &p) {
    // Compute noise cell coordinates and offsets
    int ix = math::floorToInt(p.x),
        iy = math::floorToInt(p.y),
        iz = math::floorToInt(p.z);

    Float dx = p.x - ix,
          dy = p.y - iy,
          dz = p.z - iz;

    // Compute gradient weights
    ix &= (NOISE_PERM_SIZE-1);
    iy &= (NOISE_PERM_SIZE-1);
    iz &= (NOISE_PERM_SIZE-1);

    Float w000 = grad(ix,   iy,   iz,   dx,   dy,   dz);
    Float w100 = grad(ix+1, iy,   iz,   dx-1, dy,   dz);
    Float w010 = grad(ix,   iy+1, iz,   dx,   dy-1, dz);
    Float w110 = grad(ix+1, iy+1, iz,   dx-1, dy-1, dz);
    Float w001 = grad(ix,   iy,   iz+1, dx,   dy,   dz-1);
    Float w101 = grad(ix+1, iy,   iz+1, dx-1, dy,   dz-1);
    Float w011 = grad(ix,   iy+1, iz+1, dx,   dy-1, dz-1);
    Float w111 = grad(ix+1, iy+1, iz+1, dx-1, dy-1, dz-1);

    // Compute trilinear interpolation of weights
    Float wx = noiseWeight(dx),
          wy = noiseWeight(dy),
          wz = noiseWeight(dz);

    Float x00 = math::lerp(wx, w000, w100),
          x10 = math::lerp(wx, w010, w110),
          x01 = math::lerp(wx, w001, w101),
          x11 = math::lerp(wx, w011, w111),
          y0 = math::lerp(wy, x00, x10),
          y1 = math::lerp(wy, x01, x11);

    return math::lerp(wz, y0, y1);
}

Float Noise::fbm(const Point &p, const Vector &dpdx,
        const Vector &dpdy, Float omega, int maxOctaves) {
    // Compute number of octaves for antialiased FBm
    Float s2 = std::max(dpdx.lengthSquared(), dpdy.lengthSquared());
    Float foctaves = std::min((Float) maxOctaves, 1.f - .5f * math::log2(s2));
    int octaves = (int) foctaves;

    // Compute sum of octaves of noise for FBm
    Float sum = 0., lambda = 1., o = 1.;
    for (int i = 0; i < octaves; ++i) {
        sum += o * perlinNoise(lambda * p);
        lambda *= 1.99f;
        o *= omega;
    }
    Float partialOctave = foctaves - octaves;
    sum += o * math::smoothStep((Float) .3f, (Float) .7f, partialOctave)
        * perlinNoise(lambda * p);
    return sum;
}

Float Noise::turbulence(const Point &p, const Vector &dpdx,
        const Vector &dpdy, Float omega, int maxOctaves) {
    // Compute number of octaves for antialiased FBm
    Float s2 = std::max(dpdx.lengthSquared(), dpdy.lengthSquared());
    Float foctaves = std::min((Float) maxOctaves, 1.f - .5f * math::log2(s2));
    int octaves = (int) foctaves;

    // Compute sum of octaves of noise for turbulence
    Float sum = 0., lambda = 1., o = 1.;
    for (int i = 0; i < octaves; ++i) {
        sum += o * std::abs(perlinNoise(lambda * p));
        lambda *= 1.99f;
        o *= omega;
    }
    Float partialOctave = foctaves - octaves;
    sum += o * math::smoothStep((Float) .3f, (Float) .7f, partialOctave)
        * std::abs(perlinNoise(lambda * p));
    return sum;
}

MTS_NAMESPACE_END

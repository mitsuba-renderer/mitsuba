# Surface scattering models
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

macro(add_bsdf)
  add_mts_plugin(${ARGN} MTS_HW TYPE bsdf)
endmacro()

# Basic library of smooth and rough materials
add_bsdf(diffuse         diffuse.cpp)
add_bsdf(dielectric      dielectric.cpp ior.h)
add_bsdf(conductor       conductor.cpp)
add_bsdf(plastic         plastic.cpp ior.h)
add_bsdf(roughdiffuse    roughdiffuse.cpp)
add_bsdf(roughdielectric roughdielectric.cpp microfacet.h ior.h)
add_bsdf(roughconductor  roughconductor.cpp microfacet.h)
add_bsdf(roughplastic    roughplastic.cpp microfacet.h ior.h)

# Materials that act as modifiers
add_bsdf(bumpmap      bumpmap.cpp)
add_bsdf(normalmap    normalmap.cpp)
add_bsdf(twosided     twosided.cpp)
add_bsdf(mask         mask.cpp)
add_bsdf(mixturebsdf  mixturebsdf.cpp)
add_bsdf(blendbsdf    blendbsdf.cpp)
add_bsdf(coating      coating.cpp ior.h)
add_bsdf(roughcoating roughcoating.cpp ior.h)

# Other materials
add_bsdf(ward       ward.cpp)
add_bsdf(phong      phong.cpp)
add_bsdf(difftrans  difftrans.cpp)
add_bsdf(hk         hk.cpp)
add_bsdf(null       null.cpp)
add_bsdf(thindielectric thindielectric.cpp)


if (BOOST_SPIRIT_WORKS)
  add_bsdf(irawan   irawan.h irawan.cpp)
endif ()

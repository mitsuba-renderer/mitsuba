# Phase functions
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

macro(add_phase)
  add_mts_plugin(${ARGN} TYPE phase)
endmacro()

add_phase(isotropic    isotropic.cpp)
add_phase(hg           hg.cpp)
add_phase(rayleigh     rayleigh.cpp)
add_phase(kkay         kkay.cpp)
add_phase(microflake   microflake.cpp microflake_fiber.h)
add_phase(mixturephase mixturephase)

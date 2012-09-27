# Sub-surface integrators
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

macro(add_subsurface)
  add_mts_plugin(${ARGN} TYPE subsurface)
endmacro()

add_subsurface(dipole    dipole.cpp
  irrproc.h irrproc.cpp irrtree.h irrtree.cpp bluenoise.h bluenoise.cpp)

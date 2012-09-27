# Participating media
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

macro(add_medium)
  add_mts_plugin(${ARGN} TYPE medium)
endmacro()

add_medium(homogeneous   homogeneous.cpp maxexp.h)
add_medium(heterogeneous heterogeneous.cpp)

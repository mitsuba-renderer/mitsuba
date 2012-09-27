# Reconstruction filters
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

macro(add_rfilter)
  add_mts_plugin(${ARGN} TYPE rfilter)
endmacro()

add_rfilter(box        box.cpp)
add_rfilter(tent       tent.cpp)
add_rfilter(lanczos    lanczos.cpp)
add_rfilter(mitchell   mitchell.cpp)
add_rfilter(catmullrom catmullrom.cpp)
add_rfilter(gaussian   gaussian.cpp)

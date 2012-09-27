# Volumetric data sources
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

macro(add_volume)
  add_mts_plugin(${ARGN} TYPE volume)
endmacro()

add_volume(constvolume constvolume.cpp)
add_volume(gridvolume  gridvolume.cpp)
add_volume(hgridvolume hgridvolume.cpp)
add_volume(volcache    volcache.cpp)

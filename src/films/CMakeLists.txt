# Film implementations
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

macro(add_film)
  add_mts_plugin(${ARGN} TYPE film)
endmacro()

add_film(mfilm   mfilm.cpp cnpy.h cnpy.cpp)
add_film(ldrfilm ldrfilm.cpp annotations.h banner.h MTS_HW)
add_film(hdrfilm hdrfilm.cpp annotations.h banner.h MTS_HW)

if (OPENEXR_FOUND)
  include_directories(${ILMBASE_INCLUDE_DIRS} ${OPENEXR_INCLUDE_DIRS})
  add_film(tiledhdrfilm tiledhdrfilm.cpp)
endif()

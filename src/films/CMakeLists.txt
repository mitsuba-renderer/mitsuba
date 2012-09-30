# Film implementations
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

macro(add_film)
  add_mts_plugin(${ARGN} TYPE film)
endmacro()

add_film(mfilm   mfilm.cpp)
add_film(ldrfilm ldrfilm.cpp banner.h)
add_film(hdrfilm hdrfilm.cpp banner.h)

if (OPENEXR_FOUND)
  include_directories(${ILMBASE_INCLUDE_DIRS} ${OPENEXR_INCLUDE_DIRS})
  add_film(tiledhdrfilm tiledhdrfilm.cpp)
endif()

# Light sources
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

macro(add_emitter)
  add_mts_plugin(${ARGN} TYPE emitter)
endmacro()

add_emitter(area        area.cpp)
add_emitter(collimated  collimated.cpp)
add_emitter(directional directional.cpp)
add_emitter(point       point.cpp)
add_emitter(constant    constant.cpp)
add_emitter(envmap      envmap.cpp MTS_HW)
add_emitter(sky         sky.cpp sunsky/skymodel.cpp sunsky/skymodel.h sunsky/skymodeldata.h sunsky/sunmodel.h)
add_emitter(sun         sun.cpp sunsky/sunmodel.h)
add_emitter(sunsky      sunsky.cpp sunsky/sunmodel.h)
add_emitter(spot        spot.cpp   MTS_HW)

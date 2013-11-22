# Cameras
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

macro(add_sensor)
  add_mts_plugin(${ARGN} TYPE sensor)
endmacro()

add_sensor(radiancemeter      radiancemeter.cpp)
add_sensor(irradiancemeter    irradiancemeter.cpp)
add_sensor(fluencemeter       fluencemeter.cpp)
add_sensor(perspective        perspective.cpp)
add_sensor(perspective_rdist  perspective_rdist.cpp)
add_sensor(thinlens           thinlens.cpp)
add_sensor(orthographic       orthographic.cpp)
add_sensor(telecentric        telecentric.cpp)
add_sensor(spherical          spherical.cpp)

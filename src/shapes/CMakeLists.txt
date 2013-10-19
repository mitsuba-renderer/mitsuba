# Intersection shapes
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

macro(add_shape)
  add_mts_plugin(${ARGN} TYPE shape)
endmacro()

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

add_shape(obj        obj.cpp MTS_HW)
add_shape(serialized serialized.cpp)
add_shape(rectangle  rectangle.cpp)
add_shape(disk       disk.cpp)
add_shape(sphere     sphere.cpp)
add_shape(cylinder   cylinder.cpp)
add_shape(cube       cube.cpp)
add_shape(hair       hair.h hair.cpp)
add_shape(shapegroup shapegroup.h shapegroup.cpp)
add_shape(instance   instance.h instance.cpp)
add_shape(heightfield heightfield.cpp)
#add_shape(deformable deformable.cpp)
add_shape(ply ply.cpp ply/ply_parser.cpp 
  ply/byte_order.hpp ply/config.hpp ply/io_operators.hpp
  ply/ply.hpp ply/ply_parser.hpp)

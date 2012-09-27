# Texture types
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

macro(add_texture)
  add_mts_plugin(${ARGN} TYPE texture)
endmacro()

include_directories(${ILMBASE_INCLUDE_DIRS})

add_texture(bitmap       bitmap.cpp MTS_HW)
add_texture(scale        scale.cpp  MTS_HW)
add_texture(gridtexture  gridtexture.cpp)
add_texture(checkerboard checkerboard.cpp)
add_texture(vertexcolors vertexcolors.cpp)
add_texture(wireframe    wireframe.cpp)
add_texture(curvature    curvature.cpp)

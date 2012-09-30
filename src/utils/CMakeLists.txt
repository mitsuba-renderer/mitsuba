# Utilities
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

macro(add_utility)
  add_mts_plugin(${ARGN} TYPE utility)
endmacro()

include_directories(${ILMBASE_INCLUDE_DIRS})

add_utility(addimages      addimages.cpp)
add_utility(joinrgb        joinrgb.cpp)
add_utility(cylclip        cylclip.cpp MTS_HW)
add_utility(kdbench        kdbench.cpp)
add_utility(tonemap        tonemap.cpp)
#add_utility(rdielprec      rdielprec.cpp)

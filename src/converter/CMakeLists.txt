# COLLADA/OBJ importer (mtsimport)
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

include_directories(${COLLADA_INCLUDE_DIRS} ${XERCES_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR})

set(SRCS
  collada.cpp
  converter.h
  converter.cpp
  obj.cpp
)

# COLLADA-DOM 2.4 introduces namespaces, which are not source-compatible
if (COLLADA_NAMESPACE)
  add_definitions(${COLLADA_DEFINITIONS} -DCOLLADA_DOM_2_4)
endif()

# Create a static library which is also used by mtsgui. Adding the dependent
# libraries will allow transitive linking.
add_library(mtsconverter_lib STATIC ${SRCS})
target_link_libraries(mtsconverter_lib
	${COLLADA_LIBRARIES} ${XERCES_LIBRARIES} ${OPENGL_glu_LIBRARY})

add_mts_exe(mtsimport mtsimport.cpp ${SRCS}
  ${MTS_DARWIN_STUB} ${MTS_WINDOWS_STUB} MTS_HW LINK_LIBRARIES mtsconverter_lib
  RES_DESCRIPTION "Mitsuba COLLADA 1.4 & Wavefront OBJ importer")
